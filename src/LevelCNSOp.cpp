#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file LevelCNSOp.cpp
 *
 * \brief Member functions for LevelCNSOp
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "AMRLevel.H"
#include "MultiBlockCoordSys.H"
#include "LevelGridMetrics.H"
#include "UnitNormalsF_F.H"
#include "SPMD.H"
#include "LevelFluxRegister.H"
#include "MOLUtilities.H"
#include "FourthOrderUtil.H"
#include "BlockRegister.H"
#include "CHMatrixOps.H"
#include "NonLinearSolverLapack.H"

//----- Internal -----//

#include "LevelCNSOp.H"
#include "LevelMappedFunc.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CRDutil.H"
#include "DataTemp.H"
#include "DebugOut.H"
#include "LGintegrator.H"
#include "ARKUtil.H"
#ifdef CH_CTHR
#include "ThreadTeamArchitectChombo.H"
#endif


/*******************************************************************************
 *
 * Class LevelCNSOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_boxes The disjoint box layout
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_U     Conservative state in physical space.  This is
 *                      often set in this class.  Think of this more
 *                      as data space rather than in/out.
 *  \param[in]  a_UExchangeCopier
 *                      Copier for \<U\>
 *  \param[in]  a_JUExchangeCopier
 *                      Copier for \<JU\>
 *  \param[in]  a_timeInterpolator
 *                      Used to set coarsened-fine data at an RK4
 *                      stage
 *  \param[in]  a_hasCoarserGrid
 *                      T - A coarser level exists and is in use
 *  \param[in]  a_hasFinerGrid
 *                      T - A finer level exists and is in use
 *//*-----------------------------------------------------------------*/

LevelCNSOp::LevelCNSOp(const DisjointBoxLayout& a_boxes,
                       LevelGridMetrics&        a_levelGridMetrics,
                       LevelData<FArrayBox>&    a_U,
                       const Copier&            a_UExchangeCopier,
                       const Copier&            a_JUExchangeCopier,
                       const bool&              a_hasCoarserGrid,
                       const bool&              a_hasFinerGrid,
                       BoxPartCache&            a_boxPartCache)
  :
  m_boxes(a_boxes),
  m_levelGridMetrics(a_levelGridMetrics),
  m_U(a_U),
  m_UExchangeCopier(a_UExchangeCopier),
  m_JUExchangeCopier(a_JUExchangeCopier),
  m_timeInterpolator(nullptr),
  m_hasCoarserGrid(a_hasCoarserGrid),
  m_hasFinerGrid(a_hasFinerGrid),
  m_boxPartCache(a_boxPartCache),
  m_patchOp(a_levelGridMetrics),
  m_MMBSingleLevelOp(a_levelGridMetrics),
  m_minConvDt(1.E9),
  m_minDiffDt(1.E9),
  m_minChemDt(1.E9),
  m_minLocalDt(1.E9),
  m_minConvDtCell(IntVect_zero),
  m_minDiffDtCell(IntVect_zero),
  m_minChemDtCell(IntVect_zero),
  m_lockMinLocalDt(ATOMIC_FLAG_INIT),
  m_firstDt(true),
  m_prevDt(0.),
  m_unitNormalsDefined(false)
{
  m_nlsolver = std::make_unique<NonLinearSolverLapack>();
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::define(const int       a_level,
                   const RealVect& a_dx)
{
  CRD::msg << CRD::fv3 << "LevelCNSOp::define (level: "
           << a_level << ")" << CRD::end;
  m_dx = a_dx;
  m_level = a_level;

  // Define intermediate data for a stage (there is an argument for recomputing
  // some of this data but for now, we store all that is required)
  {
    m_WcellPntLvl.define(
      m_boxes,
      CRDparam::g_CRDPhysics->numPrimitive(),
      CRDparam::queryGhosts(CRDparam::NumGhostWcellPnt)*IntVect_unit);
    m_WcellAvgLvl.define(
      m_boxes,
      CRDparam::g_CRDPhysics->numPrimitive(),
      CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg)*IntVect_unit);
    // Where a Riemann solve is required
    //**FIXME Remove TngBdry when new data structure is available
    const int numGhostWpFace =
      std::max(
        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
        std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
                 CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan)));
    m_WfaceAvgLvl.define(
      m_boxes,
      CRDparam::g_CRDPhysics->numPrimitive(),
      numGhostWpFace*IntVect_unit);
  }

  // Define m_flattening.  It is always defined with at least 1 component
  {
    const int numGhost = std::max(
      CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
      std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
               CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan)));
    int numComp;
    if (CRDparam::g_useFlattening)
      {
        numComp = 1 + c_numDebugPlotVar;
      }
    else
      {
        numComp = std::max(1, (int)c_numDebugPlotVar);
      }
    m_flattening.define(m_boxes, numComp, numGhost*IntVect::Unit);
  }

  // Define m_faceAvgTimeAvgPlotData
  //**NOTE: It would be good to move this to a plotting class to clean up
  //        the code and help in-situ processing a little bit more
  //**NOTE: For now, m_faceAvgTimeAvgPlotData and m_faceAvgPlotData have the
  //        same number of components. So, both sets will plot the combination
  //        of both inputs (this can be fixed later)
  {
    const int numGhost = 0;
    int numComp = 0;
    int plotVariables = CRDparam::g_plotLoFaceAvgComps;
    plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
    if (plotVariables & CRDparam::PlotPrimitive)
      {
        numComp += CRDparam::g_CRDPhysics->numPrimitive();
      }
    if (plotVariables & CRDparam::PlotFluxes)
      {
        if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
          {
            numComp += CRDparam::g_CRDPhysics->numFluxes();
          }
        if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
          {
            numComp += CRDparam::g_CRDPhysics->numFluxes();
          }
      }
    if (plotVariables & CRDparam::PlotTurbulentComps)
      {
        // This list currently includes
        // (a) <u_i*u_j> = SpaceDim*SpaceDim
        // (b) <p*p> = 1
        // (c) <T*T> = 1
        // (d) <rho*rho> = 1
        // (e) <tau_sgs_ij> = SpaceDim*SpaceDim
        // (f) <mag_tau_w> = 1
        // (g) <tau_w_vector> (wall-normal space) = SpaceDim
        // (h) <mag_tau_w_sqrd> = 1
        // This leads to 5 + SpaceDim
        // + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim) components
        numComp += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
      }
    numComp = std::max(numComp, 1); // Always needs at least 1 component
    m_faceAvgTimeAvgPlotData.define(m_boxes, numComp, numGhost*IntVect_unit);
  }
  // Define m_faceAvgPlotData
  {
    const int numGhost = 0;
    int numComp = 0;
    int plotVariables = CRDparam::g_plotLoFaceAvgComps;
    plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
    if (plotVariables & CRDparam::PlotPrimitive)
      {
        numComp += CRDparam::g_CRDPhysics->numPrimitive();
      }
    if (plotVariables & CRDparam::PlotFluxes)
      {
        if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
          {
            numComp += CRDparam::g_CRDPhysics->numFluxes();
          }
        if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
          {
            numComp += CRDparam::g_CRDPhysics->numFluxes();
          }
      }
    if (plotVariables & CRDparam::PlotTurbulentComps)
      {
        // This list currently includes
        // (a) <u_i*u_j> = SpaceDim*SpaceDim
        // (b) <p*p> = 1
        // (c) <T*T> = 1
        // (d) <rho*rho> = 1
        // (e) <tau_sgs_ij> = SpaceDim*SpaceDim
        // (f) <mag_tau_w> = 1
        // (g) <tau_w_vector> (wall-normal space) = SpaceDim
        // (h) <mag_tau_w_sqrd> = 1
        // This leads to 5 + SpaceDim
        // + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim) components
        numComp += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
      }
    numComp = std::max(numComp, 1); // Always needs at least 1 component
    m_faceAvgPlotData.define(m_boxes, numComp, numGhost*IntVect_unit);
  }
  // Define m_timeAvgData
  {
    const int numGhost = 0;
    int numComp = 1; // always need at least 1 component
    if (CRDparam::g_plotTimeAvgTurb)
      {
        // If used, shear-stress from extrapolation goes in first comp above
        if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
          {
            numComp += 1; // Wall shear-stress from LES model
          }
        // Time-averaged velocity on the low-face of each cell in all directions
        // These are face-averaged values
        numComp += SpaceDim*SpaceDim;
        // Time-averaged velocity on the low-face of each cell in all directions
        // These are face-centered values
        numComp += SpaceDim*SpaceDim;
        // Time-avg velocity stress tensor on low-face of each cell in all dirs
        // These are face-averaged values
        numComp += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        // Time-avg velocity stress tensor on low-face of each cell in all dirs
        // These are face-centered values
        numComp += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);

        // Time-avg density and temperature in the cells
        numComp += 1 + 1;

        // New components
        numComp += 2; // rho*u and rho*u*u
        numComp += 1; // <tau_{wall-model}^2>
        // Time-avg <SGS> stress tensor on x-faces
        numComp += (SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        numComp += (SpaceDim-1); // Actual wall shear-stress
                                 // Momentum flux on wall face -- only use the
                                 // wall-tangential components
      }
    m_timeAvgData.define(m_boxes, numComp, numGhost*IntVect_unit);
  }
  // Define total-time after initialization of sharp time-averaging filter
  m_totalTimeAfterFilterInit = 0;

  // Define m_sgske. This plots the current, instantaneous SGS KE estimate
  m_sgske.define(m_boxes, 1, IntVect::Zero);
  m_Jsgske.define(m_boxes, 1, IntVect::Zero);

  // Define data structures for coarsened SGS KE calculation
  // NOTE: this is only for single-level single-block problems
  defineSGSKE();

  const int numGhostJU = 1; //**FIXME: this should probably be fine
  if (m_hasCoarserGrid && (CRDparam::g_useSGSCoarsening) &&
      (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex))
    {
      m_currModelJU.define(m_boxes,
                           CRDparam::g_CRDPhysics->numConservative(),
                           numGhostJU*IntVect::Unit);
    }

  // Define special operators as required
#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
  if (CRDparam::g_explicitFilterType == CRDparam::ExplicitFilterSpectral)
    {
      // This expects Cartesian grids.  Special considerations are needed for
      // multiblock.
      CH_assert(!m_levelGridMetrics.isMultiBlock());
      const ProblemDomain& problemDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(0);
      CRD::msg << CRD::fh1 << CRD::verb(SpectralFilter::c_verbosity)
               << "Initializing FFTW" << CRD::end;
      m_spectralFilter.define(problemDomain,
                              m_boxes,
                              CRDparam::g_CRDPhysics->numConservative(),
                              CRDparam::g_spectralFilterProfile);
      m_spectralFilter.setFilterParam(CRDparam::g_spectralFilterParam);
    }
  if (m_level == 0 &&
      CRDparam::g_turbForcingType == CRDparam::TurbForcingSpectral)
    {
      // This expects Cartesian grids.  Special considerations are needed for
      // multiblock.
      CH_assert(!m_levelGridMetrics.isMultiBlock());
      const ProblemDomain& problemDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(0);
      CRD::msg << CRD::fh1 << CRD::verb(SpectralForcing::c_verbosity)
               << "Initializing FFTW for spectral forcing" << CRD::end;
      m_spectralForcing.define(problemDomain, m_boxes);
    }
#endif
  // Define m_bndryCellData and m_bndryNtJ
  m_bndryCellData.define(m_boxes);
  m_bndryNtJ.define(m_boxes);
  int numComps = 2*(CRDparam::g_CRDPhysics->numPrimitive());
  // We only need values stored right along boundaries
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      Box ghostBox = disjointBox;
      ghostBox.grow(
        std::max(
          std::max(CRDparam::queryGhosts(CRDparam::NumGhostWcellExtTngBdry),
                   CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry)),
          std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan),
                   CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm))));
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box bndryFaceBox = disjointBox;
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          // Low side info
          CRDparam::g_CNSIBC->getBoundaryFaces(
            bndryFaceBox, ghostBox, blockDomain, dir, Side::Lo);
          if (!bndryFaceBox.isEmpty())
            {
              Box bndryBox = bndryFaceBox; // Store for bndryNtJ
              bndryFaceBox.grow(dir, 1); // Grow on both low and high sides
              bndryFaceBox.enclosedCells(dir); // Make this a cell box
              int indx = MultiBlockRegions::FaceTag::indexFace(dir,Side::Lo);
              m_bndryCellData[dit][indx].define(bndryFaceBox, numComps);
              // Shift into the domain for wall-model evaluation use
              bndryBox.shift(dir, 1);
              m_bndryNtJ[dit][indx].define(bndryBox, SpaceDim*SpaceDim);
            }
          // High side info
          CRDparam::g_CNSIBC->getBoundaryFaces(
            bndryFaceBox, ghostBox, blockDomain, dir, Side::Hi);
          if (!bndryFaceBox.isEmpty())
            {
              Box bndryBox = bndryFaceBox; // Store for bndryNtJ
              bndryFaceBox.grow(dir, 1); // Grow on both low and high sides
              bndryFaceBox.enclosedCells(dir); // Make this a cell box
              int indx = MultiBlockRegions::FaceTag::indexFace(dir,Side::Hi);
              m_bndryCellData[dit][indx].define(bndryFaceBox, numComps);
              // Shift into the domain for wall-model evaluation use
              bndryBox.shift(dir, -1);
              m_bndryNtJ[dit][indx].define(bndryBox, SpaceDim*SpaceDim);
            }
        }
    }
  // Reset m_prevDt if we've called all of these operators
  m_prevDt = 0.;
  m_totalTimeAfterFilterInit = 0;
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction) in reverse traversal
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::reverseDefine()
{
  CRD::msg << CRD::fv3 << "LevelCNSOp::reverseDefine (level: "
           << m_level << ")" << CRD::end;

  //**FIXME Unit normals aren't always used, might be beneficial to add
  //        a condition here to do this define or not
  //**FIXME Wish there wasn't a separate define for this.
  // Define the unit normals required for the Riemann solver
  defineUnitNormals(m_levelGridMetrics.m_N);

  // Single-level MMB operations
  if (m_levelGridMetrics.isMultiBlock())
    {
      m_MMBSingleLevelOp.define(
        std::max(CRDparam::g_CRDPhysics->numFluxes(),
                 CRDparam::g_CRDPhysics->numPrimitive()));
    }

  // Allocate data structures for recirculating-inlet flat-plate TBL
  CRDparam::g_CNSIBC->initializeInletDataStructures(
    m_U, m_level, m_boxes, m_levelGridMetrics, m_hasFinerGrid);
}

//**FIXME not following code notation conventions
/*
 * Sets the time interpolator for this level.
 *
 * Since we don't know at construction which time integrator we are using,
 * (i.e. whether RK4 or ARK4), we need to specify at a later time what the
 * interpolation method is. All interpolation methods are currently derived
 * from TimeInterpolatorRK4, so we just accept a pointer of that type and
 * dynamically cast it before calling some ::eval function.
 */
void
LevelCNSOp::setTimeInterpolatorPtr(
  const TimeInterpolatorRK4* a_timeInterpolator)
{
  CH_TIME("LevelCNSOp::setTimeInterpolatorPtr");
  CH_assert(a_timeInterpolator != nullptr);
  m_timeInterpolator = a_timeInterpolator;
  CH_assert(m_timeInterpolator != nullptr); // Make sure it set correctly...
}

/*--------------------------------------------------------------------*/
//  Find unit normals for applying the Riemann problem on mapped grids
/** \param[in]  a_NLev  Metrics N on the faces
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::defineUnitNormals(LevelData<FluxBox>& a_NLev)
{
  CH_TIME("LevelCNSOp::defineUnitNormals");
  CRD::msg << CRD::fv3 << "LevelCNSOp::defineUnitNormals (level: "
           << m_level << ")" << CRD::end;

  m_unitNormals.define(m_boxes);
  const int numInteriorGhosts =
    std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
             CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan));
  const int numGhosts =
    std::max(numInteriorGhosts,
             CRDparam::queryGhosts(CRDparam::NumGhostWcellExtTngBdry));
  CH_assert(numGhosts <= m_levelGridMetrics.getNReqGhostVect()[0]);

  // Define m_facePntDeltaC: DEBUG ONLY -- just for code speedup
  m_facePntDeltaC.define(m_boxes);
  // Define m_faceCoord: DEBUG ONLY -- just for code speedup
  m_faceCoord.define(m_boxes);

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      // The unit normals are required anywhere a Riemann problem is solved
      // and on boundary faces.  The Riemann problem is solved on
      // NumGhostWpFaceAvg ghosts and the boundary faces are needed out to
      // NumGhostWcellTngBdry ghosts
      Box box = grow(disjointBox, numGhosts);
      if (blockDomain.contains(box))
        {
          box = grow(disjointBox, numInteriorGhosts);
        }
      box &= blockDomain;

      // Only for m_facePntDeltaC: DEBUG ONLY -- just for code speedup
      Box box1Dom = grow(disjointBox, 1);
      box1Dom &= blockDomain;
      FluxBox& deltaCFxb = m_facePntDeltaC[dit];
      deltaCFxb.define(box1Dom, 1);
      FLUXBOXSTACKTEMP(xiFxb, box1Dom, SpaceDim);
      const BlockCoordSys& blockCoordSys =
        *(m_levelGridMetrics.getCoordSys(disjointBox));
      const RealVect dxVect = m_levelGridMetrics.getCoordSys(disjointBox)->dx();
      const Real cellVol = dxVect.product();

      // Only for m_faceCoord: DEBUG ONLY -- just for code speedup
      FluxBox& faceCoord = m_faceCoord[dit];
      Box box2Dom = grow(disjointBox, 2);
      box2Dom &= blockDomain;
      Box faceCoordBoxDom =
        grow(disjointBox, (1+CRDparam::g_sgskeCrsLevFilterRatio));
      faceCoordBoxDom &= blockDomain;
      faceCoord.define(faceCoordBoxDom, SpaceDim);
      FLUXBOXSTACKTEMP(xiFxb2, faceCoordBoxDom, SpaceDim);

      FluxBox& unitNormalFxb = m_unitNormals[dit];
      unitNormalFxb.define(box, SpaceDim*SpaceDim);
      const FluxBox &NFxb = a_NLev[dit];
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const IntVect& metricTermComponentsDir =
            BlockCoordSys::metricsCompDir(dir);
          FArrayBox& unitNormalFab = unitNormalFxb[dir];
          const FArrayBox& NFab = NFxb[dir];
          CH_assert(NFab.contains(unitNormalFab.box()));
          FORT_GETUNITNORMALS(CHF_FRA(unitNormalFab),
                              CHF_CONST_FRA(NFab),
                              CHF_CONST_INTVECT(metricTermComponentsDir),
                              CHF_CONST_INT(dir),
                              CHF_BOX(unitNormalFab.box()));

          // Only for m_facePntDeltaC: DEBUG ONLY -- just for code speedup
          Box faceBox = box1Dom;
          faceBox.surroundingNodes(dir);
          CRDparam::g_CNSIBC->getFaceCompCoordinates(faceBox,
                                                     xiFxb[dir],
                                                     dir,
                                                     blockCoordSys);
          const FArrayBox& arrXi = xiFxb[dir];
          const FArrayBox& arrDC = deltaCFxb[dir];
          const Real invDim = 1./SpaceDim;
          MD_BOXLOOP(faceBox, i)
            {
              RealVect xiLoc(D_DECL(arrXi[MD_IX(i,0)],
                                    arrXi[MD_IX(i,1)],
                                    arrXi[MD_IX(i,2)]));
              Real faceJ =
                m_levelGridMetrics.getCoordSys(disjointBox)->pointwiseJ(xiLoc);
              arrDC[MD_IX(i,0)] = std::pow((faceJ*cellVol), invDim);
            }
          Box faceBox2 = faceCoordBoxDom;
          faceBox2.surroundingNodes(dir);
          // Only for m_faceCoord: DEBUG ONLY -- just for code speedup
          CRDparam::g_CNSIBC->getFaceCoordinates(
            faceBox2, xiFxb2[dir], faceCoord[dir], dir, blockCoordSys);
        }
      const int sgsRefineRatio = CRDparam::g_sgskeCrsLevFilterRatio;
      if (!m_hasCoarserGrid && (sgsRefineRatio > 1) &&
          CRDparam::g_useSGSCoarsening && !m_levelGridMetrics.isMultiBlock())
        {
          // Fill m_cellPntJ
          Box boxJDom = grow(disjointBox, 1);
          boxJDom &= blockDomain;
          Box pntJBox = boxJDom;
          FArrayBox& cellPntJFab = m_cellPntJ[dit];
          cellPntJFab.setVal(0.);
          MD_BOXLOOP(pntJBox, i)
            {
              IntVect iv(D_DECL(i0,i1,i2));
              // Get mapped location
              RealVect xiLoc = dxVect*(iv + 0.5*IntVect_unit);
              // Get cell-centered J
              Real cellPntJ =
                m_levelGridMetrics.getCoordSys(disjointBox)->pointwiseJ(xiLoc);
              cellPntJFab[MD_IX(i, 0)] = cellPntJ;
            }
          // Fill m_crsCellPntJ
          const Box& crsDisjointBox = m_crsCellPntJ.getBoxes()[dit];
          Box crsJBox = grow(crsDisjointBox, 2);
          // Coarsen the problem domain
          const ProblemDomain& problemDomain =
            m_levelGridMetrics.getCoordSys().problemDomain(0);
          ProblemDomain coarseProbDom = problemDomain;
          coarseProbDom.coarsen(sgsRefineRatio);
          // Restrict crsJBox to the coarsened domain
          crsJBox &= coarseProbDom;
          FArrayBox& crsCellPntJ = m_crsCellPntJ[dit];
          FArrayBox& XFab = m_crsXLvl[dit];
          MD_BOXLOOP(crsJBox, i)
            {
              IntVect iv(D_DECL(i0,i1,i2));
              // Get mapped location
              RealVect xiLoc = sgsRefineRatio*dxVect*(iv + 0.5*IntVect_unit);
              // Get cell-centered J
              Real cellPntJ =
                m_levelGridMetrics.getCoordSys(disjointBox)->pointwiseJ(xiLoc);
              crsCellPntJ[MD_IX(i, 0)] = cellPntJ;
              // Get cell-centered coordinates
              RealVect xLoc =
                m_levelGridMetrics.getCoordSys(disjointBox)->realCoord(xiLoc);
              D_TERM(XFab[MD_IX(i, 0)] = xLoc[0];,
                     XFab[MD_IX(i, 1)] = xLoc[1];,
                     XFab[MD_IX(i, 2)] = xLoc[2];);
            }
          // Compute m_crsDeltaC
          FArrayBox& deltaCFab = m_crsDeltaC[dit];
          Box crsDeltaCBox = grow(crsDisjointBox, 1);
          crsDeltaCBox &= coarseProbDom;
          RealVect crsDxVect = sgsRefineRatio*dxVect;
          const Real crsCellVol = crsDxVect.product();
          MD_BOXLOOP(crsDeltaCBox, i)
            {
              Real cellPntJ = crsCellPntJ[MD_IX(i, 0)];
              deltaCFab[MD_IX(i,0)] = std::pow((cellPntJ*crsCellVol), 1./3.);
            }
        }
    }
  // Compute face-centered NtJ on near-bndry faces for LES wall-model
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box bndryFaceBox = disjointBox;
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          // Low side info
          CRDparam::g_CNSIBC->getBoundaryFaces(
            bndryFaceBox, disjointBox, blockDomain, dir, Side::Lo);
          if (!bndryFaceBox.isEmpty())
            {
              FArrayBox& NtJ = m_bndryNtJ[dit][
                MultiBlockRegions::FaceTag::indexFace(dir,Side::Lo)];
              RealVect offsetVect = 0.5*RealVect_unit;
              offsetVect[dir] = 0.; // Set the offset vector to face centers
              m_levelGridMetrics.solveNtJ(NtJ, disjointBox, offsetVect);
            }
          // High side info
          CRDparam::g_CNSIBC->getBoundaryFaces(
            bndryFaceBox, disjointBox, blockDomain, dir, Side::Hi);
          if (!bndryFaceBox.isEmpty())
            {
              FArrayBox& NtJ = m_bndryNtJ[dit][
                MultiBlockRegions::FaceTag::indexFace(dir,Side::Hi)];
              RealVect offsetVect = 0.5*RealVect_unit;
              offsetVect[dir] = 0.; // Set the offset vector to face centers
              m_levelGridMetrics.solveNtJ(NtJ, disjointBox, offsetVect);
            }
        }
    }
  m_unitNormalsDefined = true;
}

/*--------------------------------------------------------------------*/
//  Define the data structures on the coarsest level needed to compute
//  the SGS kinetic energy estimate
/** \param[in]  a_null  Nothing needed for this function yet
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::defineSGSKE()
{
  CH_TIME("LevelCNSOp::defineSGSKE");
  CRD::msg << CRD::fv3 << "LevelCNSOp::defineSGSKE (level: "
           << m_level << ")" << CRD::end;

  const int sgsRefineRatio = CRDparam::g_sgskeCrsLevFilterRatio;
  if (m_hasCoarserGrid || sgsRefineRatio == 1 ||
      !CRDparam::g_useSGSCoarsening || m_levelGridMetrics.isMultiBlock())
    {
      return;
    }
  const int tensorComps = SpaceDim*SpaceDim;
  // 1)  coarsen this DBL
  DisjointBoxLayout coarseBoxesSGS;
  coarsen(coarseBoxesSGS, m_boxes, sgsRefineRatio);
  // 2)  check the DBL to see if it is good with the number of ghost cells
  //     required by the SGS KE calculation on the coarse grid
  int numCrsJWGhosts    = 3;
  int numCrsSGSKEGhosts = 2;
  Vector<Box> coarseBoxes = coarseBoxesSGS.boxArray();
  for (int i = 0; i != coarseBoxes.size(); ++i)
    {
      if ((coarseBoxes[i]).shortside() < numCrsJWGhosts ||
          (coarseBoxes[i]).shortside() < numCrsSGSKEGhosts)
        {
          CRD::msg << "Base mesh does not support specified coarsening!"
                   << CRD::abort;
        }
    }
  // 3)  create a current-level LD for <JW>
  int numCurrPartJWGhost = 1;
  int numJWcellAvgComps = 1 + SpaceDim;
  m_partJWcellAvgLvl.define(m_boxes,
                            numJWcellAvgComps,
                            numCurrPartJWGhost*IntVect::Unit);
  // 4)  create a coarse LD for <JW> with ghosts needed for computing <SGSKE>
  m_crsCellAvgJW.define(
    coarseBoxesSGS, numJWcellAvgComps, numCrsJWGhosts*IntVect::Unit);
  // 5)  create an averaging operator for use with coarse and fine <JW>
  m_cellAvgJWAvgOp.define(
    m_boxes, coarseBoxesSGS, numJWcellAvgComps, sgsRefineRatio);
  // 6)  create an exchange operator for coarse <JW>
  m_crsCellAvgJWExchange.exchangeDefine(coarseBoxesSGS,
                                        numCrsJWGhosts*IntVect::Unit,
                                        true);
  // 7)  create a coarse LD for cell-centered J
  int numCrsJGhosts = 2;
  m_crsCellPntJ.define(coarseBoxesSGS, 1, numCrsJGhosts*IntVect::Unit);
  // 8)  create a current-level LD for <J grad W N^t/J>
  m_cellAvgJWGradNtJLvl.define(m_boxes, tensorComps, IntVect::Zero);
  // 9)  create a coarse LD for <J grad W N^t/J>
  int numCrsJGradGhosts = 2;
  m_crsCellAvgJPhysGrad.define(coarseBoxesSGS,
                               tensorComps,
                               numCrsJGradGhosts*IntVect::Unit);
  // 10) create an averaging operator for use with fine <J grad W N^t/J>
  m_cellAvgJPhysGradWAvgOp.define(m_boxes,
                                  coarseBoxesSGS,
                                  tensorComps,
                                  sgsRefineRatio);
  // 11) create an exchange operator for coarse <J grad W N^t/J>
  m_crsCellAvgJGradWExchange.exchangeDefine(coarseBoxesSGS,
                                            numCrsJGradGhosts*IntVect::Unit,
                                            true);
  // 12) create a current-level LD for cell-centered J
  int numCellPntJGhosts = 1;
  m_cellPntJ.define(m_boxes, 1, numCellPntJGhosts*IntVect::Unit);
  // 13) create a coarse LD for cell-centered physical-space coordinates
  int numCellPntCoordGhosts = 2;
  m_crsXLvl.define(coarseBoxesSGS,SpaceDim,numCellPntCoordGhosts*IntVect::Unit);
  // 14) create a coarse LD for grid-scale cell-cutoff-length
  int numDeltaCGhosts = 1;
  m_crsDeltaC.define(coarseBoxesSGS, 1, numDeltaCGhosts*IntVect::Unit);
  // 15) create a coarse LD for <J-SGSKE> with needed ghosts for interpolation
  m_crsCellAvgJSGSKE.define(coarseBoxesSGS, 1, numCrsSGSKEGhosts*IntVect::Unit);
  // 16) create an exchange operator for coarse <SGSKE>
  m_crsCellAvgSGSKEExchange.exchangeDefine(coarseBoxesSGS,
                                           numCrsSGSKEGhosts*IntVect::Unit,
                                           true);
  // 17) create a fine LD for <SGSKE> with same number of ghosts as fine m_U
  int numFnCellAvgUGhosts = CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg);
  m_fnCellAvgSGSKE.define(m_boxes, 1, numFnCellAvgUGhosts*IntVect::Unit);
  // 18) create a fine LD for <J-SGSKE>
  int numFnCellAvgJSGSKE = CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg);
  m_fnCellAvgJSGSKE.define(m_boxes, 1, numFnCellAvgJSGSKE*IntVect::Unit);
  // 19) create an interpolator operator for coarse to fine <SGSKE>
  int numInterpComp = 1;
  const ProblemDomain& problemDomain =
    m_levelGridMetrics.getCoordSys().problemDomain(0);
  m_crsToFnSGSKEInterpolator.define(
    m_boxes, numInterpComp, sgsRefineRatio, problemDomain);
  // 20) create an exchange operator for fine <SGSKE>
  m_fnCellAvgSGSKEExchange.exchangeDefine(m_boxes,
                                          numFnCellAvgUGhosts*IntVect::Unit,
                                          true);
  // 21) create an exchange operator for fine <J-SGSKE>
  m_fnCellAvgJSGSKEExchange.exchangeDefine(m_boxes,
                                           numFnCellAvgJSGSKE*IntVect::Unit,
                                           true);
}

/*--------------------------------------------------------------------*/
//  Compute new timestep from latest estimate of the solution
/** \param[in]  a_U     Average solution state from which to compute
 *                      the new time step.  Must be defined on 1 ghost
 *                      cell everywhere except across domain
 *                      boundaries
 *  \param[in]  a_time  Current solution time
 *  \return             The new time step from all physics'
 *                      constraints
 *//*-----------------------------------------------------------------*/

Real
LevelCNSOp::computeNewDt(const LevelData<FArrayBox>& a_U,
                         const Real&                 a_time,
                         const Real&                 a_currentDt,
                         const Real&                 a_prevStepDt) const
{
  CH_TIME("LevelCNSOp::computeNewDt");
  CRD::msg << CRD::fv2 << "LevelCNSOp::computeNewDt (level: " << m_level
           << ")" << CRD::end;

  Real dtNewLocal = m_minLocalDt;
  if (m_firstDt)
    {
      Real dtNewLocalBox = computeFirstDt(a_U, a_time);
      dtNewLocal = std::min(dtNewLocal, dtNewLocalBox);
    }
  CRD::msg << CRD::fv2 << "Local min time step: " << dtNewLocal << CRD::end;

  Real deltat_pid = 1.E+100;
  if(m_UminusUhat.size() == 3
     && CRDparam::g_additiveRK
     && CRDparam::g_ARKUsePIDControl)
    {
      /// Step size control constants
      Real kappa = 0.9;
      Real epsilon = CRDparam::g_ARKPIDControlEps;
      Real orderOfEmbed = 3.;
      Real omega_n = a_currentDt/a_prevStepDt;
      Real ki = 0.25;
      Real kp = 0.14;
      Real kD = 0.1;
      Real alpha = (ki + kp + (2.*omega_n/(1.+omega_n)*kD))/orderOfEmbed;
      Real beta = (kp + 2.*omega_n*kD)/orderOfEmbed;
      Real gamma = ((2.*omega_n*omega_n)/(1.+omega_n))*kD/orderOfEmbed;

      deltat_pid = kappa*a_currentDt
        *std::pow((epsilon/m_UminusUhat[2]),alpha)
        *std::pow((m_UminusUhat[1]/epsilon),beta)
        *std::pow((epsilon/m_UminusUhat[0]),gamma);
    }

//--Global

  struct ReduceLoc
  {
    Real dt;
    int rank;
  };
  ReduceLoc global[4];
  global[0].dt = dtNewLocal;
  global[0].rank = procID();
  global[1].dt = m_minConvDt;
  global[1].rank = procID();
  global[2].dt = m_minDiffDt;
  global[2].rank = procID();
  global[3].dt = m_minChemDt;
  global[3].rank = procID();

  // Real dtNewGlobal = dtNewLocal;
  // Real dtNewConv = m_minConvDt;
  // Real dtNewDiff = m_minDiffDt;
  // Real dtNewChem = m_minChemDt;

#ifdef CH_MPI
  ReduceLoc local[4];
  local[0].dt = dtNewLocal;
  local[0].rank = procID();
  local[1].dt = m_minConvDt;
  local[1].rank = procID();
  local[2].dt = m_minDiffDt;
  local[2].rank = procID();
  local[3].dt = m_minChemDt;
  local[3].rank = procID();
  MPI_Allreduce(local,
                global,
                4,
                MPI_CH_REAL_INT,
                MPI_MINLOC,
                Chombo_MPI::comm);
  MPI_Bcast(m_minConvDtCell.dataPtr(), SpaceDim, MPI_INT, global[1].rank,
            Chombo_MPI::comm);
  MPI_Bcast(m_minDiffDtCell.dataPtr(), SpaceDim, MPI_INT, global[2].rank,
            Chombo_MPI::comm);
  MPI_Bcast(m_minChemDtCell.dataPtr(), SpaceDim, MPI_INT, global[3].rank,
            Chombo_MPI::comm);
#endif
  Real dtNewGlobal = global[0].dt;
  Real dtNewConv   = global[1].dt;
  Real dtNewDiff   = global[2].dt;
  Real dtNewChem   = global[3].dt;
  Real dtNewPid    = deltat_pid;

  // Use the smaller of the global PID time step size and
  // the global physics time step size
  if(CRDparam::g_additiveRK == true)
    {
      dtNewGlobal = std::min(dtNewGlobal, std::max(dtNewPid, dtNewChem));
    }

  m_firstDt = false;
  if (dtNewGlobal > CRDparam::g_maxDt && CRDparam::g_maxDt > 0.)
    {
      dtNewGlobal = CRDparam::g_maxDt;
    }
  if (CRDparam::g_verboseDt)
    {
      CRD::msg << "Time step on level " << m_level << CRD::h2;
      CRD::msg << "Global total time step\n" << dtNewGlobal << CRD::var;
      if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
        {
          CRD::msg << "Global inertial time step\n" << dtNewConv << " at "
                   << m_minConvDtCell << " rank " << global[1].rank << CRD::var;
        }
      if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
        {
          CRD::msg << "Global viscous time step\n" << dtNewDiff << " at "
                   << m_minDiffDtCell << " rank " << global[2].rank
                   << CRD::var;
        }
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          CRD::msg << "Global chemical step\n" << dtNewChem << " at "
                   << m_minChemDtCell << " rank " << global[3].rank
                   << CRD::var;
        }
      if (CRDparam::g_ARKUsePIDControl)
        {
          CRD::msg << "Global PID step\n" << dtNewPid << " at "
                   << m_minPidDtCell << " rank " << 0 << CRD::var;
        }
      //     CRD::msg << CRD::fv3 << "MinimumSpeciesIdx: " << m_cnlim
      //              << " specName: " << CRDparam::g_speciesNames[m_cnlim]
      //              << CRD::end;
    }
  // Reset member data for min dt
  m_minLocalDt = 1.E9;
  m_minConvDt = 1.E9;
  m_minDiffDt = 1.E9;
  m_minChemDt = 1.E9;
  return dtNewGlobal;
}

/*--------------------------------------------------------------------*/
//  Add artificial viscosity to a_JUnew
/** Flux registers are also updated
 *  \param[in]  a_JUnew Solution in conserved state.
 *  \param[out] a_JUnew Contribution of artificial viscosity added.
 *  \param[in]  a_Uold  Solution at start of the time step
 *  \param[out] a_fnFluxRegister
 *                      Register with the next finer level
 *  \param[out] a_crFluxRegister
 *                      Register with the next coarser level
 *  \param[in]  a_dt    Weight for divergence operator
 *  \param[in]  a_time  Current time
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::addArtificialViscosity(
  LevelData<FArrayBox>&       a_JUnew,
  LevelFluxRegister&          a_fnFluxRegister,
  LevelFluxRegister&          a_crFluxRegister,
  const LevelData<FArrayBox>& a_Uold,
  const LevelData<FArrayBox>& a_WOld,
  const Real                  a_dt,
  const Real                  a_time) const
{
  CH_TIME("LevelCNSOp::addArtificialViscosity");
  CRD::msg << CRD::fv3 << "LevelCNSOp::addArtificialViscosity (level: "
           << m_level << ")" << CRD::end;
  CH_assert(CRDparam::g_useArtificialViscosity);

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);

      FArrayBox& JUnewFab           = a_JUnew[dit];
      const FArrayBox& UoldFab      = a_Uold[dit];
      const FArrayBox& WoldFab      = a_WOld[dit];
      const FluxBox& Nfxb           = m_levelGridMetrics.m_N[dit];
      const FArrayBox& Jfab         = m_levelGridMetrics.m_J[dit];
      const FluxBox& unitNormalsFxb = m_unitNormals[dit];

      FLUXBOXSTACKTEMP(artViscNtF, box, CRDparam::g_CRDPhysics->numFluxes());

      // Compute the flux and update JUnew
      m_patchOp.addArtificialViscosity(box,
                                       blockDomain,
                                       JUnewFab,
                                       artViscNtF,
                                       UoldFab,
                                       WoldFab,
                                       Nfxb,
                                       Jfab,
                                       unitNormalsFxb,
                                       m_dx,
                                       a_dt,
                                       a_time,
                                       m_level);

      // Update the flux registers
      m_patchOp.updateFluxRegisters(a_crFluxRegister,
                                    a_fnFluxRegister,
                                    artViscNtF,
                                    box,
                                    dit(),
                                    a_dt,
                                    m_hasCoarserGrid,
                                    m_hasFinerGrid);
    }
}

/*--------------------------------------------------------------------*/
//  Perform species correction
/** FIXME: This might not work with AMR or flux register stuff
 *  \param[in]  a_JUnew Solution in conserved state.
 *  \param[out] a_JUnew Solution with species correction implemented
 *//*-----------------------------------------------------------------*/

//**FIXME I don't think application of this in a JU context makes much sense.
//**      I'm not convinced averages need sum to 1 (or density).  Point values,
//**      yes, averages, no.  On the other hand, if density is the true value,
//**      then these are just fractions of density and redistribution provides
//**      a consistent state without any adverse effect.
//**
//**      Revisit later code (i.e., consToPrim).  Is the trusted density the
//**      density or the sum of rhocn?

void
LevelCNSOp::speciesCorrection(LevelData<FArrayBox>& a_JUnew) const
{
  if (!(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf) ||
      !CRDparam::g_useSpeciesCorrection)
    {
      return;
    }
  CH_TIME("LevelCNSOp::speciesCorrection");
  CRD::msg << CRD::fv3 << "LevelCNSOp::speciesCorrection (level: "
           << m_level << ")" << CRD::end;

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);

      FArrayBox& JUFab = a_JUnew[dit];
      // Compute the flux and update JUnew
      m_patchOp.speciesCorrection(box,
                                  blockDomain,
                                  JUFab,
                                  m_dx);
    }
}

/*--------------------------------------------------------------------*/
//  Evaluate a_RHS (which is d(JU)/dt) at the current time based on
//  a_JU for a stage of Runge-Kutta
/** \param[out]  a_RHS  d(JU)/dt
 *  \param[in]   a_JU   Solution at beginning of stage
 *  \param[in]   a_stage
 *                      Stage index [0-3] for RK4
 *  \param[in]   a_stageTime
 *                      Time at start of stage
 *  \param[in]   a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]   a_fnFluxRegister
 *                      Flux register between this and a finer level
 *  \param[in]   a_crFluxRegister
 *                      Flux register between this and a coarser level
 *  \param[in]   a_dt   Time step on this level.  This is for the
 *                      complete time step, not just the stage.
 *  \param[in]   a_timeOld
 *                      Time at start of time step on this level
 *  \param[in]   a_crTimeOld
 *                      Time at start of time step on coarser level
 *  \param[in]   a_crTimeNew
 *                      Time at end of time step on coarser level
 *  \param[in]   a_subcycleParams
 *                      Parameters from subcyling algorithm (see
 *                      AMRLevel.H for more info)
 *  \param[in]   a_termFlags
 *                      Flags that indicate which terms of the NSE
 *                      to use. See the LevelCNS::Terms enum, stiff
 *                      means reacting source term, nonstiff is the rest
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::evalRHS(RHS&               a_RHS,
                    SOLN&              a_JU,
                    const int          a_stage,
                    const Real         a_stageTime,
                    const Real         a_stageWeight,
                    LevelFluxRegister& a_fnFluxRegister,
                    LevelFluxRegister& a_crFluxRegister,
                    const Real         a_dt,
                    const Real         a_timeOld,
                    const Real         a_crTimeOld,
                    const Real         a_crTimeNew,
                    SubcycleParams     a_subcycleParams,
                    const int          a_termFlags,
                    SOLN&              a_WOld)
{
  CH_assert(a_termFlags != Terms::None); // If none, why call this function?
  // Ensure we have at least some valid value in the flags
  CH_assert(a_termFlags & Terms::NonStiff || a_termFlags & Terms::Stiff);
  CH_TIME("LevelCNSOp::evalRHS");
  CRD::msg << CRD::fv2 << "LevelCNSOp::evalRHS (level: " << m_level
           << ", stage: " << a_stage << ')';
  if (CRDparam::g_additiveRK)
    {
      CRD::msg << CRD::fv2 << " stiff: " << a_termFlags;
      CNSIBC::s_lastRKStage = (a_stage == 5);
    }
  else
    {
      CNSIBC::s_lastRKStage = (a_stage == 3);
    }
  CRD::msg << CRD::fv2 << CRD::end;
  CNSIBC::s_firstRKStage = (a_stage == 0);

  // Find <U> in all cells and ghosts and store in m_U
  if ((a_termFlags & Terms::NonStiff) || (CRDparam::g_reactionOrder == 4))
    {
      fillGhostsRK4AndComputeU(a_JU,
                               a_stage,
                               a_dt,
                               a_timeOld,
                               a_crTimeOld,
                               a_crTimeNew,
                               a_subcycleParams);

      // Enforce positivity of species after exchange (also renormalize them)
      // speciesCorrection(m_U);
    }

  // LES subgrid-scale calculation
  const int sgsRefineRatio = CRDparam::g_sgskeCrsLevFilterRatio;
  if (((!m_hasCoarserGrid && m_levelGridMetrics.isMultiBlock()) ||
       !CRDparam::g_useSGSCoarsening ||
       ((sgsRefineRatio == 1) && !m_hasCoarserGrid)) &&
      (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex) &&
      (a_termFlags & Terms::NonStiff) &&
      !(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
    {
      // Compute SGSKE component in m_U
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = m_boxes[dit];
          CRDparam::g_CRDPhysics->cellAvgSGSKineticEnergy(
            m_U[dit],
            a_JU[dit],
            m_levelGridMetrics,
            dit(),
            false, //**FIXME: false = 2nd-order, true = 4th-order
            disjointBox);
        }

      // Find <U> in all cells and ghosts again
      fillGhostsRK4AndComputeU(a_JU,
                               a_stage,
                               a_dt,
                               a_timeOld,
                               a_crTimeOld,
                               a_crTimeNew,
                               a_subcycleParams);

      // Enforce positivity of species after exchange (also renormalize them)
      speciesCorrection(m_U);
    }

  // LES subgrid-scale kinetic energy calculation
  if ((a_termFlags & Terms::NonStiff) &&
      !(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
    {
      CRDparam::g_CRDPhysics->storeSGSKineticEnergy(m_sgske, m_U, true);
      CRDparam::g_CRDPhysics->storeSGSKineticEnergy(m_Jsgske, a_JU, true);
    }

  // For recirculating-inlet flat-plate, copy interior to face and rescale
  CRDparam::g_CNSIBC->copyInteriorToInletAndRescale(
    m_U, m_level, a_timeOld, a_dt, a_stage, m_boxes, m_levelGridMetrics);

  const bool setFlattening = (CRDparam::g_useFlattening && (a_stage == 0));
  // Set number of species transport equations
  const int numSpecies = CRDparam::g_numSpecies;
  const int numReactions = CRDparam::g_numReactions;
  // Face-averaged plot variables
  int plotVariables = CRDparam::g_plotLoFaceAvgComps;
  plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;

/*=================================================================*
 * Embedded lambda tEvalRHS1
 *=================================================================*/

  auto tEvalRHS1 =
    [this, &a_RHS, &a_WOld, setFlattening, numSpecies, numReactions,
     plotVariables, a_stage, a_stageTime, a_stageWeight,
     &a_fnFluxRegister, &a_crFluxRegister, a_dt, a_timeOld,
     a_crTimeOld, a_crTimeNew, a_subcycleParams, a_termFlags,
     &a_JU]
    (DataIndex didx, Box workingBox, Box disjointBox)
    {
      CH_assert(disjointBox.contains(workingBox));

      FArrayBox& WOld = a_WOld[didx];

      // Create a box that represents the working box and ghost cells
      // that would not be owned by other working boxes.
      // This box can be used to safely write to shared BaseFabs
      Box threadsafeGrowedWorkingBox = workingBox;
      // Find the directions from the workingBox that are not contained
      // by the disjointBox
      std::vector<std::pair<int, Side::LoHiSide>> growDirPair;
      for (const int dir : EachDir)
        {
          if (!disjointBox.contains(adjCellLo(workingBox, dir, 1)))
            {
              growDirPair.emplace_back(dir, Side::LoHiSide::Lo);
            }
          if (!disjointBox.contains(adjCellHi(workingBox, dir, 1)))
            {
              growDirPair.emplace_back(dir, Side::LoHiSide::Hi);
            }
        }
      int maxNumGhosts = std::max(std::max(std::max(
        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
        std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
                  CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan))),
        CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg)),
        CRDparam::queryGhosts(CRDparam::NumGhostWcellPnt));
      for (const std::pair<int, Side::LoHiSide>& growDir : growDirPair)
        {
          threadsafeGrowedWorkingBox.growDir(growDir.first, growDir.second, maxNumGhosts);
        }

      // The problem or block domain for the box should be computed once here
      // and passed to all functions.  This is for safety since it must be
      // computed from a box that is <= disjoint.
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);

      //**FIXME MMB needs special shown below
      // const Box& FBox1inDomain = m_grow1inDomainLayout[dit];
      // FluxBox FfaceAvg(FBox1inDomain, SpaceDim*m_numFluxes);
      Box box1Dom = grow(workingBox, 1);
      box1Dom &= blockDomain;
      Box box2Dom = grow(workingBox, 2);
      box2Dom &= blockDomain;

//--The following are defined here to allow use in both inertial and viscous
//--calculations.  The inertial routine will set the boundary conditions
//--sufficient for inertial and viscous.

      // The point conservative state U in the cells.
      Box boxUcellPntDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostUcellPnt));
      boxUcellPntDom &= blockDomain;
      FABSTACKTEMP(UcellPntFab,
                   boxUcellPntDom,
                   CRDparam::g_CRDPhysics->numConservative());

      // The average primitive state W in the cells.
      Box boxWcellAvgDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
      boxWcellAvgDom &= grow(blockDomain, 2);
      FABSTACKTEMP(WcellAvgFab,
                   boxWcellAvgDom,
                   CRDparam::g_CRDPhysics->numPrimitive());
      // FArrayBox& WcellAvgFab = m_WcellAvgLvl[didx];
      CH_assert(WcellAvgFab.box().contains(boxWcellAvgDom));
      // We use setVal because some operations may be applied to unset data,
      // later to be overwritten.  Sometimes it helps to comment this out to
      // view the unset data.
      WcellAvgFab.setVal(0.);

      // The point primitive state W in the cells.
      Box boxWcellPntDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostWcellPnt));
      boxWcellPntDom &= blockDomain;
      FABSTACKTEMP(WcellPntFab,
                   boxWcellPntDom,
                   CRDparam::g_CRDPhysics->numPrimitive());
      // FArrayBox& WcellPntFab = m_WcellPntLvl[didx];
      CH_assert(WcellPntFab.box().contains(boxWcellPntDom));

      // This is the average primitive state W on the faces
      Box boxWfaceAvgDom =
        grow(workingBox, std::max(
               CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTngBdry),
               std::max(CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgNrm),
                        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTan))));
      boxWfaceAvgDom &= blockDomain;
      // If the above includes <T>, then the native primitive state needs to
      // be computed on a slightly larger box.  Arrays are allocated to this
      // size, but the full primitive state is only available on the above
      // box for viscous calculations.
      Box boxWpFaceAvgDom = grow(
        workingBox,
        std::max(
          CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
          std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
                   CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan))));
      boxWpFaceAvgDom &= blockDomain;
      FLUXBOXSTACKTEMP(WfaceAvgFxb,
                       boxWpFaceAvgDom,
                       CRDparam::g_CRDPhysics->numPrimitive());
      // FluxBox& WfaceAvgFxb = m_WfaceAvgLvl[didx];
      CH_assert(WfaceAvgFxb.box().contains(boxWpFaceAvgDom));

//==INERTIAL (HYPERBOLIC) FLUX

      // Compute primitive values
      Box boxWfromUavgDom =
        grow(workingBox, CRDparam::queryGhosts(
               CRDparam::NumGhostWfromUcellAvg));
      boxWfromUavgDom &= blockDomain;
      const int numWcomp  = CRDparam::g_CRDPhysics->numPrimitive();
      // FAB that holds primitive state based on average of U
      FABSTACKTEMP(WfromUavgFab, boxWfromUavgDom, numWcomp);

      WfromUavgFab.setVal(0.0);
      // Now convert from cons to prim of cell data
      m_patchOp.computeWcell(UcellPntFab,
                             WcellPntFab,
                             WcellAvgFab,
                             WfromUavgFab,
                             didx,
                             workingBox,
                             blockDomain,
                             m_U[didx],
                             WOld);

      Interval WAvgIval = WcellAvgFab.interval();
      Interval WOldIval = WOld.interval();
      CH_assert(WOld.box().contains(boxWfromUavgDom));
      Box boxWfromUavgDomCopyBack = boxWfaceAvgDom;
      boxWfromUavgDomCopyBack &= threadsafeGrowedWorkingBox;
      WOld.copy(WfromUavgFab,
                boxWfromUavgDomCopyBack,
                WAvgIval.begin(),
                boxWfromUavgDomCopyBack,
                WOldIval.begin(),
                WOldIval.size());

//--Compute <W> on faces using primitive values

      m_patchOp.computeWfaceAvg(disjointBox,
                                workingBox,
                                blockDomain,
                                WfaceAvgFxb,
                                WcellAvgFab,
                                WfromUavgFab,
                                m_flattening[didx],
                                m_unitNormals[didx],
                                didx,
                                a_stageTime,
                                a_dt,
                                a_stageWeight,
                                setFlattening,
                                m_level);

//--Copy temporary data to member data

      Box boxWcellAvgDomCopyBack = boxWcellAvgDom;
      Box boxWpFaceAvgDomCopyBack = boxWpFaceAvgDom;
      Box boxWcellPntDomCopyBack = boxWcellPntDom;

      boxWcellAvgDomCopyBack &= threadsafeGrowedWorkingBox;
      boxWpFaceAvgDomCopyBack &= threadsafeGrowedWorkingBox;
      boxWcellPntDomCopyBack &= threadsafeGrowedWorkingBox;

      m_WcellAvgLvl[didx].copy(WcellAvgFab, boxWcellAvgDomCopyBack);
      m_WfaceAvgLvl[didx].copy(WfaceAvgFxb, boxWpFaceAvgDomCopyBack);
      m_WcellPntLvl[didx].copy(WcellPntFab, boxWcellPntDomCopyBack);

//--Update the block registers which solve a Riemann problem at block interfaces

      //**FIXME Should single blocks use this if periodic?
      if (m_levelGridMetrics.isMultiBlock())
        {
          const Interval primIntv(
            0, CRDparam::g_CRDPhysics->numPrimitive() - 1);
          m_MMBSingleLevelOp.setFlux(workingBox, WfaceAvgFxb, didx, primIntv);
        }
    };

/*=================================================================*
 * Embedded lambda tEvalRHS2
 *=================================================================*/

  auto tEvalRHS2 =
    [this, &a_RHS, &a_WOld, setFlattening, numSpecies, numReactions,
     plotVariables, a_stage, a_stageTime, a_stageWeight,
     &a_fnFluxRegister, &a_crFluxRegister, a_dt, a_timeOld,
     a_crTimeOld, a_crTimeNew, a_subcycleParams, a_termFlags,
     &a_JU]
    (DataIndex didx, Box workingBox, Box disjointBox)
    {
      CH_assert(disjointBox.contains(workingBox));

      FArrayBox& RHSfab = a_RHS[didx];
      FABSTACKTEMP(RHSchk, RHSfab.box(), RHSfab.nComp());
      FArrayBox& WOld = a_WOld[didx];

      // The problem or block domain for the box should be computed once here
      // and passed to all functions.  This is for safety since it must be
      // computed from a box that is <= disjoint.
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);

      //**FIXME MMB needs special shown below
      // const Box& FBox1inDomain = m_grow1inDomainLayout[dit];
      // FluxBox FfaceAvg(FBox1inDomain, SpaceDim*m_numFluxes);
      Box box1Dom = grow(workingBox, 1);
      box1Dom &= blockDomain;
      Box box2Dom = grow(workingBox, 2);
      box2Dom &= blockDomain;

      // The final flux
      FLUXBOXSTACKTEMP(NtFlux,
                       workingBox,
                       CRDparam::g_CRDPhysics->numFluxes());
      NtFlux.setVal(0.);

//--The following are defined here to allow use in both inertial and viscous
//--calculations.  The inertial routine will set the boundary conditions
//--sufficient for inertial and viscous.

      // The point conservative state U in the cells.
      Box boxUcellPntDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostUcellPnt));
      boxUcellPntDom &= blockDomain;
      FABSTACKTEMP(UcellPntFab,
                   boxUcellPntDom,
                   CRDparam::g_CRDPhysics->numConservative());

      // The average primitive state W in the cells.
      Box boxWcellAvgDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
      boxWcellAvgDom &= grow(blockDomain, 2);
      // FABSTACKTEMP(WcellAvgFab,
      //              boxWcellAvgDom,
      //              CRDparam::g_CRDPhysics->numPrimitive());
      // FArrayBox& WcellAvgFab = m_WcellAvgLvl[didx];
      FABSTACKTEMP(WcellAvgFab,
                   boxWcellAvgDom,
                   CRDparam::g_CRDPhysics->numPrimitive());
      WcellAvgFab.copy(m_WcellAvgLvl[didx], boxWcellAvgDom);

      CH_assert(WcellAvgFab.box().contains(boxWcellAvgDom));

      // The point primitive state W in the cells.
      Box boxWcellPntDom =
        grow(workingBox, CRDparam::queryGhosts(CRDparam::NumGhostWcellPnt));
      boxWcellPntDom &= blockDomain;
      // FABSTACKTEMP(WcellPntFab,
      //              boxWcellPntDom,
      //              CRDparam::g_CRDPhysics->numPrimitive());
      // FArrayBox& WcellPntFab = m_WcellPntLvl[didx];
      FABSTACKTEMP(WcellPntFab,
                   boxWcellPntDom,
                   CRDparam::g_CRDPhysics->numPrimitive());
      WcellPntFab.copy(m_WcellPntLvl[didx], boxWcellPntDom);

      CH_assert(WcellPntFab.box().contains(boxWcellPntDom));

      // This is the average primitive state W on the faces
      Box boxWfaceAvgDom =
        grow(workingBox, std::max(
               CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTngBdry),
               std::max(CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgNrm),
                        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTan))));
      boxWfaceAvgDom &= blockDomain;
      // If the above includes <T>, then the native primitive state needs to
      // be computed on a slightly larger box.  Arrays are allocated to this
      // size, but the full primitive state is only available on the above
      // box for viscous calculations.
      Box boxWpFaceAvgDom = grow(
        workingBox,
        std::max(
          CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
          std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
                   CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan))));
      boxWpFaceAvgDom &= blockDomain;
      // FLUXBOXSTACKTEMP(WfaceAvgFxb,
      //                  boxWpFaceAvgDom,
      //                  CRDparam::g_CRDPhysics->numPrimitive());
      //FluxBox& WfaceAvgFxb = m_WfaceAvgLvl[didx];
      FLUXBOXSTACKTEMP(WfaceAvgFxb,
                   boxWpFaceAvgDom,
                   CRDparam::g_CRDPhysics->numPrimitive());
      WfaceAvgFxb.copy(m_WfaceAvgLvl[didx], boxWpFaceAvgDom);

      CH_assert(WfaceAvgFxb.box().contains(boxWpFaceAvgDom));

      // This is the point values of primitive state W on the faces
      Box boxWfacePntDom =
        grow(
          workingBox,
          std::max(
            CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTngBdry),
            std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntNrm),
                     CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntTan))));
      boxWfacePntDom &= blockDomain;
      FLUXBOXSTACKTEMP(WfacePntFxb,
                       boxWfacePntDom,
                       CRDparam::g_CRDPhysics->numPrimitive());

      // Mapped cell-averaged turbulent source terms
      FABSTACKTEMP(turbSourceAvgFab, workingBox,
                   CRDparam::g_CRDPhysics->numFluxes());
      turbSourceAvgFab.setVal(0.);

      // Physical-space velocity gradients on faces (for LES wall-model)
      FLUXBOXSTACKTEMP(facePntVelGradFxb, box1Dom, SpaceDim*SpaceDim);
      // Combined stress fluxes for LES wall-model (SGS and viscous)
      FLUXBOXSTACKTEMP(stressFluxFxb, box1Dom, SpaceDim*SpaceDim);
      stressFluxFxb.setVal(0.);

      // Inverse of the time steps, set to 1 larger than workingBox for
      // source terms but only use workingBox box for time steps
      FABSTACKTEMP(invDtFab, box1Dom, 1);
      invDtFab.setVal(0.);

//==INERTIAL (HYPERBOLIC) FLUX

      if (a_termFlags & Terms::NonStiff)
        {

//--Compute <F> using primitive values

          // This is the full fourth-order <F>, stored as (dir, var) with var
          // contiguous
          Box boxFfaceAvgDom =
            grow(workingBox,
                 CRDparam::queryGhosts(CRDparam::NumGhostFfromWpFaceAvg));
          boxFfaceAvgDom &= blockDomain;
          FLUXBOXSTACKTEMP(flux,
                           boxFfaceAvgDom,
                           SpaceDim*CRDparam::g_CRDPhysics->numFluxes());
          // This is a second-order <F> computed from WAvgFace, and used for
          // gradients.  This is only need if NumGhostInertialFfaceAvg == 0.
          Box boxFfromWpAvgDom;  // Empty
          if (CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg) == 0)
            {
              boxFfromWpAvgDom =
                grow(workingBox,
                     CRDparam::queryGhosts(CRDparam::NumGhostFfromWpFaceAvg));
              boxFfromWpAvgDom &= blockDomain;
            }
          FLUXBOXSTACKTEMP(fluxFromWpAvg,
                           boxFfromWpAvgDom,
                           SpaceDim*CRDparam::g_CRDPhysics->numFluxes());

          m_patchOp.addInertialFlux(disjointBox,
                                    workingBox,
                                    blockDomain,
                                    flux,
                                    fluxFromWpAvg,
                                    WfaceAvgFxb,
                                    WfacePntFxb,
                                    WcellAvgFab,
                                    WcellPntFab,
                                    m_faceAvgPlotData[didx],
                                    m_flattening[didx],
                                    m_bndryCellData[didx],
                                    m_bndryNtJ[didx],
                                    m_U[didx],
                                    m_unitNormals[didx],
                                    didx,
                                    a_stageTime,
                                    a_dt,
                                    m_prevDt,
                                    a_stageWeight,
                                    setFlattening,
                                    m_level);

//--Convert <F> to <NtF>

          // In case only using viscous model, avoid this code since 'flux'
          // is not set.  Inertial update only sets BC conditions in
          // WcellAvgFab.
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              const FluxBox& Nfxb = m_levelGridMetrics.m_N[didx];
              if (CRDparam::g_cartesian)
                {
                  CH_TIME("LevelCNS::evalRHS::CartCopy");
                  const int numFluxes = CRDparam::g_CRDPhysics->numFluxes();
                  for (int dir = 0; dir != SpaceDim; ++dir)
                    {
                      Box faceBox = workingBox;
                      faceBox.surroundingNodes(dir);
                      NtFlux[dir].copy(flux[dir],
                                       faceBox,
                                       dir*numFluxes,
                                       faceBox,
                                       0,
                                       numFluxes);
                    }
                }
              else
                {
                  // Depending on number of ghosts available, may use flux
                  // itself or fluxFromWAvg to compute gradients of fluxes.
                  FluxBox* fluxGrad =
                    (CRDparam::queryGhosts(
                      CRDparam::NumGhostInertialFfaceAvg) == 0)
                    ? &fluxFromWpAvg : &flux;
                  const Interval fluxIntv(
                    0, CRDparam::g_CRDPhysics->numFluxes() - 1);
                  //**FIXME CS should come from AMRLevel
                  const BlockCoordSys* blockCoordSys =
                    m_levelGridMetrics.getCoordSys(workingBox);
                  bool fourthOrder = false;
                  // Only if both convolve and deconvolve are set to 2nd do we
                  // perform the product rule at 2nd order
                  if (CRDparam::g_faceConvolveFlatten < 2 ||
                      CRDparam::g_faceDeconvolveFlatten < 2)
                    {
                      fourthOrder = true;
                    }
                  blockCoordSys->
                    computeMetricTermProductAverage(NtFlux,
                                                    flux,
                                                    Nfxb,
                                                    SpaceDim,
                                                    *fluxGrad,
                                                    workingBox,
                                                    fourthOrder,
                                                    fluxIntv,
                                                    fluxIntv,
                                                    0,
                                                    &blockDomain);
                }

//--Copy data to LevelData for plotting (RK stage for plotting can be changed)

              if ((plotVariables & CRDparam::PlotPrimitive) && (a_stage == 0))
                {
                  m_faceAvgPlotData[didx].copy(
                    WfaceAvgFxb, 0, 0, CRDparam::g_CRDPhysics->numPrimitive());
                }

//--Copy data to LevelData for plotting (RK stage for plotting can be changed)

              if ((plotVariables & CRDparam::PlotFluxes) && (a_stage == 0))
                {
                  int cStart = 0;
                  if (plotVariables & CRDparam::PlotPrimitive)
                    {
                      cStart = CRDparam::g_CRDPhysics->numPrimitive();
                    }
                  m_faceAvgPlotData[didx].copy(
                    NtFlux, 0, cStart, CRDparam::g_CRDPhysics->numFluxes());
                }


              Real minConvDt = std::numeric_limits<Real>::max();
              IntVect minConvDtCell;

              // See "High-Order, Finite-Volume Methods in Mapped
              // Coordinates", Colella et al., JCP 2011 for derivation of
              // this constraint.
              const Real stabilityConstraint = 1.3925;
              CRDparam::g_CRDPhysics->
                getMaxWaveSpeedEvalRHS(blockDomain,
                                       workingBox,
                                       disjointBox,
                                       invDtFab,
                                       WfacePntFxb,
                                       Nfxb,
                                       m_levelGridMetrics.m_J[didx],
                                       m_levelGridMetrics,
                                       stabilityConstraint,
                                       m_dx,
                                       minConvDt,
                                       minConvDtCell);

              while (m_lockMinLocalDt.test_and_set(std::memory_order_acquire));
              if (minConvDt < m_minConvDt)
              {
                m_minConvDt = minConvDt;
                m_minConvDtCell = minConvDtCell;
              }
              m_lockMinLocalDt.clear(std::memory_order_release);

            }

          // if (workingBox.contains(IntVect{ 659, 2, 23 }))
          //   {
          //     RHSchk.setVal(0.);
          //     m_patchOp.fluxDivergence(workingBox, RHSchk, NtFlux, m_dx);
          //     CRD::msg.setMaxPrecFloatSN();
          //     CRD::msg << "!!-NAME:     "
          //              << CRDparam::g_CRDPhysics->consvStateName(9)
          //              << CRD::end;
          //     CRD::msg << "!!-Inertial: " << RHSchk(IntVect{ 659, 2, 23 }, 9)
          //              << CRD::end;
          //     CRD::msg.setFloatDefault();
          //   }
          }

      if (a_termFlags & Terms::NonStiff)
        {

//==VISCOUS (ELLIPTIC) FLUX

          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              // The final flux
              FLUXBOXSTACKTEMP(NtViscFlux,
                               workingBox,
                               CRDparam::g_CRDPhysics->numFluxes());
              NtViscFlux.setVal(0.);
              // Face-averaged flux G
              FLUXBOXSTACKTEMP(fluxViscousfaceAvg,
                               workingBox,
                               SpaceDim*CRDparam::g_CRDPhysics->numFluxes());
              // Face-centered flux G
              FLUXBOXSTACKTEMP(fluxViscousfacePnt,
                               box1Dom,
                               SpaceDim*CRDparam::g_CRDPhysics->numFluxes());

              Real minDiffDt = std::numeric_limits<Real>::max();
              IntVect minDiffDtCell;

              m_patchOp.addViscousFlux(workingBox,
                                       blockDomain,
                                       fluxViscousfaceAvg,
                                       fluxViscousfacePnt,
                                       invDtFab,
                                       turbSourceAvgFab,
                                       WcellAvgFab,
                                       WcellPntFab,
                                       m_timeAvgData[didx],
                                       m_faceAvgPlotData[didx],
                                       WfaceAvgFxb,
                                       WfacePntFxb,
                                       facePntVelGradFxb,
                                       stressFluxFxb,
                                       m_facePntDeltaC[didx],
                                       m_faceCoord[didx],
                                       m_unitNormals[didx],
                                       didx,
                                       a_dt,
                                       a_stageTime,
                                       m_totalTimeAfterFilterInit,
                                       m_level,
                                       minDiffDt,
                                       minDiffDtCell);


              while (m_lockMinLocalDt.test_and_set(std::memory_order_acquire));
              if (minDiffDt < m_minDiffDt)
              {
                m_minDiffDt = minDiffDt;
                m_minDiffDtCell = minDiffDtCell;
              }
              m_lockMinLocalDt.clear(std::memory_order_release);

              if (CRDparam::g_cartesian)
                {
                  CH_TIME("LevelCNS::evalRHS::CartCopy");
                  const int numFluxes = CRDparam::g_CRDPhysics->numFluxes();
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      Box faceBox = workingBox;
                      faceBox.surroundingNodes(faceDir);
                      // Note: the viscous flux dyad is stored with unit-stride
                      // in space.  So the state components of the flux for a
                      // given direction are separated by SpaceDim in the flux
                      // dyad.  Starting component depends on direction.
                      int idxFluxDyad = faceDir;
                      // This is a loop over the state vector
                      for (int idxFlux = 0; idxFlux != numFluxes; ++idxFlux)
                        {
                          NtFlux[faceDir].minus(fluxViscousfaceAvg[faceDir],
                                                faceBox,
                                                idxFluxDyad,
                                                idxFlux,
                                                1);
                          idxFluxDyad += SpaceDim;
                        }
                    }
                }
              else
                {
                  const BlockCoordSys* blockCoordSys =
                    m_levelGridMetrics.getCoordSys(workingBox);
                  const FluxBox& Nfxb = m_levelGridMetrics.m_N[didx];
                  const Interval fluxIntv(
                    0, CRDparam::g_CRDPhysics->numFluxes() - 1);
                  bool fourthOrder = false;
                  // Only if both convolve and deconvolve are set to 2nd do we
                  // perform the product rule at 2nd order
                  if (CRDparam::g_faceConvolveFlatten < 2 ||
                      CRDparam::g_faceDeconvolveFlatten < 2)
                    {
                      fourthOrder = true;
                    }
                  blockCoordSys->
                    computeMetricTermProductAverage(NtViscFlux,
                                                    fluxViscousfaceAvg,
                                                    Nfxb,
                                                    SpaceDim,
                                                    fluxViscousfacePnt,
                                                    workingBox,
                                                    fourthOrder,
                                                    fluxIntv,
                                                    fluxIntv,
                                                    1,
                                                    &blockDomain);


                  NtFlux.minus(NtViscFlux, workingBox, 0, 0,
                               CRDparam::g_CRDPhysics->numFluxes());

//--Update the block registers which average single-level fluxes at block
//--interfaces for viscous flux only

                  //**FIXME Should single blocks use this if periodic?
                  if (m_levelGridMetrics.isMultiBlock())
                    {
                      for (const int dir : EachDir)
                        {
                          NtViscFlux[dir].negate();
                        }
                      const Interval fluxIntv(
                        0, CRDparam::g_CRDPhysics->numFluxes() - 1);
                      m_MMBSingleLevelOp.setFlux(workingBox,
                                                 NtViscFlux,
                                                 didx,
                                                 fluxIntv);
                    }

                  // if (workingBox.contains(IntVect{ 659, 2, 23 }))
                  //   {
                  //     RHSchk.setVal(0.);
                  //     m_patchOp.fluxDivergence(
                  //       workingBox, RHSchk, NtViscFlux, m_dx);
                  //     CRD::msg << "!!-Viscous: "
                  //              << -RHSchk(IntVect{ 659, 2, 23 }, 9)
                  //              << CRD::end;
                  //   }
                }

//--Copy data to LevelData for plotting (RK stage for plotting can be changed)
              if ((plotVariables & CRDparam::PlotFluxes) && (a_stage == 0))
                {
                  int cStart = 0;
                  if (plotVariables & CRDparam::PlotPrimitive)
                    {
                      cStart = CRDparam::g_CRDPhysics->numPrimitive();
                    }
                  if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
                    {
                      cStart += CRDparam::g_CRDPhysics->numFluxes();
                    }
                  m_faceAvgPlotData[didx].copy(
                    NtViscFlux, 0, cStart,
                    CRDparam::g_CRDPhysics->numFluxes());
                }
            }

//--Update the flux registers

          m_patchOp.updateFluxRegisters(a_crFluxRegister,
                                        a_fnFluxRegister,
                                        NtFlux,
                                        workingBox,
                                        didx,
                                        a_stageWeight,
                                        m_hasCoarserGrid,
                                        m_hasFinerGrid);

//--Evaluate the divergence to get the RHS

          m_patchOp.fluxDivergence(workingBox, RHSfab, NtFlux, m_dx);

//--Add any turbulent source terms (for WMLES -- source physics not required) 

          if (CRDparam::g_turbModelType)
            {
              // Before adding source terms, compute wall-model terms using
              // <U> update (RHSfab)
              // NOTE: doesn't account for source terms added to <U> RHSfab
              //       However, it's easy to move this below addSourceTerm
              //       BE CAREFUL if this is done. It may be that addSourceTerm
              //       eventually has a dependency on turbSourceAvgFab
              CRDparam::g_CRDPhysics->updateEtaZero(
                turbSourceAvgFab, WcellAvgFab, a_JU[didx], facePntVelGradFxb,
                WfaceAvgFxb, NtFlux, stressFluxFxb, RHSfab, m_unitNormals[didx],
                blockDomain, m_levelGridMetrics, didx, disjointBox,
                workingBox, a_dt);
              const int cTurbBegin =
                CRDparam::g_CRDPhysics->turbConsInterval().begin();
              const int numTurbComp =
                CRDparam::g_CRDPhysics->turbConsInterval().size();
              RHSfab.plus(turbSourceAvgFab,cTurbBegin,cTurbBegin,numTurbComp);
            }

//--Add any source terms

          if (CRDparam::g_physicsModels & CRDparam::PhysicsSource)
            {
              // FAB to add the source term
              FABSTACKTEMP(ScellPntFab,
                           box1Dom,
                           CRDparam::g_CRDPhysics->numFluxes());
              ScellPntFab.setVal(0.);
              // Add source term
              m_patchOp.addSourceTerm(ScellPntFab,
                                      invDtFab,
                                      WcellPntFab,
                                      m_U[didx],
                                      WfaceAvgFxb,
                                      a_stageTime,
                                      a_stageWeight,
                                      m_level,
                                      disjointBox,
                                      blockDomain,
                                      box1Dom,
                                      didx,
                                      m_globalKE,
                                      m_globalHelicity);
              FABSTACKTEMP(mappedScellAvgFab,
                           workingBox,
                           CRDparam::g_CRDPhysics->numFluxes());
              // Solve for <JS> from (S)
              m_patchOp.solveJSfromS(mappedScellAvgFab,
                                     ScellPntFab,
                                     workingBox,
                                     box1Dom,
                                     didx,
                                     blockDomain,
                                     true); // fourth-order
              // Add to the RHS
              RHSfab.plus(mappedScellAvgFab,
                          workingBox,
                          0,
                          0,
                          mappedScellAvgFab.nComp());
            }

//--Perform some time-averaging of turbulent quantities
          if (CNSIBC::s_firstRKStage &&
              (a_stageTime >= CRDparam::g_startTimeAvgTime))
            {
              PatchMappedFunc::timeAvgData(m_faceAvgTimeAvgPlotData[didx],
                                           m_faceAvgPlotData[didx], blockDomain,
                                           workingBox, a_dt, a_stageTime,
                                           m_totalTimeAfterFilterInit);
            }
        }  // Solve non-stiff

      if (a_termFlags & Terms::Stiff)
        {

//--Add reaction source terms

          // Now convert from cons to prim of cell data
          const Box boxWcellPntFabRxnDom =
            (CRDparam::g_reactionOrder == 2) ? workingBox : box1Dom;
          if (!(a_termFlags & Terms::NonStiff))
            {
              const bool fourthOrderWpntState =
                (CRDparam::g_reactionOrder == 2) ? false : true;
              m_patchOp.computeWstate(WcellPntFab,
                                      m_U[didx],
                                      WOld,
                                      boxWcellPntFabRxnDom,
                                      blockDomain,
                                      true,
                                      fourthOrderWpntState);
            }
          // Box on which to compute cell-centered reactions
          const Box boxRxnCellPntFabBox =
            (CRDparam::g_reactionOrder == 2) ? workingBox : box1Dom;

          // Reactions can only occur if we have at least one reaction set
          if ((CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf) &&
              (numReactions > 0) &&
              (a_stageTime > CRDparam::g_reactionStartTime) &&
              (a_stageTime < CRDparam::g_reactionEndTime))
            {
              // First index for species transport equations
              const int UcompStart = CRDparam::g_CRDPhysics->
                speciesConsInterval().begin();
              // Number of components that have reactions
              const int numRCT = numSpecies;
              // FAB to solve the reaction source term
              FABSTACKTEMP(RCTcellPntFab, boxRxnCellPntFabBox, numRCT);
              // Solve for rho*omega, return species associated with smallest
              // dt
  //             if (a_stage == 3)
  //               {
  // if (boxRxnCellPntFabBox.contains(IntVect{ 284, 11 }))
  //   {
  //     const int cPrimSpecBeg =
  //       CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  //     const int cPres = CRDparam::g_CRDPhysics->pressureIndex();
  //     const int cTemp = CRDparam::g_CRDPhysics->temperatureIndex();
  //     CRD::msg << "Intro frac: " << WcellPntFab(IntVect{ 284, 11 }, cPrimSpecBeg + 4) << CRD::end;
  //     CRD::msg << "      pres: " << WcellPntFab(IntVect{ 284, 11 }, cPres) << CRD::end;
  //     CRD::msg << "      temp: " << WcellPntFab(IntVect{ 284, 11 }, cTemp) << CRD::end;
  //   }

                  Real minChemDt = std::numeric_limits<Real>::max();
                  IntVect minChemDtCell;
                  int cnlim;

                  cnlim = CRDparam::g_CRDPhysics->addReactionSource(
                    boxRxnCellPntFabBox,
                    RCTcellPntFab,
                    invDtFab,
                    WcellPntFab,
                    a_stageTime,
                    m_level,
                    minChemDt,
                    minChemDtCell);

                  while (m_lockMinLocalDt.test_and_set(std::memory_order_acquire));
                  if (minChemDt < m_minChemDt)
                  {
                    m_cnlim = cnlim;
                    m_minChemDt = minChemDt;
                    m_minChemDtCell = minChemDtCell;
                  }
                  m_lockMinLocalDt.clear(std::memory_order_release);

  // if (boxRxnCellPntFabBox.contains(IntVect{ 284, 11 }))
  //   {
  //     CRD::msg << "End rate: " << RCTcellPntFab(IntVect{ 284, 11 }, 4)
  //              << CRD::end;
  //   }
  //               }
  //             else
                // {
                //   // Real advanceDt = (a_stage < 2) ? 0.5*a_dt : a_dt;
                //   Real advanceDt = a_dt;
                //   cnlim = CRDparam::g_CRDPhysics->ARSwithDiagonalFluxCorrection(
                //     boxRxnCellPntFabBox,
                //     RCTcellPntFab,
                //     invDtFab,
                //     WcellPntFab,
                //     UcellPntFab,
                //     advanceDt,
                //     m_level,
                //     m_minChemDt,
                //     m_minChemDtCell);
                // }

              FABSTACKTEMP(mappedSource,
                           workingBox,
                           numRCT);
              // Solve for <J rho omega> from (rho omega)
              m_patchOp.solveJSfromS(mappedSource,
                                     RCTcellPntFab,
                                     workingBox,
                                     boxRxnCellPntFabBox,
                                     didx,
                                     blockDomain,
                                     CRDparam::g_reactionOrder == 4);

              // Add the reaction source term to the species transport
              // equations

        // if (workingBox.contains(IntVect{ 284, 11 }))
        //   {
        //     CRD::msg << "!!-Reaction: "
        //              << mappedSource(IntVect{ 284, 11 }, 4)
        //              << CRD::end;
        //   }
              const int srcComp = 0;
              RHSfab.plus(mappedSource,
                          workingBox,
                          srcComp,
                          UcompStart,
                          numRCT);
              if (cnlim >= 0)
                {
                  CRD::msg << CRD::fv3 << "Chemical limiting species: "
                           << minChemDt << " "
                           << CRDparam::g_speciesNames[cnlim] << CRD::end;
                }
            }
        }

      // Invert sum of time steps
      Real minLocalDt = 1.E9;
      MD_ARRAY_RESTRICT(arrDtFab, invDtFab);
      MD_BOXLOOP(workingBox, i)
        {
          Real invdt = arrDtFab[MD_IX(i,0)];
          if (1./invdt > 1.E-15 && invdt > 1.E-15 && invdt == invdt)
            {
              minLocalDt = std::min(minLocalDt, 1./invdt);
            }
        }
      while (m_lockMinLocalDt.test_and_set(std::memory_order_acquire));
      m_minLocalDt = std::min(m_minLocalDt, minLocalDt);
      m_lockMinLocalDt.clear(std::memory_order_release);
  };

//--Loop over the boxes ("patches") and evaluate the fluxes.

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      DataIndex didx = dit();
      FArrayBox& RHSfab = a_RHS[didx];
      RHSfab.setVal(0.);

#ifdef CH_CTHR
      if (CRDparam::g_threads)
        {
          std::vector<Box> patchBoxes;
          std::pair<IntVect, IntVect> boxPart =
            cachedBoxPart(disjointBox,
                          CRDparam::g_threads->numThreadsPerTeam());
          const IntVect splitDim = boxPart.second - 1;
          stc::nestedLoop(IntVect_zero, splitDim,
                          [&, partDim = boxPart.first]
                          (const IntVect& iv)
                            {
                              const IntVect lo = disjointBox.smallEnd() +
                                partDim*iv;
                              const IntVect hi = disjointBox.smallEnd() +
                                partDim*(iv + IntVect_unit) - IntVect_unit;
                              patchBoxes.emplace_back(lo, hi);
                            });
          std::vector<std::function<void(DataIndex, Box, Box)>>
            subTasks{ patchBoxes.size(),
                      std::function<void(DataIndex, Box, Box)>(tEvalRHS1) };
          std::vector<ThreadTools::ChomboSubTaskArgs_t> args;
          for (const auto& workingBox : patchBoxes)
            {
              args.emplace_back(didx, workingBox, disjointBox);
            }
          CRDparam::g_threads->addNewTask(subTasks, args);
        }
      else
#endif
        {
          tEvalRHS1(didx, disjointBox, disjointBox);
        }
    }  // Loop over boxes

#ifdef CH_CTHR
  if (CRDparam::g_threads)
    {
      CRDparam::g_threads->wait();
    }
#endif

  // LES SGS KE calculation for coarsest level -- single-block only
  estimateCrsLevSingleBlockSGSKE(a_JU);

//--Do inter-block Riemann solve here

  // At interfaces between block boundaries, set flux via a Riemann solve
  if (m_levelGridMetrics.isMultiBlock() && (a_termFlags & Terms::NonStiff))
    {
      int numWcomp = CRDparam::g_CRDPhysics->numNativePrimitive();
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
        }
      const Interval primIntv(0, numWcomp - 1);
      m_MMBSingleLevelOp.exchange(false);  // False for no sign considerations
      m_MMBSingleLevelOp.applyAtBoundary(
        m_boxes,                    // Boxes
        [this, numWcomp, primIntv]  // Start of lambda
        (const Box&           a_cells,
         const int            a_dir,
         const Side::LoHiSide a_side,
         const DataIndex&     a_didx,
         const FArrayBox&     a_WfaceAvgLR)
        {
          const Interval primIntv(0, numWcomp - 1);
          Box riemannBox(a_cells);  // Cells at interior
          riemannBox.shiftHalf(a_dir, Side::sign(a_side));  // Faces on boundary
          // The fabs are shifted left and right so we cannot use the same one
          // for both the left and right state.  Here we just make a copy of
          // both
          Box LRbox(a_WfaceAvgLR.box());
          FABSTACKTEMP(Wminus, LRbox, numWcomp);
          Wminus.copy(a_WfaceAvgLR, LRbox, 0, LRbox, 0, numWcomp);
          FABSTACKTEMP(Wplus,  LRbox, numWcomp);
          Wplus .copy(a_WfaceAvgLR, LRbox, 0, LRbox, 0, numWcomp);
          const FluxBox& unitNormalFxb = m_unitNormals[a_didx];
          PatchCNSOp::preRiemann(Wplus,
                                 Wminus,
                                 unitNormalFxb,
                                 a_dir,
                                 riemannBox);
          FArrayBox& WfaceAvgStar = m_WfaceAvgLvl[a_didx][a_dir];
          CRDparam::g_CRDPhysics->riemann(
            WfaceAvgStar,  // On dir-FACEs
            Wplus,         // On cells to left of idir-FACEs
            Wminus,        // On cells to right of idir-FACEs
            a_dir,
            riemannBox);
          PatchCNSOp::postRiemann(WfaceAvgStar,
                                  unitNormalFxb,
                                  a_dir,
                                  riemannBox);
        }); 
    }

//--In-situ processing of domain-sums and point-values (case specific)
  computeTurbForcingSums(m_WcellPntLvl, m_globalKE, m_globalHelicity);

//--Loop over the boxes ("patches") and evaluate the fluxes.

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      DataIndex didx = dit();
      FArrayBox& RHSfab = a_RHS[didx];
      RHSfab.setVal(0.);

#ifdef CH_CTHR
      if (CRDparam::g_threads)
        {
          std::vector<Box> patchBoxes;
          std::pair<IntVect, IntVect> boxPart =
            cachedBoxPart(disjointBox,
                          CRDparam::g_threads->numThreadsPerTeam());
          const IntVect splitDim = boxPart.second - 1;
          stc::nestedLoop(IntVect_zero, splitDim,
                          [&, partDim = boxPart.first]
                          (const IntVect& iv)
                            {
                              const IntVect lo = disjointBox.smallEnd() +
                                partDim*iv;
                              const IntVect hi = disjointBox.smallEnd() +
                                partDim*(iv + IntVect_unit) - IntVect_unit;
                              patchBoxes.emplace_back(lo, hi);
                            });
          std::vector<std::function<void(DataIndex, Box, Box)>>
            subTasks{ patchBoxes.size(),
                      std::function<void(DataIndex, Box, Box)>(tEvalRHS2) };
          std::vector<ThreadTools::ChomboSubTaskArgs_t> args;
          for (const auto& workingBox : patchBoxes)
            {
              args.emplace_back(didx, workingBox, disjointBox);
            }
          CRDparam::g_threads->addNewTask(subTasks, args);
        }
      else
#endif
        {
          tEvalRHS2(didx, disjointBox, disjointBox);
        }
    }  // Loop over boxes

  // For flat-plate, spatially average m_timeAvgData (problem specific)
  CRDparam::g_CNSIBC->spatiallyAverageData(m_timeAvgData,
                                           a_stage, a_stageTime);

#ifdef CH_CTHR
  if (CRDparam::g_threads)
    {
      CRDparam::g_threads->wait();
    }
#endif

  // At interfaces between block boundaries, reflux with the average flux
  // computed from both sides.
  if (m_levelGridMetrics.isMultiBlock() && (a_termFlags & Terms::NonStiff) &&
      (CRDparam::g_physicsModels & CRDparam::PhysicsViscous))
    {
      const Interval fluxIntv(
        0, CRDparam::g_CRDPhysics->numFluxes() - 1);
      m_MMBSingleLevelOp.exchange(false);  // False for no sign considerations
      m_MMBSingleLevelOp.refluxAverage(a_RHS, m_dx, fluxIntv);
    }

//--Do in-situ processing of domain-sums and point-values (case specific)
  // Placed here so that finalized RHS can be used if necessary
  CRDparam::g_CNSIBC->inSituSumPntProcessing(m_faceAvgPlotData,
                                             m_WfaceAvgLvl,
                                             m_WcellAvgLvl,
                                             m_WcellPntLvl,
                                             m_unitNormals,
                                             m_levelGridMetrics,
                                             a_stageTime,
                                             a_stage,
                                             m_level);

  m_firstDt = false;
  if (CNSIBC::s_lastRKStage)
    {
      m_prevDt = a_dt;
      if (a_timeOld >= CRDparam::g_startTimeAvgTime)
        {
          m_totalTimeAfterFilterInit += a_dt;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Evaluate the NonStiff terms of a_RHS (which is d(JU)/dt) at the
//  current time based on a_JU for a stage of Runge-Kutta.
/** \param[out]  a_RHS  d(JU)/dt
 *  \param[in]   a_JU   Solution at beginning of stage
 *  \param[in]   a_stage
 *                      Stage index [0-3] for RK4
 *  \param[in]   a_stageTime
 *                      Time at start of stage
 *  \param[in]   a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]   a_fnFluxRegister
 *                      Flux register between this and a finer level
 *  \param[in]   a_crFluxRegister
 *                      Flux register between this and a coarser level
 *  \param[in]   a_dt   Time step on this level.  This is for the
 *                      complete time step, not just the stage.
 *  \param[in]   a_timeOld
 *                      Time at start of time step on this level
 *  \param[in]   a_crTimeOld
 *                      Time at start of time step on coarser level
 *  \param[in]   a_crTimeNew
 *                      Time at end of time step on coarser level
 *  \param[in]   a_subcycleParams
 *                      Parameters from subcyling algorithm (see
 *                      AMRLevel.H for more info)
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::evalNonStiff(RHS&               a_RHS,
                         SOLN&              a_JU,
                         const int          a_stage,
                         const Real         a_stageTime,
                         const Real         a_stageWeight,
                         LevelFluxRegister& a_fnFluxRegister,
                         LevelFluxRegister& a_crFluxRegister,
                         const Real         a_dt,
                         const Real         a_timeOld,
                         const Real         a_crTimeOld,
                         const Real         a_crTimeNew,
                         SubcycleParams     a_subcycleParams,
                         SOLN&              a_WOld)
{
  CH_TIME("LevelCNSOp::evalNonStiff");
  CRD::msg << CRD::fv3 << "LevelCNSOp::evalNonStiff (level: " << m_level << ")"
           << CRD::end;
  int termFlags = LevelCNSOp::Terms::NonStiff;
  this->evalRHS(a_RHS,
                a_JU,
                a_stage,
                a_stageTime,
                a_stageWeight,
                a_fnFluxRegister,
                a_crFluxRegister,
                a_dt,
                a_timeOld,
                a_crTimeOld,
                a_crTimeNew,
                a_subcycleParams,
                termFlags,
                a_WOld);
}

/*--------------------------------------------------------------------*/
//  Evaluate the stiff terms of a_RHS (which is d(JU)/dt) at the
//  current time based on a_JU for a stage of Runge-Kutta.
/** \param[out]  a_RHS  d(JU)/dt
 *  \param[in]   a_JU   Solution at beginning of stage
 *  \param[in]   a_stage
 *                      Stage index [0-3] for RK4
 *  \param[in]   a_stageTime
 *                      Time at start of stage
 *  \param[in]   a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]   a_fnFluxRegister
 *                      Flux register between this and a finer level
 *  \param[in]   a_crFluxRegister
 *                      Flux register between this and a coarser level
 *  \param[in]   a_dt   Time step on this level.  This is for the
 *                      complete time step, not just the stage.
 *  \param[in]   a_timeOld
 *                      Time at start of time step on this level
 *  \param[in]   a_crTimeOld
 *                      Time at start of time step on coarser level
 *  \param[in]   a_crTimeNew
 *                      Time at end of time step on coarser level
 *  \param[in]   a_subcycleParams
 *                      Parameters from subcyling algorithm (see
 *                      AMRLevel.H for more info)
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::evalStiff(RHS&               a_RHS,
                      SOLN&              a_JU,
                      const int          a_stage,
                      const Real         a_stageTime,
                      const Real         a_stageWeight,
                      LevelFluxRegister& a_fnFluxRegister,
                      LevelFluxRegister& a_crFluxRegister,
                      const Real         a_dt,
                      const Real         a_timeOld,
                      const Real         a_crTimeOld,
                      const Real         a_crTimeNew,
                      SubcycleParams     a_subcycleParams,
                      SOLN&              a_WOld)
{
  CH_TIME("LevelCNSOp::evalStiff");
  CRD::msg << CRD::fv4 << "LevelCNSOp::evalStiff (level: " << m_level << ")"
           << CRD::end;
  int termFlags = LevelCNSOp::Terms::Stiff;
  this->evalRHS(a_RHS,
                a_JU,
                a_stage,
                a_stageTime,
                a_stageWeight,
                a_fnFluxRegister,
                a_crFluxRegister,
                a_dt,
                a_timeOld,
                a_crTimeOld,
                a_crTimeNew,
                a_subcycleParams,
                termFlags,
                a_WOld);
}

/*--------------------------------------------------------------------*/
//  Evaluate the stiff terms of a_RHS (which is d(JU)/dt) at the
//  current time based on a_JU for a stage of Runge-Kutta.
//  Works on a box instead of a LevelData. Is not 4th order.
/** \param[out]  a_RHS  d(JU)/dt
 *  \param[in]   a_JU   Solution at beginning of stage
 *  \param[in]   a_box  Box to compute the stiff term over
 *  \param[in]   a_stageTime
 *                      Time at start of stage
 *  \param[in]   a_dit  Data iterator associated with a_RHS and a_JU
 * *//*-----------------------------------------------------------------*/

void
LevelCNSOp::evalStiff(FArrayBox&         a_RHS,
                      const FArrayBox&   a_JU,
                      const Box&         a_box,
                      const Real         a_stageTime,
                      const DataIterator a_dit,
                      FArrayBox&         a_WOld)
{
  CH_TIME("LevelCNSOp::evalStiff");
  CRD::msg << CRD::fv4 << "LevelCNSOp::evalStiff (level: " << m_level << ")"
           << CRD::end;

  CH_assert(a_RHS.contains(a_box));
  CH_assert(a_JU.contains(a_box));
  //  Set number of species transport equations
  const int numSpecies = CRDparam::g_numSpecies;
  const int numReactions = CRDparam::g_numReactions;

  int numConservative = CRDparam::g_CRDPhysics->numConservative();

  a_RHS.setVal(0.); // Probably not needed

  if((CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf) &&
     (numReactions > 0) && (a_stageTime > CRDparam::g_reactionStartTime)
     && (a_stageTime < CRDparam::g_reactionEndTime))
    {
      // Convert form JU to U for this box with no stencil, not 4th order accurate
      FABSTACKTEMP(UFab, a_box, numConservative);
      CH_assert(a_JU.contains(a_box));
      m_levelGridMetrics.computeValidULowOrder(UFab, a_JU, a_box, a_dit);

      // Compute primitive
      FABSTACKTEMP(WFab, a_box, CRDparam::g_CRDPhysics->numPrimitive());
      WFab.setVal(0.);
      CRDparam::g_CRDPhysics->consToPrim(WFab, UFab, a_box, a_WOld);
      CRDparam::g_CRDPhysics->extraPrimitiveState(WFab, a_box);

      FABSTACKTEMP(invDtFab, a_box, 1);
      FABSTACKTEMP(reactingSource, a_box, numSpecies);

      // This function gives us dS_dU
      CRDparam::g_CRDPhysics->addReactionSource(a_box,
                                                reactingSource,
                                                invDtFab,
                                                WFab,
                                                a_stageTime,
                                                m_level,
                                                m_minChemDt,
                                                m_minChemDtCell);

      FABSTACKTEMP(mappedSource, a_box, numSpecies);

      // Do low order conversion back to JS
      m_levelGridMetrics.computeJULowOrder(mappedSource,
                                           reactingSource,
                                           a_box,
                                           a_dit);

      const int USpecStart =
        CRDparam::g_CRDPhysics->speciesConsInterval().begin();
      const int srcComp = mappedSource.interval().begin();
      a_RHS.copy(mappedSource, a_box, srcComp, a_box, USpecStart, numSpecies);
    }
}

/*--------------------------------------------------------------------*/
// Solves the nonlinear problem from ARK using the solver selected during
// initialization/definition
/** \param[out]  a_newSoln   Solution of the nonlinear problem (JU^i)
 *  \param[in]   a_prevStageSoln
 *                           Solution of the previous ARK stage (JU^{i-1})
 *  \param[in]   a_prevTimeSoln
 *                           Previous time step solution (JU^n)
 *  \param[in]   a_rhs       The NonStiff portion of the RHS of the nonlinear
 *                           problem (X^i)
 *  \param[in]   a_stage     Current stage in the ARK method
 *  \param[in]   a_time      Current solution time
 *  \param[in]   a_stageWeight
 *                           Contribution of this stage to the final flux
 *                           over the time step
 *  \param[in]   a_finerFluxRegister
 *                           Flux register between this and a finer level
 *  \param[in]   a_coarserFluxRegister
 *                          Flux register between this and a coarser level
 *  \param[in]   a_dt       Time step on this level.  This is for the
 *                          complete time step, not just the stage.
 *  \param[in]   a_timeOld
 *                          Time at start of time step on this level
 *  \param[in]   a_crTimeOld
 *                          Time at start of time step on coarser level
 *  \param[in]   a_crTimeNew
 *                          Time at end of time step on coarser level
 *  \param[in]   a_subcycleParams
 *                          Parameters from subcyling algorithm (see
 *                          AMRLevel.H for more info)
 *//*-----------------------------------------------------------------*/
void
LevelCNSOp::solve(SOLN&              a_newSoln,
                  const SOLN&        a_prevStageSoln,
                  const SOLN&        a_prevTimeSoln,
                  const RHS&         a_rhs,
                  int                a_stage,
                  Real               a_time,
                  Real               a_stageweight,
                  LevelFluxRegister& a_finerFluxRegister,
                  LevelFluxRegister& a_coarserFluxRegister,
                  const Real         a_dt,
                  const Real         a_timeOld,
                  const Real         a_timeCoarseOld,
                  const Real         a_timeCoarseNew,
                  SubcycleParams     a_subcycleParams,
                  SOLN&              a_WOld)
{
  CH_TIME("LevelCNSOp::solve");
  CRD::msg << CRD::fv3 << "LevelCNSOp::solve (level: " << m_level
           << ")" << CRD::end;
  m_nlsolver->solve(a_newSoln,
                    a_prevStageSoln,
                    a_prevTimeSoln,
                    a_rhs,
                    a_stage,
                    a_time,
                    a_stageweight,
                    a_finerFluxRegister,
                    a_coarserFluxRegister,
                    a_dt,
                    a_timeOld,
                    a_timeCoarseOld,
                    a_timeCoarseNew,
                    a_subcycleParams,
                    this,
                    a_WOld);
}

/*--------------------------------------------------------------------*/
// Calculates the reaction source term Jacobian dS/dU for a single cell
/** \param[out]  a_reactionJacobian The resultant Jacobian
 *  \param[in]   a_JU        The state to use for Jacobian calcs
 *  \param[in]   a_box       Box to compute the reaction Jacobian over
 *  \param[in]   a_stageTime
 *                           Contribution of this stage to the final flux
 *                           over the time step
 *  \param[in]   a_dt        Time step on this level.  This is for the
 *                           complete time step, not just the stage.
 *  \param[in]   a_dit       Data iterator associated with a_RHS and a_JU
 *//*-----------------------------------------------------------------*/
void
LevelCNSOp::calcRxnJacobian(FArrayBox&          a_reactionJacobian,
                            const FArrayBox&    a_JU,
                            const Box&          a_box,
                            const Real          a_stageTime,
                            const Real          a_dt,
                            const DataIterator& a_dit,
                            FArrayBox&          a_WOld)
{
  CH_TIME("LevelCNSOp::calcRxnJacobian");
  CRD::msg << CRD::fv4 << "LevelCNSOp::calcRxnJacobian (level: " << m_level
           << ")" << CRD::end;

  CH_assert(a_reactionJacobian.contains(a_box));
  CH_assert(a_JU.contains(a_box));
  //  Set number of species transport equations
  const int numReactions = CRDparam::g_numReactions;
  const int numPrimitive = CRDparam::g_CRDPhysics->numPrimitive();
  const int numConservative = CRDparam::g_CRDPhysics->numConservative();
  const Interval jacIval = a_reactionJacobian.interval();
  a_reactionJacobian.setVal(0., a_box, jacIval.begin(), jacIval.size());

  if((CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf) &&
     (numReactions > 0) && (a_stageTime > CRDparam::g_reactionStartTime)
     && (a_stageTime < CRDparam::g_reactionEndTime))
    {
      // Convert form JU to U for this box with no stencil, not 4th order accurate
      FArrayBox UFab(a_box, numConservative);
      CH_assert(a_JU.contains(a_box));
      m_levelGridMetrics.computeValidULowOrder(UFab, a_JU, a_box, a_dit);

      // Compute primitive
      FArrayBox WFab(a_box, numPrimitive);
      WFab.setVal(0.);
      CRDparam::g_CRDPhysics->consToPrim(WFab, UFab, a_box, a_WOld);
      CRDparam::g_CRDPhysics->extraPrimitiveState(WFab, a_box);

      // This function gives us dS_dU
      CRDparam::g_CRDPhysics->computeReactionJacobian(a_box,
                                                      a_reactionJacobian,
                                                      WFab,
                                                      a_dt);
    }
}

/*--------------------------------------------------------------------*/
//  Increments a_lhs solution by a_rhs right hand side multiplied by a_scale
//  That is, does a_lhs += a_scale*a_rhs
/** \param[out]  a_lhs   Solution that will be incremented
 *  \param[in]   a_rhs   The data that will be added to a_lhs
 *  \param[in]   a_scale A scalar that scales the rhs
 *//*-----------------------------------------------------------------*/
void
LevelCNSOp::increment(SOLN&       a_lhs,
                      const RHS&  a_rhs,
                      const Real& a_scale)
{
  CH_TIME("LevelCNSOp::increment");
  const int nComp = a_lhs.nComp();
  CH_assert(a_rhs.nComp() == nComp);
  CH_assert(a_lhs.getBoxes() == a_rhs.getBoxes());
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& solnFab = a_lhs[dit];
      const FArrayBox& dataFab = a_rhs[dit];
      const Box& box = dataFab.box();
      CH_assert(solnFab.box().contains(box));
      solnFab.plus(dataFab,
                   box,
                   box,
                   a_scale,
                   0, // srccomp
                   0, // destcomp
                   nComp);
      // if (box.contains(IntVect{ 659, 2, 23 }))
      //   {
      //     CRD::msg << "!!-ARK add   "
      //              << dataFab(IntVect{ 659, 2, 23 }, 9)*a_scale << CRD::end;
      //     CRD::msg << "!!-ARK toget " << solnFab(IntVect{ 659, 2, 23 }, 9)
      //              << CRD::end;
      //   }
    }
}

/*--------------------------------------------------------------------*/
//  Update the solution (a_JU += dt*a_RHS)
/** This is used to update the RK intermediate and the new solution
 *  \param[out] a_JU    Update this solution
 *  \param[in]  a_RHS   RHS
 *  \param[in]  a_dt    time step for update
 *  \param[in]  a_stage RK stage (for debugging)
 *  \param[in]  a_toNewSoln
 *                      T: update to solution JU^(n+1) at end of RK
 *                         procedure
 *                      F: update to RK intermediate from JU^n
 *  \param[in]  a_oldSoln
 *                      JU^(n)
 *  The last 3 parameters are only for debugging and may not be set by
 *  all RK methods.  If not set, the defaults are -1, true, nullptr.
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::updateODE(SOLN&       a_JU,
                      const RHS&  a_RHS,
                      Real        a_dt,
                      const int   a_stage,
                      bool        a_toNewSoln,
                      const SOLN* a_oldSoln) const
{
  CH_TIME("LevelCNSOp::updateODE");
  CRD::msg << CRD::fv3 << "LevelCNSOp::updateODE (level: " << m_level << ")"
           << CRD::end;

  const int nComp = a_JU.nComp();
  CH_assert(nComp == a_RHS.nComp());

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_boxes[dit];
      const int idxBlk = m_boxes.blockIndex(dit);
      CH_assert(idxBlk != -1);
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(idxBlk);
      FArrayBox& JUfab = a_JU[dit];
      const FArrayBox& RHSfab = a_RHS[dit];

      JUfab.plus(RHSfab, disjointBox, disjointBox, a_dt, 0, 0, nComp);
      // Do not want...?  No because this is performed at the end of a full time
      // step and that should be sufficient
      if (CRDparam::g_numSpecies != 0)
        {
          m_patchOp.speciesCorrection(disjointBox,
                                      blockDomain,
                                      JUfab,
                                      m_dx);
        }
      // if (a_stage == 3 && a_toNewSoln
      //     && disjointBox.contains(IntVect{ 284, 11 }))
      //   {
      //     CRD::msg << "!!-Total: "
      //              << (JUfab(IntVect{ 284, 11 }, CRDparam::g_CRDPhysics
      //                       ->speciesConsInterval().begin() + 4) -
      //       (*a_oldSoln)[dit](IntVect{ 284, 11 }, CRDparam::g_CRDPhysics
      //                         ->speciesConsInterval().begin() + 4))/(6*a_dt)
      //              << CRD::end;
      //   }
    }
}

/*--------------------------------------------------------------------*/
//  Define a_newJU to match a_JU, including ghost cells
/**
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::defineSolnData(SOLN&       a_newJU,
                           const SOLN& a_JU) const
{
  CH_TIME("LevelCNSOp::defineSolnData");
  a_newJU.define(a_JU.getBoxes(), a_JU.nComp(), a_JU.ghostVect());
}

/*--------------------------------------------------------------------*/
//  Define type a_newRHS based on a_JU, including required ghost cells
/** \note
 *  <ul>
 *    <li> No ghost cells are defined in a_newRHS
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::defineRHSData(RHS&        a_newRHS,
                          const SOLN& a_JU) const
{
  CH_TIME("LevelCNSOp::defineRHSData");
  a_newRHS.define(a_JU.getBoxes(), a_JU.nComp(), IntVect::Zero);
  for(DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rhsFab = a_newRHS[dit];
      rhsFab.setVal(0.);
    }
}


/*--------------------------------------------------------------------*/
//  Copy a_srcJU to a_dstJU
/**
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::copySolnData(SOLN&       a_dstJU,
                         const SOLN& a_srcJU) const
{
  CH_TIME("LevelCNSOp::copySolnData");
  CRD::msg << CRD::fv3 << "LevelCNSOp::copySolnData (level: "
           << m_level << ")" << CRD::end;
  a_srcJU.copyTo(a_dstJU);
}


/*--------------------------------------------------------------------*/
//  Fill in ghost cells of a_U if doing a restart
/**
 *
 *  \param[in]  a_U     Conservative state in physical space
 *  \param[out] a_U     1 layer of ghost cells filled.
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::restartExchange(LevelData<FArrayBox>& a_U) const
{
  CRD::msg << CRD::fv3 << "LevelCNSOp::restartExchange (level: "
           << m_level << ")" << CRD::end;
  a_U.exchange(m_UExchangeCopier);
}

/*--------------------------------------------------------------------*/
//  Compute global sums for turbulence forcing
/** \param[in]  a_WcellAvgLvl     
 *                      Cell-centered primitive state from which to compute
 *                      the global sums
 *  \param[out] a_globalKE
 *                      The domain-summed kinetic energy
 *  \param[out] a_globalHelicity
 *                      The domain-summed helicity
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::computeTurbForcingSums(
  const LevelData<FArrayBox>& a_WcellPntLvl,
  Real&                       a_globalKE,
  Real&                       a_globalHelicity) const
{
  CH_TIME("LevelCNSOp::computeTurbForcingSums");
  if (!CRDparam::g_useTurbForce) { return; }
  Real locMeanKE = 0.;
  Real locMeanHelicity = 0.;
  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();
  const Real vol = m_dx.product()/(CRDparam::g_domainLength.product());
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      
      const Box& box = a_WcellPntLvl.disjointBoxLayout()[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;
      const FArrayBox& WcellPntFab = a_WcellPntLvl[dit];
      // (1) Compute and volume-weighted-sum the cell-averaged kinetic energy
      FABSTACKTEMP(cellPntKE, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const Real rho = WcellPntFab[MD_IX(i, cRho)];
          const RealVect vel(D_DECL(WcellPntFab[MD_IX(i, cVel)],
                                    WcellPntFab[MD_IX(i, cVel+1)],
                                    WcellPntFab[MD_IX(i, cVel+2)]));
          cellPntKE[MD_IX(i, 0)] = 0.5*rho*(vel.radSquared());
        }
      FABSTACKTEMP(cellAvgKE, box, 1);
      // If higher-order, convolve the kinetic energy
      int order = ((CRDparam::g_cellDeconvolveFlatten < 2) &&
                   (CRDparam::g_cellConvolveFlatten < 2)) ? 4 : 2;
      int convolveSign = -1;
      CRDutil::deconvolve(cellAvgKE, cellPntKE, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum the cell-averaged kinetic energy
      MD_BOXLOOP(box, i)
        {
          locMeanKE += vol*cellAvgKE[MD_IX(i, 0)];
        }

      // (2) Compute the cell-centered velocity gradient
      const bool fourthOrder = (order == 4) ? true : false;
      FABSTACKTEMP(cellPntVelGrad, box1Dom, SpaceDim*SpaceDim);
      for (int row = 0; row != SpaceDim; ++row)
        {
          for (int col = 0; col != SpaceDim; ++col)
            {
              const int valComp = cVel + row;
              const int derivComp = row + col*SpaceDim;
              PatchMappedFunc::cellAvgDerivFromCellAvgCS(
                cellPntVelGrad, WcellPntFab, valComp, derivComp, blockDomain,
                box1Dom, m_levelGridMetrics.dxVect(), col, fourthOrder, false);
            }
        }

      // Map the cell-centered velocity gradient
      FABSTACKTEMP(cellPntPhysVelGrad, box1Dom, SpaceDim*SpaceDim);
      PatchMappedFunc::gradientCStoPS(
        box1Dom, cellPntPhysVelGrad, cellPntVelGrad,
        m_levelGridMetrics.m_cellNtJ[dit]);

      // Compute the cell-centered vorticity
      const int numVorticityComp = PatchMappedFunc::m_numCurlComps;
      FABSTACKTEMP(cellPntVorticity, box1Dom, numVorticityComp);
      PatchMappedFunc::curlPS(box1Dom, cellPntVorticity, cellPntPhysVelGrad);

      // (3) Compute the cell-averaged helicity
      FABSTACKTEMP(cellPntHelicity, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const RealVect vorticity(D_DECL(cellPntVorticity[MD_IX(i, 0)],
                                          cellPntVorticity[MD_IX(i, 1)],
                                          cellPntVorticity[MD_IX(i, 2)]));
          const RealVect vel(D_DECL(WcellPntFab[MD_IX(i, cVel)],
                                    WcellPntFab[MD_IX(i, cVel+1)],
                                    WcellPntFab[MD_IX(i, cVel+2)]));
          cellPntHelicity[MD_IX(i, 0)] = vorticity.dotProduct(vel);
        }
      // Convolve the cell-centered helicity
      FABSTACKTEMP(cellAvgHelicity, box1Dom, 1);
      CRDutil::deconvolve(cellAvgHelicity, cellPntHelicity, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum the cell-averaged helicity
      MD_BOXLOOP(box, i)
        {
          locMeanHelicity += vol*cellAvgHelicity[MD_IX(i, 0)];
        }
    }

  Real meanGlobalKE = locMeanKE;
  Real meanGlobalHelicity = locMeanHelicity;
#ifdef CH_MPI
  MPI_Allreduce(&locMeanKE, &meanGlobalKE, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanHelicity, &meanGlobalHelicity, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  a_globalKE = meanGlobalKE;
  a_globalHelicity = meanGlobalHelicity;
}

/*--------------------------------------------------------------------*/
//  Set wall-model data from finer level
/** \param[in]  a_JU    \<\JU\> state on this level
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::setWallModelFromFiner(LevelData<FArrayBox>& a_JU)
{
  CH_TIME("LevelCNSOp::setWallModelFromFiner");
  if (!(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex) ||
      !(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
    {
      return;
    }
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);
      FArrayBox& JUFab = a_JU[dit];
      CRDparam::g_CRDPhysics->setWallModelFromFiner(
        JUFab, blockDomain, box, m_levelGridMetrics, dit());
    }
}

/*--------------------------------------------------------------------*/
//  Compute coarsest level LES SGS KE for single-block problems
/** \param[in]  a_null  Nothing input yet
 *  \return             Nothing yet
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::estimateCrsLevSingleBlockSGSKE(LevelData<FArrayBox>& a_JU)
{
  CH_TIME("LevelCNSOp::estimateCrsLevSingleBlockSGSKE");
  const int sgsRefineRatio = CRDparam::g_sgskeCrsLevFilterRatio;
  if (m_levelGridMetrics.isMultiBlock() || m_hasCoarserGrid ||
      !CRDparam::g_useSGSCoarsening || (sgsRefineRatio == 1) ||
      !(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex))
    {
      return;
    }
  // Coarsen the problem domain
  const ProblemDomain& problemDomain =
    m_levelGridMetrics.getCoordSys().problemDomain(0);
  ProblemDomain coarseProbDom = problemDomain;
  coarseProbDom.coarsen(sgsRefineRatio);

  // Loop over the boxes and compute mapped gradients
  const int tensorComps = SpaceDim*SpaceDim;
  const int sgsKEcomp = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);

      // Set up some boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      const FArrayBox& JFab = m_levelGridMetrics.m_J[dit];
      const FArrayBox& cellAvgW = m_WcellAvgLvl[dit];
      FArrayBox& cellAvgPartW = m_partJWcellAvgLvl[dit];

      // Fill cellAvgPartW with density and velocity components
      cellAvgPartW.copy(cellAvgW, 0, 0, SpaceDim+1);

      // Get the dx vector on this block
      RealVect dxVect = m_levelGridMetrics.getCoordSys(box)->dx();

      // 1) compute <grad W>_i from <W>_i+1/2
      FABSTACKTEMP(cellAvgVelGrad, box2Dom, tensorComps);
      PatchMappedFunc::cellAvgGradFromFaceAvgCS(
        cellAvgVelGrad, m_WfaceAvgLvl[dit],
        CRDparam::g_CRDPhysics->velocityInterval(), 0, box2Dom, dxVect);

      // 2) compute grad W_i from <grad W>_i
      FABSTACKTEMP(cellPntVelGrad, box1Dom, tensorComps);
      const Interval velGradIntv = cellAvgVelGrad.interval();
      CRDutil::deconvolve(
        cellPntVelGrad, cellAvgVelGrad, box1Dom, blockDomain, velGradIntv,
        4, 1, false);

      // 3) compute grad W N^t/J
      // Note that for the sake of lessening overhead, we reuse cellAvgVelGrad
      // as space in which to compute cellPntVelPhysGrad
      const FArrayBox& NtJFab = m_levelGridMetrics.m_cellNtJ[dit];
      PatchMappedFunc::gradientCStoPS(
        box1Dom, cellAvgVelGrad, cellPntVelGrad, NtJFab);

      // 4) compute J*grad W N^t/J
      const FArrayBox& JpntFab = m_cellPntJ[dit];
      for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
        {
          MD_BOXLOOP(box1Dom, i)
            {
              cellAvgVelGrad[MD_IX(i, comp)] *= JpntFab[MD_IX(i, 0)];
            }
        }

      // 5) compute <J grad W N^t/J>_i
      // Note that for the sake of lessening overhead, we reuse cellPntVelGrad
      // as space in which to compute cellAvgJWGradNtJ
      CRDutil::deconvolve(
        cellPntVelGrad, cellAvgVelGrad, box, blockDomain, velGradIntv, 4,
        -1, false);

      // 6) fill LevelData with <J grad W N^t/J>
      m_cellAvgJWGradNtJLvl[dit].copy(cellPntVelGrad);

      // 7) compute <JW> from <W> and <J> (density and velocity)
      const int cellAvgJWComps = 1 + SpaceDim;
      FABSTACKTEMP(cellAvgJW, box, cellAvgJWComps);
      fourthOrderCellProd(
        cellAvgJW, cellAvgPartW, JFab, box, blockDomain, true);

      // 8) fill LevelData with <JW>
      m_partJWcellAvgLvl[dit].copy(cellAvgJW);
    }

  // 9)  average-down <JW> from the base-mesh to the coarser mesh
  m_cellAvgJWAvgOp.averageToCoarse(m_crsCellAvgJW, m_partJWcellAvgLvl);

  // 10) average-down <J grad W N^t/J> from the base-mesh to the coarser mesh
  m_cellAvgJPhysGradWAvgOp.averageToCoarse(m_crsCellAvgJPhysGrad,
                                           m_cellAvgJWGradNtJLvl);

  // 11) exchange <JW>
  m_crsCellAvgJW.exchange(m_crsCellAvgJWExchange);

  // 12) exchange <J grad W N^t/J>
  m_crsCellAvgJPhysGrad.exchange(m_crsCellAvgJGradWExchange);

  // 13) compute SGSKE
  CRDparam::g_CRDPhysics->crsLevCellAvgSGSKineticEnergy(m_crsCellAvgJSGSKE,
                                                        m_crsCellAvgJW,
                                                        m_crsCellAvgJPhysGrad,
                                                        m_crsCellPntJ,
                                                        m_crsXLvl,
                                                        m_crsDeltaC,
                                                        coarseProbDom);

  // 14) exchange coarse <J-SGSKE>
  m_crsCellAvgJSGSKE.exchange(m_crsCellAvgSGSKEExchange);

  // 15) interpolate <J-SGSKE> from coarse to fine
  m_crsToFnSGSKEInterpolator.interpToFine(m_fnCellAvgJSGSKE,m_crsCellAvgJSGSKE);

  // 16) exchange fine <J-SGSKE>
  m_fnCellAvgJSGSKE.exchange(m_fnCellAvgJSGSKEExchange);

  // 17) compute <SGSKE> on current level ("fine" mesh)
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);

      // Set up some boxes
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Deconvolve <J-SGSKE>
      const FArrayBox& JSGSKE = m_fnCellAvgJSGSKE[dit];
      FABSTACKTEMP(cellPntSGSKE, box1Dom, 1);
      CRDutil::deconvolve(
        cellPntSGSKE, JSGSKE, box1Dom, blockDomain, Interval(0,0), 4, 1, false);

      // Divide J-SGSKE by J
      const FArrayBox& JpntFab = m_cellPntJ[dit];
      MD_BOXLOOP(box1Dom, i)
        {
          cellPntSGSKE[MD_IX(i, 0)] /= JpntFab[MD_IX(i, 0)];
        }

      // Convolve SGSKE
      FArrayBox& SGSKE = m_fnCellAvgSGSKE[dit];
      CRDutil::deconvolve(
        SGSKE, cellPntSGSKE, box, blockDomain, Interval(0,0), 4, -1, false);
    }

  // 18) exchange fine <SGSKE>
  m_fnCellAvgSGSKE.exchange(m_fnCellAvgSGSKEExchange);

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(box);

      // Set up some boxes
      Box box3Dom = grow(box, 3);
      box3Dom &= blockDomain;

      // 19) interpolate fine <SGSKE> to faces
      const FArrayBox& cellAvgSGSKE = m_fnCellAvgSGSKE[dit];
      FLUXBOXSTACKTEMP(faceAvgSGSKE, box3Dom, 1);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          FArrayBox& faceAvgSGSKEFab = faceAvgSGSKE[dir];
          Box faceBox = box3Dom;
          faceBox.surroundingNodes(dir);
          PatchMappedFunc::faceAvgValFromCellAvgCS(
            faceAvgSGSKEFab, cellAvgSGSKE, 0, 0, blockDomain, faceBox,
            dir, 4, false);
        }
      // 20) copy <SGSKE>_i+1/2 to m_WfaceAvgLvl
      m_WfaceAvgLvl[dit].copy(faceAvgSGSKE, 0, sgsKEcomp, 1);
      // 21) copy <SGSKE>_i to m_WcellAvgLvl
      m_WcellAvgLvl[dit].copy(cellAvgSGSKE, 0, sgsKEcomp, 1);
    }

  // 22) copy to fine <JU>
  m_fnCellAvgJSGSKE.copyTo(Interval(0,0), a_JU, Interval(sgsKEcomp,sgsKEcomp));
}

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
/*-------------------------------------------------------------------------*//**
 *  \brief Add spectral forcing field to solution momentum field on all levels
 *  \param[inout] a_U       Solution LevelData
 *  \param[in]    a_Source  Spectral forcing field generated by
 *                          SpectralForcing::calcSpectralForce() and stored as
 *                          AMRLevelCNS member LevelData.
 *  \param[in]    a_dt      time step of current (already completed) level
 *                          advance
 *  \return                 The total injected kinetic energy
 *//*-------------------------------------------------------------------------*/

Real
LevelCNSOp::addSpectralForce(LevelData<FArrayBox>&    a_U,
                             LevelData<FArrayBox>&    a_Source,
                             Real                     a_dt)
{
  CH_TIME("LevelCNSOp::addSpectralForce");
  CRD::msg << CRD::fv2 << "LevelCNSOp::addSpectralForce" << CRD::end;

  // level-specific grid size (if entire domain was covered by this level!)
  const Real inv_ncells = m_dx.product() / CRDparam::g_domainLength.product();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int velIndx = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();

  // -----------------------------------------------------------------
  // Subtract the mean momentum perturbation from a_Source
  Real localMean[4] = {0., 0., 0., 0.};
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      FArrayBox& UFab = a_U[dit];
      FArrayBox& SFab = a_Source[dit];
      MD_ARRAY_RESTRICT(arrU, UFab);
      MD_ARRAY_RESTRICT(arrS, SFab);
      MD_BOXLOOP(box, i)
        {
          const Real sqrtrho = std::sqrt(arrU[MD_IX(i, rhoIndx)]);
          localMean[0] += sqrtrho * arrS[MD_IX(i, 0)];
          localMean[1] += sqrtrho * arrS[MD_IX(i, 1)];
          localMean[2] += sqrtrho * arrS[MD_IX(i, 2)];
          localMean[3] += sqrtrho;
        }
    }

  Real globalMean[4];
  MPI_Allreduce(&localMean, &globalMean, 4,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  globalMean[0] /= globalMean[3];
  globalMean[1] /= globalMean[3];
  globalMean[2] /= globalMean[3];

  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      FArrayBox& SFab = a_Source[dit];
      MD_ARRAY_RESTRICT(arrS, SFab);
      MD_BOXLOOP(box, i)
        {
          arrS[MD_IX(i, 0)] -= globalMean[0];
          arrS[MD_IX(i, 1)] -= globalMean[1];
          arrS[MD_IX(i, 2)] -= globalMean[2];
        }
    }

  // -----------------------------------------------------------------
  // Calculate the unscaled kinetic energy of the perturbations

  localMean[0] = 0.;
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      FArrayBox& UFab = a_U[dit];
      FArrayBox& SFab = a_Source[dit];
      MD_ARRAY_RESTRICT(arrU, UFab);
      MD_ARRAY_RESTRICT(arrS, SFab);
      MD_BOXLOOP(box, i)
        {
          localMean[0] += (  arrU[MD_IX(i, velIndx)]   * arrS[MD_IX(i, 0)]
                           + arrU[MD_IX(i, velIndx+1)] * arrS[MD_IX(i, 1)]
                           + arrU[MD_IX(i, velIndx+2)] * arrS[MD_IX(i, 2)])
                          / std::sqrt(arrU[MD_IX(i, rhoIndx)]);
        }
    }

  MPI_Allreduce(&localMean[0], &globalMean[0], 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  globalMean[0] *= inv_ncells;

  // -----------------------------------------------------------------
  // Apply the perturbations with the correct KE injection rate

  // note that globalMean is a sum, and not an average, this is important!
  const Real scale = (CRDparam::g_spectralForcingEps * a_dt) / globalMean[0];
  CRD::msg << CRD::fv2 << "Scaling factor \n" << scale << CRD::var;

  localMean[0] = 0.;
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = m_boxes[dit];
      FArrayBox& UFab = a_U[dit];
      FArrayBox& SFab = a_Source[dit];
      MD_ARRAY_RESTRICT(arrU, UFab);
      MD_ARRAY_RESTRICT(arrS, SFab);
      MD_BOXLOOP(box, i)
        {
          Real& rho = arrU[MD_IX(i, rhoIndx)];
          Real& M0 = arrU[MD_IX(i, velIndx)];
          Real& M1 = arrU[MD_IX(i, velIndx+1)];
          Real& M2 = arrU[MD_IX(i, velIndx+2)];

          const Real EkinOld = (M0*M0 + M1*M1 + M2*M2) / rho;

          const Real sqrtrho = std::sqrt(rho);
          M0 += scale * sqrtrho * arrS[MD_IX(i, 0)];
          M1 += scale * sqrtrho * arrS[MD_IX(i, 1)];
          M2 += scale * sqrtrho * arrS[MD_IX(i, 2)];

          const Real EkinNew = (M0*M0 + M1*M1 + M2*M2) / rho;

          localMean[0] += EkinNew - EkinOld; // accumulation of injected KE
        }
    }

  globalMean[0] = 0.;
  MPI_Allreduce(&localMean[0], &globalMean[0], 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  const Real prod_rate = 0.5 * globalMean[0] * inv_ncells / a_dt;
  CRD::msg << CRD::fv2 << "KE inj. rate [J/(m^3 s)] \n" << prod_rate << CRD::var;

  const Real injectedKE = 0.5 * globalMean[0] * m_dx.product();
  CRD::msg << CRD::fv2 << "injected KE [J]\n" << injectedKE << CRD::var;
  return injectedKE;
}
#endif  /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */

/*--------------------------------------------------------------------*/
//  Computes the max norm: ||JU-\hat{JU}||_inf
/** Computes the max norm of JU - \hat{JU} for use in PID step size control
 *  calculations. Computes the max norm on this processor, the global
 *  max norm is found in computeNewDt.
 *
 *  Stores the result in member data so that it may be used by other
 *  functions of this class.
 *
 *  \param[in] a_JU     Conservative state in computational space
 *  \param[in] a_JUhat  Conservative state in computational space that
 *                      is associated with the embedded ARK4 method
 *//*-----------------------------------------------------------------*/
Real LevelCNSOp::computeErrorEstimate(const LevelData<FArrayBox>& a_JU,
                                      const LevelData<FArrayBox>& a_JUhat) const
{
  CH_TIME("LevelCNSOp::computeUminusUhatNorm");
  CRD::msg << CRD::fv2 << "LevelCNSOp::computeUminusUhatNorm (level: " << m_level
           << ")" << CRD::end;

  const Interval& interval = CRDparam::g_CRDPhysics->speciesConsInterval();
  Real localMaxNorm = -1.E+99;
  // CRD::msg << "  Name   FabNorm    MaxNorm" << CRD::end;
  for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& JUdata = a_JU[dit];
      const FArrayBox& JUHatdata = a_JUhat[dit];
      const Box& disjointBox = m_boxes[dit];
      Real fabMaxNorm = m_patchOp.computeUminusUhatMaxNorm(JUdata,
                                                           JUHatdata,
                                                           disjointBox,
                                                           interval);
      localMaxNorm = std::max(localMaxNorm, fabMaxNorm);
    }

  // We need to find the global localMaxNorm
  struct ReduceLoc
  {
    Real maxNorm;
    int rank;
  };
  ReduceLoc global[1];
  global[0].maxNorm = localMaxNorm;
  global[0].rank = procID();

#ifdef CH_MPI
  ReduceLoc local[1];
  local[0].maxNorm = localMaxNorm;
  local[0].rank = procID();
  MPI_Allreduce(local,
                global,
                1,
                MPI_CH_REAL_INT,
                MPI_MAXLOC,
                Chombo_MPI::comm);
#endif

  Real maxNormGlobal = global[0].maxNorm;
  m_UminusUhat.push(maxNormGlobal);

  return maxNormGlobal;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Fill in ghost cells of m_U to RK4 intermediates
/** The fourth-order routine requires 5 layers of ghosts of \<U\>.
 *  This routine fills valid ghosts, invalid ghosts, and extra-block
 *  ghosts.  Notably, ghosts outside the problem domain are NOT set by
 *  this routine.  I.e., there is no application of BC.  These are the
 *  steps to fill the ghost cells:
 *  <ol>
 *    <li> Fill invalid ghost cells across coarse-fine interfaces.
 *         This fill sets \<U\> in the ghost cells and \<JU\> in one
 *         layer of cells
 *    <li> Fill one layer of \<JU\> in valid ghost cells.
 *    <li> Compute \<U\> in all valid cells
 *    <li> Fill the valid ghosts of \<U\> by exchange
 *    <li> Fill the extra-block ghosts of \<U\> using multiblock
 *         exchange
 *  </ol>
 *  U itself is member data for the class LevelCNSOp.  This may change
 *  when a switch to DenseRK4Output is made.
 *
 *  \param[in]  a_JU    Conservative state in computational space
 *  \param[out] a_JU    1 layer of ghost cells filled.
 *  \param[in]  a_stage RK4 stage index [0-3]
 *  \param[in]  a_timeOld
 *                      Time at beginning of time-step on this level
 *                      (not the stage time)
 *  \param[in]  a_crTimeOld
 *                      Time at beginning of coarser time step
 *  \param[in]  a_crTimeNew
 *                      Time at end of coarser time step
 *  \param[in]  a_subcycleParams
 *                      Parameters from subcyling algorithm (see
 *                      AMRLevel.H for more info)
 *//*-----------------------------------------------------------------*/

void
LevelCNSOp::fillGhostsRK4AndComputeU(
  LevelData<FArrayBox>& a_JU,
  const int             a_stage,
  const Real            a_dt,
  const Real            a_timeOld,
  const Real            a_crTimeOld,
  const Real            a_crTimeNew,
  SubcycleParams        a_subcycleParams)
{
  CH_TIME("LevelCNSOp::fillGhostsRK4AndComputeU");
  CRD::msg << CRD::fv3 << "LevelCNSOp::fillGhostsRK4andComputeU (level: "
           << m_level << ")" << CRD::end;

//**FIXME: I think what happens for MB is that all block ghosts are
//**       interpolated and then possibly overwritten by an exchange if there
//**       is data on the other side.  We may want to change this once
//**       everything is working correctly.

#if 0
  // set ghost cells to bad values so it is apparent exchanges happen correct
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      MD_BOXLOOP(m_U[dit].box(), i)
        {
          if (!m_boxes[dit].contains(MD_GETIV(i)))
            {
              m_U[dit][MD_IX(i, 0)] = 123456789e300;
            }
        }
    }
#endif

//--Fill invalid ghosts of <U> (and 1 layer of <JU>) if there is a coarser
//--level

  if (m_hasCoarserGrid)
    {
      // If this is stage 1, then the coarsened-fine data was set in
      // AMRLevelCNS::advance.
      if (a_stage > 0)
        {
          // Start of time-step (*NOT* stage) as time interpolation
          // coefficient in range [0:1] of coarse time-step
          Real crDt = a_crTimeNew - a_crTimeOld;
          Real alpha = (a_timeOld - a_crTimeOld)/crDt;
          if (a_subcycleParams.isFirstSubcycle)
            {
              alpha = 0.;
            }
          CH_assert((alpha >= 0.) && (alpha < 1.));

          // This will internally store the coarsened-fine data used to fill
          // ghosts cells for this level.
          const Real dtRatio = a_dt/crDt;
          m_levelGridMetrics.timeIntermediate(
            *m_timeInterpolator,
            alpha,
            dtRatio,
            a_stage,
            CRDparam::g_CRDPhysics->velocityInterval(),
            a_JU,
            true);
        }

      // LES subgrid-scale calculation
      if ((CRDparam::g_useSGSCoarsening) &&
          m_levelGridMetrics.isMultiBlock() &&
	  !(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
        {
          LevelData<FArrayBox>& crLevU = m_levelGridMetrics.presetCrLevU();
          LevelData<FArrayBox>& crLevJU =
            m_levelGridMetrics.presetCr2ThisInterpolatorCrFnLevJU();
          const LevelGridMetrics& coarserLGM =
                *(m_levelGridMetrics.getCoarserLevelGridMetrics());
          for (DataIterator dit = crLevU.getBoxes().dataIterator();
               dit.ok(); ++dit)
            {
              const Box disjointBox = crLevU.getBoxes()[dit];
              CRDparam::g_CRDPhysics->cellAvgSGSKineticEnergy(
                crLevU[dit],
                crLevJU[dit],
                coarserLGM,
                dit(),
                false,
                disjointBox);
            }
          //**NOTE: The following is a hack necessary for SGS KE coarsening
          //        with AMR
          // Exchange crLevU again
          crLevU.exchange();
          // Fill all extra-block ghosts of crLevU again
          coarserLGM.multiblockExchangeU(
            crLevU, CRDparam::g_CRDPhysics->velocityInterval());
          // Copy crLevU to crFnLevU again
          crLevU.copyTo(m_levelGridMetrics.presetCr2ThisInterpolatorCrFnLevU(),
                        m_levelGridMetrics.crFnLevUCopier());
        }

      // We can now use one routine which will find <U> and fill
      // the ghosts.  Note that this also fills 1 layer of ghosts in <JU>.
      m_levelGridMetrics.fillFineGhostCells(m_U, a_JU);

      if ((CRDparam::g_useSGSCoarsening) &&
          m_levelGridMetrics.isMultiBlock() &&
	  !(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
        {
          // now we have to interpolate the entire KE coarse grid to the fine
          const Interval turbIntv = CRDparam::g_CRDPhysics->turbConsInterval();
          const int sgsKEComp = turbIntv.begin();
          const Interval sgsKEIntv(sgsKEComp,  sgsKEComp);
          m_levelGridMetrics.completeInterpCr2FnComp(
            m_currModelJU,
            m_levelGridMetrics.presetCrLevU(),
            m_levelGridMetrics.presetCr2ThisInterpolatorCrFnLevJU(),
            CRDparam::g_CRDPhysics->velocityInterval(),
            sgsKEIntv);
          m_currModelJU.copyTo(Interval(sgsKEComp, sgsKEComp), a_JU,
                               Interval(sgsKEComp, sgsKEComp));
        }
    }

//--We just filled the invalid ghosts and now have to compute <U> in valid
//--cells, valid ghosts, and extra-block ghosts.  Note that these are the
//--same operations usually performed in MappedLevelData

  // Compute <U> in all valid cells.  For this, we need to exchange 1 layer of
  // ghosts in <JU>
  if (CRDparam::g_cartesian)
    {
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          m_U[dit].copy(a_JU[dit], m_boxes[dit]);
        }
     }
  else
    {
      a_JU.exchange(m_JUExchangeCopier);
      // Only if both convolve and deconvolve are set to 2nd do we perform
      // the inverse product rule at 2nd order
      if (CRDparam::g_cellConvolveFlatten < 2 ||
          CRDparam::g_cellDeconvolveFlatten < 2)
        {
          // Note: block boundaries are treated the same as domain boundaries in
          // that 2nd-order one-sided derivatives are used there.
          m_levelGridMetrics.computeValidU(m_U, a_JU);
        }
      else
        {
          for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
            {
              m_levelGridMetrics.computeValidULowOrder(m_U[dit],
                                                       a_JU[dit],
                                                       m_boxes[dit],
                                                       dit);
            }
        }
    }

//--Fill all the valid ghosts of <U> by exchange

  m_U.exchange(m_UExchangeCopier);

//--Fill the extra-block ghosts of <U>

  m_levelGridMetrics.multiblockExchangeU(
    m_U, CRDparam::g_CRDPhysics->velocityInterval());
}


/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Compute first time step
/** \param[in]  a_U     Average solution state from which to compute
 *                      the new time step.  Must be defined on 1 ghost
 *                      cell everywhere except across domain
 *                      boundaries
 *  \param[in]  a_time  Current solution time
 *  \return             The new time step from all physics'
 *                      constraints
 *//*-----------------------------------------------------------------*/

Real
LevelCNSOp::computeFirstDt(const LevelData<FArrayBox>& a_U,
                           const Real&                 a_time) const
{
  Real dtNewLocal = 1.E9;
  for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = m_boxes[dit];
      const BlockDomain& blockDomain =
        m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box box1 = grow(disjointBox, 1);
      Box box1Dom = box1 & blockDomain;
      FABSTACKTEMP(WFab, box1, CRDparam::g_CRDPhysics->numPrimitive());
      CH_assert(a_U[dit].box().contains(box1Dom));
      CRDparam::g_CRDPhysics->consToPrim(WFab, a_U[dit], box1Dom,  WFab);
      CRDparam::g_CRDPhysics->extraPrimitiveState(WFab, box1Dom);
      // FAB of 1/dt_conv + 1/dt_visc + 1/dt_chem
      FABSTACKTEMP(invDtFab, disjointBox, 1);
      invDtFab.setVal(0.);

//--Inertial

      if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
        {
          // See "High-Order, Finite-Volume Methods in Mapped Coordinates",
          //Colella et al., JCP 2011 for derivation of this constraint.
          const Real stabilityConstraint = 1.3925;
          CRDparam::g_CRDPhysics->getMaxWaveSpeed(blockDomain,
                                                  disjointBox, // Solve box
                                                  disjointBox, // Disjoint box
                                                  invDtFab,
                                                  WFab,
                                                  m_levelGridMetrics.m_N[dit],
                                                  m_levelGridMetrics.m_J[dit],
                                                  m_levelGridMetrics,
                                                  stabilityConstraint,
                                                  m_dx,
                                                  m_minConvDt,
                                                  m_minConvDtCell);
        }

//--Viscous **FIXME: possibly account for thermal conductivity coefficient

      if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
        {
          Real lambdaDMax = 16./3.;
          CRDparam::g_CRDPhysics->ellipticDt(blockDomain,
                                             disjointBox, // Solve box
                                             disjointBox, // Disjoint box
                                             invDtFab,
                                             WFab,
                                             lambdaDMax,
                                             m_levelGridMetrics.m_NtJ[dit],
                                             m_levelGridMetrics,
                                             m_dx,
                                             m_minDiffDt,
                                             m_minDiffDtCell);
        }

//--Chemical

      if ((CRDparam::g_numReactions > 0) &&
          (a_time > CRDparam::g_reactionStartTime) &&
          (a_time < CRDparam::g_reactionEndTime))
        {
          // Number of components that have reactions
          const int numRCT = CRDparam::g_numSpecies;
          // Dummy fab to store the rates
          FABSTACKTEMP(dummyFAB, disjointBox, numRCT);
          // Solves for rho*omega
          CRDparam::g_CRDPhysics->addReactionSource(disjointBox,
                                                    dummyFAB,
                                                    invDtFab,
                                                    WFab,
                                                    a_time,
                                                    m_level,
                                                    m_minChemDt,
                                                    m_minChemDtCell);
        }
      // Invert sum of time steps
      Real minLocalDt = 1.E9;
      MD_ARRAY_RESTRICT(arrDtFab, invDtFab);
      MD_BOXLOOP(disjointBox, i)
        {
          Real invdt = arrDtFab[MD_IX(i,0)];
          if (1./invdt > 1.E-15 && invdt > 1.E-15 && invdt == invdt)
            {
              minLocalDt = std::min(minLocalDt, 1./invdt);
            }
        }
      dtNewLocal = std::min(dtNewLocal, minLocalDt);
    }
  m_firstDt = false;
  return dtNewLocal;
}
