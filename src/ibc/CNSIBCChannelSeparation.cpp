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
 * \file CNSIBCChannelSeparation.cpp
 *
 * \brief Member functions for CNSIBCChannelSeparation
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "UnitNormalsF_F.H"
#include "RootSolver.H"

#include "AMRIO.H" // DEBUG ONLY

//----- Internal -----//

#include "CNSIBCChannelSeparation.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"
#include "PatchCNSOp.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "CRDutil.H"
#include "LoHiCenter.H"

/*******************************************************************************
 *
 * Class CNSIBCRecirculatingInletTFP: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCChannelSeparation::CNSIBCChannelSeparation()
  :
  CNSIBCGeneralized(),
  m_delta(0.),
  m_deltaOut(0.),
  m_phi1(0.5),
  m_phi2(100.),
  m_psi1(50.),
  m_psi2(200.),
  m_lambda(0.),
  m_samplePlaneLoc(0.75),
  m_yPlusMaxGuess(1.e10),
  m_numLevels(0.),
  m_numCellsInletBlock(IntVect_zero),
  m_inletBlockDxVect(RealVect_zero),
  m_baseMeshSize(IntVect_zero),
  m_onlyInlet(0),
  m_useWallModel(0),
  m_newGeometry(0),
  m_deltaRecyPercent(1.)
{
  readBCInfo();

  for (int i = 0; i != m_numLevels; ++i)
    {
      std::unique_ptr<LevelData<FArrayBox>> singleBoxSampledPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxSampledPlane.push_back(std::move(singleBoxSampledPlane));
      std::unique_ptr<LevelData<FArrayBox>> singleBoxInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxInletPlane.push_back(std::move(singleBoxInletPlane));
      std::unique_ptr<LevelData<FArrayBox>> multiBoxInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_multiBoxInletPlane.push_back(std::move(multiBoxInletPlane));
      std::unique_ptr<Copier> singleBoxSampledPlaneCopier =
        std::make_unique<Copier>();
      m_singleBoxSampledPlaneCopier.push_back(
        std::move(singleBoxSampledPlaneCopier));
      std::unique_ptr<Copier> multiBoxInletPlaneCopier =
        std::make_unique<Copier>();
      m_multiBoxInletPlaneCopier.push_back(std::move(multiBoxInletPlaneCopier));
      std::unique_ptr<Copier> multiBoxInletExchangeCopier =
        std::make_unique<Copier>();
      m_multiBoxInletExchangeCopier.push_back(
        std::move(multiBoxInletExchangeCopier));
      std::vector<Real> domainVector;
      m_uVelSTAvg.push_back(domainVector);
      m_uVelSAvg.push_back(domainVector);
      m_vVelSTAvg.push_back(domainVector);
      m_vVelSAvg.push_back(domainVector);
      m_wVelSAvg.push_back(domainVector);
      m_rhoSTAvg.push_back(domainVector);
      m_rhoSAvg.push_back(domainVector);
      m_tempSTAvg.push_back(domainVector);
      m_tempSAvg.push_back(domainVector);
      m_yLoc.push_back(domainVector);
      // DEBUG ONLY
      m_yLocMapped.push_back(domainVector);
      // END DEBUG ONLY
      m_zLoc.push_back(domainVector);
      m_inletMeanInitialized.push_back(0);
      m_levelDefined.push_back(0);
      m_deltaRecy.push_back(0.);
      m_timeCounter.push_back(0);
      m_stageCounter.push_back(0);
      m_currDt.push_back(0.);
      m_prevDt.push_back(0.);
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCChannelSeparation::~CNSIBCChannelSeparation()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/

const char *const
CNSIBCChannelSeparation::IBCName() const
{
  return "Smooth-ramp channel separation";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\>.
 *  \param[out] a_U     State on the level to be initialized in this
 *                      routine
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_time  Initial time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::initialize(
  LevelData<FArrayBox>&      a_U,
  LevelGridMetrics&          a_gridMetrics,
  const LayoutData<FluxBox>& a_unitNormals,
  const Real                 a_time,
  const int                  a_level) const
{
  CH_TIME("CNSIBCChannelSeparation::initialize");
#if (CH_SPACEDIM != 1)
  const int numWVar   = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx  = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx   = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx  = CRDparam::g_CRDPhysics->temperatureIndex();

  // Constant freestream state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf   = state.velocity()[0]; // freestream velocity
  const Real rho_inf = state.density();     // freestream density
  const Real T_inf   = state.temperature(); // freestream temperature
  const Real p_inf   = state.pressure();    // freestream pressure
  const Real mu      = CRDparam::g_mu;
  const Real gamma   = CRDparam::g_gamma;
  const Real Rval    = CRDparam::g_R;

  const Real Pr_t = 0.89; // Turbulent Prandtl number (Urbin & Knight 2001)
  const Real M_inf_sq = u_inf*u_inf/(gamma*Rval*T_inf);
  const Real d_0 = 0.5*(gamma - 1.)*M_inf_sq*Pr_t;
  // Using Crocco-Busemann boundary layer approximation (compressible)
  const Real T_w = T_inf*(1. + d_0);
  const Real T_w_star = T_w/T_inf;

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box4Dom = grow(box, 4);
      box4Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box4Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box4Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box4Dom,XiFab,XFab,blockCoordSys);

      // Get physical coordinates for y-face centers
      Box faceBox4Dom = box4Dom;
      faceBox4Dom.surroundingNodes(1); // Make this a face-box in y-direction
      FABSTACKTEMP(XiNodeFab, faceBox4Dom, SpaceDim); // Cartesian coordinates
      FABSTACKTEMP(XNodeFab, faceBox4Dom, SpaceDim); // Physical coordinates
      this->CNSIBC::getFaceCoordinates(
        faceBox4Dom, XiNodeFab, XNodeFab, 1, blockCoordSys);

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box2Dom));

      // Pointwise values of W
      FABSTACKTEMP(Wc, box4Dom, numWVar);
      Wc.setVal(0.);

      // Constant factors for use in approximating the IC TBL thickness
      // NOTE: this is very approximate and is probably incorrect, but good
      //       enough for the IC
      const Real c_0 = m_delta/0.37;
      const Real c_1 = 5./4.;
      const Real c_2 = rho_inf*u_inf/mu;

      // Approximate how far along a plate the boundary layer would be in order
      // to have the specified inlet boundary layer thickness
      const Real x_star = std::pow(c_0, c_1)*std::pow(c_2, 0.25);

      // Create a box representing an interior cell of the inlet block
      const Box inletCellBox(IntVect_zero, IntVect_unit);
      // Get the low end of the inlet block
      const RealVect lowXiVal = RealVect_zero;
      const RealVect lowXVal =
        a_gridMetrics.getCoordSys(inletCellBox)->realCoord(lowXiVal);
      const RealVect dx = a_gridMetrics.dxVect();
      RealVect faceOffsetVect = 0.5*dx;
      faceOffsetVect[1] = 0.;

      FABSTACKTEMP(deltaFab, box4Dom, 1); // Bndry-layer-thickness fab
      FABSTACKTEMP(modXFab, box4Dom, SpaceDim);
      // Approximate the boundary-layer thickness for the initial condition
      MD_BOXLOOP(box4Dom, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)],
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          // Get the y-height at the plate at this location
          IntVect currIV = MD_GETIV(i);
          currIV[1] = 0; // set to the lowest cell
          RealVect wallXiLoc = dx*currIV + faceOffsetVect;
          const RealVect wallXLoc =
            a_gridMetrics.getCoordSys(box)->realCoord(wallXiLoc);
          const Real xOffset = x_star - lowXVal[0];
          Real x = loc[0] + xOffset; // x-location wrt to hypothetical plate tip
          Real Re_x = rho_inf*u_inf*x/mu; // location-based Reynolds number
          deltaFab[MD_IX(i, 0)] = 0.37*x/std::pow(Re_x, 0.2);
          D_TERM(modXFab[MD_IX(i, 0)] = x;,
                 modXFab[MD_IX(i, 1)] = loc[1] - wallXLoc[1];,
                 modXFab[MD_IX(i, 2)] = loc[2];);
          XFab[MD_IX(i, 1)] -= wallXLoc[1];
          XNodeFab[MD_IX(i, 1)] -= wallXLoc[1];
        }
      FABSTACKTEMP(meanFab, box4Dom, numWVar);
      FABSTACKTEMP(etaFab, box4Dom, 1);
      muskerIC(meanFab, etaFab, deltaFab, modXFab, XNodeFab, box4Dom, state);

      const Real hiTol = 1.e100;
      const Real loTol = -1.e100;
      // Now we need to add perturbations to streamwise and spanwise velocity
      RealVect physicalLength = m_numCellsInletBlock*m_inletBlockDxVect;
      MD_BOXLOOP(box4Dom, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)],
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          RealVect numPerturbs = m_perturbFreq;
          RealVect arg = PI*numPerturbs*loc/(physicalLength);
          // Add perturbations to x-vel and y-vel together
          Real c_7 = std::max(0., loc[1] - m_delta);
          Real c_6 = 2.*PI/(2.*m_delta);
          Real e4 = std::exp(-c_7*c_6);
          RealVect sinVal(
            D_DECL(std::sin(arg[0]),std::sin(arg[1]-PI/2.),std::sin(arg[2])));
          RealVect cosVal(
            D_DECL(std::cos(arg[0]),std::cos(arg[1]-PI/2.),std::cos(arg[2])));

          // Add perturbations to x-vel and z-vel together (decays wall-normal)
          D_TERM(
            Real u_perturb3 = e4*m_velPerturb*D_TERM(
              sinVal[0],*sinVal[1],*cosVal[2]);
            if (SpaceDim < 3)
              {
                u_perturb3 = 0.;
              },,
            Real w_perturb3 = -e4*m_velPerturb*D_TERM(
              cosVal[0],*sinVal[1],*sinVal[2]););

          RealVect velMean(D_DECL(meanFab[MD_IX(i,WvelIndx)],
                                  meanFab[MD_IX(i,WvelIndx+1)],
                                  meanFab[MD_IX(i,WvelIndx+2)]));
          RealVect velFluc(D_DECL(u_perturb3, 0., w_perturb3));
          RealVect velTotal = velMean;
          if (SpaceDim == 3)
            {
              velTotal += velFluc;
            }
          Real velTotalMag = velTotal.vectorLength();

          // Update temperature based on perturbations
          const Real T_mean =
            T_inf*(T_w_star - d_0*velTotalMag*velTotalMag*(1./(u_inf*u_inf)));
          meanFab[MD_IX(i,tempIndx)] = T_mean;
          // Update density based on temperature and pressure
          const Real rho_mean = p_inf/(T_mean*Rval);
          meanFab[MD_IX(i,rhoIndx)] = rho_mean;

          // Update velocity
          D_TERM(meanFab[MD_IX(i,WvelIndx)]   = velTotal[0];
                 CH_assert((velTotal[0] < hiTol) && (velTotal[0] > loTol));,
                 meanFab[MD_IX(i,WvelIndx+1)] = velTotal[1];
                 CH_assert((velTotal[1] < hiTol) && (velTotal[1] > loTol));,
                 meanFab[MD_IX(i,WvelIndx+2)] = velTotal[2];
                 CH_assert((velTotal[2] < hiTol) && (velTotal[2] > loTol)););

          CH_assert((rho_mean < hiTol) && (rho_mean > loTol));
          CH_assert((T_mean < hiTol) && (T_mean > loTol));
        }

      // Fill Wc with data
      MD_BOXLOOP(box4Dom, i)
        {
          // Assign the entire state
          D_TERM(Wc[MD_IX(i,WvelIndx)]   = meanFab[MD_IX(i,WvelIndx)];,
                 Wc[MD_IX(i,WvelIndx+1)] = meanFab[MD_IX(i,WvelIndx+1)];,
                 Wc[MD_IX(i,WvelIndx+2)] = meanFab[MD_IX(i,WvelIndx+2)];);
          Wc[MD_IX(i, rhoIndx)]  = meanFab[MD_IX(i,rhoIndx)];
          Wc[MD_IX(i, tempIndx)] = meanFab[MD_IX(i,tempIndx)];
        }

      // Initialize the values in UFab
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Initialize the inlet region data structures
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg
 *                      Cell-averaged conservative variables
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Current level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_GlobalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::initializeInletDataStructures(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics,
  const bool                  a_hasFinerGrid) const
{
  CH_TIME("CNSIBCChannelSeparation::initializeInletDataStructures");
  CRD::msg << CRD::fv3
           << "CNSIBCChannelSeparation::initializeInletDataStructures"
           << " (level: " << a_level << ")" << CRD::end;
  if (a_hasFinerGrid || (a_level < (m_numLevels - 1)))
    {
      if (!m_levelDefined[a_level+1])
        {
          return;
        }
    }

  //**NOTE: cp = Copy Plane, sp = Sample Plane, ip = Inlet Plane
  // Create a disjointBoxLayout that covers the plane to copy from with one box

  //**NOTE: We assume that the origin of the domain is IntVect_zero in
  //        computational space
  // 1) Get the current physical domain box and number of cells on this level
  const IntVect domainSize =
    a_disjointBoxLayout.physDomain().domainBox().size();
  const int refRatio = domainSize[0]/m_baseMeshSize[0];
  // 2) Get the number of cells between inlet and sampling plane
  const int spCellLoc = refRatio*m_numCellsInletBlock[0]*m_samplePlaneLoc;
  // 3) Get the domain box for this level covering just the inlet block
  Box domainBox = a_disjointBoxLayout.physDomain().domainBox();
  domainBox.setSmall(IntVect_zero);
  IntVect highSide = refRatio*m_numCellsInletBlock - IntVect_unit;
  domainBox.setBig(highSide); // domainBox is now the box covering the current
                              // level inlet-block
  Box inletBlockBox = domainBox;
  // 4) Get the sampling plane box covering the full possible domain space
  Box fullSamplingBox = domainBox;
  fullSamplingBox.setSmall(0, spCellLoc);
  fullSamplingBox.setBig(0, spCellLoc);
  // 5) Get the real sampling plane box at this point
  int yHiTemp = domainBox.smallEnd()[1];
  domainBox.setBig(1, yHiTemp);
  Vector<Box> tempBoxes = a_disjointBoxLayout.boxArray();
  for (int boxComp = 0; boxComp != tempBoxes.size(); ++boxComp)
    {
      Box tempBox = tempBoxes[boxComp];
      // Get the boxes at the sampling plane location
      if (tempBox.intersects(fullSamplingBox))
        {
          domainBox.minBox(tempBox);
        }
    }
  // 6) Set up one box to cover the copy plane on this level
  Box cpSingleBox = domainBox;
  cpSingleBox.setSmall(0, domainBox.smallEnd()[0] + spCellLoc);
  cpSingleBox.setBig(0, domainBox.smallEnd()[0] + spCellLoc);
  Vector<Box> singlePlaneBoxes;
  singlePlaneBoxes.push_back(cpSingleBox);
  Vector<int> copyPlaneProcList;
  copyPlaneProcList.push_back(0);
  DisjointBoxLayout copyPlaneDBL(singlePlaneBoxes, copyPlaneProcList);

  // 7) Create a disjointBoxLayout that covers the inlet plane with one box
  // Box domainInlet = inletBlockBox;
  // domainInlet.setBig(0, inletBlockBox.smallEnd()[0]);
  // Box compositeBox = domainInlet;
  // compositeBox.setBig(1, inletBlockBox.smallEnd()[1]);
  // for (int boxComp = 0; boxComp != tempBoxes.size(); ++boxComp)
  //   {
  //     Box tempBox = tempBoxes[boxComp];
  //     // Get the boxes at the sampling plane location
  //     if (tempBox.intersects(domainInlet))
  //       {
  //         compositeBox.minBox(tempBox);
  //       }
  //   }
  // Box ipSingleBox = compositeBox;
  // ipSingleBox.setBig(0, inletBlockBox.smallEnd()[0]);
  // Vector<Box> inletPlaneBoxes;
  // inletPlaneBoxes.push_back(ipSingleBox);
  // Vector<int> inletPlaneProcList;
  // inletPlaneProcList.push_back(0);
  // DisjointBoxLayout ipSingleBoxDBL(inletPlaneBoxes, inletPlaneProcList);
  // DEBUG ONLY
  Box ipSingleBox = domainBox;
  ipSingleBox.setBig(0, domainBox.smallEnd()[0]);
  Vector<Box> inletPlaneBoxes;
  inletPlaneBoxes.push_back(ipSingleBox);
  Vector<int> inletPlaneProcList;
  inletPlaneProcList.push_back(0);
  DisjointBoxLayout ipSingleBoxDBL(inletPlaneBoxes, inletPlaneProcList);
  // END DEBUG ONLY

  // 8) Create a DBL that covers the inlet plane with multiple boxes
  Vector<Box> ipMultiBoxes;
  Vector<int> ipMultiBoxProcList;
  Vector<Box> globalBoxes = a_disjointBoxLayout.boxArray();
  Vector<int> globalProcList = a_disjointBoxLayout.procIDs();
  for (int boxComp = 0; boxComp != globalBoxes.size(); ++boxComp)
    {
      Box disjointBox = globalBoxes[boxComp];
      int boxProcID = globalProcList[boxComp];
      if (disjointBox.intersects(ipSingleBox))
        {
          Box tempInletPlaneMultiBox = ipSingleBox;
          tempInletPlaneMultiBox &= disjointBox;
          ipMultiBoxes.push_back(tempInletPlaneMultiBox);
          ipMultiBoxProcList.push_back(boxProcID);
        }
    }

  // We need a sense of problemDomain and periodicity here
  ProblemDomain ipProbDom(ipSingleBox);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir >= 2)
        {
          ipProbDom.setPeriodic(dir, true);
        }
    }
  DisjointBoxLayout ipMultiBoxDBL(ipMultiBoxes, ipMultiBoxProcList, ipProbDom);

  const int numConsComp = CRDparam::g_CRDPhysics->numConservative();
  const int numPrimComp = CRDparam::g_CRDPhysics->numPrimitive();
  const int numGhost = 6;
  IntVect ghostVect = numGhost*IntVect_unit;
  ghostVect[0] = 0;
  m_singleBoxSampledPlane[a_level]->define(
    copyPlaneDBL, numConsComp, ghostVect);
  m_singleBoxInletPlane[a_level]->define(
    ipSingleBoxDBL, numPrimComp, ghostVect);
  m_multiBoxInletPlane[a_level]->define(
    ipMultiBoxDBL, numPrimComp, ghostVect);

  // Define copiers
  m_singleBoxSampledPlaneCopier[a_level]->define(a_disjointBoxLayout,
                                                 copyPlaneDBL);
  m_multiBoxInletPlaneCopier[a_level]->ghostDefine(
    m_singleBoxInletPlane[a_level]->disjointBoxLayout(),
    m_multiBoxInletPlane[a_level]->disjointBoxLayout(),
    ipProbDom, IntVect_zero, ghostVect);
  const IntVect exchangeGhostVect = ghostVect;
  m_multiBoxInletExchangeCopier[a_level]->exchangeDefine(
    m_multiBoxInletPlane[a_level]->disjointBoxLayout(),
    exchangeGhostVect);

  for (int i = 0; i != domainBox.size()[1]; ++i)
    {
      m_uVelSAvg[a_level].push_back(-1.e100);
      m_uVelSTAvg[a_level].push_back(-1.e100);
      m_vVelSAvg[a_level].push_back(-1.e100);
      m_vVelSTAvg[a_level].push_back(-1.e100);
      m_wVelSAvg[a_level].push_back(-1.e100);
      m_rhoSAvg[a_level].push_back(-1.e100);
      m_rhoSTAvg[a_level].push_back(-1.e100);
      m_tempSAvg[a_level].push_back(-1.e100);
      m_tempSTAvg[a_level].push_back(-1.e100);
      m_yLoc[a_level].push_back(0.);
      m_yLocMapped[a_level].push_back(0.);
    }
  D_TERM(
    ,,
    for (int i = 0; i != domainBox.size()[2]; ++i)
      {
        m_zLoc[a_level].push_back(0.);
      });

  Vector<Box> tempBoxes2 = a_disjointBoxLayout.boxArray();
  Box djBox = tempBoxes2[0];
  for (int boxComp = 0; boxComp != tempBoxes2.size(); ++boxComp)
    {
      if (ipSingleBox.contains(tempBoxes2[boxComp]))
        {
          djBox = tempBoxes2[boxComp];
        }
    }
  const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(djBox));
  // Get physical coordinates
  IntVect loCorner = ipSingleBox.smallEnd();
  IntVect hiCorner(D_DECL(loCorner[0],ipSingleBox.bigEnd()[1],loCorner[2]));
  Box yLocBox(loCorner, hiCorner);
  FABSTACKTEMP(XiFab, yLocBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, yLocBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(yLocBox, XiFab, XFab, blockCoordSys);
  std::vector<Real>& yLoc = m_yLoc[a_level];
  MD_BOXLOOP(yLocBox, i)
    {
      IntVect localIndex = MD_GETIV(i);
      yLoc[localIndex[1] - domainBox.smallEnd()[1]]
        = XFab[MD_IX(i, 1)] - 0.22; // Adjust for the plate height
    }
// Get physical z-coordinates
#if CH_SPACEDIM == 3
  IntVect loCornerZ = ipSingleBox.smallEnd();
  IntVect hiCornerZ(
    D_DECL(loCorner[0],loCorner[1],ipSingleBox.bigEnd()[2]));
  Box zLocBox(loCornerZ, hiCornerZ);
  FABSTACKTEMP(XiFabTwo, zLocBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFabTwo, zLocBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(zLocBox,XiFabTwo,XFabTwo,blockCoordSys);
  std::vector<Real>& zLoc = m_zLoc[a_level];
  MD_BOXLOOP(zLocBox, i)
    {
      IntVect localIndex = MD_GETIV(i);
      D_TERM(,,zLoc[localIndex[2] - domainBox.smallEnd()[2]] =
             XFabTwo[MD_IX(i, 2)];);
    }
#endif

  // Set the time-counter to negative one -- we'll increment at the start
  // of every stage
  m_timeCounter[a_level] = -1;
  m_stageCounter[a_level] = -1;
  m_levelDefined[a_level] = 1;
}

/*--------------------------------------------------------------------*/
//  Copy the interior plane to the inlet plane and rescale it
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg
 *                      Cell-averaged conservative variables
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Current level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_GlobalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::copyInteriorToInletAndRescale(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const Real                  a_t,
  const Real                  a_dt,
  const int                   a_stage,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics) const
{
  CH_TIME("CNSIBCChannelSeparation::copyInteriorToInletAndRescale");
  CRD::msg << CRD::fv3
           << "CNSIBCChannelSeparation::copyInteriorToInletAndRescale"
           << " (level: " << a_level << ")" << CRD::end;
  // Increment the stage-counter
  m_stageCounter[a_level] += 1;
  if (a_stage) {return;}

  // Increment the time-counter
  m_timeCounter[a_level] += 1;

  // Set the old time-step size
  if (m_stageCounter[a_level] == 0)
    {
      m_currDt[a_level] = a_dt;
    }
  m_prevDt[a_level] = m_currDt[a_level];
  // Set the current time-step size
  m_currDt[a_level] = a_dt;

  //**NOTE: Steps for computing the scaled inlet condition are numbered.
  //        General comments are un-numbered.
  const int numWComp = CRDparam::g_CRDPhysics->numPrimitive();

  // fill m_singleBoxSampledPlane with data from the interior
  DisjointBoxLayout spDBL = m_singleBoxSampledPlane[a_level]->getBoxes();
  for (DataIterator ditSP = spDBL.dataIterator(); ditSP.ok(); ++ditSP)
    {
      DataIterator ditIP =
        m_singleBoxInletPlane[a_level]->getBoxes().dataIterator();
      FArrayBox& sampledPlaneU = (*(m_singleBoxSampledPlane[a_level]))[ditSP];
      FArrayBox& inletPlaneW = (*(m_singleBoxInletPlane[a_level]))[ditIP];
      sampledPlaneU.setVal(0.);
      inletPlaneW.setVal(0.);
    }
  a_U.copyTo((*(m_singleBoxSampledPlane[a_level])),
             (*(m_singleBoxSampledPlaneCopier[a_level])));
  for (DataIterator ditSP = spDBL.dataIterator(); ditSP.ok(); ++ditSP)
    {
      DataIterator ditIP =
        m_singleBoxInletPlane[a_level]->getBoxes().dataIterator();
      FArrayBox& sampledPlaneU = (*(m_singleBoxSampledPlane[a_level]))[ditSP];
      FArrayBox& inletPlaneW = (*(m_singleBoxInletPlane[a_level]))[ditIP];
      Box spBox = spDBL[ditSP];
      Box ipBox = m_singleBoxInletPlane[a_level]->getBoxes()[ditIP];

      FABSTACKTEMP(WcellPntFab, spBox, numWComp);
      PatchCNSOp::computeWpntCell(WcellPntFab,
                                  sampledPlaneU,
                                  sampledPlaneU, // dummy input that's not used
                                  spBox,
                                  a_disjointBoxLayout.physDomain(),
                                  true,
                                  false);

     CNSIBCChannelSeparation::spaceTimeAverageVariables(
       WcellPntFab, spBox, a_t, a_dt, a_level);

     // Initialize y-mapping
     {
       // const int numWComp = CRDparam::g_CRDPhysics->numPrimitive();

       // Constant freestream state
       const CRDState& state = CRDState::get(m_idxStateInit);
       const Real u_inf   = state.velocity()[0]; // freestream velocity
       const Real rho_inf = state.density();     // freestream density
       const Real mu      = CRDparam::g_mu;

       // Get some box from the current block to get the blockCoordSys
       Box box = a_U.disjointBoxLayout()[a_U.dataIterator()];
       // Get coordinate system and domain for the block
       const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
       Box faceBox = ipBox;
       faceBox.surroundingNodes(1); // Make this a face-box in y-direction
       // Get physical coordinates for y-face centers
       FABSTACKTEMP(XiNodeFab, faceBox, SpaceDim); // Cartesian coordinates
       FABSTACKTEMP(XNodeFab, faceBox, SpaceDim); // Physical coordinates
       this->CNSIBC::getFaceCoordinates(
         faceBox, XiNodeFab, XNodeFab, 1, blockCoordSys);

       // Approximate the virtual x-location of the start of the inlet
       // Essentially, how far are we from the "real" start of the bndry-layer
       const Real c_0 = m_delta/0.37;
       const Real c_1 = 5./4.;
       const Real c_2 = rho_inf*u_inf/mu;
       const Real x_star = std::pow(c_0, c_1)*std::pow(c_2, 0.25);
       FABSTACKTEMP(inletCoord, spBox, SpaceDim); // Inlet coordinate fab
       FABSTACKTEMP(deltaFab, spBox, 1);
       std::vector<Real>& linearYLoc = m_yLoc[a_level];
       std::vector<Real>& mappedYLoc = m_yLocMapped[a_level];
       // Get height of wall-adjacent cell for setting up inlet mean profile
       MD_BOXLOOP(spBox, i)
         {
           const int yIndx = MD_GETIV(i)[1];
           const Real yLoc = linearYLoc[yIndx];
           D_TERM(inletCoord[MD_IX(i, 0)] = x_star;,
                  inletCoord[MD_IX(i, 1)] = yLoc;,
                  inletCoord[MD_IX(i, 2)] = 0.;); // Just set z = 0 everywhere
           deltaFab[MD_IX(i, 0)] = m_delta;

           // Set up the mapping between the inlet and the recycle plane
           // DEBUG ONLY
           // Test matching profiles based on matching mean velocity profiles
           int zIndx = 0;
           D_TERM(,,zIndx = MD_GETIV(i)[2];);
           if (zIndx == 0)
             {
               const Real R = CRDparam::g_R;
               const Real gamma = CRDparam::g_gamma;
               const Real Pr_t = 0.89; // Turbulent Prandtl number
               const Real T_inf = state.temperature();
               const Real M_inf_sq = u_inf*u_inf/(gamma*R*T_inf);
               const Real a_0 = 0.5*(gamma - 1.)*M_inf_sq*Pr_t;
               const Real alpha = std::sqrt(a_0/(1 + a_0));

               // Using Crocco-Busemann boundary layer approximation
               // Approximate density at wall using temperature at wall
               const Real rho_w = rho_inf/(1. + a_0);
               const Real nu_w  = mu/rho_w;

               // Estimate the friction velocity
               const Real Re_theta = 2.E4; // Momentum thickness Re
               //**FIXME: should be user defined
               const Real Z = Re_theta/1000.; // For Coles' wake param
               // Coles' wake parameter for
               // estimating friction velocity
               // (Wenzel 2018 p.453)
               const Real pi1 =
                 0.66*(1. - std::exp(-0.4*std::sqrt(Z) - 0.48*Z));
               const Real u_vd_inf = (u_inf/alpha)*std::asin(alpha);

               // Solve for delta_100 at the inlet
               const Real maxDelta = 1.5*m_delta;
               const Real minDelta = m_delta;
               int iterBrent2 = 0;
               int errorBrent2 = 0;
               const deltaFunc& g =
                 deltaFunc(m_delta,nu_w,u_vd_inf,pi1,m_yPlusMaxGuess);
               Real delta_100_inlet = RootSolver::BrentER(
                 iterBrent2, errorBrent2, g, minDelta, maxDelta);
               if (errorBrent2 != 0 || delta_100_inlet != delta_100_inlet)
                 {
                   CRD::msg << "Bad delta value: " << delta_100_inlet
                            << CRD::error;
                 }

               // Solve for u_tau at the inlet
               Real uTauMin = 0.;
               Real uTauMax = m_yPlusMaxGuess*nu_w/yLoc;
               int iterBrent  = 0;
               int errorBrent = 0;
               const UTauMuskerFunc& f =
                 UTauMuskerFunc(delta_100_inlet,nu_w,u_vd_inf,pi1);
               Real u_tau_inlet = RootSolver::BrentER(
                 iterBrent, errorBrent, f, uTauMin, uTauMax);
               if (errorBrent != 0 || u_tau_inlet != u_tau_inlet)
                 {
                   CRD::msg << "Bad uTau value: " << u_tau_inlet
                            << CRD::error;
                 }

               // Estimate delta_99 at the recycling plane
               RealVect physicalLength =
                 m_numCellsInletBlock*m_inletBlockDxVect;
               const Real deltaX =
                 physicalLength[0]*m_samplePlaneLoc;
               const Real b_0 = m_delta/0.37;
               const Real b_1 = 5./4.;
               const Real b_2 = rho_inf*u_inf/mu;
               const Real x_star = std::pow(b_0, b_1)*std::pow(b_2, 0.25);
               const Real x = x_star + deltaX;
               const Real Re_x = rho_inf*u_inf*x/mu;
               const Real delta_guess = 0.37*x/std::pow(Re_x, 0.2);
               const Real delta_99_recy = m_deltaRecyPercent*delta_guess;

               // Solve for delta_100 at the recycling plane
               const Real maxDelta_recy = 1.5*delta_99_recy;
               const Real minDelta_recy = delta_99_recy;
               int iterBrent3 = 0;
               int errorBrent3 = 0;
               const deltaFunc& g2 =
                 deltaFunc(delta_99_recy,nu_w,u_vd_inf,pi1,m_yPlusMaxGuess);
               Real delta_100_recy = RootSolver::BrentER(
                 iterBrent3, errorBrent3, g2, minDelta_recy, maxDelta_recy);
               if (errorBrent3 != 0 || delta_100_recy != delta_100_recy)
                 {
                   CRD::msg << "Bad delta value: " << delta_100_recy
                            << CRD::error;
                 }

               // Solve for u_tau at the recycling plane
               int iterBrent4  = 0;
               int errorBrent4 = 0;
               const UTauMuskerFunc& f2 =
                 UTauMuskerFunc(delta_100_recy,nu_w,u_vd_inf,pi1);
               Real u_tau_recy = RootSolver::BrentER(
                 iterBrent4, errorBrent4, f2, uTauMin, uTauMax);
               if (errorBrent4 != 0 || u_tau_recy != u_tau_recy)
                 {
                   CRD::msg << "Bad uTau value: " << u_tau_recy << CRD::error;
                 }

               // Solve for y_recy using current y_inlet
               if (yLoc <= m_delta)
                 {
                   int iterBrent5  = 0;
                   int errorBrent5 = 0;
                   const Real maxYVal = delta_100_recy;
                   const Real minYVal = yLoc;

                   const yMatchingFunc& g3 = yMatchingFunc(
                     delta_100_inlet, delta_100_recy, yLoc, u_tau_inlet,
                     u_tau_recy, nu_w, pi1);
                   Real yRecy = RootSolver::BrentER(
                     iterBrent5, errorBrent5, g3, minYVal, maxYVal);
                   if (errorBrent5 != 0 || yRecy != yRecy)
                     {
                       CRD::msg << "Bad y-value: " << yRecy
                                << CRD::error;
                     }
                   mappedYLoc[yIndx] = yRecy;
                   // pout() << "yLoc <= d: yIndx = " << yIndx
                   //        << "; yLoc = " << yLoc << "; yRecy = " << yRecy
                   //        << "; yRatio = " << (yLoc/yRecy) << std::endl;
                 }
               else
                 {
                   // Set this to the yLoc*(delta_recy/m_delta) value
                   // This assumes linear scaling above boundary layer edge
                   mappedYLoc[yIndx] = yLoc*(delta_99_recy/m_delta);
                   // pout() << "yLoc > d: yIndx = " << yIndx
                   //        << "; yLoc = " << yLoc << "; yRecy = "
                   //        << mappedYLoc[yIndx] << "; yRatio = "
                   //        << (yLoc/mappedYLoc[yIndx]) << std::endl;
                 }
             }
           // END DEBUG ONLY
         }
     }

     CNSIBCChannelSeparation::recycleFluctuations(
       inletPlaneW, WcellPntFab, spBox, ipBox, a_level, a_t);

     // // DEBUG ONLY
     // char fileName[64];
     // const int currTimeStep = m_timeCounter[a_level];
     // sprintf(fileName, "inletFluctuations.level-%d.%06d.%dd.hdf5",
     //         a_level, currTimeStep, SpaceDim);
     // Vector<string> names(numWComp);
     // char compNameString[64];
     // for (int compName = 0; compName != numWComp; ++compName)
     //   {
     //     sprintf(compNameString, "inst-%s-target",
     //             CRDparam::g_CRDPhysics->primStateName(compName));
     //     names[compName] = compNameString;
     //   }
     // writeFABname(&inletPlaneW, fileName, names);
    }

  // Set m_multiBoxInletPlane to zero so that we don't introduce invalid values
  // at the recycling plane in the invalid ghost cells
  DisjointBoxLayout mbIPDBL = m_multiBoxInletPlane[a_level]->getBoxes();
  for (DataIterator ditIP = mbIPDBL.dataIterator(); ditIP.ok(); ++ditIP)
    {
      FArrayBox& multiBoxInletFab = (*(m_multiBoxInletPlane[a_level]))[ditIP];
      multiBoxInletFab.setVal(0.);
    }
  // First, copy m_singleBoxInletPlane to m_multiBoxInletPlane
  const Interval copyIntv(0, numWComp-1);
  m_singleBoxInletPlane[a_level]->copyTo(
    copyIntv, (*(m_multiBoxInletPlane[a_level])),
    copyIntv, (*(m_multiBoxInletPlaneCopier[a_level])));
  // Then fill periodic ghost cells of m_multiBoxInletPlane
  m_multiBoxInletPlane[a_level]->exchange(
    (*(m_multiBoxInletExchangeCopier[a_level])));
  // // Finally, fill multiBoxInletPlane invalid ghosts from singleBoxInletPlane
  // m_singleBoxInletPlane[a_level]->copyTo(
  //   copyIntv, (*(m_multiBoxInletPlane[a_level])),
  //   copyIntv, (*(m_multiBoxInletPlaneCopier[a_level])));

  // Additionally, spatially average some data here and print it to the pout
  // files. Compute the separation and reattachment locations as well
#if (CH_SPACEDIM == 3)
  if (!m_useWallModel || m_onlyInlet) { return; }
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int turbCompBegin = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int etaIndx = turbCompBegin + 1 + 1;
  const Real mu = CRDparam::g_mu;
  // 1) Create a near wall box over which to average
  const Box domBox = a_U.getBoxes().physDomain().domainBox();
  Box avgDomBox = a_U.getBoxes().physDomain().domainBox();
  D_TERM(,avgDomBox.setBig(1, avgDomBox.smallEnd()[1]);,
         avgDomBox.setBig(2, avgDomBox.smallEnd()[2]););
  FABSTACKTEMP(domainTauFab, avgDomBox, 2);
  domainTauFab.setVal(0.);
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_U.disjointBoxLayout()[dit];

      // Get the block coordinate system
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(disjointBox));

      // Take the first layer of cells on the low z-side of disjointBox
      Box avgBox = disjointBox;
      avgBox.setBig(1, avgBox.smallEnd(1));
      avgBox.setBig(2, avgBox.smallEnd(2));
      // Assume the entire dataset needs to be spatially averaged
      const FArrayBox& UFab = a_U[dit()];
      // Get the physical space coordinates on the wall face
      FABSTACKTEMP(faceXFab, avgBox, SpaceDim);
      FABSTACKTEMP(faceXiFab, avgBox, SpaceDim);
      CRDparam::g_CNSIBC->getFaceCoordinates(
        avgBox, faceXiFab, faceXFab, 1, blockCoordSys);
      // Average the data in the z-direction
      if (avgBox.smallEnd(1) == avgDomBox.smallEnd(1))
        {
          MD_BOXLOOP(avgBox, i)
            {
              // Create the z-direction box to average over
              Box zAvgBox(MD_GETIV(i),
                          MD_GETIV(i) + BASISV(2)*(disjointBox.size(2) - 1));
              // Sum up the values
              Real avgVal = 0.;
              Real avgX = 0.;
              MD_BOXLOOP(zAvgBox, j)
                {
                  Real rho = UFab[MD_IX(j, rhoIndx)];
                  Real rhoU = UFab[MD_IX(j, WvelIndx)];
                  Real signRhoU = 1.;
                  if (rhoU < 0.)
                    {
                      signRhoU = -1.;
                    }
                  Real nu = mu/rho;
                  Real eta = std::max(0., UFab[MD_IX(j, etaIndx)]);
                  avgVal += signRhoU*eta*nu;
                  avgX += faceXFab[MD_IX(i, 0)];
                }
              // Assign the average to domainFab
              IntVect domFabIV = MD_GETIV(i);
              domFabIV[2] = avgDomBox.smallEnd(2);
              domainTauFab[MD_IV(domFabIV, 0)] += avgVal;
              domainTauFab[MD_IV(domFabIV, 1)] += avgX;
            }
        }
    }
  // Now for the fun part -- we loop over the domainFab and MPI sum everything
  // Basically, we want to do a linear-in, linear-out with a buffer
  // This should reduce the number of MPI_Allreduce calls and the runtime
  int bufferSize = (avgDomBox.numPts())*2;
  std::vector<Real> local(bufferSize, 0.);
  std::vector<Real> global(bufferSize, 0.);
  // Fill the buffer
  int linearIndex = 0;
  for (int comp = 0; comp != 2; ++comp)
    {
      MD_BOXLOOP(avgDomBox, i)
        {
          local[linearIndex] = domainTauFab[MD_IX(i, comp)];
          global[linearIndex] = local[linearIndex];
          ++linearIndex;
        }
    }
  // Sum the buffer
#ifdef CH_MPI
  MPI_Allreduce(local.data(),
                global.data(),
                bufferSize,
                MPI_CH_REAL,
                MPI_SUM,
                Chombo_MPI::comm);
#endif
  // Read out the buffer
  linearIndex = 0;
  for (int comp = 0; comp != 2; ++comp)
    {
      MD_BOXLOOP(avgDomBox, i)
        {
          domainTauFab[MD_IX(i, comp)] = (global[linearIndex])/domBox.size(2);
          ++linearIndex;
        }
    }

  // Compute the first and last sign change
  // Real currSign = 1.;
  // Real prevSign = 1.;
  std::vector<Real> sepLoc;
  Box sepBox = avgDomBox;
  sepBox.setSmall(0, (avgDomBox.smallEnd()[0])+1);
  const int MD_ID(o, 0); // Offset in x-direction
  MD_BOXLOOP(avgDomBox, i)
    {
      const int xIndx = MD_GETIV(i)[0];
      if (xIndx > 0)
        {
          Real tauR = domainTauFab[MD_IX(i, 0)];
          Real tauL = domainTauFab[MD_OFFSETIX(i,-,o,0)];
          Real tauSign = tauR*tauL;
          if (tauSign < 0)
            {
              sepLoc.push_back(domainTauFab[MD_IX(i, 1)]);
            }
        }
    }
  if ((procID() == 0) && !(sepLoc.empty()))
    {
      CRD::msg << "Tau: Level: " << a_level << "; Time: " << a_t
               << "; sepLoc1 = " << sepLoc.front() << "; sepLoc2 = "
               << sepLoc.back() << CRD::var;
    }
  else if (procID() == 0)
    {
      Real frontLoc = -1.5;
      CRD::msg << "Tau: Level: " << a_level << "; Time: " << a_t
               << "; sepLoc1 = " << frontLoc << "; sepLoc2 = "
               << frontLoc << CRD::var;
    }

  // Now print the spanwise averages of velocity at 9 unique points
  if (a_level == 1) // DEBUG ONLY -- 1 for refRatio = 4; 2 for refRatio = 2
    {
      // Create a vector of boxes representing the points of interest
      std::vector<Box> boxList;
      IntVect ivHi = domBox.bigEnd();
      IntVect ivLo = domBox.smallEnd();

      if (!m_newGeometry)
        {
          // Box 1
          ivHi[0] = 167;
          ivLo[0] = 167;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box1(ivLo, ivHi);
          boxList.push_back(box1);

          // Box 2
          ivHi[0] = 354;
          ivLo[0] = 354;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box2(ivLo, ivHi);
          boxList.push_back(box2);

          // Box 3
          ivHi[0] = 384;
          ivLo[0] = 384;
          ivHi[1] = 8;
          ivLo[1] = 8;
          Box box3(ivLo, ivHi);
          boxList.push_back(box3);

          // Box 4
          ivHi[0] = 430;
          ivLo[0] = 430;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box4(ivLo, ivHi);
          boxList.push_back(box4);

          // Box 5
          ivHi[0] = 480;
          ivLo[0] = 480;
          ivHi[1] = 24;
          ivLo[1] = 24;
          Box box5(ivLo, ivHi);
          boxList.push_back(box5);

          // Box 6
          ivHi[0] = 529;
          ivLo[0] = 529;
          ivHi[1] = 44;
          ivLo[1] = 44;
          Box box6(ivLo, ivHi);
          boxList.push_back(box6);

          // Box 7
          ivHi[0] = 576;
          ivLo[0] = 576;
          ivHi[1] = 44;
          ivLo[1] = 44;
          Box box7(ivLo, ivHi);
          boxList.push_back(box7);

          // Box 8
          ivHi[0] = 729;
          ivLo[0] = 729;
          ivHi[1] = 32;
          ivLo[1] = 32;
          Box box8(ivLo, ivHi);
          boxList.push_back(box8);

          // Box 9
          ivHi[0] = 816;
          ivLo[0] = 816;
          ivHi[1] = 19;
          ivLo[1] = 19;
          Box box9(ivLo, ivHi);
          boxList.push_back(box9);
        }
      else
        {
          // Box 1 -- physical location (-0.62, 0.238)
          ivHi[0] = 169;
          ivLo[0] = 169;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box1(ivLo, ivHi);
          boxList.push_back(box1);

          // Box 2 -- physical location (-0.15, 0.238)
          ivHi[0] = 355;
          ivLo[0] = 355;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box2(ivLo, ivHi);
          boxList.push_back(box2);

          // Box 3 -- physical location (0, 0.237)
          ivHi[0] = 384;
          ivLo[0] = 384;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box3(ivLo, ivHi);
          boxList.push_back(box3);

          // Box 4 -- physical location (0.25, 0.215)
          ivHi[0] = 432;
          ivLo[0] = 432;
          ivHi[1] = 9;
          ivLo[1] = 9;
          Box box4(ivLo, ivHi);
          boxList.push_back(box4);

          // Box 5 -- physical location (0.5, 0.16)
          ivHi[0] = 480;
          ivLo[0] = 480;
          ivHi[1] = 25;
          ivLo[1] = 25;
          Box box5(ivLo, ivHi);
          boxList.push_back(box5);

          // Box 6 -- physical location (0.75, 0.13)
          ivHi[0] = 528;
          ivLo[0] = 528;
          ivHi[1] = 52;
          ivLo[1] = 52;
          Box box6(ivLo, ivHi);
          boxList.push_back(box6);

          // Box 7 -- physical location (1, 0.108)
          ivHi[0] = 576;
          ivLo[0] = 576;
          ivHi[1] = 52;
          ivLo[1] = 52;
          Box box7(ivLo, ivHi);
          boxList.push_back(box7);

          // Box 8 -- physical location (1.3, 0.07)
          ivHi[0] = 729;
          ivLo[0] = 729;
          ivHi[1] = 35;
          ivLo[1] = 35;
          Box box8(ivLo, ivHi);
          boxList.push_back(box8);

          // Box 9 -- physical location (1.75, 0.039)
          ivHi[0] = 816;
          ivLo[0] = 816;
          ivHi[1] = 20;
          ivLo[1] = 20;
          Box box9(ivLo, ivHi);
          boxList.push_back(box9);
        }

      for (int boxComp = 0; boxComp != boxList.size(); ++boxComp)
        {
          RealVect meanVel = RealVect_zero;
          Box avgBox = boxList[boxComp];
          for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
            {
              const Box disjointBox = a_U.disjointBoxLayout()[dit];
              const FArrayBox& UFab = a_U[dit()];

              // Intersect this box with the box in the list
              if (disjointBox.intersects(avgBox))
                {
                  Box localBox = avgBox;
                  localBox &= disjointBox;
                  MD_BOXLOOP(localBox, i)
                    {
                      Real rho = UFab[MD_IX(i, rhoIndx)];
                      D_TERM(meanVel[0] += UFab[MD_IX(i, WvelIndx)]/rho;,
                             meanVel[1] += UFab[MD_IX(i, WvelIndx+1)]/rho;,
                             meanVel[2] += UFab[MD_IX(i, WvelIndx+2)]/rho;);
                    }
                }
            }
          int bufferSize = 3;
          std::vector<Real> local(bufferSize, 0.);
          std::vector<Real> global(bufferSize, 0.);
          // Fill the buffer
          int linearIndex = 0;
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              local[linearIndex] = meanVel[comp];
              global[linearIndex] = local[linearIndex];
              ++linearIndex;
            }
          // Sum the buffer
#ifdef CH_MPI
          MPI_Allreduce(local.data(),
                        global.data(),
                        bufferSize,
                        MPI_CH_REAL,
                        MPI_SUM,
                        Chombo_MPI::comm);
#endif
          // Read out the buffer
          linearIndex = 0;
          RealVect avgVel = RealVect_zero;
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              avgVel[comp] = (global[linearIndex])/domBox.size(2);
              ++linearIndex;
            }
          if (procID() == 0)
            {
              CRD::msg << "Point " << boxComp << "; u = " << avgVel[0]
                       << "; v = " << avgVel[1] << "; w = " << avgVel[2]
                       << CRD::var;
            }
        }
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCChannelSeparation::haveExactSol() const
{
  return false;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the imposed (exterior or farfield) primitive state at flow BCC
/** State is set according to reference conditions
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *  \param[in]  a_disjointBox
 *                      Complete current interior disjoint box
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_bcInfo) const
{
  CH_TIME("CNSIBCChannelSeparation::setImposedBCprimState");
  // Set all the component variables and intervals
  const int numWComp = CRDparam::g_CRDPhysics->numPrimitive();
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    {
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    case CRDparam::DomainBCTypeFarfield:
    {
      break;
    }
    case CRDparam::DomainBCTypeOutflow:
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      // const Real alpha = a_bcInfo.m_relaxCBCStateParam; // 1/0 = non/reflect
      // Create a mean boundary layer profile for the outlet state
      const Real outlet_delta = m_deltaOut;
      FABSTACKTEMP(etaInletFab, a_boundaryFaceBox, 1);
      FABSTACKTEMP(meanInletFab, a_boundaryFaceBox, numWComp);
      // Set meanInletFab to a_Wface so all unchanged components remain
      for (int comp = 0; comp != numWComp; ++comp)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              meanInletFab[MD_IX(i, comp)] = state(comp);
            }
        }
      FABSTACKTEMP(deltaFab, a_boundaryFaceBox, 1);
      deltaFab.setVal(outlet_delta); // Set to the predicted BL thickness
      // Get coordinate system and domain for the block
      FABSTACKTEMP(XFab, a_boundaryFaceBox, SpaceDim);
      FABSTACKTEMP(XiFab, a_boundaryFaceBox, SpaceDim);
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(a_disjointBox));
      this->CNSIBC::getFaceCoordinates(
        a_boundaryFaceBox, XiFab, XFab, 0, blockCoordSys);
      // Get y-face coordinates
      Box faceBox = a_boundaryFaceBox;
      faceBox.growLo(0, 1); // Grow this to cover one interior cell
      faceBox.enclosedCells(0); // Make this a cell box
      faceBox.shift(0, 1); // Shift this to align x-index with a_boundaryFaceBox
      faceBox.surroundingNodes(1); // Make this a face box in y-direction
      FABSTACKTEMP(XNodeFab, faceBox, SpaceDim);
      FABSTACKTEMP(XiNodeFab, faceBox, SpaceDim);
      this->CNSIBC::getFaceCoordinates(
        faceBox, XiNodeFab, XNodeFab, 1, blockCoordSys);
      if (m_onlyInlet)
        {
          // Shift the wall-normal direction down by 0.22
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              XFab[MD_IX(i, 1)] -= 0.22;
            }
          MD_BOXLOOP(faceBox, i)
            {
              XNodeFab[MD_IX(i, 1)] -= 0.22;
            }
        }
      muskerIC(meanInletFab, etaInletFab, deltaFab, XFab, XNodeFab,
               a_boundaryFaceBox, state);
      for (int c = 0; c != numWComp; ++c)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              // const Real interior = a_Wcell[MD_IX(i, c)];
              const Real exterior = meanInletFab[MD_IX(i, c)];
              // a_Wface[MD_IX(i, c)] = alpha*interior + (1. - alpha)*exterior;
              a_Wface[MD_IX(i, c)] = exterior;
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeInflow:
    case CRDparam::DomainBCTypeRelaxedCBCIn:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      const Real u_inf = state.velocity()[0];
      // We need to use the previously set boundary data
      // 1) loop over the entire DBL for m_multiBoxInletPlane
      // 2) if the box is contained in a_disjointBox, use this dataIterator
      // 3) get the correct FArrayBox using the correct dataIterator
      DisjointBoxLayout mbipDBL = m_multiBoxInletPlane[a_level]->getBoxes();
      DataIterator ditIPCurrent = mbipDBL.dataIterator();
      for (DataIterator ditIP = mbipDBL.dataIterator(); ditIP.ok(); ++ditIP)
        {
          const Box box = mbipDBL[ditIP];
          if (a_disjointBox.contains(box))
            {
              ditIPCurrent = ditIP;
            }
        }
      FArrayBox& scaledFab =
        (*(m_multiBoxInletPlane[a_level]))[ditIPCurrent];
      scaledFab.shiftHalf(0, -1);

      FABSTACKTEMP(etaInletFab, a_boundaryFaceBox, 1);
      FABSTACKTEMP(meanInletFab, a_boundaryFaceBox, numWComp);
      // Set meanInletFab to a_Wface so all unchanged components remain
      for (int comp = 0; comp != numWComp; ++comp)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              meanInletFab[MD_IX(i, comp)] = a_Wface[MD_IX(i, comp)];
            }
        }
      FABSTACKTEMP(deltaFab, a_boundaryFaceBox, 1);
      deltaFab.setVal(m_delta); // Set to the incoming BL thickness
      // Get coordinate system and domain for the block
      FABSTACKTEMP(XFab, a_boundaryFaceBox, SpaceDim);
      FABSTACKTEMP(XiFab, a_boundaryFaceBox, SpaceDim);
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(a_disjointBox));
      this->CNSIBC::getFaceCoordinates(
        a_boundaryFaceBox, XiFab, XFab, 0, blockCoordSys);
      // Get y-face coordinates
      Box faceBox = a_boundaryFaceBox;
      faceBox.growHi(0, 1); // Grow this to cover one interior cell
      faceBox.enclosedCells(0); // Make this a cell box
      faceBox.surroundingNodes(1); // Make this a face box in y-direction
      FABSTACKTEMP(XNodeFab, faceBox, SpaceDim);
      FABSTACKTEMP(XiNodeFab, faceBox, SpaceDim);
      this->CNSIBC::getFaceCoordinates(
        faceBox, XiNodeFab, XNodeFab, 1, blockCoordSys);
      // Shift the wall-normal direction down by 0.22
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          XFab[MD_IX(i, 1)] -= 0.22;
        }
      MD_BOXLOOP(faceBox, i)
        {
          XNodeFab[MD_IX(i, 1)] -= 0.22;
        }
      muskerIC(meanInletFab, etaInletFab, deltaFab, XFab, XNodeFab,
               a_boundaryFaceBox, state);
      // Add a small perturbation to prevent laminarization (0%)
      const Real perturbMag = 0.0*u_inf;
      const RealVect physicalLength = m_numCellsInletBlock*m_inletBlockDxVect;
      const Real timeShift = u_inf*a_time;
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)] - timeShift,
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          RealVect numPerturbs = m_perturbFreq;
          RealVect arg = PI*numPerturbs*loc/(physicalLength);
          // Add perturbations to x-vel and y-vel together
          Real c_7 = std::max(0., loc[1] - m_delta);
          Real c_6 = 2.*PI/(2.*m_delta);
          Real e4 = std::exp(-c_7*c_6);
          RealVect sinVal(
            D_DECL(std::sin(arg[0]),std::sin(arg[1]-PI/2.),std::sin(arg[2])));
          RealVect cosVal(
            D_DECL(std::cos(arg[0]),std::cos(arg[1]-PI/2.),std::cos(arg[2])));

          // Add perturbations to x-vel and z-vel together (decays wall-normal)
          D_TERM(
            Real u_perturb3 = e4*perturbMag*D_TERM(
              sinVal[0],*sinVal[1],*cosVal[2]);
            if (SpaceDim < 3)
              {
                u_perturb3 = 0.;
              },,
            Real w_perturb3 = -e4*perturbMag*D_TERM(
              cosVal[0],*sinVal[1],*sinVal[2]););
          if (loc[1] <= m_delta)
            {
              D_TERM(meanInletFab[MD_IX(i, cVel)] += u_perturb3;,
                     meanInletFab[MD_IX(i, cVel+1)] += 0.;,
                     meanInletFab[MD_IX(i, cVel+2)] += w_perturb3;);
            }
        }
      // Fill a_Wface with the changed components from meanInletFab
      for (int comp = 0; comp != numWComp; ++comp)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              a_Wface[MD_IX(i, comp)] = meanInletFab[MD_IX(i, comp)]
                + scaledFab[MD_IX(i, comp)];
            }
        }
      scaledFab.shiftHalf(0, 1);
      break;
    }
    default:
      CH_assert(false);
    }
}

/*--------------------------------------------------------------------*/
//  Output pnt-values and domain-sums specific to simulation
/** \param[in] a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in] a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in] a_WcellPntFab
 *                      Cell-centered primitive state
 *  \param[in] a_unitNormalsFxb
 *                      Unit-normal basis for faces
 *  \param[in] a_gridMetrics
 *                      Grid metrics for current level
 *  \param[in] a_time   Current time on current level
 *  \param[in] a_level  Current level
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::inSituSumPntProcessing(
  const LevelData<FluxBox>&   a_faceAvgPlotData,
  const LevelData<FluxBox>&   a_WfaceAvgFxb,
  const LevelData<FArrayBox>& a_WcellAvgFab,
  const LevelData<FArrayBox>& a_WcellPntFab,
  const LayoutData<FluxBox>&  a_unitNormalsFxb,
  const LevelGridMetrics&     a_gridMetrics,
  const Real                  a_time,
  const int                   a_stage,
  const int                   a_level) const
{
  CH_TIME("CNSIBCChannelSeparation::inSituSumPntProcessing");

  int plotVar = CRDparam::g_plotLoFaceAvgComps;
  plotVar |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
  int plotPrimAndTurb = ((plotVar & CRDparam::PlotPrimitive) &&
                         (plotVar & CRDparam::PlotTurbulentComps)) ? 1 : 0;

  if (!m_onlyInlet || !plotPrimAndTurb || a_stage ||
      (a_level != (m_numLevels-1))) { return; }

  const int cVel = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();

  // Indexing for a_faceAvgPlotData (need to abstract indexing away from here)
  int cStart = CRDparam::g_CRDPhysics->numPrimitive();
  if ((plotVar & CRDparam::PlotFluxes) &&
      (CRDparam::g_physicsModels & CRDparam::PhysicsInertial))
    {
      cStart += CRDparam::g_CRDPhysics->numFluxes();
    }
  if ((plotVar & CRDparam::PlotFluxes) &&
      (CRDparam::g_physicsModels & CRDparam::PhysicsViscous))
    {
      cStart += CRDparam::g_CRDPhysics->numFluxes();
    }
  const int cResolvedUU = cStart;
  const int cSGSUU = cStart + 3 + (SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
  const int cTauWStreamwise = cSGSUU + (SpaceDim*(SpaceDim-1.)/2.+SpaceDim) + 1;

  // Print <u_i*u_j>, <u_i>, <tau_w_streamwise>
  //       6        + 3    + 1
  const int numStressComp = (SpaceDim*(SpaceDim-1.)/2.+SpaceDim);
  const int numPrintComp = 10;

  // 1) Create a near wall box over which to average
  const Box domBox = a_WcellAvgFab.getBoxes().physDomain().domainBox();
  Box avgDom1Box = domBox;
  D_TERM(,avgDom1Box.setBig(1, domBox.smallEnd()[1]);,
         avgDom1Box.setBig(2, domBox.smallEnd()[2]););
  FABSTACKTEMP(domainAvg1Fab, avgDom1Box, numPrintComp);
  domainAvg1Fab.setVal(0.);
  // 2) Create a box halfway through the boundary layer
  const int halfDeltaYIndx =
    0.03125*(domBox.bigEnd()[1] - domBox.smallEnd()[1]);
  Box avgDom2Box = domBox;
  D_TERM(,avgDom2Box.setSmall(1, halfDeltaYIndx);
         avgDom2Box.setBig(1, halfDeltaYIndx);,
         avgDom2Box.setBig(2, domBox.smallEnd()[2]););
  FABSTACKTEMP(domainAvg2Fab, avgDom2Box, numPrintComp);
  domainAvg2Fab.setVal(0.);

  for (DataIterator dit = a_WcellAvgFab.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_WcellAvgFab.disjointBoxLayout()[dit];
      // const BlockDomain& blockDomain =
      //   a_gridMetrics.getCoordSys().problemDomain(box);
      const FArrayBox& WfaceAvgYFab = a_WfaceAvgFxb[dit][1];
      const FArrayBox& faceAvgPlotYFab = a_faceAvgPlotData[dit][1];

      // Plot wall-surface point values at 0.2, 0.4, 0.6, 0.8*L_x
      // and plot point values at y = 0.03125*L_y, x = 0.2, 0.4, 0.6, 0.8*L_x
      const int domXLength = domBox.bigEnd()[0] - domBox.smallEnd()[0];
      const int domYLength = domBox.bigEnd()[1] - domBox.smallEnd()[1];
      std::vector<IntVect> intVectList;
      IntVect pnt1 = IntVect_zero;
      pnt1[0] = 0.2*domXLength;
      intVectList.push_back(pnt1);
      IntVect pnt2 = IntVect_zero;
      pnt2[0] = 0.4*domXLength;
      intVectList.push_back(pnt2);
      IntVect pnt3 = IntVect_zero;
      pnt3[0] = 0.6*domXLength;
      intVectList.push_back(pnt3);
      IntVect pnt4 = IntVect_zero;
      pnt4[0] = 0.8*domXLength;
      intVectList.push_back(pnt4);

      pnt1[1] = 0.03125*domYLength;
      intVectList.push_back(pnt1);
      pnt2[1] = 0.03125*domYLength;
      intVectList.push_back(pnt2);
      pnt3[1] = 0.03125*domYLength;
      intVectList.push_back(pnt3);
      pnt4[1] = 0.03125*domYLength;
      intVectList.push_back(pnt4);

      for (int cVect = 0; cVect != intVectList.size(); ++cVect)
        {
          IntVect idxVect = intVectList[cVect];
          if (box.contains(intVectList[cVect]))
            {
              Real uu = faceAvgPlotYFab(idxVect, cResolvedUU)
                + faceAvgPlotYFab(idxVect, cSGSUU);
              Real uv = faceAvgPlotYFab(idxVect, cResolvedUU+1)
                + faceAvgPlotYFab(idxVect, cSGSUU+1);
              Real uw = faceAvgPlotYFab(idxVect, cResolvedUU+2)
                + faceAvgPlotYFab(idxVect, cSGSUU+2);
              Real vv = faceAvgPlotYFab(idxVect, cResolvedUU+3)
                + faceAvgPlotYFab(idxVect, cSGSUU+3);
              Real vw = faceAvgPlotYFab(idxVect, cResolvedUU+4)
                + faceAvgPlotYFab(idxVect, cSGSUU+4);
              Real ww = faceAvgPlotYFab(idxVect, cResolvedUU+5)
                + faceAvgPlotYFab(idxVect, cSGSUU+5);
              Real tau_w = faceAvgPlotYFab(idxVect, cTauWStreamwise);
              pout() << "Pnt: Time: " << a_time << ", Point: " << cVect
                     << ", u: " << WfaceAvgYFab(idxVect, cVel)
                     << ", v: " << WfaceAvgYFab(idxVect, cVel+1)
                     << ", w: " << WfaceAvgYFab(idxVect, cVel+2)
                     << ", uu: " << uu << ", uv: " << uv << ", uw: " << uw
                     << ", vv: " << vv << ", vw: " << vw << ", ww: " << ww
                     << ", tau_w: " << tau_w
                     << std::endl;
              pout() << "" << std::endl;
            }
        }

      // Take the first layer of cells on the low z-side of disjointBox
      Box avgBox1 = box; // Set to the entire domain box
      avgBox1.setSmall(1, domBox.smallEnd(1)); // First set small side
      avgBox1.setBig(1, domBox.smallEnd(1));   // Now set large side
      avgBox1.setBig(2, box.smallEnd(2));      // Project onto back plane
      Box avgBox2 = box;
      avgBox2.setSmall(1, domBox.smallEnd(1)); // Avoid making empty box
      avgBox2.setBig(1, domBox.bigEnd(1));     // Avoid making empty box
      avgBox2.setSmall(1, halfDeltaYIndx);     // Now set to correct index
      avgBox2.setBig(1, halfDeltaYIndx);       // Now set to correct index
      avgBox2.setBig(2, box.smallEnd(2));      // Project onto back plane

      // Plot wall-surface average values at 0.2, 0.4, 0.6, 0.8*L_x
      if (box.contains(avgBox1))
        {
          int cVal = 0;
          // Avg <u_i>
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              MD_BOXLOOP(avgBox1, i)
                {
                  // Create the z-direction box to average over
                  Box zAvgBox(MD_GETIV(i),
                              MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
                  // Sum up the values
                  Real avgVal = 0.;
                  MD_BOXLOOP(zAvgBox, j)
                    {
                      Real uVal = WfaceAvgYFab[MD_IX(j, cVel+comp)];
                      avgVal += uVal;
                    }
                  // Assign the average to domainFab
                  IntVect domFabIV = MD_GETIV(i);
                  domFabIV[2] = domBox.smallEnd(2);
                  domainAvg1Fab[MD_IV(domFabIV, cVal)] += avgVal;
                }
              ++cVal;
            }
          // Avg <u_i*u_j>
          for (int comp = 0; comp != numStressComp; ++comp)
            {
              MD_BOXLOOP(avgBox1, i)
                {
                  // Create the z-direction box to average over
                  Box zAvgBox(MD_GETIV(i),
                              MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
                  // Sum up the values
                  Real avgVal = 0.;
                  MD_BOXLOOP(zAvgBox, j)
                    {
                      Real uuVal = faceAvgPlotYFab[MD_IX(j, cResolvedUU+comp)];
                      Real sgsVal = faceAvgPlotYFab[MD_IX(j, cSGSUU+comp)];
                      avgVal += uuVal + sgsVal;
                    }
                  // Assign the average to domainFab
                  IntVect domFabIV = MD_GETIV(i);
                  domFabIV[2] = domBox.smallEnd(2);
                  domainAvg1Fab[MD_IV(domFabIV, cVal)] += avgVal;
                }
              ++cVal;
            }
          // Avg wall shear stress
          MD_BOXLOOP(avgBox1, i)
            {
              // Create the z-direction box to average over
              Box zAvgBox(MD_GETIV(i),
                          MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
              // Sum up the values
              Real avgVal = 0.;
              MD_BOXLOOP(zAvgBox, j)
                {
                  Real tauVal = faceAvgPlotYFab[MD_IX(j, cTauWStreamwise)];
                  avgVal += tauVal;
                }
              // Assign the average to domainFab
              IntVect domFabIV = MD_GETIV(i);
              domFabIV[2] = domBox.smallEnd(2);
              domainAvg1Fab[MD_IV(domFabIV, cVal)] += avgVal;
            }
        }

      // Plot average values at y = 0.03125*L_y, x = 0.2, 0.4, 0.6, 0.8*L_x
      if (box.contains(avgBox2))
        {
          int cVal = 0;
          // Avg <u_i>
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              MD_BOXLOOP(avgBox2, i)
                {
                  // Create the z-direction box to average over
                  Box zAvgBox(MD_GETIV(i),
                              MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
                  // Sum up the values
                  Real avgVal = 0.;
                  MD_BOXLOOP(zAvgBox, j)
                    {
                      Real uVal = WfaceAvgYFab[MD_IX(j, cVel+comp)];
                      avgVal += uVal;
                    }
                  // Assign the average to domainFab
                  IntVect domFabIV = MD_GETIV(i);
                  domFabIV[2] = domBox.smallEnd(2);
                  domainAvg2Fab[MD_IV(domFabIV, cVal)] += avgVal;
                }
              ++cVal;
            }
          // Avg <u_i*u_j>
          for (int comp = 0; comp != numStressComp; ++comp)
            {
              MD_BOXLOOP(avgBox2, i)
                {
                  // Create the z-direction box to average over
                  Box zAvgBox(MD_GETIV(i),
                              MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
                  // Sum up the values
                  Real avgVal = 0.;
                  MD_BOXLOOP(zAvgBox, j)
                    {
                      Real uuVal = faceAvgPlotYFab[MD_IX(j, cResolvedUU+comp)];
                      Real sgsVal = faceAvgPlotYFab[MD_IX(j, cSGSUU+comp)];
                      avgVal += uuVal + sgsVal;
                    }
                  // Assign the average to domainFab
                  IntVect domFabIV = MD_GETIV(i);
                  domFabIV[2] = domBox.smallEnd(2);
                  domainAvg2Fab[MD_IV(domFabIV, cVal)] += avgVal;
                }
              ++cVal;
            }
          // Avg wall shear stress
          MD_BOXLOOP(avgBox2, i)
            {
              // Create the z-direction box to average over
              Box zAvgBox(MD_GETIV(i),
                          MD_GETIV(i)+BASISV(2)*(box.size(2)-1));
              // Sum up the values
              Real avgVal = 0.;
              MD_BOXLOOP(zAvgBox, j)
                {
                  Real tauVal = faceAvgPlotYFab[MD_IX(j, cTauWStreamwise)];
                  avgVal += tauVal;
                }
              // Assign the average to domainFab
              IntVect domFabIV = MD_GETIV(i);
              domFabIV[2] = domBox.smallEnd(2);
              domainAvg2Fab[MD_IV(domFabIV, cVal)] += avgVal;
            }
        }
    }

  // Linear-in, linear-out with a buffer
  int buffer1Size = (avgDom1Box.numPts())*numPrintComp;
  std::vector<Real> local1(buffer1Size, 0.);
  std::vector<Real> global1(buffer1Size, 0.);
  // Fill the buffer
  int linearIndex = 0;
  for (int comp = 0; comp != numPrintComp; ++comp)
    {
      MD_BOXLOOP(avgDom1Box, i)
        {
          local1[linearIndex] = domainAvg1Fab[MD_IX(i, comp)];
          global1[linearIndex] = local1[linearIndex];
          ++linearIndex;
        }
    }
  // Sum the buffer
#ifdef CH_MPI
  MPI_Allreduce(local1.data(),
                global1.data(),
                buffer1Size,
                MPI_CH_REAL,
                MPI_SUM,
                Chombo_MPI::comm);
#endif
  // Read out the buffer
  linearIndex = 0;
  for (int comp = 0; comp != numPrintComp; ++comp)
    {
      int zLength = domBox.size(2);
      MD_BOXLOOP(avgDom1Box, i)
        {
          domainAvg1Fab[MD_IX(i, comp)] = (global1[linearIndex])/zLength;
          ++linearIndex;
        }
    }

  // Linear-in, linear-out with a buffer
  int buffer2Size = (avgDom2Box.numPts())*numPrintComp;
  std::vector<Real> local2(buffer2Size, 0.);
  std::vector<Real> global2(buffer2Size, 0.);
  // Fill the buffer
  linearIndex = 0;
  for (int comp = 0; comp != numPrintComp; ++comp)
    {
      MD_BOXLOOP(avgDom2Box, i)
        {
          local2[linearIndex] = domainAvg2Fab[MD_IX(i, comp)];
          global2[linearIndex] = local2[linearIndex];
          ++linearIndex;
        }
    }
  // Sum the buffer
#ifdef CH_MPI
  MPI_Allreduce(local2.data(),
                global2.data(),
                buffer2Size,
                MPI_CH_REAL,
                MPI_SUM,
                Chombo_MPI::comm);
#endif
  // Read out the buffer
  linearIndex = 0;
  for (int comp = 0; comp != numPrintComp; ++comp)
    {
      int zLength = domBox.size(2);
      MD_BOXLOOP(avgDom2Box, i)
        {
          domainAvg2Fab[MD_IX(i, comp)] = (global2[linearIndex])/zLength;
          ++linearIndex;
        }
    }

  for (DataIterator dit = a_WcellAvgFab.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_WcellAvgFab.disjointBoxLayout()[dit];

      // Plot wall-surface point values at 0.2, 0.4, 0.6, 0.8*L_x
      // and plot point values at y = 0.03125*L_y, x = 0.2, 0.4, 0.6, 0.8*L_x
      const int domXLength = domBox.bigEnd()[0] - domBox.smallEnd()[0];
      const int domYLength = domBox.bigEnd()[1] - domBox.smallEnd()[1];
      std::vector<IntVect> intVectList1;
      IntVect pnt1 = IntVect_zero;
      pnt1[0] = 0.2*domXLength;
      intVectList1.push_back(pnt1);
      IntVect pnt2 = IntVect_zero;
      pnt2[0] = 0.4*domXLength;
      intVectList1.push_back(pnt2);
      IntVect pnt3 = IntVect_zero;
      pnt3[0] = 0.6*domXLength;
      intVectList1.push_back(pnt3);
      IntVect pnt4 = IntVect_zero;
      pnt4[0] = 0.8*domXLength;
      intVectList1.push_back(pnt4);

      std::vector<IntVect> intVectList2;
      pnt1[1] = 0.03125*domYLength;
      intVectList2.push_back(pnt1);
      pnt2[1] = 0.03125*domYLength;
      intVectList2.push_back(pnt2);
      pnt3[1] = 0.03125*domYLength;
      intVectList2.push_back(pnt3);
      pnt4[1] = 0.03125*domYLength;
      intVectList2.push_back(pnt4);

      for (int cVect = 0; cVect != intVectList1.size(); ++cVect)
        {
          IntVect idxVect = intVectList1[cVect];
          if (box.contains(intVectList1[cVect]))
            {
              Real uu = domainAvg1Fab(idxVect, SpaceDim);
              Real uv = domainAvg1Fab(idxVect, SpaceDim+1);
              Real uw = domainAvg1Fab(idxVect, SpaceDim+2);
              Real vv = domainAvg1Fab(idxVect, SpaceDim+3);
              Real vw = domainAvg1Fab(idxVect, SpaceDim+4);
              Real ww = domainAvg1Fab(idxVect, SpaceDim+5);
              Real tau_w = domainAvg1Fab(idxVect, SpaceDim+6);
              pout() << "Avg: Time: " << a_time << ", Point: " << cVect
                     << ", u: " << domainAvg1Fab(idxVect, cVel)
                     << ", v: " << domainAvg1Fab(idxVect, cVel+1)
                     << ", w: " << domainAvg1Fab(idxVect, cVel+2)
                     << ", uu: " << uu << ", uv: " << uv << ", uw: " << uw
                     << ", vv: " << vv << ", vw: " << vw << ", ww: " << ww
                     << ", tau_w: " << tau_w << std::endl;
              pout() << "" << std::endl;
            }
        }
      for (int cVect = 0; cVect != intVectList2.size(); ++cVect)
        {
          IntVect idxVect = intVectList2[cVect];
          if (box.contains(intVectList2[cVect]))
            {
              Real uu = domainAvg2Fab(idxVect, SpaceDim);
              Real uv = domainAvg2Fab(idxVect, SpaceDim+1);
              Real uw = domainAvg2Fab(idxVect, SpaceDim+2);
              Real vv = domainAvg2Fab(idxVect, SpaceDim+3);
              Real vw = domainAvg2Fab(idxVect, SpaceDim+4);
              Real ww = domainAvg2Fab(idxVect, SpaceDim+5);
              Real tau_w = domainAvg2Fab(idxVect, SpaceDim+6);
              pout() << "Avg: Time: " << a_time << ", Point: " << (cVect+4)
                     << ", u: " << domainAvg2Fab(idxVect, cVel)
                     << ", v: " << domainAvg2Fab(idxVect, cVel+1)
                     << ", w: " << domainAvg2Fab(idxVect, cVel+2)
                     << ", uu: " << uu << ", uv: " << uv << ", uw: " << uw
                     << ", vv: " << vv << ", vw: " << vw << ", ww: " << ww
                     << ", tau_w: " << tau_w << std::endl;
              pout() << "" << std::endl;
            }
        }
    }
}

/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Read any information related to the IBC from input
/** \note
 *  <ul>
 *    <li> No output should be printed aside from errors
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::readBCInfo()
{
  ParmParse ppIBC("ibc");

  ppIBC.query("inlet_boundary_layer_thickness", m_delta);
  ppIBC.query("outlet_boundary_layer_thickness", m_deltaOut);
  ppIBC.query("recycling_plane_location", m_samplePlaneLoc);
  ppIBC.query("delta_recy_scaling", m_deltaRecyPercent);
  ppIBC.query("y_plus_max_guess", m_yPlusMaxGuess);
  ppIBC.query("vel_perturb_magnitude", m_velPerturb);
  std::vector<Real> perturbFreq(SpaceDim);
  ppIBC.queryarr("perturb_freq", perturbFreq, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_perturbFreq.dataPtr(),
                                           &perturbFreq.front());

  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf = state.velocity()[0]; // freestream velocity
  m_lambda = m_delta/u_inf;

  // We need to know how many levels on which we could possibly require
  // sampling-plane recycling
  ParmParse ppAMR("amr");

  int maxLevel = 0;
  ppAMR.query("max_level", maxLevel);
  m_numLevels = maxLevel + 1;

  // Read number of cells in the inlet block base mesh
  ParmParse ppGRID("grid");
  D_TERM(ppGRID.query("l_length", m_numCellsInletBlock[0]);,
         ppGRID.query("r_length", m_numCellsInletBlock[1]);,
         ppGRID.query("z_length", m_numCellsInletBlock[2]););
  // Read mesh spacing in the inlet block base mesh
  D_TERM(ppGRID.query("l_dx", m_inletBlockDxVect[0]);,
         ppGRID.query("r_dx", m_inletBlockDxVect[1]);,
         ppGRID.query("z_dx", m_inletBlockDxVect[2]););

  // Construct the base problem domain
  const int numGhostCells = 8;
  const int gapSize = 3*numGhostCells;
  D_TERM(m_baseMeshSize[0] = m_numCellsInletBlock[0] + gapSize + 48 + gapSize
         + m_numCellsInletBlock[1];,
         m_baseMeshSize[1] = 32;,
         m_baseMeshSize[2] = m_numCellsInletBlock[2];);

  // If this is the new geometry, we need to change a few things
  ParmParse ppCOORDSYS("coordsys");
  m_newGeometry = 0;
  ppCOORDSYS.query("new_geometry", m_newGeometry);
  if (m_newGeometry)
    {
      std::vector<int> numCells(SpaceDim);
      ppGRID.queryarr("num_cells", numCells, 0, SpaceDim);
      SpaceDimArray<int, int>::loadFromArray(m_numCellsInletBlock.dataPtr(),
                                             &numCells.front());
      D_TERM(ppCOORDSYS.query("num_x_cells_inlet_block",
                              m_numCellsInletBlock[0]);,,);
      int numCellsRampBlock = 1;
      ppCOORDSYS.query("num_x_cells_total_ramp_block", numCellsRampBlock);
      int numCellsOutletBlock = 1;
      ppCOORDSYS.query("num_x_cells_outlet_block", numCellsOutletBlock);
      m_baseMeshSize = m_numCellsInletBlock;
      m_baseMeshSize[0] = m_numCellsInletBlock[0] + gapSize + numCellsRampBlock
        + gapSize + numCellsOutletBlock;
      int numCellsRampBlockPreRamp = 1;
      ppCOORDSYS.query("num_x_cells_ramp_block_before_ramp",
                       numCellsRampBlockPreRamp);
      int numCellsRampBlockPostRamp = 1;
      ppCOORDSYS.query("num_x_cells_ramp_block_after_ramp",
                       numCellsRampBlockPostRamp);
      // The ramp is unit length and must always have a set number of cells
      // in the input file
      const Real physDx = 1./(numCellsRampBlock - numCellsRampBlockPreRamp
                              - numCellsRampBlockPostRamp);
      const Real yHeight = 0.74 - 0.22; // This is fixed for the inlet block
      Real zWidth = 0.128;
      ppCOORDSYS.query("spanwise_domain_width", zWidth);
      D_TERM(,Real physDy = yHeight/m_numCellsInletBlock[1];,
             Real physDz = zWidth/m_numCellsInletBlock[2];);
      D_TERM(m_inletBlockDxVect[0] = physDx;,
             m_inletBlockDxVect[1] = physDy;,
             m_inletBlockDxVect[2] = physDz;);
    }

  ppCOORDSYS.query("only_simulate_channel_inlet", m_onlyInlet);
  if (m_onlyInlet)
    {
      m_baseMeshSize[0] = m_numCellsInletBlock[0];
    }

  ParmParse ppTurb("turb");
  ppTurb.query("use_wall_model", m_useWallModel);

  m_readInput = true;
}

/*--------------------------------------------------------------------*/
//  Space-time average the flow variables
/** \note
 *  <ul>
 *    <li> No output should be printed aside from errors
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::spaceTimeAverageVariables(
  const FArrayBox& a_WcellPntFab,
  const Box        a_samplePlaneBox,
  const Real       a_t,
  const Real       a_dt,
  const int        a_level) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::spaceTimeAverageVariables");
  // NOTE: space-time averaging doesn't work too well for computing the
  //       instantaneous fluctuations to superimpose at the inlet.
  //       Instead, we essentially use spatial-averaging.
  //       However, we need a long-time-averaged velocity for the
  //       boundary-layer thickness calculation at the recycle plane.
  //       So, we store both the short-time-averaged and long-time-averaged
  //       velocity at the recycling plane.

  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const IntVect nCells = a_samplePlaneBox.size();
  int nCellsZ = 1;
  D_TERM(,,nCellsZ = nCells[2];);

  // Compute spanwise-averaged u-vel (streamwise) and v-vel (normal)
  //    a) {u}^n = \Delta t/T\avg{u}^n + (1 - \Delta t/T)\avg{u}^n-1
  //    b) t_1 = m_lambda*m_psi1
  //    c) if t_current < t_1, T = T_1 = m_phi1*m_lambda
  //    d) t_2 = m_lambda*(m_psi1 + m_psi2)
  //    e) if t_current >= t_1 && t_current < t_2, T = T_2 = m_phi2*m_lambda
  //    f) if t_current >= t_2, T = T_3 = T_2 + t_current - t_2

  // Set up time-averaging window sizes
  Real T = m_phi2*m_lambda + a_t - m_lambda*(m_psi1 + m_psi2);
  if (a_t < (m_lambda*m_psi1))
    {
      T = m_phi1*m_lambda;
    }
  else if (a_t < (m_lambda*(m_psi1 + m_psi2)))
    {
      T = m_phi2*m_lambda;
    }

  // Compute the spanwise-averaged primitive state -- once per time-step
  for (int yLoc = 0; yLoc != nCells[1]; ++yLoc)
    {
      // set up the one-dimensional averaging box at the current y-index
      IntVect loIntVect = a_samplePlaneBox.smallEnd();
      IntVect hiIntVect = a_samplePlaneBox.bigEnd();
      loIntVect[1] = yLoc;
      hiIntVect[1] = yLoc;
      Box averageBox(loIntVect, hiIntVect);

      // average the variables over the 1D box
      Real rhoAvg = 0.;
      D_TERM(Real uVelAvg = 0.;, Real vVelAvg = 0.;, Real wVelAvg = 0.;);
      Real TAvg = 0.;
      MD_BOXLOOP(averageBox, i)
        {
          rhoAvg += a_WcellPntFab[MD_IX(i, rhoIndx)];
          D_TERM(
            uVelAvg += a_WcellPntFab[MD_IX(i, WvelIndx)];,
            vVelAvg += a_WcellPntFab[MD_IX(i, WvelIndx+1)];,
            wVelAvg += a_WcellPntFab[MD_IX(i, WvelIndx+2)];);
          TAvg += a_WcellPntFab[MD_IX(i, tempIndx)];
        }

      // divide by the number of cells to finish off the spatial averaging
      rhoAvg /= nCellsZ;
      D_TERM(uVelAvg /= nCellsZ;, vVelAvg /= nCellsZ;, wVelAvg /= nCellsZ;);
      TAvg /= nCellsZ;

      // Save the instantaneous spatial averages
      D_TERM(m_uVelSAvg[a_level][yLoc] = uVelAvg;,
             m_vVelSAvg[a_level][yLoc] = vVelAvg;,
             m_wVelSAvg[a_level][yLoc] = wVelAvg;);
      m_rhoSAvg[a_level][yLoc] = rhoAvg;
      m_tempSAvg[a_level][yLoc] = TAvg;

      // Temporally average the spatial averages
      const Real rhoSTAvg = m_rhoSTAvg[a_level][yLoc];
      const Real uVelSTAvg = m_uVelSTAvg[a_level][yLoc];
      const Real vVelSTAvg = m_vVelSTAvg[a_level][yLoc];
      const Real tempSTAvg = m_tempSTAvg[a_level][yLoc];
      Real c_0 = a_dt/T;
      // we place these checks here for restarts so that only the current
      // value of the average is used for this time step -- while we may
      // be able to guess a close-enough value for the previous space-time
      // average, it still won't be as good as just using the current
      // spatially-averaged instantaneous value
      if ((uVelSTAvg <= -1.e100) ||
          (vVelSTAvg <= -1.e100) ||
          (rhoSTAvg  <= -1.e100) ||
          (tempSTAvg <= -1.e100))
        {
          c_0 = 1.;
        }
      D_TERM(m_uVelSTAvg[a_level][yLoc] = c_0*uVelAvg + (1. - c_0)*uVelSTAvg;,
             m_vVelSTAvg[a_level][yLoc] = c_0*vVelAvg + (1. - c_0)*vVelSTAvg;,);
      m_rhoSTAvg[a_level][yLoc] = c_0*rhoAvg + (1. - c_0)*rhoSTAvg;
      m_tempSTAvg[a_level][yLoc] = c_0*TAvg + (1. - c_0)*tempSTAvg;
    }
}

/*--------------------------------------------------------------------*/
//  Compute the inlet recycled-state consisting of fluctuations
/** \param[out] a_inletPlaneW
 *                      Recycled inlet state
 *  \param[in]  a_spBox Box covering the sampling plane
 *  \param[in]  a_stage Box covering the inlet plane
 *  \param[in]  a_t     Current time
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::recycleFluctuations(
  FArrayBox&       a_inletPlaneW,
  const FArrayBox& a_WcellPntFab,
  const Box&       a_spBox,
  const Box&       a_ipBox,
  const int        a_level,
  const Real       a_t) const
{
  CH_TIME("CNSIBCChannelSeparation::recycleFluctuations");
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int nCellsY = m_yLoc[a_level].size();
  D_TERM(,,const int nCellsZ = m_zLoc[a_level].size(););

  // Constant freestream state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf = state.velocity()[0]; // freestream velocity
  const Real rho_inf = state.density(); // freestream density

  // 1) Specify delta at the recycle plane
  Real physLengthY = m_inletBlockDxVect[1]*m_numCellsInletBlock[1];
  Real delta_recy = physLengthY;
  if (!a_level)
    {
      for (int yIndx = 0; yIndx != nCellsY; ++yIndx)
        {
          Real yLoc = m_yLoc[a_level][nCellsY - yIndx - 1];
          //**NOTE: we iterate from the top so we can catch the last point that
          //        has an avgUVel >= 0.99 the target freestream velocity
          Real avgUVel = m_uVelSTAvg[a_level][nCellsY - yIndx - 1];
          if (avgUVel >= 0.99*u_inf)
            {
              delta_recy = yLoc;
            }
        }
      // We know that, initiallly, delta_recy will be wrong due to transients
      // So, we force it to be some value until some pre-specified time
      Real physLengthX = m_inletBlockDxVect[0]*m_numCellsInletBlock[0];
      Real deltaX = physLengthX*m_samplePlaneLoc;
      Real c_0 = deltaX/m_delta;
      Real c_1 = std::pow(0.27, 6./5.);
      Real Re_delta = rho_inf*u_inf*m_delta/(CRDparam::g_mu);
      Real c_2 = std::pow(Re_delta, -1./5.);
      Real c_3 = 1. + c_0*c_1*c_2;
      Real c_4 = std::pow(c_3, 5./6.);
      Real deltaGuess = c_4*m_delta;
      Real beta = 200*m_lambda; // 200 BL-thickness times
      Real zeta = 2000*m_lambda; // 2000 BL-thickness times
      Real alpha = 0.5 + 0.5*std::tanh((1./beta)*(a_t - zeta));
      delta_recy = (1. - alpha)*deltaGuess + alpha*delta_recy;
      delta_recy = m_deltaRecyPercent*deltaGuess;
      for (int i = 0; i != m_numLevels; ++i)
        {
          m_deltaRecy[i] = delta_recy;
        }
    }
  delta_recy = m_deltaRecy[0];

  // We're adding a temporally varying z-direction shifting to the fluctuations
  D_TERM(,,
         Real physLengthZ = m_inletBlockDxVect[2]*m_numCellsInletBlock[2];
         Real z_offset = 0.5*physLengthZ;
         Real alpha_t_1 = 0.8*m_delta;
         Real omega_1 = 1.;
         Real tL1 = 80.*m_lambda;
         Real alpha_t_2 = 0.4*m_delta;
         Real omega_2 = 2.;
         Real tL2 = 40.*m_lambda;
         Real alpha_t_3 = 0.2*m_delta;
         Real omega_3 = 4.;
         Real tL3 = 20.*m_lambda;
         Real z_shift = z_offset + alpha_t_1*std::sin(2.*PI*omega_1*a_t/tL1)
         + alpha_t_2*std::sin(2.*PI*omega_2*a_t/tL2)
         + alpha_t_3*std::sin(2.*PI*omega_3*a_t/tL3);
         // Always perform shifting
         // if (a_t > (240.*m_lambda)) // ~6.6 flow-throughs
         //   {
         //     // Bring spanwise shifting to a halt at a periodic point so that
         //     // the shifting location is at least continuous in time
         //     z_shift = 0.;
         //   }
    );

  // Compute the fluctuation rescaling
  FABSTACKTEMP(rescaledFluc, a_spBox, CRDparam::g_CRDPhysics->numPrimitive());
  rescaledFluc.setVal(0.);
  MD_BOXLOOP(a_spBox, i)
    {
      IntVect inletLoc = MD_GETIV(i);
      inletLoc[0] = a_ipBox.smallEnd()[0]; // set this to the inlet plane
      // Get the indices for the current wall-normal location
      const Real yIdx = i1 - a_spBox.smallEnd()[1];
      // const Real yLoc = m_yLoc[a_level][yIdx];
      // fluctuations for inlet
      //const Real y_s = (delta_recy/m_delta)*yLoc;
      const Real y_s = m_yLocMapped[a_level][yIdx];
      // find the y-indices that bound y_s
      Real y_hi = m_yLoc[a_level][nCellsY - 1];
      Real y_lo = m_yLoc[a_level][0];
      int  hiIndx = nCellsY - 1;
      int  loIndx = 0;
      for (int yIndx = 0; yIndx != nCellsY; ++yIndx)
        {
          Real hiYLoc = m_yLoc[a_level][nCellsY - 1 - yIndx];
          Real loYLoc = m_yLoc[a_level][yIndx];
          if (hiYLoc >= y_s)
            {
              y_hi = hiYLoc;
              hiIndx = nCellsY - 1 - yIndx;
            }
          if (loYLoc <= y_s)
            {
              y_lo = loYLoc;
              loIndx = yIndx;
            }
        }
      IntVect hiLoc = MD_GETIV(i);
      IntVect loLoc = MD_GETIV(i);
      hiLoc[1] = hiIndx; // set this to the lower bounding cell on the
                         // recycling plane (next cell below y_s loc)
      loLoc[1] = loIndx; // set this to the higher bounding cell on the
                         // recycling plane (next cell above y_s loc)
      Real rho_mean_hi = m_rhoSAvg[a_level][hiIndx];
      Real rho_mean_lo = m_rhoSAvg[a_level][loIndx];
      D_TERM(
        Real u_mean_hi = m_uVelSAvg[a_level][hiIndx];
        Real u_mean_lo = m_uVelSAvg[a_level][loIndx];,
        Real v_mean_hi = m_vVelSAvg[a_level][hiIndx];
        Real v_mean_lo = m_vVelSAvg[a_level][loIndx];,
        Real w_mean_hi = m_wVelSAvg[a_level][hiIndx];
        Real w_mean_lo = m_wVelSAvg[a_level][loIndx];);
      Real T_mean_hi = m_tempSAvg[a_level][hiIndx];
      Real T_mean_lo = m_tempSAvg[a_level][loIndx];

      Real rho_fluc_h = a_WcellPntFab[MD_IV(hiLoc, rhoIndx)] - rho_mean_hi;
      Real rho_fluc_l = a_WcellPntFab[MD_IV(loLoc, rhoIndx)] - rho_mean_lo;
      D_TERM(
        Real u_fluc_h = a_WcellPntFab[MD_IV(hiLoc, WvelIndx)]-u_mean_hi;
        Real u_fluc_l = a_WcellPntFab[MD_IV(loLoc, WvelIndx)]-u_mean_lo;,
        Real v_fluc_h = a_WcellPntFab[MD_IV(hiLoc, WvelIndx+1)]-v_mean_hi;
        Real v_fluc_l = a_WcellPntFab[MD_IV(loLoc, WvelIndx+1)]-v_mean_lo;,
        Real w_fluc_h = a_WcellPntFab[MD_IV(hiLoc, WvelIndx+2)]-w_mean_hi;
        Real w_fluc_l = a_WcellPntFab[MD_IV(loLoc, WvelIndx+2)]-w_mean_lo;);
      Real T_fluc_h = a_WcellPntFab[MD_IV(hiLoc, tempIndx)] - T_mean_hi;
      Real T_fluc_l = a_WcellPntFab[MD_IV(loLoc, tempIndx)] - T_mean_lo;

      Real rho_fluc = rho_fluc_l;
      D_TERM(Real u_fluc = u_fluc_l;,
             Real v_fluc = v_fluc_l;,
             Real w_fluc = w_fluc_l;);
      Real T_fluc = T_fluc_l;
      if (hiIndx != loIndx)
        {
          rho_fluc = rho_fluc_l
            + (y_s - y_lo)*(rho_fluc_h - rho_fluc_l)/(y_hi - y_lo);
          D_TERM(
            u_fluc = u_fluc_l+(y_s-y_lo)*(u_fluc_h-u_fluc_l)/(y_hi-y_lo);,
            v_fluc = v_fluc_l+(y_s-y_lo)*(v_fluc_h-v_fluc_l)/(y_hi-y_lo);,
            w_fluc = w_fluc_l+(y_s-y_lo)*(w_fluc_h-w_fluc_l)/(y_hi-y_lo););
          T_fluc = T_fluc_l + (y_s-y_lo)*(T_fluc_h-T_fluc_l)/(y_hi-y_lo);
        }
      if (y_s >= 1.25*delta_recy)
        {
          rho_fluc = 0.;
          D_TERM(u_fluc = 0.;,
                 v_fluc = 0.;,
                 w_fluc = 0.;);
          T_fluc = 0.;
        }

      rescaledFluc[MD_IX(i, rhoIndx)] = rho_fluc;
      D_TERM(rescaledFluc[MD_IX(i, WvelIndx)] = u_fluc;,
             rescaledFluc[MD_IX(i, WvelIndx+1)] = v_fluc;,
             rescaledFluc[MD_IX(i, WvelIndx+2)] = w_fluc;);
      rescaledFluc[MD_IX(i, tempIndx)] = T_fluc;
    }

  // Add in the spanwise shifting in time
  MD_BOXLOOP(a_spBox, i)
    {
      IntVect inletLoc = MD_GETIV(i);
      inletLoc[0] = a_ipBox.smallEnd()[0]; // set this to the inlet plane

#if (CH_SPACEDIM == 3)
      int hiIndx = nCellsZ - 1;
      int loIndx = 0;
      Real zIdx = 0.;
      Real zLoc = 0.;
      Real zLocTrue = 0.;
      D_TERM(
        ,,zIdx = i2 - a_spBox.smallEnd()[2];
        zLoc = m_zLoc[a_level][zIdx];
        zLocTrue = zLoc + z_shift;
        if (zLocTrue < 0.)
          {
            zLocTrue += physLengthZ;
          }
        else if (zLocTrue > physLengthZ)
          {
            zLocTrue -= physLengthZ;
          });

      // First, identify the low z index
      Real z_hi = m_zLoc[a_level][nCellsZ - 1];
      Real z_lo = m_zLoc[a_level][0];
      bool exactCheck = false;
      // Check if zLocTrue is between z_hi and z_lo
      if ((zLocTrue < z_lo) || (zLocTrue > z_hi))
        {
          loIndx = nCellsZ - 1;
          hiIndx = 0;
          D_TERM(,,z_hi = m_zLoc[a_level][0] + physLengthZ;
                 z_lo = m_zLoc[a_level][nCellsZ - 1];);
          if (zLocTrue < z_lo) // Make sure zLocTrue is on the high side
            {
              D_TERM(,,zLocTrue += physLengthZ;);
            }
        }
      else // If it is between z_hi and z_lo, iterate through and set the values
        {
          for (int zIndx = 0; zIndx != nCellsZ; ++zIndx)
            {
              Real loZLoc = m_zLoc[a_level][zIndx];
              // Check to see if zLocTrue equals any of the 
              if (loZLoc == zLocTrue)
                {
                  exactCheck = true;
                  z_lo = loZLoc;
                  z_hi = loZLoc;
                  loIndx = zIndx;
                  hiIndx = zIndx;
                }
              else if (loZLoc < zLocTrue) // Set z_lo to last value < zLocTrue
                {
                  z_lo = loZLoc;
                  loIndx = zIndx;
                  hiIndx = zIndx + 1;
                }
            }
          z_hi = m_zLoc[a_level][hiIndx];
          z_lo = m_zLoc[a_level][loIndx];
        }
#endif

      IntVect hiLoc = MD_GETIV(i);
      IntVect loLoc = MD_GETIV(i);
      D_TERM(,,hiLoc[2] = hiIndx;
               loLoc[2] = loIndx;);

      // Left fluctuations
      Real rho_fluc_lo = rescaledFluc[MD_IV(loLoc, rhoIndx)];
      D_TERM(
        Real u_fluc_lo = rescaledFluc[MD_IV(loLoc, WvelIndx)];,
        Real v_fluc_lo = rescaledFluc[MD_IV(loLoc, WvelIndx+1)];,
        Real w_fluc_lo = rescaledFluc[MD_IV(loLoc, WvelIndx+2)];);
      Real T_fluc_lo = rescaledFluc[MD_IV(loLoc, tempIndx)];

      Real rho_fluc = rho_fluc_lo;
      D_TERM(
        Real u_fluc = u_fluc_lo;,
        Real v_fluc = v_fluc_lo;,
        Real w_fluc = w_fluc_lo;);
      Real T_fluc = T_fluc_lo;
#if (CH_SPACEDIM == 3)

      // Right fluctuations
      Real rho_fluc_hi = rescaledFluc[MD_IV(hiLoc, rhoIndx)];
      D_TERM(
        Real u_fluc_hi = rescaledFluc[MD_IV(hiLoc, WvelIndx)];,
        Real v_fluc_hi = rescaledFluc[MD_IV(hiLoc, WvelIndx+1)];,
        Real w_fluc_hi = rescaledFluc[MD_IV(hiLoc, WvelIndx+2)];);
      Real T_fluc_hi = rescaledFluc[MD_IV(hiLoc, tempIndx)];

      if (!exactCheck)
        {
          Real DeltaZ = z_hi - z_lo;
          Real diffZ = zLocTrue - z_lo;
          rho_fluc = rho_fluc_lo + diffZ*(rho_fluc_hi - rho_fluc_lo)/DeltaZ;
          D_TERM(
            u_fluc = (u_fluc_lo + diffZ*(u_fluc_hi - u_fluc_lo)/DeltaZ);,
            v_fluc = (v_fluc_lo + diffZ*(v_fluc_hi - v_fluc_lo)/DeltaZ);,
            w_fluc = (w_fluc_lo + diffZ*(w_fluc_hi - w_fluc_lo)/DeltaZ););
          T_fluc = T_fluc_lo + diffZ*(T_fluc_hi - T_fluc_lo)/DeltaZ;
        }
#endif

      D_TERM(
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx)]   = u_fluc;,
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx+1)] = v_fluc;,
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx+2)] = w_fluc;);
      a_inletPlaneW[MD_IV(inletLoc, tempIndx)] = T_fluc;

      // Version 2
      a_inletPlaneW[MD_IV(inletLoc, rhoIndx)] = rho_fluc;
      a_inletPlaneW[MD_IV(inletLoc, presIndx)] = 0.;
    }
}

/*--------------------------------------------------------------------*/
//  Compute an approximate initial condition based on Musker 1979
/** \param[out] a_meanFab
 *                      Approximate mean initial condition
 *  \param[in]  a_etaFab
 *                      Approximate eta_0 initial condition
 *  \param[in]  a_XFab  Physical space coordinates
 *  \param[in]  a_box   Box over which to initialize
 *//*-----------------------------------------------------------------*/

void
CNSIBCChannelSeparation::muskerIC(FArrayBox&       a_meanFab,
                                  FArrayBox&       a_etaFab,
                                  const FArrayBox& a_deltaFab,
                                  const FArrayBox& a_XFab,
                                  const FArrayBox& a_XNodeFab,
                                  const Box&       a_box,
                                  const CRDState&  a_state) const
{
  CH_TIME("CNSIBCChannelSeparation::muskerIC");

  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const Real R = CRDparam::g_R;
  const Real gamma = CRDparam::g_gamma;
  const Real mu = CRDparam::g_mu;
  const Real Pr_t = 0.89; // Turbulent Prandtl number (Urbin & Knight 2001)

  // Constant freestream state
  const Real u_inf = a_state.velocity()[0]; // freestream velocity
  const Real T_inf = a_state.temperature(); // freestream temperature
  const Real p_inf = a_state.pressure();    // freestream pressure
  const Real rho_inf = a_state.density();   // freestream density
  CH_assert(rho_inf > 0.);
  CH_assert(p_inf > 0.);
  CH_assert(T_inf > 0.);

  const Real M_inf_sq = u_inf*u_inf/(gamma*R*T_inf);
  const Real c_0 = 0.5*(gamma - 1.)*M_inf_sq*Pr_t;
  const Real alpha = std::sqrt(c_0/(1 + c_0));

  // Using Crocco-Busemann boundary layer approximation (compressible)
  const Real T_w = T_inf*(1. + c_0);
  const Real T_w_star = T_w/T_inf;
  // Approximate density at the wall using temperature at the wall
  const Real rho_w = rho_inf*(T_inf/T_w);
  const Real nu_w  = mu/rho_w;

  // Estimate the friction velocity
  const Real Re_theta = 2.E4; // Momentum thickness Reynolds number
                              //**FIXME: should be user defined (or estimated)
  const Real Z = Re_theta/1000.; // Constant for estimating Coles' wake param
  // Coles' wake parameter for estimating friction velocity (Wenzel 2018 p.453)
  const Real pi1 = 0.66*(1. - std::exp(-0.4*std::sqrt(Z) - 0.48*Z));

  const Real hiTol = 1.e100;
  const Real loTol = -1.e100;
  FABSTACKTEMP(uTau, a_box, 1);
  const int MD_ID(o, 1);
  MD_BOXLOOP(a_box, i)
    {
      // Approximate the boundary-layer thickness for the initial condition
      RealVect loc(D_DECL(a_XFab[MD_IX(i, 0)],
                          a_XFab[MD_IX(i, 1)],
                          a_XFab[MD_IX(i, 2)]));
      Real delta = a_deltaFab[MD_IX(i, 0)];
      // Friction velocity (based on Musker 1979)
      Real uTauMin = 0.;
      Real uTauMax = m_yPlusMaxGuess*nu_w/loc[1];
      int iterBrent  = 0;
      int errorBrent = 0;

      const Real k = 0.41;
      const Real u_vd_inf = (u_inf/alpha)*std::asin(alpha);

      const Real maxDelta = 1.5*delta;
      const Real minDelta = delta;
      int iterBrent2 = 0;
      int errorBrent2 = 0;
      const deltaFunc& g = deltaFunc(delta,nu_w,u_vd_inf,pi1,m_yPlusMaxGuess);
      Real delta_100 =
        RootSolver::BrentER(iterBrent2, errorBrent2, g, minDelta, maxDelta);
      if (errorBrent2 != 0 || delta_100 != delta_100)
        {
          CRD::msg << "Bad delta value: " << delta_100 << CRD::error;
        }

      const UTauMuskerFunc& f = UTauMuskerFunc(delta_100,nu_w,u_vd_inf,pi1);
      Real u_tau =
        RootSolver::BrentER(iterBrent, errorBrent, f, uTauMin, uTauMax);
      if (errorBrent != 0 || u_tau != u_tau)
        {
          CRD::msg << "Bad uTau value: " << u_tau << CRD::error;
        }
      uTau[MD_IX(i, 0)] = u_tau;

      // Now, we can use a subcell resolution to average the variables
      // within a cell and reduce the profile error
      const int numSubCells = 20;
      const Real loYNode = a_XNodeFab[MD_IX(i, 1)];
      const Real hiYNode = a_XNodeFab[MD_OFFSETIX(i,+,o,1)];
      const Real subCellDy = (hiYNode - loYNode)/numSubCells;
      // Initialize values to zero
      Real T_mean_sum = 0.;
      Real rho_mean_sum = 0.;
      Real u_vel_sum = 0.;
      Real y_vel_sum = 0.;
      for (int cell = 0; cell != numSubCells; ++cell)
        {
          // Modify the yLoc value to match the center of the current subcell
          loc[1] = loYNode + subCellDy*(0.5 + cell);

          // Estimate the streamwise velocity (Musker 1979)
          const Real y_p = loc[1]*u_tau/nu_w;
          const Real b_0 = 5.424*std::atan2((2.*y_p - 8.15), 16.7);
          const Real b_1 = std::pow((y_p + 10.6), 9.6);
          const Real b_2 = (y_p*y_p - 8.15*y_p + 86);
          const Real b_3 = std::log10(b_1/(b_2*b_2));
          const Real b_4 = loc[1]/delta_100;
          const Real b_5 = b_4*b_4;
          const Real b_6 = b_5*b_4;
          const Real b_7 = 2.44*(pi1*(6.*b_5 - 4.*b_6) + b_5*(1. - b_4));
          const Real u_vd = u_tau*(b_0 + b_3 - 3.52 + b_7);
          Real u_vel = std::sin(u_vd*alpha/u_inf)*u_inf/alpha;
          if (loc[1] > delta_100)
            {
              u_vel = u_inf;
            }

          // Estimate the temperature using Crocco-Busemann
          const Real T_mean =
            T_inf*(T_w_star - c_0*u_vel*u_vel*(1./(u_inf*u_inf)));
          // Compute mean density
          const Real rho_mean = p_inf/(T_mean*R);
          // Estimate wall-normal velocity through du/dx and then integrating
          // -dv/dy and applying v=0 at the wall
          Real yLocMod = loc[1];
          if (loc[1] > delta_100) // If > 1 BL thickness, y-vel is constant
            {
              yLocMod = delta_100;
            }
          const Real e_0 = -(yLocMod*4.*u_tau)/(5.*k*k*loc[0]);
          const Real e_1 = (1. - std::log(yLocMod*u_tau/nu_w));
          const Real e_2 = (1./k) + (u_inf*std::asin(alpha)/(alpha*u_tau));
          const Real y_vel = e_0*e_1*(1./e_2);

          u_vel_sum += u_vel/numSubCells;
          y_vel_sum += y_vel/numSubCells;
          rho_mean_sum += rho_mean/numSubCells;
          T_mean_sum += T_mean/numSubCells;
        }
      const Real eta_0 = u_tau*u_tau/nu_w;
      a_etaFab[MD_IX(i, 0)] = eta_0;

      CH_assert(rho_mean_sum > 0.);
      CH_assert(T_mean_sum > 0.);
      CH_assert(u_vel_sum < hiTol && u_vel_sum > loTol);
      a_meanFab[MD_IX(i, rhoIndx)] = rho_mean_sum;
      D_TERM(
        a_meanFab[MD_IX(i, WvelIndx)] = u_vel_sum;,
        a_meanFab[MD_IX(i, WvelIndx+1)] = y_vel_sum;,
        a_meanFab[MD_IX(i, WvelIndx+2)] = 0.;); // z-velocity must be zero
      a_meanFab[MD_IX(i, tempIndx)] = T_mean_sum;
      a_meanFab[MD_IX(i, presIndx)] = p_inf;

      // // Estimate the streamwise velocity (Musker 1979)
      // const Real y_p = loc[1]*u_tau/nu_w;
      // const Real b_0 = 5.424*std::atan2((2.*y_p - 8.15), 16.7);
      // const Real b_1 = std::pow((y_p + 10.6), 9.6);
      // const Real b_2 = (y_p*y_p - 8.15*y_p + 86);
      // const Real b_3 = std::log10(b_1/(b_2*b_2));
      // const Real b_4 = loc[1]/delta;
      // const Real b_5 = b_4*b_4;
      // const Real b_6 = b_5*b_4;
      // const Real b_7 = 2.44*(pi1*(6.*b_5 - 4.*b_6) + b_5*(1. - b_4));
      // const Real u_vd = u_tau*(b_0 + b_3 - 3.52 + b_7);
      // Real u_vel = std::sin(u_vd*alpha/u_inf)*u_inf/alpha;
      // if (loc[1] > delta)
      //   {
      //     u_vel = u_inf;
      //   }
      // const Real eta_0 = u_tau*u_tau/nu_w;
      // a_etaFab[MD_IX(i, 0)] = eta_0;

      // // Estimate the temperature using Crocco-Busemann
      // const Real T_mean = T_inf*(T_w_star - c_0*u_vel*u_vel*(1./(u_inf*u_inf)));
      // // Compute mean density
      // const Real rho_mean = p_inf/(T_mean*R);

      // CH_assert(rho_mean > 0.);
      // CH_assert(T_mean > 0.);
      // CH_assert(u_vel < hiTol && u_vel > loTol);
      // a_meanFab[MD_IX(i, rhoIndx)] = rho_mean;
      // D_TERM(
      //   a_meanFab[MD_IX(i, WvelIndx)] = u_vel;,
      //   a_meanFab[MD_IX(i, WvelIndx+1)] = 0.;, // we don't really know y-vel
      //   a_meanFab[MD_IX(i, WvelIndx+2)] = 0.;); // z-velocity must be zero
      // a_meanFab[MD_IX(i, tempIndx)] = T_mean;
      // a_meanFab[MD_IX(i, presIndx)] = p_inf;
    }
}
