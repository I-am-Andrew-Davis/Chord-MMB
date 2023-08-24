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
 * \file CNSIBCRecirculatingInletTFP.cpp
 *
 * \brief Member functions for CNSIBCRecirculatingInletTFP
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

#include "CNSIBCRecirculatingInletTFP.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"
#include "PatchCNSOp.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodVMS.H"
#include "CRDutil.H"
#include "ViscousTensor4thOrderOpF_F.H"
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

CNSIBCRecirculatingInletTFP::CNSIBCRecirculatingInletTFP()
  :
  CNSIBCGeneralized(),
  m_delta(0.),
  m_alpha(4.),
  m_beta(0.2),
  m_phi1(0.5),
  m_phi2(40.),
  m_psi1(120.),
  m_psi2(200.),
  m_lambda(0.),
  m_samplePlaneLoc(0.75),
  m_yPlusMaxGuess(1.e10),
  m_velPerturb(0.),
  m_numLevels(0.),
  m_useAvgYVel(0),
  m_inflowMethod(-1),
  m_outflowMethod(-1),
  m_useWallModel(0)
{
  readBCInfo();

  for (int i = 0; i != m_numLevels; ++i)
    {
      std::unique_ptr<LevelData<FArrayBox>> singleBoxSampledPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxSampledPlane.push_back(std::move(singleBoxSampledPlane));
      std::unique_ptr<LevelData<FArrayBox>> singleBoxSampledInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxSampledInletPlane.push_back(
        std::move(singleBoxSampledInletPlane));
      std::unique_ptr<LevelData<FArrayBox>> singleBoxInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxInletPlane.push_back(std::move(singleBoxInletPlane));
      // DEBUG ONLY
      std::unique_ptr<LevelData<FArrayBox>> singleBoxInletPlaneRelaxed =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxInletPlaneRelaxed.push_back(
        std::move(singleBoxInletPlaneRelaxed));
      // END DEBUG ONLY
      // DEBUG ONLY
      std::unique_ptr<LevelData<FArrayBox>> singleBoxInletPlaneInst =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxInletPlaneInst.push_back(std::move(singleBoxInletPlaneInst));
      // END DEBUG ONLY
      std::unique_ptr<LevelData<FArrayBox>> singleBoxMeanInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_singleBoxMeanInletPlane.push_back(std::move(singleBoxMeanInletPlane));
      std::unique_ptr<LevelData<FArrayBox>> multiBoxInletPlane =
        std::make_unique<LevelData<FArrayBox>>();
      m_multiBoxInletPlane.push_back(std::move(multiBoxInletPlane));
      // DEBUG ONLY
      std::unique_ptr<LevelData<FArrayBox>> multiBoxInletPlaneRelaxed =
        std::make_unique<LevelData<FArrayBox>>();
      m_multiBoxInletPlaneRelaxed.push_back(
        std::move(multiBoxInletPlaneRelaxed));
      // END DEBUG ONLY
      // DEBUG ONLY
      std::unique_ptr<LevelData<FArrayBox>> multiBoxInletPlaneInst =
        std::make_unique<LevelData<FArrayBox>>();
      m_multiBoxInletPlaneInst.push_back(std::move(multiBoxInletPlaneInst));
      // END DEBUG ONLY
      std::unique_ptr<Copier> singleBoxSampledPlaneCopier =
        std::make_unique<Copier>();
      m_singleBoxSampledPlaneCopier.push_back(
        std::move(singleBoxSampledPlaneCopier));
      std::unique_ptr<Copier> singleBoxSampledInletPlaneCopier =
        std::make_unique<Copier>();
      m_singleBoxSampledInletPlaneCopier.push_back(
        std::move(singleBoxSampledInletPlaneCopier));
      std::unique_ptr<Copier> multiBoxInletPlaneCopier =
        std::make_unique<Copier>();
      m_multiBoxInletPlaneCopier.push_back(std::move(multiBoxInletPlaneCopier));
      // DEBUG ONLY
      std::unique_ptr<Copier> multiBoxInletPlaneRelaxedCopier =
        std::make_unique<Copier>();
      m_multiBoxInletPlaneRelaxedCopier.push_back(
        std::move(multiBoxInletPlaneRelaxedCopier));
      // END DEBUG ONLY
      // DEBUG ONLY
      std::unique_ptr<Copier> multiBoxInletPlaneInstCopier =
        std::make_unique<Copier>();
      m_multiBoxInletPlaneInstCopier.push_back(
        std::move(multiBoxInletPlaneInstCopier));
      // END DEBUG ONLY
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
      m_timeAfterInit.push_back(0.);
      m_stageCounter.push_back(0);
      m_currDt.push_back(0.);
      m_prevDt.push_back(0.);
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCRecirculatingInletTFP::~CNSIBCRecirculatingInletTFP()
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
CNSIBCRecirculatingInletTFP::IBCName() const
{
  return "Recycled-inlet turbulent flat-plate";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCRecirculatingInletTFP::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** For flat-plate, don't tag at inlet, outlet, or periodic boundaries
 *  \param[in]  a_tagBufferSize
 *                      Requested tag buffer size (should be
 *                      respected).
 *
 *  \note
 *  <ul>
 *    <li> Allocate levels and methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory*
CNSIBCRecirculatingInletTFP::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  //**FIXME: this would most likely be a problem with multiblock
  //**FIXME: this isn't entirely correct even for single-block
  //         We assume the lower end of the domain starts at (0,0,0)
  const Box domainBox(IntVect::Zero, CRDparam::g_domainBaseSize);
  const IntVect domLoVect = domainBox.smallEnd();
  const IntVect domHiVect = domainBox.bigEnd();
  // Inlet non-refined region
  Box inletBox = domainBox;
  inletBox.setBig(0, domainBox.smallEnd()[0] + 8);
  // Outlet non-refined region
  Box outletBox = domainBox;
  outletBox.setSmall(0, domainBox.bigEnd()[0] - 8);
  // Periodic non-refined regions
  Box leftPeriodicBox = domainBox;
  Box rightPeriodicBox = domainBox;
  D_TERM(,,leftPeriodicBox.setBig(2, domainBox.smallEnd()[2] + 4);
           rightPeriodicBox.setSmall(2, domainBox.bigEnd()[2] - 4););
  // Assign the boxes
  int numRestrictBoxes = (SpaceDim - 1)*2;
  std::vector<Box> restrictBoxes(numRestrictBoxes);
  restrictBoxes[0].define(inletBox);
  restrictBoxes[1].define(outletBox);
  if (SpaceDim > 2)
    {
      restrictBoxes[2].define(leftPeriodicBox);
      restrictBoxes[3].define(rightPeriodicBox);
    }
  CNSIBCCombustionReference::setTagMethodLevel(a_tagBufferSize, tagLevel,
                                               restrictBoxes);
  // Return the factory
  return new TagLevelFactory(tagLevel);
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
CNSIBCRecirculatingInletTFP::initialize(
  LevelData<FArrayBox>&      a_U,
  LevelGridMetrics&          a_gridMetrics,
  const LayoutData<FluxBox>& a_unitNormals,
  const Real                 a_time,
  const int                  a_level) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::initialize");
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

      FABSTACKTEMP(deltaFab, box4Dom, 1); // Bndry-layer-thickness fab
      FABSTACKTEMP(modXFab, box4Dom, SpaceDim);
      // Approximate the boundary-layer thickness for the initial condition
      MD_BOXLOOP(box4Dom, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)],
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          Real x = loc[0] + x_star; // x-location wrt to hypothetical plate tip
          Real Re_x = rho_inf*u_inf*x/mu; // location-based Reynolds number
          deltaFab[MD_IX(i, 0)] = 0.37*x/std::pow(Re_x, 0.2);
          D_TERM(modXFab[MD_IX(i, 0)] = x;,
                 modXFab[MD_IX(i, 1)] = loc[1];,
                 modXFab[MD_IX(i, 2)] = loc[2];);
        }
      FABSTACKTEMP(meanFab, box4Dom, numWVar);
      FABSTACKTEMP(etaFab, box4Dom, 1);
      // For the IC, we don't really need to offset the velocity profile
      Real virtualWallDy = 0.0;
      muskerIC(meanFab,etaFab,deltaFab,modXFab,XNodeFab,virtualWallDy,box4Dom);

      // Now we need to add perturbations to streamwise and spanwise velocity
      MD_BOXLOOP(box4Dom, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)],
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          RealVect physicalLength = CRDparam::g_physicalLength;
          RealVect numPerturbs = m_perturbFreq;
          RealVect arg = PI*numPerturbs*loc/(physicalLength);
          // Add perturbations to x-vel and y-vel together
          Real c_7 = std::max(0., loc[1] - 0.5*m_delta);
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
          D_TERM(meanFab[MD_IX(i,WvelIndx)]   = velTotal[0];,
                 meanFab[MD_IX(i,WvelIndx+1)] = velTotal[1];,
                 meanFab[MD_IX(i,WvelIndx+2)] = velTotal[2];);
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
      // NOTE: We overwrite the default SSV eta_0 initialization because we
      //       have an analytical estimate for the flat-plate initial eta_0
      if (CRDparam::g_turbModelType && m_useWallModel)
        {
          // Note that we only care about the y direction in this case
          int yDir = 1;
          const int turbCompBegin =
            CRDparam::g_CRDPhysics->turbConsInterval().begin();
          const int etaIndx = turbCompBegin + yDir + 1;
          Box boundaryFaceBoxLo;
          this->CNSIBC::getBoundaryFaces(
            boundaryFaceBoxLo, box4Dom, blockDomain, yDir, Side::Lo);
          if (!boundaryFaceBoxLo.isEmpty())
            {
              // Check boundary condition
              BoundaryIndex bcIdx;
              bcIdx.define(
                a_gridMetrics.getCoordSys().whichBlock(box), yDir, Side::Lo);
              BCInfo domBC = this->CNSIBC::getDomainBC(bcIdx);
              // We only support the use of walls in this function
              if ((CRDparam::DomainBCTypeAllWall & domBC.m_type) &&
                  !(CRDparam::DomainBCTypeMixed & domBC.m_type))
                {
                  MD_BOXLOOP(box2Dom, i)
                    {
                      Real eta_0 = etaFab[MD_IX(i, 0)];
                      UFab[MD_IX(i, etaIndx)] = eta_0;
                    }
                }
            }
        }
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
CNSIBCRecirculatingInletTFP::initializeInletDataStructures(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics,
  const bool                  a_hasFinerGrid) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::initializeInletDataStructures");
  CRD::msg << CRD::fv3
           << "CNSIBCRecirculatingInletTFP::initializeInletDataStructures"
           << " (level: " << a_level << ")" << CRD::end;
  if (a_hasFinerGrid || (a_level < (m_numLevels - 1)))
    {
      if (!m_levelDefined[a_level+1])
        {
          return;
        }
    }
  // If using multiBlock, we need to test all of the following setup very
  // diligently and carefully
  if (a_gridMetrics.isMultiBlock())
    {
      CRD::msg << "CNSIBCRecirculatingInletTFP::initializeInletDataStructures: "
               << "Be careful with multiBlock!" << CRD::error;
    }

  //**NOTE: cp = Copy Plane, sp = Sample Plane, ip = Inlet Plane
  // Create a disjointBoxLayout that covers the plane to copy from with one box

  //**NOTE: We assume that the origin of the domain is IntVect_zero in
  //        computational space
  // 1) Get the number of cells on the base mesh (this assumes single-block)
  const IntVect numCellsBaseMesh = CRDparam::g_domainBaseSize;
  // 2) Get the current physical domain box and number of cells on this level
  const IntVect domainSize =
    a_disjointBoxLayout.physDomain().domainBox().size();
  // 3) Get the number of cells between inlet and sampling plane
  const int spCellLoc = (domainSize[0] - 1)*m_samplePlaneLoc;
  // 4) Get the actual box at the sampling plane
  Box domainBox = a_disjointBoxLayout.physDomain().domainBox();
  int yHiTemp = domainBox.smallEnd()[1];
  domainBox.setBig(1, yHiTemp);
  Vector<Box> tempBoxes = a_disjointBoxLayout.boxArray();
  for (int boxComp = 0; boxComp != tempBoxes.size(); ++boxComp)
    {
      Box tempBox = tempBoxes[boxComp];
      // Get the boxes at the sampling plane location
      domainBox.minBox(tempBox);
    }
  // 4) Set up one box to cover the copy plane on this level
  Box cpSingleBox = domainBox;
  cpSingleBox.setSmall(0, domainBox.smallEnd()[0] + spCellLoc);
  cpSingleBox.setBig(0, domainBox.smallEnd()[0] + spCellLoc);
  Vector<Box> singlePlaneBoxes;
  singlePlaneBoxes.push_back(cpSingleBox);
  Vector<int> copyPlaneProcList;
  copyPlaneProcList.push_back(0);
  DisjointBoxLayout copyPlaneDBL(singlePlaneBoxes, copyPlaneProcList);

  // 5) Create a disjointBoxLayout that covers the inlet plane with one box
  Box ipSingleBox = domainBox;
  ipSingleBox.setBig(0, domainBox.smallEnd()[0]);
  Vector<Box> inletPlaneBoxes;
  inletPlaneBoxes.push_back(ipSingleBox);
  Vector<int> inletPlaneProcList;
  inletPlaneProcList.push_back(0);
  DisjointBoxLayout ipSingleBoxDBL(inletPlaneBoxes, inletPlaneProcList);

  // 6) Create a DBL that covers the inlet plane with multiple boxes
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
  m_singleBoxSampledInletPlane[a_level]->define(
    ipSingleBoxDBL, numConsComp, ghostVect);
  m_singleBoxInletPlane[a_level]->define(
    ipSingleBoxDBL, numPrimComp, ghostVect);
  // DEBUG ONLY
  m_singleBoxInletPlaneRelaxed[a_level]->define(
    ipSingleBoxDBL, numPrimComp, ghostVect);
  // END DEBUG ONLY
  // DEBUG ONLY
  m_singleBoxInletPlaneInst[a_level]->define(
    ipSingleBoxDBL, numPrimComp, ghostVect);
  // END DEBUG ONLY
  m_singleBoxMeanInletPlane[a_level]->define(
    copyPlaneDBL, numPrimComp, IntVect_zero);
  m_multiBoxInletPlane[a_level]->define(
    ipMultiBoxDBL, numPrimComp, ghostVect);
  // DEBUG ONLY
  m_multiBoxInletPlaneRelaxed[a_level]->define(
    ipMultiBoxDBL, numPrimComp, ghostVect);
  // END DEBUG ONLY
  // DEBUG ONLY
  m_multiBoxInletPlaneInst[a_level]->define(
    ipMultiBoxDBL, numPrimComp, ghostVect);
  // END DEBUG ONLY

  // Define copiers
  m_singleBoxSampledPlaneCopier[a_level]->define(a_disjointBoxLayout,
                                                 copyPlaneDBL);
  m_singleBoxSampledInletPlaneCopier[a_level]->define(a_disjointBoxLayout,
                                                      ipSingleBoxDBL);
  m_multiBoxInletPlaneCopier[a_level]->ghostDefine(
    m_singleBoxInletPlane[a_level]->disjointBoxLayout(),
    m_multiBoxInletPlane[a_level]->disjointBoxLayout(),
    ipProbDom, IntVect_zero, ghostVect);
  // DEBUG ONLY
  m_multiBoxInletPlaneRelaxedCopier[a_level]->define(
    m_multiBoxInletPlaneRelaxed[a_level]->disjointBoxLayout(),
    m_singleBoxInletPlaneRelaxed[a_level]->disjointBoxLayout());
  // END DEBUG ONLY
  // DEBUG ONLY
  m_multiBoxInletPlaneInstCopier[a_level]->define(
    m_multiBoxInletPlaneInst[a_level]->disjointBoxLayout(),
    m_singleBoxInletPlaneInst[a_level]->disjointBoxLayout());
  // END DEBUG ONLY
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
        = XFab[MD_IX(i, 1)];
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
CNSIBCRecirculatingInletTFP::copyInteriorToInletAndRescale(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const Real                  a_t,
  const Real                  a_dt,
  const int                   a_stage,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::copyInteriorToInletAndRescale");
  CRD::msg << CRD::fv3
           << "CNSIBCRecirculatingInletTFP::copyInteriorToInletAndRescale"
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

  // Constant freestream state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf   = state.velocity()[0]; // freestream velocity
  const Real rho_inf = state.density();     // freestream density
  const Real mu      = CRDparam::g_mu;

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
  // fill m_singleBoxSampledInletPlane with data from the interior inlet
  a_U.copyTo((*(m_singleBoxSampledInletPlane[a_level])),
             (*(m_singleBoxSampledInletPlaneCopier[a_level])));
  for (DataIterator ditSP = spDBL.dataIterator(); ditSP.ok(); ++ditSP)
    {
      DataIterator ditIP =
        m_singleBoxInletPlane[a_level]->getBoxes().dataIterator();
      FArrayBox& sampledPlaneU = (*(m_singleBoxSampledPlane[a_level]))[ditSP];
      FArrayBox& inletPlaneW = (*(m_singleBoxInletPlane[a_level]))[ditIP];
      Box spBox = spDBL[ditSP];
      Box ipBox = m_singleBoxInletPlane[a_level]->getBoxes()[ditIP];

      // convert the conservative variables to primitive variables
      FABSTACKTEMP(WcellPntFab, spBox, numWComp);
      PatchCNSOp::computeWpntCell(WcellPntFab,
                                  sampledPlaneU,
                                  sampledPlaneU, // dummy input that's not used
                                  spBox,
                                  a_disjointBoxLayout.physDomain(),
                                  true,
                                  false);
      // convert the conservative variables to primitive variables
      FArrayBox& inletPlaneInteriorU =
        (*(m_singleBoxSampledInletPlane[a_level]))[ditIP];
      Box spIBox = inletPlaneInteriorU.box();
      FABSTACKTEMP(WcellInletSampledPntFab, spIBox, numWComp);
      PatchCNSOp::computeWpntCell(WcellInletSampledPntFab,
                                  inletPlaneInteriorU,
                                  inletPlaneInteriorU, // dummy input, not used
                                  spIBox,
                                  a_disjointBoxLayout.physDomain(),
                                  true,
                                  false);

     CNSIBCRecirculatingInletTFP::spaceTimeAverageVariables(
       WcellPntFab, spBox, a_t, a_dt, a_level);

     // Compute the inlet mean state if not already computed
     FArrayBox& meanInletFab = (*(m_singleBoxMeanInletPlane[a_level]))[ditSP];
     int meanInitialized = m_inletMeanInitialized[a_level];
     if (!meanInitialized)
       {
         // Get some box from the current block to get the blockCoordSys
         Box box = a_U.disjointBoxLayout()[a_U.dataIterator()];
         // Get coordinate system and domain for the block
         const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
         Box faceBox = meanInletFab.box();
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
         Real virtualWallDy = 0.0;
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
                 const Real deltaX =
                   CRDparam::g_physicalLength[0]*m_samplePlaneLoc;
                 const Real b_0 = m_delta/0.37;
                 const Real b_1 = 5./4.;
                 const Real b_2 = rho_inf*u_inf/mu;
                 const Real x_star = std::pow(b_0, b_1)*std::pow(b_2, 0.25);
                 const Real x = x_star + deltaX;
                 const Real Re_x = rho_inf*u_inf*x/mu;
                 const Real delta_guess = 0.37*x/std::pow(Re_x, 0.2);
                 const Real delta_99_recy = 1.1*delta_guess;

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
         FABSTACKTEMP(etaInletFab, spBox, 1);
         muskerIC(meanInletFab, etaInletFab, deltaFab, inletCoord, XNodeFab,
                  virtualWallDy, spBox);
         m_inletMeanInitialized[a_level] = 1;
       }

     CNSIBCRecirculatingInletTFP::recycleFluctuations(
       inletPlaneW, meanInletFab, WcellPntFab, WcellInletSampledPntFab,
       spBox, ipBox, a_level, a_t);
    }
  // // DEBUG ONLY
  // // Plot m_singleBoxInletPlane if we are at the plot frequency
  // int plotFreq = 10;
  // int plotTest = m_timeCounter[a_level] % plotFreq;
  // if (!plotTest && !(procID()))
  //   {
  //     // Plot the instantaneous values
  //     const int currTimeStep = m_timeCounter[a_level];
  //     DataIterator ditIP =
  //       m_singleBoxInletPlane[a_level]->getBoxes().dataIterator();
  //     FArrayBox& inletPlaneW = (*(m_singleBoxInletPlane[a_level]))[ditIP];
  //     Box ipBox = m_singleBoxInletPlane[a_level]->getBoxes()[ditIP];
  //     FABSTACKTEMP(inletTempFab, ipBox, numWComp);
  //     inletTempFab.copy(inletPlaneW);
  //     Vector<string> names(numWComp);
  //     char compNameString[64];
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "inst-%s-target",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     char fileName[64];
  //     sprintf(fileName, "inletTargetValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletTempFab, fileName, names);

  //     // Plot the spanwise averaged values
  //     FABSTACKTEMP(inletAvgFab, ipBox, numWComp);
  //     inletAvgFab.setVal(0.);
  //     const IntVect nCells = ipBox.size();
  //     int nCellsZ = 1;
  //     D_TERM(,,nCellsZ = nCells[2];);
  //     const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  //     const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  //     const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  //     // Compute the spanwise-averaged primitive state -- once per time-step
  //     for (int yLoc = 0; yLoc != nCells[1]; ++yLoc)
  //       {
  //         // set up the one-dimensional averaging box at the current y-index
  //         IntVect loIntVect = ipBox.smallEnd();
  //         IntVect hiIntVect = ipBox.bigEnd();
  //         loIntVect[1] = yLoc;
  //         hiIntVect[1] = yLoc;
  //         Box averageBox(loIntVect, hiIntVect);

  //         // average the variables over the 1D box
  //         Real rhoAvg = 0.;
  //         D_TERM(Real uVelAvg = 0.;, Real vVelAvg = 0.;, Real wVelAvg = 0.;);
  //         Real TAvg = 0.;
  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             rhoAvg += inletTempFab[MD_IX(i, rhoIndx)];
  //             D_TERM(
  //               uVelAvg += inletTempFab[MD_IX(i, WvelIndx)];,
  //               vVelAvg += inletTempFab[MD_IX(i, WvelIndx+1)];,
  //               wVelAvg += inletTempFab[MD_IX(i, WvelIndx+2)];);
  //             TAvg += inletTempFab[MD_IX(i, tempIndx)];
  //           }

  //         // divide by the number of cells to finish off the spatial averaging
  //         rhoAvg /= nCellsZ;
  //         D_TERM(uVelAvg /= nCellsZ;, vVelAvg /= nCellsZ;, wVelAvg /= nCellsZ;);
  //         TAvg /= nCellsZ;

  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             inletAvgFab[MD_IX(i, rhoIndx)] = rhoAvg;
  //             D_TERM(inletAvgFab[MD_IX(i, WvelIndx)] = uVelAvg;,
  //                    inletAvgFab[MD_IX(i, WvelIndx+1)] = vVelAvg;,
  //                    inletAvgFab[MD_IX(i, WvelIndx+2)] = wVelAvg;);
  //             inletAvgFab[MD_IX(i, tempIndx)] = TAvg;
  //           }
  //       }
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "spanwiseAvg-%s-target",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     sprintf(fileName, "inletSpanwiseAvgTargetValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletAvgFab, fileName, names);
  //   }
  // // END DEBUG ONLY

  // First, copy m_singleBoxInletPlane to m_multiBoxInletPlane
  const Interval copyIntv(0, numWComp-1);
  m_singleBoxInletPlane[a_level]->copyTo(
    copyIntv, (*(m_multiBoxInletPlane[a_level])),
    copyIntv, (*(m_multiBoxInletPlaneCopier[a_level])));
  // Then fill periodic ghost cells of m_multiBoxInletPlane
  m_multiBoxInletPlane[a_level]->exchange(
    (*(m_multiBoxInletExchangeCopier[a_level])));
  // Finally, fill multiBoxInletPlane invalid ghosts from singleBoxInletPlane
  m_singleBoxInletPlane[a_level]->copyTo(
    copyIntv, (*(m_multiBoxInletPlane[a_level])),
    copyIntv, (*(m_multiBoxInletPlaneCopier[a_level])));

  // Update the current timeAfterInit
  m_timeAfterInit[a_level] += a_dt;
}

/*--------------------------------------------------------------------*/
//  Add body force
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg
 *                      Cell-averaged conservative variables
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
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
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCRecirculatingInletTFP::addSourceTerm(
  FArrayBox&           a_sourceFab,
  FArrayBox&           a_invDtFab,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_UcellAvg,
  const FluxBox&       a_WfaceAvgFxb,
  const ProblemDomain& a_domain,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const Real           a_stageWeight,
  const int            a_level,
  const Box&           a_disjointBox,
  const Box&           a_solveBox,
  const DataIndex&     a_didx,
  const Real           a_GlobalKE,
  const Real           a_globalHelicity) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::addSourceTerm");
  // Boxes on which to work
  Box box1Dom = grow(a_solveBox, 1);
  box1Dom &= a_domain;

  // Index variables
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cPres = CRDparam::g_CRDPhysics->pressureIndex();
  const int energyIndx = CRDparam::g_CRDPhysics->energyFluxIndex();

  // Freestream state variables
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real rho_inf = state.density();
  const Real u_inf = state.velocity()[0];
  const Real p_inf = state.pressure();
  const Real gamma = CRDparam::g_gamma;

  // ----------------------------- //
  //  Sponge Region Forcing Terms  //
  // ----------------------------- //

  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_solveBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_solveBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(a_solveBox, XiFab, XFab, blockCoordSys);
  MD_BOXLOOP(a_solveBox, i)
    {
      RealVect loc(
        D_DECL(XFab[MD_IX(i,0)],XFab[MD_IX(i,1)],XFab[MD_IX(i,2)]));
      Real c_0 = 0.;
      if ((loc[1] >= 0.085852) && (loc[1] <= 0.128778))
        {
          Real x1 = (-loc[1] + 0.128778)/(0.128778 - 0.085852);
          Real x3 = x1*x1*x1;
          Real x4 = x3*x1;
          Real x5 = x4*x1;
          c_0 = 1. - 6.*x5 + 15.*x4 - 10.*x3;
        }
      else if (loc[1] > 0.128778)
        {
          c_0 = 1.;
        }
      // Only apply sponge to pressure and density
      Real c_1 = 100.0;
      Real rho = a_Wcell[MD_IX(i, rhoIndx)];
      Real u_vel = a_Wcell[MD_IX(i, velIndx)];
      Real press = a_Wcell[MD_IX(i, cPres)];
      a_sourceFab[MD_IX(i, rhoIndx)] = c_1*c_0*(rho_inf - rho);
      a_sourceFab[MD_IX(i, velIndx)] = c_1*c_0*(rho_inf*u_inf - rho*u_vel);
      a_sourceFab[MD_IX(i, energyIndx)] = c_1*c_0*(p_inf - press)/(gamma - 1.);
    }

  // // Compute cell-averaged Cartesian gradients of velocity
  // FABSTACKTEMP(cellAvgGradVel, box1Dom, SpaceDim*SpaceDim);
  // PatchMappedFunc::gradient2OCS(
  //   box1Dom, cellAvgGradVel, a_Wcell, a_domain,
  //   CRDparam::g_CRDPhysics->velocityInterval(), a_gridMetrics.dxVect());
  // // Interpolate cell-averaged gradients of velocity to faces
  // const int interpOrder = 2;
  // FLUXBOXSTACKTEMP(faceAvgGradVel, a_solveBox, SpaceDim*SpaceDim);
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     Box faceBox = surroundingNodes(a_solveBox, dir);
  //     for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
  //       {
  //         PatchMappedFunc::faceAvgValFromCellAvgCS(
  //           faceAvgGradVel[dir], cellAvgGradVel, comp, comp,
  //           a_domain, faceBox, dir, interpOrder, false);
  //       }
  //   }
  // // Map the velocity gradients to physical space
  // FLUXBOXSTACKTEMP(faceAvgPhysGradVel, a_solveBox, SpaceDim*SpaceDim);
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     Box faceBox = surroundingNodes(a_solveBox, dir);
  //     PatchMappedFunc::gradientCStoPS(
  //       faceBox, faceAvgPhysGradVel[dir],
  //       faceAvgGradVel[dir], a_gridMetrics.m_NtJ[a_didx][dir]);
  //   }
  // // Get coordinate system and domain for the block
  // const BlockCoordSys& blockCoordSys =
  //   *(a_gridMetrics.getCoordSys(a_disjointBox));
  // // Compute the distance from the inlet
  // FLUXBOXSTACKTEMP(XFxb, a_solveBox, SpaceDim);
  // FLUXBOXSTACKTEMP(XiFxb, a_solveBox, SpaceDim);
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     Box faceBox = surroundingNodes(a_solveBox, dir);
  //     this->CNSIBC::getFaceCoordinates(
  //       faceBox, XiFxb[dir], XFxb[dir], dir, blockCoordSys);
  //   }
  // // Variables for distance-based attenuation coefficient
  // const Real numDeltas = 10.; // should be user defined
  // const Real forceDist = numDeltas*m_delta;
  // // Compute the divergence of velocity
  // FLUXBOXSTACKTEMP(faceAvgPhysVelDiv, a_solveBox, 1);
  // const IntVect divIndx(D_DECL(PatchMappedFunc::getGradientComp(0,0),
  //                              PatchMappedFunc::getGradientComp(1,1),
  //                              PatchMappedFunc::getGradientComp(2,2)));
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     Box faceBox = a_solveBox;
  //     faceBox.surroundingNodes(dir);
  //     FArrayBox& gradVel = faceAvgPhysGradVel[dir];
  //     FArrayBox& velDiv = faceAvgPhysVelDiv[dir];
  //     FArrayBox& XFab = XFxb[dir];
  //     MD_BOXLOOP(faceBox, i)
  //       {
  //         Real c0 = 0.0001; // diffusion coefficient -- should be user defined
  //         // Attenuation coefficient based on physical location
  //         Real c1 = std::max((1. - XFab[MD_IX(i,0)])/forceDist, 0.);
  //         velDiv[MD_IX(i, 0)] = c0*c1*c1*(
  //           D_TERM(gradVel[MD_IX(i, divIndx[0])],
  //                + gradVel[MD_IX(i, divIndx[1])],
  //                + gradVel[MD_IX(i, divIndx[2])]));
  //       }
  //   }
  // // Set up the fluxes for momentum and energy
  // //**NOTE: SpaceDim*(SpaceDim + 1) is for momentum (SpaceDim) and energy (1)
  // FLUXBOXSTACKTEMP(fluxSourceAvg, a_solveBox, SpaceDim*(SpaceDim + 1));
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     // First, we zero the flux values
  //     FArrayBox& fluxSourceAvgFab = fluxSourceAvg[dir];
  //     fluxSourceAvgFab.setVal(0.);
  //     // Set up the momentum flux here -- loop over the velocity components
  //     const FArrayBox& velDiv = faceAvgPhysVelDiv[dir];
  //     for (int comp = 0; comp != SpaceDim; ++comp)
  //       {
  //         for (int gradDir = 0; gradDir != SpaceDim; ++gradDir)
  //           {
  //             // Dissipation only contributes to the velocity-direction flux
  //             if (comp == gradDir)
  //               {
  //                 int compIndx = gradDir + comp*SpaceDim;
  //                 MD_BOXLOOP(a_solveBox, i)
  //                   {
  //                     fluxSourceAvgFab[MD_IX(i, compIndx)] =
  //                       velDiv[MD_IX(i, 0)];
  //                   }
  //               }
  //           }
  //       }

  //     // Set up the energy flux here
  //     for (int comp = 0; comp != SpaceDim; ++comp)
  //       {
  //         int compIndx = SpaceDim*SpaceDim + comp;
  //         const FArrayBox& faceAvgWFab = a_WfaceAvgFxb[dir];
  //         MD_BOXLOOP(a_solveBox, i)
  //           {
  //             fluxSourceAvgFab[MD_IX(i, compIndx)] =
  //               faceAvgWFab[MD_IX(i, velIndx + comp)]*velDiv[MD_IX(i, 0)];
  //           }
  //       }
  //   }

  // // Transform the fluxes for mapped grids
  // const int numFluxes = SpaceDim + 1;
  // FLUXBOXSTACKTEMP(NtSourceFlux, a_solveBox, numFluxes);
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     NtSourceFlux[dir].setVal(0.);
  //   }
  // bool fourthOrder = false;
  // const FluxBox& Nfxb = a_gridMetrics.m_N[a_didx];
  // const Interval fluxIntv(0, numFluxes - 1);
  // blockCoordSys.computeMetricTermProductAverage(NtSourceFlux,
  //                                               fluxSourceAvg,
  //                                               Nfxb,
  //                                               SpaceDim,
  //                                               fluxSourceAvg,
  //                                               a_solveBox,
  //                                               fourthOrder,
  //                                               fluxIntv,
  //                                               fluxIntv,
  //                                               1,
  //                                               &a_domain);

  // // Compute the divergence of the flux terms to generate source terms
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     const Real dirDx = a_gridMetrics.dxVect()[dir];
  //     const FArrayBox& NtSourceFluxFab = NtSourceFlux[dir];
  //     const int MD_ID(o, dir);
  //     // Momentum source terms
  //     for (int comp = 0; comp != SpaceDim; ++comp)
  //       {
  //         MD_BOXLOOP(a_solveBox, i)
  //           {
  //             Real hiFlux = NtSourceFluxFab[MD_OFFSETIX(i,+,o,comp)];
  //             Real loFlux = NtSourceFluxFab[MD_IX(i,comp)];
  //             a_sourceFab[MD_IX(i, velIndx + comp)] -= (hiFlux - loFlux)/dirDx;
  //           }
  //       }
  //     const int energyComp = SpaceDim;
  //     MD_BOXLOOP(a_solveBox, i)
  //       {
  //         Real hiFlux = NtSourceFluxFab[MD_OFFSETIX(i,+,o,energyComp)];
  //         Real loFlux = NtSourceFluxFab[MD_IX(i,energyComp)];
  //         a_sourceFab[MD_IX(i, energyIndx)] -= (hiFlux - loFlux)/dirDx;
  //       }
  //   }

  // // Get physical coordinates
  // FABSTACKTEMP(XiFab, a_solveBox, SpaceDim);  // Cartesian coordinates
  // FABSTACKTEMP(XFab, a_solveBox, SpaceDim);   // Physical coordinates
  // this->CNSIBC::getCellCoordinates(a_solveBox, XiFab, XFab, blockCoordSys);
  // // Modify density just slightly
  // const Real densityCoeff = 0.001;
  // MD_BOXLOOP(a_solveBox, i)
  //   {
  //     // Attenuation coefficient based on physical location
  //     Real c1 = std::max((1. - XFab[MD_IX(i,0)])/forceDist, 0.);
  //     a_sourceFab[MD_IX(i, rhoIndx)] +=
  //       c1*c1*densityCoeff*(rho_inf - a_Wcell[MD_IX(i, rhoIndx)]);
  //   }
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCRecirculatingInletTFP::haveExactSol() const
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
CNSIBCRecirculatingInletTFP::setImposedBCprimState(
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
  CH_TIME("CNSIBCRecirculatingInletTFP::setImposedBCprimState");
  // Set all the component variables and intervals
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int cVel       = velIntv.begin();
  const int cRho       = CRDparam::g_CRDPhysics->densityIndex();
  const int cPres      = CRDparam::g_CRDPhysics->pressureIndex();
  const int cTemp      = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const Real Rgas  = CRDparam::g_R;
  const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
  // const int numWComp = CRDparam::g_CRDPhysics->numPrimitive();

  // Constant freestream state
  // const Real u_inf   = state.velocity()[0]; // freestream velocity
  // const Real T_inf   = state.temperature(); // freestream temperature
  // const Real gamma   = CRDparam::g_gamma;
  // const Real Pr_t = 0.89; // Turbulent Prandtl number (Urbin & Knight 2001)
  // const Real M_inf_sq = u_inf*u_inf/(gamma*Rgas*T_inf);
  // const Real d_0 = 0.5*(gamma - 1.)*M_inf_sq*Pr_t;

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    case CRDparam::DomainBCTypeInflow:
    {
      Box exteriorCellBox = a_boundaryFaceBox;
      exteriorCellBox.growLo(0, 1);
      exteriorCellBox.enclosedCells(0);
      if ((m_inflowMethod == 0) || (m_inflowMethod == 2))
        {
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
          // // DEBUG ONLY
          // FArrayBox& relaxedFab =
          //   (*(m_multiBoxInletPlaneRelaxed[a_level]))[ditIPCurrent];
          // relaxedFab.shiftHalf(0, -1);
          // FArrayBox& instFab =
          //   (*(m_multiBoxInletPlaneInst[a_level]))[ditIPCurrent];
          // instFab.shiftHalf(0, -1);
          // MD_BOXLOOP(a_boundaryFaceBox, i)
          //   {
          //     instFab[MD_IX(i, cRho)] = a_Wface[MD_IX(i, cRho)];
          //     D_TERM(instFab[MD_IX(i, cVel)] = a_Wface[MD_IX(i, cVel)];,
          //            instFab[MD_IX(i, cVel+1)] = a_Wface[MD_IX(i, cVel+1)];,
          //            instFab[MD_IX(i, cVel+2)] = a_Wface[MD_IX(i, cVel+2)];);
          //     instFab[MD_IX(i, cPres)] = a_Wface[MD_IX(i, cPres)];
          //     instFab[MD_IX(i, cTemp)] = a_Wface[MD_IX(i, cTemp)];
          //   }
          // // END DEBUG ONLY

          // Apply the state to the boundary face
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              RealVect velInlet(D_DECL(scaledFab[MD_IX(i, cVel)],
                                       scaledFab[MD_IX(i, cVel+1)],
                                       scaledFab[MD_IX(i, cVel+2)]));
              Real rhoInlet = scaledFab[MD_IX(i, cRho)];
              Real tempInlet = scaledFab[MD_IX(i, cPres)]/(Rgas*rhoInlet);
              for (const int idxVel : EachDir)
                {
                  a_Wface[MD_IX(i, cVel + idxVel)] =
                    velInlet[idxVel];
                }
              a_Wface[MD_IX(i, cPres)] = scaledFab[MD_IX(i, cPres)];
              // a_Wface[MD_IX(i, cPres)] = state.pressure();
              // Inflow cannot become outflow -- always set rho and temperature
              for (int j = 0; j != numSpecies; ++j)
                {
                  a_Wface[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
                }
              a_Wface[MD_IX(i, cTemp)] = tempInlet;
              // a_Wface[MD_IX(i, cRho)] =
              //   a_Wface[MD_IX(i, cPres)]/(Rgas*a_Wface[MD_IX(i, cTemp)]);
              a_Wface[MD_IX(i, cRho)] = rhoInlet;

              // // DEBUG ONLY
              // relaxedFab[MD_IX(i, cRho)] = a_Wface[MD_IX(i, cRho)];
              // D_TERM(relaxedFab[MD_IX(i,cVel)] = a_Wface[MD_IX(i, cVel)];,
              //        relaxedFab[MD_IX(i,cVel+1)] = a_Wface[MD_IX(i,cVel+1)];,
              //       relaxedFab[MD_IX(i,cVel+2)] = a_Wface[MD_IX(i,cVel+2)];);
              // relaxedFab[MD_IX(i, cTemp)] = a_Wface[MD_IX(i, cTemp)];
              // relaxedFab[MD_IX(i, cPres)] = a_Wface[MD_IX(i, cPres)];
              // // END DEBUG ONLY
            }
          scaledFab.shiftHalf(0, 1);
          // // DEBUG ONLY
          // // First try a Riemann solution between instFab and relaxedFab
          // FABSTACKTEMP(tempRiemannFab, a_boundaryFaceBox, numWComp);
          // FABSTACKTEMP(tempRelaxedFab, a_boundaryFaceBox, numWComp);
          // tempRelaxedFab.copy(relaxedFab);
          // FABSTACKTEMP(tempInstFab, a_boundaryFaceBox, numWComp);
          // tempInstFab.copy(instFab);
          // CRDparam::g_CRDPhysics->riemannBC(tempRiemannFab,
          //                                   tempRelaxedFab,
          //                                   tempInstFab,
          //                                   a_unitNormalBasisFab,
          //                                   0,
          //                                   Side::Lo,
          //                                   a_boundaryFaceBox);
          // instFab.copy(tempRiemannFab);
          // instFab.shiftHalf(0, 1);
          // relaxedFab.shiftHalf(0, 1);
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    case CRDparam::DomainBCTypeFarfield:
    {
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    case CRDparam::DomainBCTypeOutflow:
    {
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCIn:
    {
      break;
    }
    default:
      CH_assert(false);
    }
}

/*--------------------------------------------------------------------*/
//  Set exterior face-state for relaxed characteristic BCs
/** \param[out] a_WfaceAvgExterior
 *                      Exterior face-averaged value for Riemann solve
 *  \param[in]  a_bndryCellFab
 *                      Cell-averaged values exterior/interior to bndry.
 *                      These values are at current and previous
 *                      time-points
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged values (including interior cells)
 *  \param[in]  a_WfaceAvgDirFab
 *                      Face-averaged values (including interior faces)
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_bcInfo
 *                      Current boundary type
 *  \param[in]  a_disjointBox
 *                      Comment
 *  \param[in]  a_totalFaceBox
 *                      Comment
 *  \param[in]  a_gridMetrics
 *                      Grid-metrics on the current level
 *  \param[in]  a_dir   Face-normal direction
 *  \param[in]  a_side  Side of box on which bndry resides (low/high)
 *  \param[in]  a_time  Current time
 *  \param[in]  a_prevDt
 *                      Time-step size at the previous time-step
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBCRecirculatingInletTFP::setRelaxedCBCPrimState(
  FArrayBox&              a_WfaceAvgExterior,
  FArrayBox&              a_bndryCellFab,
  FArrayBox&              a_WcellAvgFab,
  const FArrayBox&        a_WfaceAvgDirFab,
  const FArrayBox&        a_unitNormalBasisFab,
  const BoundaryIndex&    a_bcIdx,
  const BCInfo&           a_bcInfo,
  const Box&              a_disjointBox,
  const Box&              a_totalFaceBox,
  const LevelGridMetrics& a_gridMetrics,
  const int               a_dir,
  const Side::LoHiSide    a_side,
  const Real              a_time,
  const Real              a_prevDt,
  const int               a_level) const
{
  CH_TIME("CNSIBC::setRelaxedCBCPrimState");
  // Set up some boxes for determining interior and exterior cells
  Box interiorCellBox = a_totalFaceBox;
  Box exteriorCellBox = a_totalFaceBox;
  Box totalCellBox = a_totalFaceBox;
  totalCellBox.grow(a_dir, 1); // Grow both the low and high sides
  totalCellBox.enclosedCells(a_dir); // Cell-box covering in/exterior cells
  if (a_side == Side::Lo)
    {
      interiorCellBox.growHi(a_dir, 1);
      interiorCellBox.enclosedCells(a_dir);
      exteriorCellBox.growLo(a_dir, 1);
      exteriorCellBox.enclosedCells(a_dir);
    }
  else
    {
      interiorCellBox.growLo(a_dir, 1);
      interiorCellBox.enclosedCells(a_dir);
      exteriorCellBox.growHi(a_dir, 1);
      exteriorCellBox.enclosedCells(a_dir);
    }
  int numPrim = CRDparam::g_CRDPhysics->numPrimitive();

  // Set all the component variables and intervals
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int cVel       = velIntv.begin();
  const int cRho       = CRDparam::g_CRDPhysics->densityIndex();
  const int cPres      = CRDparam::g_CRDPhysics->pressureIndex();
  const int cTemp      = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const Real Rgas  = CRDparam::g_R;
  const Real gamma = CRDparam::g_gamma;
  const Real f_0 = 1./m_prevDt[a_level]; // for use in computing the time-deriv

  const CRDState& state = CRDState::get(a_bcInfo.m_idxState);

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    {
      // If this is the first stage of the first step, we need to fill the
      // previous time-values with the current time value
      int numStages = 4;
      if (CRDparam::g_additiveRK)
        {
          numStages = 6;
        }
      int currStage = m_stageCounter[a_level] % numStages;
      if (m_stageCounter[a_level] == 0)
        {
          for (int comp = 0; comp != numPrim; ++comp)
            {
              const int pastIndx = comp + numPrim;
              MD_BOXLOOP(interiorCellBox, i)
                {
                  a_bndryCellFab[MD_IX(i, pastIndx)] =
                    a_WcellAvgFab[MD_IX(i, comp)];
                }
            }
        }
      // Compute the wave-amplitudes for the interior state at the current time
      // using time-data
      const int numWaves = 1 + SpaceDim + 1; // we're excluding other waves here
      FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
      MD_BOXLOOP(interiorCellBox, i)
        {
          // Primitive state at current time
          const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
          const RealVect currVel(
            D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                   a_WcellAvgFab[MD_IX(i, cVel+1)],
                   a_WcellAvgFab[MD_IX(i, cVel+2)]));
          const Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
          // Primitive state at previous time
          const Real prevRho = a_bndryCellFab[MD_IX(i, numPrim+cRho)];
          const RealVect prevVel(
            D_DECL(a_bndryCellFab[MD_IX(i, numPrim+cVel)],
                   a_bndryCellFab[MD_IX(i, numPrim+cVel+1)],
                   a_bndryCellFab[MD_IX(i, numPrim+cVel+2)]));
          const Real prevPress = a_bndryCellFab[MD_IX(i, numPrim+cPres)];
          // Time derivatives of the primitive state
          const Real rhoTimeDeriv = f_0*(currRho - prevRho);
          const RealVect velTimeDeriv = f_0*(currVel - prevVel);
          const Real pressTimeDeriv = f_0*(currPress - prevPress);

          // Current speed of sound
          const Real currC = std::sqrt(gamma*currPress/currRho);

          // Wave-amplitudes of interior state
          // (1) left-moving acoustic (downward)
          intWaveAmpCurrTime[MD_IX(i, cRho)] =
            currRho*currC*velTimeDeriv[1] - pressTimeDeriv;
          //**NOTE: Remember that y-velocity is the acoustic direction
          D_TERM(
            // (2) entropy wave in x-direction
            intWaveAmpCurrTime[MD_IX(i, cVel)] = -velTimeDeriv[0];,
            // (3) entropy wave in y-direction
            intWaveAmpCurrTime[MD_IX(i, cVel+1)] =
              pressTimeDeriv - currC*currC*rhoTimeDeriv;,
            // (4) entropy wave in z-direction
            intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
          // (5) right-moving acoustic (upward)
          intWaveAmpCurrTime[MD_IX(i, cPres)] =
            -pressTimeDeriv - currRho*currC*velTimeDeriv[1];
        }

      //**NOTE: The exterior wave amplitudes are assumed to be zero

      // Advance the previous-time exterior state using the current time
      // exterior wave-amplitudes and the previous-time interior wave-amplitudes
      FABSTACKTEMP(extState2, interiorCellBox, numPrim);
      MD_BOXLOOP(interiorCellBox, i)
        {
          // Primitive state at current time
          const Real currRho = state.density();
          const Real currPress = state.pressure();
          // Current speed of sound
          const Real currC = std::sqrt(gamma*currPress/currRho);
          // Primitive state at previous time
          const Real prevRho = currRho;
          const RealVect prevVel(D_DECL(state.velocity()[0],
                                        state.velocity()[1],
                                        state.velocity()[2]));
          const Real prevPress = currPress;
          // Interior primitive state at previous time
          const RealVect prevVelInt(
            D_DECL(a_bndryCellFab[MD_IX(i, numPrim+cVel)],
                   a_bndryCellFab[MD_IX(i, numPrim+cVel+1)],
                   a_bndryCellFab[MD_IX(i, numPrim+cVel+2)]));
          // Waves -- L1, L2, L4 come from exterior, L3, L5 come from interior
          const Real L_1 = 0.;
          const RealVect L_vel(D_DECL(0.,
                                      intWaveAmpCurrTime[MD_IX(i, cVel+1)],
                                      0.));
          const Real L_5 = intWaveAmpCurrTime[MD_IX(i, cPres)];
          // Advance the primitive state in time
          extState2[MD_IX(i, cRho)] = prevRho
            + m_prevDt[a_level]*((-1./(currC*currC))*(
                                   L_vel[1] + 0.5*(L_5 + L_1)));
          D_TERM(
            extState2[MD_IX(i, cVel)] = prevVel[0]
            + m_prevDt[a_level]*(-L_vel[0]);,
            extState2[MD_IX(i, cVel+1)] = prevVelInt[1]
            + m_prevDt[a_level]*((-1./(2.*currRho*currC))*(L_5 - L_1));,
            extState2[MD_IX(i, cVel+2)] = prevVel[2]
            + m_prevDt[a_level]*(-L_vel[2]););
          extState2[MD_IX(i, cPres)] = prevPress
            + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
        }

      // Set all variables to characteristic farfield state
      if (m_timeCounter[a_level] < 0)
        {
          // Set all variables except vertical velocity to farfield state
          for (int comp = 0; comp != numPrim; ++comp)
            {
              if (comp != cVel+1)
                {
                  const Real stateVal = state(comp);
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_WfaceAvgExterior[MD_IX(i, comp)] = stateVal;
                    }
                }
            }
        }
      else
        {
          extState2.shift(a_dir, 1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velFarfield(D_DECL(extState2[MD_IX(i, cVel)],
                                          extState2[MD_IX(i, cVel+1)],
                                          extState2[MD_IX(i, cVel+2)]));
              Real rhoFarfield = extState2[MD_IX(i, cRho)];
              Real pressFarfield = extState2[MD_IX(i, cPres)];
              Real tempFarfield = pressFarfield/(Rgas*rhoFarfield);

              // First set all components to the state values -- takes care of
              // anything not currently adjusted by characteristic formulation
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_WfaceAvgExterior[MD_IX(i, comp)] = state(comp);
                }
              // Set all values currently adjusted by characteristics
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoFarfield;
              D_TERM(
                a_WfaceAvgExterior[MD_IX(i, cVel)] = velFarfield[0];,
                a_WfaceAvgExterior[MD_IX(i, cVel+1)] = velFarfield[1];,
                a_WfaceAvgExterior[MD_IX(i, cVel+2)] = velFarfield[2];);
              a_WfaceAvgExterior[MD_IX(i, cPres)] = pressFarfield;
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempFarfield;
            }
        }
      // Move the current time-values to the previous times
      if (m_timeCounter[a_level])
        {
          if (currStage == (numStages - 1))
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    {
      if (m_outflowMethod == 0)
        {
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Compute wave-amplitudes for interior state at current time
          // using time-data
          const int numWaves = 1 + SpaceDim + 1; // excluding other waves here
          FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                                            a_WcellAvgFab[MD_IX(i, cVel+1)],
                                            a_WcellAvgFab[MD_IX(i, cVel+2)]));
              const Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                intWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          //**NOTE: The exterior wave amplitudes are assumed to be zero

          // Advance previous-time exterior state using current time exterior
          // wave-amplitudes and previous-time interior wave-amplitudes
          FABSTACKTEMP(extState2, interiorCellBox, numPrim);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = state.density();
              const Real currPress = state.pressure();
              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);
              // Primitive state at previous time
              const Real prevRho = currRho;
              const RealVect prevVel(D_DECL(state.velocity()[0],
                                            state.velocity()[1],
                                            state.velocity()[2]));
              const Real prevPress = currPress;
              // Interior primitive state at previous time
              const RealVect prevVelInt(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              // Waves -- L1 comes from exterior, all others come from interior
              const Real L_1 = 0.;
              const RealVect L_vel(
                D_DECL(intWaveAmpCurrTime[MD_IX(i, cVel)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+1)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+2)]));
              const Real L_5 = intWaveAmpCurrTime[MD_IX(i, cPres)];
              // Advance the primitive state in time
              extState2[MD_IX(i, cRho)] = prevRho
                + m_prevDt[a_level]*((-1./(currC*currC))*(
                                       L_vel[0] + 0.5*(L_5 + L_1)));
              D_TERM(
                extState2[MD_IX(i, cVel)] = prevVelInt[0]
                + m_prevDt[a_level]*((-1./(2.*currRho*currC))*(L_5 - L_1));,
                extState2[MD_IX(i, cVel+1)] = prevVelInt[1]
                + m_prevDt[a_level]*(-L_vel[1]);,
                extState2[MD_IX(i, cVel+2)] = prevVelInt[2]
                + m_prevDt[a_level]*(-L_vel[2]););
              extState2[MD_IX(i, cPres)] = prevPress
                + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
            }

          // Set all variables to characteristic outflow state
          extState2.shift(a_dir, 1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velFarfield(D_DECL(extState2[MD_IX(i, cVel)],
                                          extState2[MD_IX(i, cVel+1)],
                                          extState2[MD_IX(i, cVel+2)]));
              Real rhoFarfield = extState2[MD_IX(i, cRho)];
              Real pressFarfield = extState2[MD_IX(i, cPres)];
              Real tempFarfield = pressFarfield/(Rgas*rhoFarfield);

              // First set all components to the state values -- takes care of
              // anything not currently adjusted by characteristic formulation
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_WfaceAvgExterior[MD_IX(i, comp)] = state(comp);
                }
              // Set all values currently adjusted by characteristics
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoFarfield;
              D_TERM(
                a_WfaceAvgExterior[MD_IX(i, cVel)] = velFarfield[0];,
                a_WfaceAvgExterior[MD_IX(i, cVel+1)] = velFarfield[1];,
                a_WfaceAvgExterior[MD_IX(i, cVel+2)] = velFarfield[2];);
              a_WfaceAvgExterior[MD_IX(i, cPres)] = pressFarfield;
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempFarfield;
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
        }
      else if (m_outflowMethod == 1)
        {
          const int MD_ID(o, a_dir);
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_OFFSETIX(i,-,o,comp)];
                    }
                }
            }
          // Compute wave-amplitudes for interior state at current time
          // using time-data
          const int numWaves = 1 + SpaceDim + 1; // excluding other waves here
          FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                                            a_WcellAvgFab[MD_IX(i, cVel+1)],
                                            a_WcellAvgFab[MD_IX(i, cVel+2)]));
              const Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                intWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Compute wave-amplitudes for exterior state at current time
          // using time-data
          FABSTACKTEMP(extWaveAmpCurrTime, exteriorCellBox, numWaves);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_OFFSETIX(i,+,o,cRho)];
              const RealVect currVel(
                D_DECL(a_WcellAvgFab[MD_OFFSETIX(i,+,o,cVel)],
                       a_WcellAvgFab[MD_OFFSETIX(i,+,o,cVel+1)],
                       a_WcellAvgFab[MD_OFFSETIX(i,+,o,cVel+2)]));
              const Real currPress = state.pressure();
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                extWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Advance previous-time exterior state using current time exterior
          // wave-amplitudes and previous-time interior wave-amplitudes
          FABSTACKTEMP(extState2, exteriorCellBox, numPrim);
          intWaveAmpCurrTime.shift(a_dir, 1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              // Primitive state at previous time
              const RealVect prevVel(D_DECL(state.velocity()[0],
                                            state.velocity()[1],
                                            state.velocity()[2]));
              // Interior primitive state at previous time
              const RealVect prevVelInt(
                D_DECL(a_bndryCellFab[MD_OFFSETIX(i,-,o,cVel+numPrim)],
                       a_bndryCellFab[MD_OFFSETIX(i,-,o,cVel+1+numPrim)],
                       a_bndryCellFab[MD_OFFSETIX(i,-,o,cVel+2+numPrim)]));
              // Exterior primitive state at previous time
              const Real prevRhoExt = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const Real prevPressExt = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              const RealVect prevVelExt(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              // Waves -- L1 comes from exterior, all others come from interior
              const Real L_1 = 0.1*extWaveAmpCurrTime[MD_IX(i, cRho)];
              const RealVect L_vel(
                D_DECL(intWaveAmpCurrTime[MD_IX(i, cVel)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+1)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+2)]));
              const Real L_5 = intWaveAmpCurrTime[MD_IX(i, cPres)];
              // Exterior primitive state at current time
              const Real currRhoExt = a_bndryCellFab[MD_IX(i, cRho)];
              const Real currPressExt = a_bndryCellFab[MD_IX(i, cPres)];
              const Real currCExt = std::sqrt(gamma*currPressExt/currRhoExt);
              // Advance the primitive state in time
              extState2[MD_IX(i, cRho)] = prevRhoExt
                + m_prevDt[a_level]*((-1./(currCExt*currCExt))*(
                                       L_vel[0] + 0.5*(L_5 + L_1)));
              D_TERM(
                extState2[MD_IX(i, cVel)] = prevVelInt[0]
                + m_prevDt[a_level]*(
                  (-1./(2.*currRhoExt*currCExt))*(L_5 - L_1));,
                extState2[MD_IX(i, cVel+1)] = prevVelInt[1]
                + m_prevDt[a_level]*(-L_vel[1]);,
                extState2[MD_IX(i, cVel+2)] = prevVelInt[2]
                + m_prevDt[a_level]*(-L_vel[2]););
              extState2[MD_IX(i, cPres)] = prevPressExt
                + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
            }

          // Set all variables to characteristic outflow state
          //extState2.shift(a_dir, 1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velFarfield(D_DECL(extState2[MD_IX(i, cVel)],
                                          extState2[MD_IX(i, cVel+1)],
                                          extState2[MD_IX(i, cVel+2)]));
              Real rhoFarfield = extState2[MD_IX(i, cRho)];
              Real pressFarfield = extState2[MD_IX(i, cPres)];
              Real tempFarfield = pressFarfield/(Rgas*rhoFarfield);

              // First set all components to the state values -- takes care of
              // anything not currently adjusted by characteristic formulation
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_WfaceAvgExterior[MD_IX(i, comp)] = state(comp);
                }
              // Set all values currently adjusted by characteristics
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoFarfield;
              D_TERM(
                a_WfaceAvgExterior[MD_IX(i, cVel)] = velFarfield[0];,
                a_WfaceAvgExterior[MD_IX(i, cVel+1)] = velFarfield[1];,
                a_WfaceAvgExterior[MD_IX(i, cVel+2)] = velFarfield[2];);
              a_WfaceAvgExterior[MD_IX(i, cPres)] = pressFarfield;
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempFarfield;
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                      MD_BOXLOOP(exteriorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WfaceAvgExterior[MD_IX(i, comp)];
                        }
                    }
                }
            }
        }
      else if (m_outflowMethod == 2)
        {
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Compute wave-amplitudes for interior state at current time
          // using time-data
          const int numWaves = 1 + SpaceDim + 1; // excluding other waves here
          FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                                            a_WcellAvgFab[MD_IX(i, cVel+1)],
                                            a_WcellAvgFab[MD_IX(i, cVel+2)]));
              const Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                intWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          //**NOTE: The exterior wave amplitudes are assumed to be zero

          // Advance previous-time exterior state using current time exterior
          // wave-amplitudes and previous-time interior wave-amplitudes
          FABSTACKTEMP(extState2, interiorCellBox, numPrim);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              const Real currPress = 0.5*(state.pressure())
                + 0.5*a_WcellAvgFab[MD_IX(i, cPres)];
              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);
              // Primitive state at previous time
              // const Real prevRho = currRho;
              // const RealVect prevVel(D_DECL(state.velocity()[0],
              //                               state.velocity()[1],
              //                               state.velocity()[2]));
              // const Real prevPress = currPress;
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = 0.5*(state.pressure())
                + 0.5*a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Waves -- L1 comes from exterior, all others come from interior
              const Real L_1 = 0.;
              const RealVect L_vel(
                D_DECL(intWaveAmpCurrTime[MD_IX(i, cVel)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+1)],
                       intWaveAmpCurrTime[MD_IX(i, cVel+2)]));
              const Real L_5 = intWaveAmpCurrTime[MD_IX(i, cPres)];
              // Advance the primitive state in time
              extState2[MD_IX(i, cRho)] = prevRho
                + m_prevDt[a_level]*((-1./(currC*currC))*(
                                       L_vel[0] + 0.5*(L_5 + L_1)));
              D_TERM(
                extState2[MD_IX(i, cVel)] = prevVel[0]
                + m_prevDt[a_level]*((-1./(2.*currRho*currC))*(L_5 - L_1));,
                extState2[MD_IX(i, cVel+1)] = prevVel[1]
                + m_prevDt[a_level]*(-L_vel[1]);,
                extState2[MD_IX(i, cVel+2)] = prevVel[2]
                + m_prevDt[a_level]*(-L_vel[2]););
              extState2[MD_IX(i, cPres)] = prevPress
                + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
            }

          // Set all variables to characteristic outflow state
          extState2.shift(a_dir, 1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velFarfield(D_DECL(extState2[MD_IX(i, cVel)],
                                          extState2[MD_IX(i, cVel+1)],
                                          extState2[MD_IX(i, cVel+2)]));
              Real rhoFarfield = extState2[MD_IX(i, cRho)];
              Real pressFarfield = extState2[MD_IX(i, cPres)];
              Real tempFarfield = pressFarfield/(Rgas*rhoFarfield);

              // First set all components to the state values -- takes care of
              // anything not currently adjusted by characteristic formulation
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_WfaceAvgExterior[MD_IX(i, comp)] = state(comp);
                }
              // Set all values currently adjusted by characteristics
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoFarfield;
              D_TERM(
                a_WfaceAvgExterior[MD_IX(i, cVel)] = velFarfield[0];,
                a_WfaceAvgExterior[MD_IX(i, cVel+1)] = velFarfield[1];,
                a_WfaceAvgExterior[MD_IX(i, cVel+2)] = velFarfield[2];);
              a_WfaceAvgExterior[MD_IX(i, cPres)] = pressFarfield;
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempFarfield;
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
        }
      else if (m_outflowMethod == 3)
        {
          // Set all variables except pressure to interior state
          // Set pressure to reference state
          MD_BOXLOOP(exteriorCellBox, i)
            {
              Real rho_interior = a_WcellAvgFab[MD_IX(i,cRho)];
              a_WfaceAvgExterior[MD_IX(i,cRho)] = rho_interior;
              D_TERM(
                a_WfaceAvgExterior[MD_IX(i, cVel)] =
                a_WcellAvgFab[MD_IX(i,cVel)];,
                a_WfaceAvgExterior[MD_IX(i, cVel+1)] =
                a_WcellAvgFab[MD_IX(i,cVel+1)];,
                a_WfaceAvgExterior[MD_IX(i, cVel+2)] =
                a_WcellAvgFab[MD_IX(i,cVel+2)];);
              a_WfaceAvgExterior[MD_IX(i, cPres)] = state.pressure();
            }
        }
      else if (m_outflowMethod == 4)
        {
          const int MD_ID(o, a_dir);
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_OFFSETIX(i,-,o,comp)];
                    }
                }
            }
          MD_BOXLOOP(exteriorCellBox, i)
            {
              const Real cp = CRDparam::g_CRDPhysics->cp();
              // Interior estimate of velocity
              RealVect intVel;
              for (const int idxVel : EachDir)
                {
                  intVel[idxVel] =
                    a_WcellAvgFab[MD_OFFSETIX(i,-,o,cVel+idxVel)];
                }
              // Previous interior velocity
              RealVect intPrevVel;
              for (const int idxVel : EachDir)
                {
                  intPrevVel[idxVel] =
                    a_bndryCellFab[MD_OFFSETIX(i,-,o,cVel+numPrim+idxVel)];
                }
              const RealVect intVel1(intVel - intPrevVel);
              const Real intVel1Sq = stc::dot(intVel1, intVel1);
              Real intT;
              // Temperature state based on characteristics
              if (intVel[0]*Side::sign(a_bcIdx.m_side) < 0.)
                {
                  intT = state.temperature() - intVel1Sq/(2*cp); // Inflow
                }
              else
                {
                  intT = a_WcellAvgFab[MD_OFFSETIX(i,-,o,cTemp)]; // Outflow
                }
              const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
              // Whether inflow or outflow, this is the assigned pressure.
              // Here, assume a frozen mixture with gamma from the interior.
              const Real p =
                state.pressure()*std::pow(intc1, -gamma/(gamma - 1.0));
              a_WfaceAvgExterior[MD_IX(i, cPres)] = p;
              // Temperature and density assigned only if inflow
              if (intVel[0]*Side::sign(a_bcIdx.m_side) < 0.)  // Inflow
                {
                  // const Real T = state.temperature()/intc1;
                  // Assume density is alread set
                  // a_Wface[MD_IX(i, cRho)]  = p/(Rgas*T);
                  // Assume temperature uses existing density and pressure
                  a_WfaceAvgExterior[MD_IX(i, cTemp)] =
                    a_WfaceAvgExterior[MD_IX(i, cPres)]/(
                      Rgas*a_WcellAvgFab[MD_OFFSETIX(i,-,o,cRho)]);
                }
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCIn:
    {
      if (m_inflowMethod == 0)
        {
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
          scaledFab.shift(a_dir, -1);
          // DEBUG ONLY
          FArrayBox& relaxedFab =
            (*(m_multiBoxInletPlaneRelaxed[a_level]))[ditIPCurrent];
          relaxedFab.shift(a_dir, -1);
          FArrayBox& instFab =
            (*(m_multiBoxInletPlaneInst[a_level]))[ditIPCurrent];
          instFab.shift(a_dir, -1);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              instFab[MD_IX(i, cRho)] = a_WfaceAvgExterior[MD_IX(i, cRho)];
              D_TERM(instFab[MD_IX(i, cVel)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel)];,
                     instFab[MD_IX(i, cVel+1)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel+1)];,
                     instFab[MD_IX(i, cVel+2)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel+2)];);
              instFab[MD_IX(i, cPres)] = a_WfaceAvgExterior[MD_IX(i, cPres)];
              instFab[MD_IX(i, cTemp)] = a_WfaceAvgExterior[MD_IX(i, cTemp)];
            }
          // END DEBUG ONLY

          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        scaledFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Compute wave-amplitudes for interior state at current time
          // using time-data
          const int numWaves = 1 + SpaceDim + 1; // exclude other waves here
          FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
          const int MD_ID(o, a_dir);
          Real alpha = 0.8;
          Real beta = 0.2;
          Real omega = 0.0;
          MD_BOXLOOP(interiorCellBox, i)
            {
              Real yLoc = m_yLoc[a_level][MD_GETIV(i)[1]];
              beta = 0.5*(0.0 + 0.8) + 0.5*(0.0 - 0.8)*std::tanh(
               (2./(0.25*m_delta))*(yLoc - 1.25*m_delta));
              beta = 0.;
              alpha = 0.;
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              RealVect currVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                                      a_WcellAvgFab[MD_IX(i, cVel+1)],
                                      a_WcellAvgFab[MD_IX(i, cVel+2)]));
              currVel[0] = (1. - omega)*currVel[0]
                + omega*scaledFab[MD_OFFSETIX(i,-,o,cVel)];
              Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
              // Relaxation pressure at current time -- exterior state trying
              // to come into the domain
              const Real relaxPress = scaledFab[MD_OFFSETIX(i,-,o,cPres)];
              currPress = (1. - alpha)*currPress + alpha*relaxPress;
              // pout() << "relaxPress = " << relaxPress << "; currPress = "
              //        << currPress << std::endl;
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                intWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Compute wave-amplitudes for exterior state at current time
          // using time-data
          FABSTACKTEMP(extWaveAmpCurrTime, exteriorCellBox, numWaves);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              Real yLoc = m_yLoc[a_level][MD_GETIV(i)[1]];
              beta = 0.5*(0.0 + 0.8) + 0.5*(0.0 - 0.8)*std::tanh(
               (2./(0.25*m_delta))*(yLoc - 1.25*m_delta));
              beta = 0.;
              // Primitive state at current time
              const Real currRho = scaledFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(scaledFab[MD_IX(i, cVel)],
                                            scaledFab[MD_IX(i, cVel+1)],
                                            scaledFab[MD_IX(i, cVel+2)]));
              Real currPress = scaledFab[MD_IX(i, cPres)];
              // Relaxation pressure at current time -- interior state trying
              // to come out of the domain
              const Real relaxPress = a_WcellAvgFab[MD_OFFSETIX(i,+,o,cPres)];
              currPress = (1. - beta)*currPress + beta*relaxPress;
              // Update scaledFab pressure
              scaledFab[MD_IX(i, cPres)] = currPress;
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of exterior state
              // (1) left-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                extWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Advance previous-time exterior state using current time exterior
          // wave-amplitudes and previous-time interior wave-amplitudes
          intWaveAmpCurrTime.shift(a_dir, -1);
          FABSTACKTEMP(extState2, exteriorCellBox, numPrim);
          Real alpha_1 = 0.9;
          MD_BOXLOOP(exteriorCellBox, i)
            {
              Real yLoc = m_yLoc[a_level][MD_GETIV(i)[1]];
              alpha_1 = 0.5*(0.0 + 1.0) + 0.5*(0.0 - 1.0)*std::tanh(
               (2./(0.25*m_delta))*(yLoc - 1.25*m_delta));
              // Primitive state at current time
              const Real currRho = scaledFab[MD_IX(i, cRho)];
              const Real currPress = scaledFab[MD_IX(i, cPres)];
              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Waves -- L1 comes from interior, all others come from exterior
              const Real L_1 = alpha_1*intWaveAmpCurrTime[MD_IX(i, cRho)]
                + (1. - alpha_1)*extWaveAmpCurrTime[MD_IX(i, cRho)];
              const RealVect L_vel(
                D_DECL(extWaveAmpCurrTime[MD_IX(i, cVel)],
                       extWaveAmpCurrTime[MD_IX(i, cVel+1)],
                       extWaveAmpCurrTime[MD_IX(i, cVel+2)]));
              const Real L_5 = extWaveAmpCurrTime[MD_IX(i, cPres)];
              // Advance the primitive state in time
              extState2[MD_IX(i, cRho)] = prevRho
                + m_prevDt[a_level]*((-1./(currC*currC))*(
                                       L_vel[0] + 0.5*(L_5 + L_1)));
              D_TERM(
                extState2[MD_IX(i, cVel)] = prevVel[0]
                + m_prevDt[a_level]*((-1./(2.*currRho*currC))*(L_5 - L_1));,
                extState2[MD_IX(i, cVel+1)] = prevVel[1]
                + m_prevDt[a_level]*(-L_vel[1]);,
                extState2[MD_IX(i, cVel+2)] = prevVel[2]
                + m_prevDt[a_level]*(-L_vel[2]););
              extState2[MD_IX(i, cPres)] = prevPress
                + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
            }

          // Apply the state to the boundary face
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velInlet(D_DECL(extState2[MD_IX(i, cVel)],
                                       extState2[MD_IX(i, cVel+1)],
                                       extState2[MD_IX(i, cVel+2)]));
              Real rhoInlet = extState2[MD_IX(i, cRho)];
              Real tempInlet = extState2[MD_IX(i, cPres)]/(Rgas*rhoInlet);
              for (const int idxVel : EachDir)
                {
                  a_WfaceAvgExterior[MD_IX(i, cVel + idxVel)] =
                    velInlet[idxVel];
                }
              a_WfaceAvgExterior[MD_IX(i, cPres)] = extState2[MD_IX(i, cPres)];
              // Inflow cannot become outflow, so always set rho and temperature
              for (int j = 0; j != numSpecies; ++j)
                {
                  a_WfaceAvgExterior[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
                }
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempInlet;
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoInlet;

              // DEBUG ONLY
              relaxedFab[MD_IX(i, cRho)] = a_WfaceAvgExterior[MD_IX(i, cRho)];
              D_TERM(relaxedFab[MD_IX(i, cVel)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel)];,
                     relaxedFab[MD_IX(i, cVel+1)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel+1)];,
                     relaxedFab[MD_IX(i, cVel+2)] =
                     a_WfaceAvgExterior[MD_IX(i, cVel+2)];);
              relaxedFab[MD_IX(i, cTemp)] = a_WfaceAvgExterior[MD_IX(i, cTemp)];
              relaxedFab[MD_IX(i, cPres)] = a_WfaceAvgExterior[MD_IX(i, cPres)];
              // END DEBUG ONLY
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(exteriorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            scaledFab[MD_IX(i, comp)];
                        }
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
          scaledFab.shift(a_dir, 1);
          // DEBUG ONLY
          // First try a Riemann solution between instFab and relaxedFab
          instFab.shiftHalf(a_dir, 1);
          relaxedFab.shiftHalf(a_dir, 1);
          Box faceBox = exteriorCellBox;
          faceBox.shiftHalf(a_dir, 1);
          FABSTACKTEMP(tempRiemannFab, faceBox, numPrim);
          FABSTACKTEMP(tempRelaxedFab, faceBox, numPrim);
          tempRelaxedFab.copy(relaxedFab);
          FABSTACKTEMP(tempInstFab, faceBox, numPrim);
          tempInstFab.copy(instFab);
          CRDparam::g_CRDPhysics->riemannBC(tempRiemannFab,
                                            tempRelaxedFab,
                                            tempInstFab,
                                            a_unitNormalBasisFab,
                                            0,
                                            Side::Lo,
                                            faceBox);
          instFab.copy(tempRiemannFab);
          instFab.shiftHalf(a_dir, -1);
          relaxedFab.shiftHalf(a_dir, -1);
          // Now shift back to the interior
          instFab.shift(a_dir, 1);
          relaxedFab.shift(a_dir, 1);
        }
      else if (m_inflowMethod == 1)
        {
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
          scaledFab.shift(a_dir, -1);
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        scaledFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Compute wave-amplitudes for interior state at current time
          // using time-data
          const int numWaves = 1 + SpaceDim + 1; // exclude other waves here
          FABSTACKTEMP(intWaveAmpCurrTime, interiorCellBox, numWaves);
          MD_BOXLOOP(interiorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = a_WcellAvgFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                                            a_WcellAvgFab[MD_IX(i, cVel+1)],
                                            a_WcellAvgFab[MD_IX(i, cVel+2)]));
              const Real currPress = a_WcellAvgFab[MD_IX(i, cPres)];
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                intWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                intWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              intWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Compute wave-amplitudes for exterior state at current time
          // using time-data
          FABSTACKTEMP(extWaveAmpCurrTime, exteriorCellBox, numWaves);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = scaledFab[MD_IX(i, cRho)];
              const RealVect currVel(D_DECL(scaledFab[MD_IX(i, cVel)],
                                            scaledFab[MD_IX(i, cVel+1)],
                                            scaledFab[MD_IX(i, cVel+2)]));
              const Real currPress = scaledFab[MD_IX(i, cPres)];
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Time derivatives of the primitive state
              const Real rhoTimeDeriv = f_0*(currRho - prevRho);
              const RealVect velTimeDeriv = f_0*(currVel - prevVel);
              const Real pressTimeDeriv = f_0*(currPress - prevPress);

              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);

              // Wave-amplitudes of interior state
              // (1) left-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cRho)] =
                currRho*currC*velTimeDeriv[0] - pressTimeDeriv;
              D_TERM(
                // (2) entropy wave in x-direction
                extWaveAmpCurrTime[MD_IX(i, cVel)] =
                pressTimeDeriv - currC*currC*rhoTimeDeriv;,
                // (3) entropy wave in y-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+1)] = -velTimeDeriv[1];,
                // (4) entropy wave in z-direction
                extWaveAmpCurrTime[MD_IX(i, cVel+2)] = -velTimeDeriv[2];);
              // (5) right-moving acoustic
              extWaveAmpCurrTime[MD_IX(i, cPres)] =
                -pressTimeDeriv - currRho*currC*velTimeDeriv[0];
            }

          // Advance previous-time exterior state using current time exterior
          // wave-amplitudes and previous-time interior wave-amplitudes
          intWaveAmpCurrTime.shift(a_dir, -1);
          FABSTACKTEMP(extState2, exteriorCellBox, numPrim);
          MD_BOXLOOP(exteriorCellBox, i)
            {
              // Primitive state at current time
              const Real currRho = scaledFab[MD_IX(i, cRho)];
              const Real currPress = scaledFab[MD_IX(i, cPres)];
              // Current speed of sound
              const Real currC = std::sqrt(gamma*currPress/currRho);
              // Primitive state at previous time
              const Real prevRho = a_bndryCellFab[MD_IX(i, cRho+numPrim)];
              const RealVect prevVel(
                D_DECL(a_bndryCellFab[MD_IX(i, cVel+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+1+numPrim)],
                       a_bndryCellFab[MD_IX(i, cVel+2+numPrim)]));
              const Real prevPress = a_bndryCellFab[MD_IX(i, cPres+numPrim)];
              // Waves -- L1 comes from interior, all others come from exterior
              const Real L_1 = intWaveAmpCurrTime[MD_IX(i, cRho)];
              const RealVect L_vel(
                D_DECL(extWaveAmpCurrTime[MD_IX(i, cVel)],
                       extWaveAmpCurrTime[MD_IX(i, cVel+1)],
                       extWaveAmpCurrTime[MD_IX(i, cVel+2)]));
              const Real L_5 = extWaveAmpCurrTime[MD_IX(i, cPres)];
              // Advance the primitive state in time
              extState2[MD_IX(i, cRho)] = prevRho
                + m_prevDt[a_level]*((-1./(currC*currC))*(
                                       L_vel[0] + 0.5*(L_5 + L_1)));
              D_TERM(
                extState2[MD_IX(i, cVel)] = prevVel[0]
                + m_prevDt[a_level]*((-1./(2.*currRho*currC))*(L_5 - L_1));,
                extState2[MD_IX(i, cVel+1)] = prevVel[1]
                + m_prevDt[a_level]*(-L_vel[1]);,
                extState2[MD_IX(i, cVel+2)] = prevVel[2]
                + m_prevDt[a_level]*(-L_vel[2]););
              extState2[MD_IX(i, cPres)] = prevPress
                + m_prevDt[a_level]*(-0.5*(L_5 + L_1));
            }

          // Apply the state to the boundary face
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velInlet(D_DECL(scaledFab[MD_IX(i, cVel)],
                                       scaledFab[MD_IX(i, cVel+1)],
                                       scaledFab[MD_IX(i, cVel+2)]));
              Real rhoInlet = extState2[MD_IX(i, cRho)];
              Real tempInlet = extState2[MD_IX(i, cPres)]/(Rgas*rhoInlet);
              for (const int idxVel : EachDir)
                {
                  a_WfaceAvgExterior[MD_IX(i, cVel + idxVel)] =
                    velInlet[idxVel];
                }
              a_WfaceAvgExterior[MD_IX(i, cPres)] = extState2[MD_IX(i, cPres)];
              // Inflow cannot become outflow, so always set rho and temperature
              for (int j = 0; j != numSpecies; ++j)
                {
                  a_WfaceAvgExterior[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
                }
              a_WfaceAvgExterior[MD_IX(i, cTemp)] = tempInlet;
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoInlet;
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(exteriorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            scaledFab[MD_IX(i, comp)];
                        }
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
          scaledFab.shift(a_dir, 1);
        }
      else if (m_inflowMethod == 2)
        {
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
          scaledFab.shift(a_dir, -1);
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        scaledFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Apply the state to the boundary face
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velInlet(D_DECL(scaledFab[MD_IX(i, cVel)],
                                       scaledFab[MD_IX(i, cVel+1)],
                                       scaledFab[MD_IX(i, cVel+2)]));
              Real rhoInlet = scaledFab[MD_IX(i, cRho)];
              for (const int idxVel : EachDir)
                {
                  a_WfaceAvgExterior[MD_IX(i, cVel + idxVel)] =
                    velInlet[idxVel];
                }
              a_WfaceAvgExterior[MD_IX(i, cPres)] =
                a_WfaceAvgExterior[MD_IX(i, cPres)];
              // Inflow cannot become outflow, so always set rho and temperature
              a_WfaceAvgExterior[MD_IX(i, cTemp)] =
                a_WfaceAvgExterior[MD_IX(i, cPres)]/(Rgas*rhoInlet);
              a_WfaceAvgExterior[MD_IX(i, cRho)] = rhoInlet;
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(exteriorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            scaledFab[MD_IX(i, comp)];
                        }
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
          scaledFab.shift(a_dir, 1);
        }
      else if (m_inflowMethod == 3)
        {
          const int MD_ID(o, a_dir);
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
          scaledFab.shift(a_dir, -1);
          // If this is the first stage of the first step, we need to fill the
          // previous time-values with the current time value
          int numStages = 4;
          if (CRDparam::g_additiveRK)
            {
              numStages = 6;
            }
          int currStage = m_stageCounter[a_level] % numStages;
          if (m_stageCounter[a_level] == 0)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  MD_BOXLOOP(exteriorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        scaledFab[MD_IX(i, comp)];
                    }
                  MD_BOXLOOP(interiorCellBox, i)
                    {
                      a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                        a_WcellAvgFab[MD_IX(i, comp)];
                    }
                }
            }
          // Apply the state to the boundary face
          MD_BOXLOOP(exteriorCellBox, i)
            {
              RealVect velInlet(D_DECL(scaledFab[MD_IX(i, cVel)],
                                       scaledFab[MD_IX(i, cVel+1)],
                                       scaledFab[MD_IX(i, cVel+2)]));
              Real TInlet = scaledFab[MD_IX(i, cTemp)];

              const Real cp = CRDparam::g_CRDPhysics->cp();
              // Current exterior velocity
              RealVect currExtVel;
              for (const int idxVel : EachDir)
                {
                  currExtVel[idxVel] = scaledFab[MD_IX(i, cVel+idxVel)];
                }
              // Previous exterior velocity
              RealVect prevExtVel;
              for (const int idxVel : EachDir)
                {
                  prevExtVel[idxVel] =
                    a_bndryCellFab[MD_IX(i,cVel+numPrim+idxVel)];
                }
              // Current interior velocity
              RealVect currIntVel;
              for (const int idxVel : EachDir)
                {
                  currIntVel[idxVel] =
                    a_WcellAvgFab[MD_OFFSETIX(i,+,o,cVel+idxVel)];
                }
              // Previous interior velocity
              RealVect prevIntVel;
              for (const int idxVel : EachDir)
                {
                  prevIntVel[idxVel] =
                    a_bndryCellFab[MD_OFFSETIX(i,+,o,cVel+numPrim+idxVel)];
                }
              for (const int idxVel : EachDir)
                {
                  a_WfaceAvgExterior[MD_IX(i,cVel+idxVel)] = currExtVel[idxVel];
                }
              // Pressure is determined from the interior.  First, find a
              // stagnation pressure for the given reference frame.  Assume
              // gas is frozen at upstream conditions.
              const RealVect intVel1(currIntVel - prevIntVel);
              const Real intVel1Sq = stc::dot(intVel1, intVel1);
              const Real intT = a_WfaceAvgExterior[MD_IX(i, cTemp)];
              const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
              const Real p0 = a_WfaceAvgExterior[MD_IX(i, cPres)]*
                std::pow(intc1, gamma/(gamma - 1.));
              // Find imposed pressure from stagnation pressure
              const RealVect extVel1(currExtVel - prevExtVel);
              const Real extVel1Sq = stc::dot(extVel1, extVel1);
              const Real extc1 = 1 + extVel1Sq/(2*cp*TInlet);
              a_WfaceAvgExterior[MD_IX(i, cPres)] =
                p0*std::pow(extc1, -gamma/(gamma - 1.0));
              // Inflow cannot become outflow, so we always set density and
              // temperature
              a_WfaceAvgExterior[MD_IX(i, cRho)] = scaledFab[MD_IX(i, cRho)];
              a_WfaceAvgExterior[MD_IX(i, cTemp)] =
                a_WfaceAvgExterior[MD_IX(i, cPres)]/(
                  scaledFab[MD_IX(i, cRho)]*Rgas);
            }
          // Move the current time-values to the previous times
          if (m_timeCounter[a_level])
            {
              if (currStage == (numStages - 1))
                {
                  for (int comp = 0; comp != numPrim; ++comp)
                    {
                      MD_BOXLOOP(exteriorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            scaledFab[MD_IX(i, comp)];
                        }
                      MD_BOXLOOP(interiorCellBox, i)
                        {
                          a_bndryCellFab[MD_IX(i, comp+numPrim)] =
                            a_WcellAvgFab[MD_IX(i, comp)];
                        }
                    }
                }
            }
          scaledFab.shift(a_dir, 1);
        }
      break;
    }
    default:
      CH_assert(false);
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
CNSIBCRecirculatingInletTFP::readBCInfo()
{
  ParmParse ppIBC("ibc");

  ppIBC.query("inlet_boundary_layer_thickness", m_delta);
  ppIBC.query("recycling_plane_location", m_samplePlaneLoc);
  ppIBC.query("y_plus_max_guess", m_yPlusMaxGuess);
  ppIBC.query("vel_perturb_magnitude", m_velPerturb);
  std::vector<Real> perturbFreq(SpaceDim);
  ppIBC.queryarr("perturb_freq", perturbFreq, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_perturbFreq.dataPtr(),
                                           &perturbFreq.front());

  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf = state.velocity()[0]; // freestream velocity
  m_lambda = m_delta/u_inf;
  ppIBC.query("recycle_mean_yVel", m_useAvgYVel);
  ppIBC.query("inflow_method", m_inflowMethod);
  ppIBC.query("outflow_method", m_outflowMethod);

  // We need to know how many levels on which we could possibly require
  // sampling-plane recycling
  ParmParse ppAMR("amr");

  int maxLevel = 0;
  ppAMR.query("max_level", maxLevel);
  m_numLevels = maxLevel + 1;

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
CNSIBCRecirculatingInletTFP::spaceTimeAverageVariables(
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

  Real currTime = m_timeAfterInit[a_level];

  // Set up time-averaging window sizes
  Real T = m_phi2*m_lambda + currTime - m_lambda*(m_psi1 + m_psi2);
  if (currTime < (m_lambda*m_psi1))
    {
      T = m_phi1*m_lambda;
    }
  else if (currTime < (m_lambda*(m_psi1 + m_psi2)))
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
//  Spatially average data if necessary
/** \param[out] a_avgData
 *                      Spatially-averaged data
 *  \param[in]  a_avgData
 *                      Data to spatially-average
 *  \param[in]  a_stage Current stage of the time-marching method
 *  \param[in]  a_t     Current time
 *//*-----------------------------------------------------------------*/

void
CNSIBCRecirculatingInletTFP::spatiallyAverageData(
  LevelData<FArrayBox>& a_avgData,
  const int             a_stage,
  const Real            a_t) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::spatiallyAverageData");
  if (a_stage != 0) { return; }
  // // DEBUG ONLY
  // // Plot m_singleBoxInletPlane if we are at the plot frequency
  // int plotFreq = 10;
  // int plotTest = m_timeCounter[0] % plotFreq;
  // const int numWComp = CRDparam::g_CRDPhysics->numPrimitive();
  // if (!plotTest)
  //   {
  //     // Copy relaxed state to singleBoxInletPlaneRelaxed
  //     const Interval copyIntv(0, numWComp-1);
  //     m_multiBoxInletPlaneRelaxed[0]->copyTo(
  //       copyIntv, (*(m_singleBoxInletPlaneRelaxed[0])),
  //       copyIntv, (*(m_multiBoxInletPlaneRelaxedCopier[0])));
  //     // Copy Riemann solved state
  //     m_multiBoxInletPlaneInst[0]->copyTo(
  //       copyIntv, (*(m_singleBoxInletPlaneInst[0])),
  //       copyIntv, (*(m_multiBoxInletPlaneInstCopier[0])));
  //   }
  // if (!plotTest && !(procID()))
  //   {
  //     // Plot the instantaneous values
  //     const int currTimeStep = m_timeCounter[0];
  //     DataIterator ditIP =
  //       m_singleBoxInletPlaneRelaxed[0]->getBoxes().dataIterator();
  //     FArrayBox& inletPlaneW =
  //       (*(m_singleBoxInletPlaneRelaxed[0]))[ditIP];
  //     Box ipBox = m_singleBoxInletPlaneRelaxed[0]->getBoxes()[ditIP];
  //     FABSTACKTEMP(inletTempFab, ipBox, numWComp);
  //     inletTempFab.copy(inletPlaneW);
  //     Vector<string> names(numWComp);
  //     char compNameString[64];
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "inst-%s-target",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     char fileName[64];
  //     sprintf(fileName, "inletRelaxedValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletTempFab, fileName, names);

  //     // Plot the spanwise averaged values
  //     FABSTACKTEMP(inletAvgFab, ipBox, numWComp);
  //     inletAvgFab.setVal(0.);
  //     const IntVect nCells = ipBox.size();
  //     int nCellsZ = 1;
  //     D_TERM(,,nCellsZ = nCells[2];);
  //     const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  //     const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  //     const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  //     // Compute the spanwise-averaged primitive state -- once per time-step
  //     for (int yLoc = 0; yLoc != nCells[1]; ++yLoc)
  //       {
  //         // set up the one-dimensional averaging box at the current y-index
  //         IntVect loIntVect = ipBox.smallEnd();
  //         IntVect hiIntVect = ipBox.bigEnd();
  //         loIntVect[1] = yLoc;
  //         hiIntVect[1] = yLoc;
  //         Box averageBox(loIntVect, hiIntVect);

  //         // average the variables over the 1D box
  //         Real rhoAvg = 0.;
  //         D_TERM(Real uVelAvg = 0.;, Real vVelAvg = 0.;, Real wVelAvg = 0.;);
  //         Real TAvg = 0.;
  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             rhoAvg += inletTempFab[MD_IX(i, rhoIndx)];
  //             D_TERM(
  //               uVelAvg += inletTempFab[MD_IX(i, WvelIndx)];,
  //               vVelAvg += inletTempFab[MD_IX(i, WvelIndx+1)];,
  //               wVelAvg += inletTempFab[MD_IX(i, WvelIndx+2)];);
  //             TAvg += inletTempFab[MD_IX(i, tempIndx)];
  //           }

  //         // divide by the number of cells to finish off the spatial averaging
  //         rhoAvg /= nCellsZ;
  //         D_TERM(uVelAvg /= nCellsZ;, vVelAvg /= nCellsZ;, wVelAvg /= nCellsZ;);
  //         TAvg /= nCellsZ;

  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             inletAvgFab[MD_IX(i, rhoIndx)] = rhoAvg;
  //             D_TERM(inletAvgFab[MD_IX(i, WvelIndx)] = uVelAvg;,
  //                    inletAvgFab[MD_IX(i, WvelIndx+1)] = vVelAvg;,
  //                    inletAvgFab[MD_IX(i, WvelIndx+2)] = wVelAvg;);
  //             inletAvgFab[MD_IX(i, tempIndx)] = TAvg;
  //           }
  //       }
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "spanwiseAvg-%s-target",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     sprintf(fileName, "inletSpanwiseAvgRelaxedValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletAvgFab, fileName, names);

  //     // Plot the instantaneous Riemann values
  //     FArrayBox& inletPlaneWRiemann =
  //       (*(m_singleBoxInletPlaneInst[0]))[ditIP];
  //     FABSTACKTEMP(inletTempRiemannFab, ipBox, numWComp);
  //     inletTempRiemannFab.copy(inletPlaneWRiemann);
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "inst-%s-Riemann",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     sprintf(fileName, "inletRiemannValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletTempRiemannFab, fileName, names);

  //     // Plot the spanwise averaged values
  //     FABSTACKTEMP(inletAvgRiemannFab, ipBox, numWComp);
  //     inletAvgRiemannFab.setVal(0.);

  //     // Compute the spanwise-averaged primitive state -- once per time-step
  //     for (int yLoc = 0; yLoc != nCells[1]; ++yLoc)
  //       {
  //         // set up the one-dimensional averaging box at the current y-index
  //         IntVect loIntVect = ipBox.smallEnd();
  //         IntVect hiIntVect = ipBox.bigEnd();
  //         loIntVect[1] = yLoc;
  //         hiIntVect[1] = yLoc;
  //         Box averageBox(loIntVect, hiIntVect);

  //         // average the variables over the 1D box
  //         Real rhoAvg = 0.;
  //         D_TERM(Real uVelAvg = 0.;, Real vVelAvg = 0.;, Real wVelAvg = 0.;);
  //         Real TAvg = 0.;
  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             rhoAvg += inletTempRiemannFab[MD_IX(i, rhoIndx)];
  //             D_TERM(
  //               uVelAvg += inletTempRiemannFab[MD_IX(i, WvelIndx)];,
  //               vVelAvg += inletTempRiemannFab[MD_IX(i, WvelIndx+1)];,
  //               wVelAvg += inletTempRiemannFab[MD_IX(i, WvelIndx+2)];);
  //             TAvg += inletTempRiemannFab[MD_IX(i, tempIndx)];
  //           }

  //         // divide by the number of cells to finish off the spatial averaging
  //         rhoAvg /= nCellsZ;
  //         D_TERM(uVelAvg /= nCellsZ;, vVelAvg /= nCellsZ;, wVelAvg /= nCellsZ;);
  //         TAvg /= nCellsZ;

  //         MD_BOXLOOP(averageBox, i)
  //           {
  //             inletAvgRiemannFab[MD_IX(i, rhoIndx)] = rhoAvg;
  //             D_TERM(inletAvgRiemannFab[MD_IX(i, WvelIndx)] = uVelAvg;,
  //                    inletAvgRiemannFab[MD_IX(i, WvelIndx+1)] = vVelAvg;,
  //                    inletAvgRiemannFab[MD_IX(i, WvelIndx+2)] = wVelAvg;);
  //             inletAvgRiemannFab[MD_IX(i, tempIndx)] = TAvg;
  //           }
  //       }
  //     for (int compName = 0; compName != numWComp; ++compName)
  //       {
  //         sprintf(compNameString, "spanwiseAvg-%s-Riemann",
  //                 CRDparam::g_CRDPhysics->primStateName(compName));
  //         names[compName] = compNameString;
  //       }
  //     sprintf(fileName, "inletSpanwiseAvgRiemannValues.%06d.%dd.hdf5",
  //             currTimeStep, SpaceDim);
  //     writeFABname(&inletAvgRiemannFab, fileName, names);
  //   }
  // // END DEBUG ONLY
  if (a_t < CRDparam::g_startTimeAvgTime) { return; }
#if (CH_SPACEDIM == 3)
  // This is going to be quite a hack and it may be terribly slow
  // We can't get the domainBox the normal way, so we'll construct it manually
  Vector<Box> tempBoxes = a_avgData.getBoxes().boxArray();
  Box domainBox = a_avgData.getBoxes().physDomain().domainBox();
  int yHiTemp = domainBox.smallEnd()[1];
  domainBox.setBig(1, yHiTemp);
  for (int boxComp = 0; boxComp != tempBoxes.size(); ++boxComp)
    {
      // Get the boxes at the sampling plane location
      domainBox.minBox(tempBoxes[boxComp]);
    }

  // Make a one-cell-thick version of this
  Box domainBox2D = domainBox;
  // Make this domain box to be one cell thick in the z-direction
  domainBox2D.setBig(2, domainBox.smallEnd(2));
  // Set up a fab for this box -- we need to use this for storing averages
  FABSTACKTEMP(domainFab, domainBox2D, a_avgData.nComp());
  // Set the fab to zero to begin with
  domainFab.setVal(0.);
  // Now, we only fill the fab where the local process has data
  for (DataIterator dit = a_avgData.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_avgData.disjointBoxLayout()[dit];

      // Take the first layer of cells on the low z-side of disjointBox
      Box avgBox = disjointBox;
      avgBox.setBig(2, avgBox.smallEnd(2));
      // Assume the entire dataset needs to be spatially averaged
      FArrayBox& avgFab = a_avgData[dit()];
      // Average the data in the z-direction
      for (int comp = 0; comp != a_avgData.nComp(); ++comp)
        {
          MD_BOXLOOP(avgBox, i)
            {
              // Create the z-direction box to average over
              Box zAvgBox(MD_GETIV(i),
                          MD_GETIV(i) + BASISV(2)*(disjointBox.size(2) - 1));
              // Sum up the values
              Real avgVal = 0.;
              MD_BOXLOOP(zAvgBox, j)
                {
                  avgVal += avgFab[MD_IX(j,comp)];
                }
              // Assign the average to domainFab
              IntVect domFabIV = MD_GETIV(i);
              domFabIV[2] = domainBox.smallEnd(2);
              domainFab[MD_IV(domFabIV, comp)] += avgVal;
            }
        }
    }
  // Now for the fun part -- we loop over the domainFab and MPI sum everything
  // Basically, we want to do a linear-in, linear-out with a buffer
  // This should reduce the number of MPI_Allreduce calls and the runtime
  int bufferSize = (domainBox2D.numPts())*(a_avgData.nComp());
  std::vector<Real> local(bufferSize, 0.);
  std::vector<Real> global(bufferSize, 0.);
  // Fill the buffer
  int linearIndex = 0;
  for (int comp = 0; comp != a_avgData.nComp(); ++comp)
    {
      MD_BOXLOOP(domainBox2D, i)
        {
          local[linearIndex] = domainFab[MD_IX(i, comp)];
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
  for (int comp = 0; comp != a_avgData.nComp(); ++comp)
    {
      MD_BOXLOOP(domainBox2D, i)
        {
          domainFab[MD_IX(i, comp)] = (global[linearIndex])/domainBox.size(2);
          ++linearIndex;
        }
    }

  // At last, we have the averages and we need to fill them into a_avgData
  for (DataIterator dit = a_avgData.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_avgData.disjointBoxLayout()[dit];
      // Fill a_avgData with the averaged values
      FArrayBox& avgFab = a_avgData[dit()];
      for (int comp = 0; comp != a_avgData.nComp(); ++comp)
        {
          MD_BOXLOOP(disjointBox, i)
            {
              IntVect domFabIV = MD_GETIV(i);
              domFabIV[2] = domainBox.smallEnd(2);
              avgFab[MD_IX(i,comp)] = domainFab[MD_IV(domFabIV, comp)];
            }
        }
    }
#endif
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
CNSIBCRecirculatingInletTFP::recycleFluctuations(
  FArrayBox&       a_inletPlaneW,
  const FArrayBox& a_meanInletFab,
  const FArrayBox& a_WcellPntFab,
  const FArrayBox& a_WcellInletSampledPntFab,
  const Box&       a_spBox,
  const Box&       a_ipBox,
  const int        a_level,
  const Real       a_t) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::recycleFluctuations");
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int nCellsY = m_yLoc[a_level].size();
  D_TERM(,,const int nCellsZ = m_zLoc[a_level].size(););
  const Real R = CRDparam::g_R;

  // Constant freestream state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf = state.velocity()[0]; // freestream velocity
  const Real rho_inf = state.density(); // freestream density
  const Real pressure_inf = state.pressure(); // freestream pressure
  const Real T_inf = state.temperature(); // freestream temperature

  // 6) compute delta at the recycle plane
  Real delta_recy = CRDparam::g_physicalLength[1];
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
      Real deltaX = CRDparam::g_physicalLength[0]*m_samplePlaneLoc;
      Real c_0 = deltaX/m_delta;
      Real c_1 = std::pow(0.27, 6./5.);
      Real Re_delta = rho_inf*u_inf*m_delta/(CRDparam::g_mu);
      Real c_2 = std::pow(Re_delta, -1./5.);
      Real c_3 = 1. + c_0*c_1*c_2;
      Real c_4 = std::pow(c_3, 5./6.);
      Real deltaGuess = c_4*m_delta;
      Real beta = 36*m_lambda; // 36 BL-thickness times (1 flow-through)
      Real zeta = 100*m_lambda; // 100 BL-thicness times (~3 flow-throughs)
      Real alpha = 0.5 + 0.5*std::tanh((1./beta)*(a_t - zeta));
      delta_recy = (1. - alpha)*deltaGuess + alpha*delta_recy;
      delta_recy = 0.90*deltaGuess;
      for (int i = 0; i != m_numLevels; ++i)
        {
          m_deltaRecy[i] = delta_recy;
        }
    }
  delta_recy = m_deltaRecy[0];

  // We're adding a temporally varying z-direction shifting to the fluctuations
  D_TERM(,,
         Real z_offset = 0.; // 0.5*CRDparam::g_physicalLength[2];
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

  // Do we use both thermo flucs plus a specified mean (old method), or do
  // we just recycle instantaneous density and compute temperature assuming
  // constant farfield pressure?
  bool oldInlet = true;

  // Compute the fluctuation rescaling
  FABSTACKTEMP(rescaledFluc, a_spBox, CRDparam::g_CRDPhysics->numPrimitive());
  rescaledFluc.setVal(0.);
  MD_BOXLOOP(a_spBox, i)
    {
      // Get the indices for the current wall-normal location
      const Real yIdx = i1 - a_spBox.smallEnd()[1];
      // const Real yLoc = m_yLoc[a_level][yIdx];
      // fluctuations for inlet
      // const Real y_s = (delta_recy/m_delta)*yLoc;
      const Real y_s = m_yLocMapped[a_level][yIdx];
      // const Real y_max = m_yLoc[a_level][nCellsY - 1];
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
      Real rho_mean_hi = m_rhoSTAvg[a_level][hiIndx];
      Real rho_mean_lo = m_rhoSTAvg[a_level][loIndx];
      if (!oldInlet)
        {
          rho_mean_hi = 0.;
          rho_mean_lo = 0.;
        }
      D_TERM(
        Real u_mean_hi = m_uVelSTAvg[a_level][hiIndx];
        Real u_mean_lo = m_uVelSTAvg[a_level][loIndx];,
        Real v_mean_hi = m_vVelSAvg[a_level][hiIndx];
        Real v_mean_lo = m_vVelSAvg[a_level][loIndx];,
        Real w_mean_hi = m_wVelSAvg[a_level][hiIndx];
        Real w_mean_lo = m_wVelSAvg[a_level][loIndx];);
      Real T_mean_hi = m_tempSTAvg[a_level][hiIndx];
      Real T_mean_lo = m_tempSTAvg[a_level][loIndx];
      if (!oldInlet)
        {
          T_mean_hi = 0.;
          T_mean_lo = 0.;
        }

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
      if (y_s >= 1.5*delta_recy)
        {
          rho_fluc = 0.;
          if (!oldInlet)
            {
              rho_fluc = rho_inf;
            }
          D_TERM(u_fluc = 0.;,
                 v_fluc = 0.;,
                 w_fluc = 0.;);
          T_fluc = 0.;
          if (!oldInlet)
            {
              T_fluc = T_inf;
            }
          // v_fluc = a_WcellPntFab[MD_IX(i, WvelIndx+1)]
          //   - m_vVelSTAvg[a_level][yIdx];
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
            zLocTrue += CRDparam::g_physicalLength[2];
          }
        else if (zLocTrue > CRDparam::g_physicalLength[2])
          {
            zLocTrue -= CRDparam::g_physicalLength[2];
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
          D_TERM(,,z_hi = m_zLoc[a_level][0] + CRDparam::g_physicalLength[2];
                 z_lo = m_zLoc[a_level][nCellsZ - 1];);
          if (zLocTrue < z_lo) // Make sure zLocTrue is on the high side
            {
              D_TERM(,,zLocTrue += CRDparam::g_physicalLength[2];);
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

      Real rho_mean = a_meanInletFab[MD_IX(i, rhoIndx)];
      D_TERM(
        Real u_mean = a_meanInletFab[MD_IX(i, WvelIndx)];,
        Real v_mean = a_meanInletFab[MD_IX(i, WvelIndx+1)];,
        Real w_mean = a_meanInletFab[MD_IX(i, WvelIndx+2)];);
      Real T_mean = a_meanInletFab[MD_IX(i, tempIndx)];

      CH_assert(rho_mean > 0.);
      CH_assert(T_mean > 0.);

      const RealVect finalVelState(D_DECL(u_mean + u_fluc,
                                          v_mean + v_fluc,
                                          w_mean + w_fluc));

      Real finalTempState = T_mean + T_fluc;
      // Version 1 -- uses farfield pressure and computes density from this
      // const Real finalPressState1 = state.pressure();
      // const Real finalRhoState1 = finalPressState1/(R*finalTempState);
      // Version 2 -- uses density fluctuations and computes pressure
      Real finalRhoState2 = rho_mean + rho_fluc;
      Real finalPressState2 = finalRhoState2*finalTempState*R;
      bool method1 = true;
      if (!oldInlet)
        {
          if (method1)
            {
              finalRhoState2 = rho_fluc;
              finalPressState2 = pressure_inf;
              finalTempState = pressure_inf/(R*rho_fluc);
            }
          else
            {
              finalRhoState2 = rho_fluc;
              finalTempState = T_fluc;
              finalPressState2 = rho_fluc*R*T_fluc;
            }
        }

      // Version 3 -- uses interior pressure and computes temperature from this
      // This provides the correct behaviour in the boundary layer, but gives
      // oscillations in density/temperature in the farfield and eventually
      // breaks down. So, we recognize that we're essentially saying the
      // exterior pressure needs to be relaxed to the interior pressure.
      // Previously, we were saying that none of the pressure needed to be
      // relaxed. The relaxation coefficient that we have in place merely
      // relaxes how much of the interior acoustic is mixed into the update
      // for the exterior state. However, we now realize that we can also
      // try to relax the incoming or outgoing acoustic to represent what we're
      // really doing -- the exterior state comes into the interior and becomes
      // the interior state and, in the process, the exterior pressure must be
      // relaxed to the interior pressure somehow. For now, we'll move this to
      // the characteristic update.
      // if (m_inflowMethod != 2)
      //   {
      //     Real interiorPress =
      //       a_WcellInletSampledPntFab[MD_IV(inletLoc, presIndx)];
      //     finalRhoState2 = rho_mean + rho_fluc;
      //     finalPressState2 = interiorPress;
      //     finalTempState = finalPressState2/(R*finalRhoState2);
      //   }

      D_TERM(
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx)]   = finalVelState[0];,
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx+1)] = finalVelState[1];,
        a_inletPlaneW[MD_IV(inletLoc, WvelIndx+2)] = finalVelState[2];);
      a_inletPlaneW[MD_IV(inletLoc, tempIndx)] = finalTempState;

      // // Version 1
      // a_inletPlaneW[MD_IV(inletLoc, rhoIndx)] = finalRhoState1;
      // a_inletPlaneW[MD_IV(inletLoc, presIndx)] = finalPressState1;

      // Version 2
      a_inletPlaneW[MD_IV(inletLoc, rhoIndx)] = finalRhoState2;
      a_inletPlaneW[MD_IV(inletLoc, presIndx)] = finalPressState2;

      CH_assert(finalTempState > 0.);
      CH_assert(finalRhoState2 > 0.);
      CH_assert(finalPressState2 > 0.);
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
CNSIBCRecirculatingInletTFP::muskerIC(FArrayBox&       a_meanFab,
                                      FArrayBox&       a_etaFab,
                                      const FArrayBox& a_deltaFab,
                                      const FArrayBox& a_XFab,
                                      const FArrayBox& a_XNodeFab,
                                      const Real       a_virtualWallDy,
                                      const Box&       a_box) const
{
  CH_TIME("CNSIBCRecirculatingInletTFP::muskerIC");

  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const Real R = CRDparam::g_R;
  const Real gamma = CRDparam::g_gamma;
  const Real mu = CRDparam::g_mu;
  const Real Pr_t = 0.89; // Turbulent Prandtl number (Urbin & Knight 2001)

  // Constant freestream state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real u_inf = state.velocity()[0]; // freestream velocity
  const Real T_inf = state.temperature(); // freestream temperature
  const Real p_inf = state.pressure();    // freestream pressure
  const Real rho_inf = state.density();   // freestream density
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

      const UTauMuskerFunc& f =
        UTauMuskerFunc(delta_100,nu_w,u_vd_inf,pi1);
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
          loc[1] = loYNode + subCellDy*(0.5 + cell) + a_virtualWallDy;

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
    }
}
