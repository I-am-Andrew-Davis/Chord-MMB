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
 * \file CNSIBCBurnerTest.cpp
 *
 * \brief Member functions for CNSIBCBurnerTest
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
 
//----- Internal -----//

#include "CNSIBCBurnerTest.H"
#include "CNSIBCCombustionTestF_F.H"
#include "CNSIBCTransientPoiseuilleF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDState.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCBurnerTest: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCBurnerTest::CNSIBCBurnerTest()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_flowDir(-1),
  m_jetDir(-1),
  m_jetLoVel(RealVect::Zero),
  m_jetHiVel(RealVect::Zero),
  m_hotSpotCenter(RealVect::Zero),
  m_hotSpotTemp(1500.),
  m_hotSpotStartTime(1.),
  m_hotSpotRunTime(-1.),
  m_energyIgnition(false),
  m_hotSpotIgn(0),
  m_gravForce(0.),
  m_jetLoc(-1.),
  m_viewerTop(0.),
  m_viewerBottom(0.),
  m_jetLength(0.),
  m_jetLoBC(CRDparam::DomainBCTypeDirichlet),
  m_jetHiBC(CRDparam::DomainBCTypeDirichlet),
  m_mappedShift(RealVect::Zero)
{
  readBCInfo();
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }
  // !!! FIXME !!! no need for this hack with the new CNSIBC method
  // Use any velocities for the jet BC's but overwrite the velocity so
  // the walls are stationary. Otherwise, the walls will move
  // BoundaryIndex bcIdx;
  // bcIdx.m_block = 0;
  // bcIdx.m_dir = m_jetDir;

  // // Set the jet type to be whatever is read in from input file
  // bcIdx.m_side = Side::Lo;
  // m_jetLoBC = getDomainBC(bcIdx).m_type;
  // m_jetLoVel = getDomainBCstate(bcIdx).m_vel;
  // BCstate bcState = getDomainBCstate(bcIdx);
  // bcState.m_vel = RealVect::Zero;
  // setDomainBCstate(getDomainBC(bcIdx).m_stateIdx, bcState);
  
  // bcIdx.m_side = Side::Hi;
  // m_jetHiBC = getDomainBC(bcIdx).m_type;
  // m_jetHiVel = getDomainBCstate(bcIdx).m_vel;
  // bcState = getDomainBCstate(bcIdx);
  // bcState.m_vel = RealVect::Zero;
  // setDomainBCstate(getDomainBC(bcIdx).m_stateIdx, bcState);
  
  // // Set the side walls with jets
  // BCInfo bc;
  // bc.m_order = 4;
  // bc.m_type = CRDparam::DomainBCTypeMixed;
  // bcIdx.m_dir = m_jetDir;
  // bcIdx.m_side = Side::Lo;
  // setDomainBC(bcIdx, bc);
  // bcIdx.m_side = Side::Hi;
  // setDomainBC(bcIdx, bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCBurnerTest::~CNSIBCBurnerTest()
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
CNSIBCBurnerTest::IBCName() const
{
  return "Cook-stove burner case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCBurnerTest::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Flow direction\n" << m_flowDir << CRD::var;
  if (m_energyIgnition)
    {
      CRD::msg << "Ignition method\nenergy source" << CRD::var;
    }
  else
    {
      CRD::msg << "Ignition method\nactivation energy" << CRD::var;
    }
  CRD::msg << "Hot spot temperature\n" << m_hotSpotTemp << CRD::var;
  CRD::msg << "Hot spot start time\n" << m_hotSpotStartTime << CRD::var;
  CRD::msg << "Hot spot end time\n" <<
    (m_hotSpotStartTime + m_hotSpotRunTime) << CRD::var;
  CRD::msg << "Gravity force\n" << m_gravForce << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** For the shock tube, tags are set based on gradients of density.
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
CNSIBCBurnerTest::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  // This ensure no tagging occurs above the test section, in the
  // region near the outlet and slip walls
  std::vector<Box> restrictBoxes(1);
  restrictBoxes[0].define(m_slipWalls[0]);
  // Shrink the box a few cells to make sure tagging is completely
  // outside of the testing section
  restrictBoxes[0].growLo(m_flowDir, -1);
  CNSIBCCombustionReference::setTagMethodLevel(a_tagBufferSize, tagLevel,
                                               restrictBoxes);
  {
    // Always refine around jets to ensure jets are captured
    IntVect loLeft(m_atJets.back().smallEnd());
    loLeft[m_flowDir] -= 1;
    IntVect hiRight(m_atJets.back().bigEnd());
    hiRight[m_flowDir] += 1;
    IntVect hiLeft(IntVect::Zero);
    hiLeft[m_flowDir] = hiRight[m_flowDir];
    hiLeft[m_jetDir] += 1;
    IntVect loRight(loLeft);
    loRight[m_jetDir] = hiRight[m_jetDir] - 1;
    Box tagBoxL(loLeft, hiLeft);
    Box tagBoxR(loRight, hiRight);
    tagBoxL.coarsen(CRDparam::g_refFromBase[CRDparam::numAMRLevel()-1]);
    tagBoxR.coarsen(CRDparam::g_refFromBase[CRDparam::numAMRLevel()-1]);
    tagLevel->appendTagMethod(new TagMethodBaseBox(tagBoxL));
    tagLevel->appendTagMethod(new TagMethodBaseBox(tagBoxR));
  }
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
CNSIBCBurnerTest::initialize(LevelData<FArrayBox>&      a_U,
                             LevelGridMetrics&          a_gridMetrics,
                             const LayoutData<FluxBox>& a_unitNormals,
                             const Real                 a_time,
                             const int                  a_level) const
{
  // Set the initial values
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(-1., rComp);
      Wc.setVal(m_initP, presIndx);
      Wc.setVal(m_initT, tempIndx);
      // Set initial velocity to 0
      Wc.setVal(0., box2Dom, WvelIndx, SpaceDim);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Wc.setVal(m_initVel[dir], WvelIndx + dir);
        }
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int wComp = sp + wCompStart;
          Wc.setVal(m_initMassFraction[sp], wComp);
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
    }
}

/*--------------------------------------------------------------------*/
//  Add hot spot source term, this is achieved by giving the reactions a
//  higher temperature which in essense lowers the activation energy
//  FIXME: This is outdated and generally should not be used
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
 *  \param[in]  a_level Grid level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCBurnerTest::addSourceTerm(FArrayBox&           a_sourceFab,
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
                                const DataIndex&     a_dataIndx,
                                const Real           a_globalKE,
                                const Real           a_globalHelicity) const
{
  // Add the source terms to a_sourceFab
  int momComp =
    CRDparam::g_CRDPhysics->vectorFluxInterval().begin() + m_flowDir;
  const int engIndx = CRDparam::g_CRDPhysics->energyFluxIndex();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  if (!(m_gravForce == 0.))
    {
      FORT_ADDPOISEUILLEFORCE(CHF_FRA(a_sourceFab),
                              CHF_BOX(a_solveBox),
                              CHF_CONST_FRA(a_Wcell),
                              CHF_CONST_REAL(m_gravForce),
                              CHF_CONST_INT(rhoIndx),
                              CHF_CONST_INT(engIndx),
                              CHF_CONST_INT(momComp));
    }
  // If hot spot hasn't started or has ended, just return 0
  if (a_time < m_hotSpotStartTime ||
      a_time > m_hotSpotStartTime + m_hotSpotRunTime)
    {
      return;
    }
  const int numSpecies = CRDparam::g_numSpecies;
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  Box currentHotSpot(m_hotSpot[a_level]);
  currentHotSpot &= a_solveBox;
  // Check if hot spot exists in this box and level
  if (currentHotSpot.ok())
    {
      // If we want to add energy source
      if (m_energyIgnition)
        {
          // Percent of time passed for ignition source
          const Real currentP = (a_time - m_hotSpotStartTime)/
            (m_hotSpotStartTime + m_hotSpotRunTime);
          // Modifies source to accomodate time marching
          const Real modP = 1./a_stageWeight;
          std::vector<Real> Cn(numSpecies);
          MD_ARRAY_RESTRICT(arrSource, a_sourceFab);
          MD_ARRAY_RESTRICT(arrW, a_Wcell);
          MD_BOXLOOP(currentHotSpot, i)
            {
              Real rho = arrW[MD_IX(i, rhoIndx)];
              Real T = arrW[MD_IX(i, tempIndx)];
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  int wcomp = wCompStart + sp;
                  Cn[sp] = arrW[MD_IX(i, wcomp)];
                }
              Real hval = CRDparam::g_CRDPhysics->enthalpy(m_hotSpotTemp,
                                                           Cn.data());
              Real hvali = CRDparam::g_CRDPhysics->enthalpy(T, Cn.data());
              // Find the difference between the current enthalpy and the
              // desired raised enthalpy
              arrSource[MD_IX(i, engIndx)] += rho*(hval-hvali)*currentP*modP;
            }
        }
    }
  return;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  !!! FIXME !!! this should be redone with the new CNSIBC method
//  Set the imposed (exterior or farfield) primitive state at flow BC
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
 *  \param[in]  a_domT  Current boundary condition being set (formixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCBurnerTest::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_domT) const
{
  // Set all the component variables and intervals
  const int numSpecies = CRDparam::g_numSpecies;
  Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx = velIntv.begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  Interval specIntv = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int specIndx = specIntv.begin();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  // get the boundary state
  //auto bcState = getDomainBCstate(a_domT);
  const CRDState& bcState = CRDState::get(a_domT.m_idxState);
#if (CH_SPACEDIM != 1)
  // Get a disjoint interior box
  Box interiorBox = adjCellBox(a_boundaryFaceBox,
                               a_bcIdx.m_dir,
                               Side::flip(a_bcIdx.m_side),
                               1);
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(interiorBox));
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_boundaryFaceBox, SpaceDim);// Cartesian coordinates
  FABSTACKTEMP(XFab, a_boundaryFaceBox, SpaceDim); // Physical coordinates
  this->CNSIBC::getFaceCoordinates(a_boundaryFaceBox, 
                                   XiFab,
                                   XFab,
                                   a_bcIdx.m_dir,
                                   blockCoordSys);
  std::vector<Real> a_inletCn(numSpecies);
  // If we are at the upper outlet
  if (a_bcIdx.m_side == Side::Hi && a_bcIdx.m_dir == m_flowDir)
    {
      for (const auto dir : EachDir)
        {
          a_Wface.setVal(m_initVel[dir], a_boundaryFaceBox, velIndx+dir, 1);
        }
      a_Wface.setVal(m_initP, a_boundaryFaceBox, presIndx, 1);
      a_Wface.setVal(m_initT, a_boundaryFaceBox, tempIndx, 1);
      Real sumRn = 0.;
      for (int sp = 0; sp != CRDparam::g_numSpecies; ++sp)
        {
          const int wcomp = sp + specIndx;
          a_Wface.setVal(bcState(wcomp), a_boundaryFaceBox, wcomp, 1);
          sumRn += CRDparam::g_CRDPhysics->speciesGasConstant(sp)*
            bcState(wcomp);
        }
      Real rho = m_initP/(sumRn*m_initT);
      a_Wface.setVal(rho, a_boundaryFaceBox, rhoIndx, 1);
    }
  else // We are at one of the inlets
    {
      Real loLoc, hiLoc, inletT;
      Real Rgas = 0.;
      std::vector<Real> inletCn(numSpecies);
      RealVect inletVel;
      int tanDir;
      // In case we are at a side jet, then the BCVel is zero for the walls
      RealVect sideJetVel;
      inletVel = bcState.velocity();
      inletT = bcState.temperature();
      for (int sp = 0; sp != CRDparam::g_numSpecies; ++sp)
        {
          const int wcomp = sp + specIndx;
          Rgas += CRDparam::g_CRDPhysics->speciesGasConstant(sp)*
            bcState(wcomp);
          inletCn[sp] = bcState(wcomp);
        }
      if (a_bcIdx.m_side == Side::Lo)
        {
          sideJetVel = m_jetLoVel;
        }
      else
        {
          sideJetVel = m_jetHiVel;
        }
      if (a_bcIdx.m_dir == m_flowDir)
        {
          loLoc = 0.;
          hiLoc = CRDparam::g_domainLength[m_jetDir];
          tanDir = m_jetDir;
        }
      else
        {
          loLoc = m_jetLoLoc[a_level];
          hiLoc = m_jetHiLoc[a_level];
          tanDir = m_flowDir;
          inletVel = sideJetVel;
        }
      Real inletGamma =
        CRDparam::g_CRDPhysics->gamma(inletT, inletCn.data());
      a_Wface.setVal(0., a_boundaryFaceBox, velIndx, SpaceDim);
      setJetBCVals(a_boundaryFaceBox, a_Wface, XFab, inletT, inletGamma,
                   inletVel[a_bcIdx.m_dir], Rgas, inletCn, loLoc, hiLoc, tanDir, a_bcIdx.m_dir);
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Set the boxes that are inflow or outflow if mixed BC is being used
//  Make sure that any boxes returned are the same ixType as the input box
/** \param[in]  a_boundaryFaceBox
 *                      Box of boundary being operated on
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_boxVect
 *                      Vector of boxes
 *//*-----------------------------------------------------------------*/

int
CNSIBCBurnerTest::setMixedBC(
  const Box&           a_boundaryFaceBox,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  Vector<Box>&         a_boxVect,
  Vector<BCInfo>&      a_domainBC) const
{
  // Boundary types on mixed boundaries, wall, jet, wall, slip wall
  a_boxVect.resize(4);
  a_domainBC.resize(4);
  // boundary condition, get order and ref state
  BCInfo bc = getDomainBC(a_bcIdx);
  // Box encompasses the entire jet
  Box boundaryJet; // Jets
  Box loWall;      // No slip wall below jet
  Box hiWall;      // No slip wall above jet
  Box slipWall;    // Slip wall
  // Set boxes
  if (a_bcIdx.m_side == Side::Lo)
    {
      // Box encompasses the entire jet
      boundaryJet = adjCellLo(m_atJets[a_level], m_jetDir, 1);
      boundaryJet.shiftHalf(a_bcIdx.m_dir,1);
      // Box of the lower wall on the left
      loWall = adjCellLo(m_belowJets[a_level], m_jetDir, 1);
      loWall.shiftHalf(a_bcIdx.m_dir,1);
      // Box of the upper wall on the left
      hiWall = adjCellLo(m_aboveJets[a_level], m_jetDir, 1);
      hiWall.shiftHalf(a_bcIdx.m_dir,1);
      // Box of the upper slip wall on the left
      slipWall = adjCellLo(m_slipWalls[a_level], m_jetDir, 1);
      slipWall.shiftHalf(a_bcIdx.m_dir,1);
      // Set lower jet type
      bc.m_type = m_jetLoBC;
      a_domainBC[1] = bc;
    }
  else
    {
      // Box encompasses the entire jet
      boundaryJet = adjCellHi(m_atJets[a_level], m_jetDir, 1);
      boundaryJet.shiftHalf(a_bcIdx.m_dir,-1);
      // Box of the lower wall on the left
      loWall = adjCellHi(m_belowJets[a_level], m_jetDir, 1);
      loWall.shiftHalf(a_bcIdx.m_dir,-1);
      // Box of the upper wall on the left
      hiWall = adjCellHi(m_aboveJets[a_level], m_jetDir, 1);
      hiWall.shiftHalf(a_bcIdx.m_dir,-1);
      // Box of the upper slip wall on the right
      slipWall = adjCellHi(m_slipWalls[a_level], m_jetDir, 1);
      slipWall.shiftHalf(a_bcIdx.m_dir,-1);
      bc.m_type = m_jetHiBC;
      a_domainBC[1] = bc;
    }
  loWall &= a_boundaryFaceBox;
  hiWall &= a_boundaryFaceBox;
  // If the current level is too coarse to include the jet, make the box type
  // the same but as an empty box
  if (boundaryJet.ok())
    {
      boundaryJet &= a_boundaryFaceBox;
    }
  else
    {
      boundaryJet.convert(a_boundaryFaceBox.ixType());
    }
  slipWall &= a_boundaryFaceBox;
  // Set the boxes for the low wall
  a_boxVect[0] = loWall;
  bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  a_domainBC[0] = bc;
  // Set the box for the jet
  a_boxVect[1] = boundaryJet;
  // Set the boxes for the high wall
  a_boxVect[2] = hiWall;
  bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  a_domainBC[2] = bc;
  // Set the boxes for the upper slip wall
  a_boxVect[3] = slipWall;
  bc.m_type = CRDparam::DomainBCTypeSlipWall;
  a_domainBC[3] = bc;
  CH_assert(minBox(a_boxVect[3],
                   minBox(a_boxVect[2],
                          minBox(a_boxVect[0], a_boxVect[1])))
            == a_boundaryFaceBox);
  return 0;
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
CNSIBCBurnerTest::readBCInfo()
{
  ParmParse ppIBC("ibc");
  
  // Set the direction of the flow of the gas and the jets
  m_flowDir = 1;
  ppIBC.query("flow_dir",m_flowDir);
  m_jetDir = m_flowDir + 1;
  if (m_jetDir > SpaceDim-1)
    {
      m_jetDir = 0;
    }
  CH_assert(m_jetDir != m_flowDir);
  // Set location of the bottom of the viewing window
  m_viewerBottom = CRDparam::g_domainOrigin[m_flowDir];
  ppIBC.query("viewer_bottom", m_viewerBottom);

  // Hot spot info
  ppIBC.query("hot_spot_temperature", m_hotSpotTemp);
  m_hotSpotStartTime = CRDparam::g_reactionStartTime;
  ppIBC.query("hot_spot_start_time", m_hotSpotStartTime);
  ppIBC.query("hot_spot_run_time", m_hotSpotRunTime);
  Real hotSpotRadius = 0.;
  ppIBC.query("hot_spot_radius", hotSpotRadius);
  std::vector<Real> inputHSCenter(SpaceDim);
  ppIBC.queryarr("hot_spot_center", inputHSCenter, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_hotSpotCenter.dataPtr(),
                                           &inputHSCenter.front());
  ppIBC.query("energy_ignition", m_energyIgnition);

  ppIBC.query("gravity_force", m_gravForce);
  if (m_gravForce > 0.)
    {
      m_gravForce *= -1.;
    }

  // Set the top of the viewing window location
  m_viewerTop = CRDparam::g_domainLength[m_flowDir];
  ppIBC.query("viewer_size", m_viewerTop);
  m_viewerTop += m_viewerBottom;
  if (m_viewerTop < 0. || m_viewerTop >= CRDparam::g_domainLength[m_flowDir])
    {
      CRD::msg << "Input (BurnerTest IBC): 'viewer_size' must be "
        "between 0 and the length" << CRD::error;
    }

  // Set the size of the jets
  m_jetLength = 0.1*CRDparam::g_domainLength[m_flowDir];
  ppIBC.query("jet_length", m_jetLength);
  if (m_jetLength < 0.)
    {
      CRD::msg << "Input (BurnerTest IBC): 'jet_length' must be "
        "greater than zero" << CRD::error;
    }
  // Set the location of the air inlets relative to the domain length in the
  // flow direction
  m_jetLoc = 0.5*CRDparam::g_domainLength[m_flowDir];
  ppIBC.query("jet_location", m_jetLoc);
  m_jetLoc += m_viewerBottom;
  if (m_jetLoc < CRDparam::g_domainOrigin[m_flowDir] ||
      m_jetLoc > (CRDparam::g_domainOrigin[m_flowDir] +
                  CRDparam::g_domainLength[m_flowDir]))
    {
      CRD::msg << "Input (BurnerTest IBC): 'jet_location' must be "
               << "greater than the origin and less than the size of the domain"
               << CRD::error;
    }

  // If log stretch mapping is used, we must take that into account
  ParmParse ppCOORD("coordsys");
  std::string COORDSYSname;
  ppCOORD.query("type", COORDSYSname);
  if (COORDSYSname == "logstretch" && ppCOORD.contains("shift"))
    {
      std::vector<Real> inputShift(SpaceDim);
      ppCOORD.getarr("shift", inputShift, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_mappedShift.dataPtr(),
                                               &inputShift.front());
    }
  // Set the boxes
  int numLevels = CRDparam::numAMRLevel();
  m_belowJets.resize(numLevels);
  m_atJets.resize(numLevels);
  m_aboveJets.resize(numLevels);
  m_slipWalls.resize(numLevels);
  m_hotSpot.resize(numLevels);
  m_jetLoLoc.resize(numLevels);
  m_jetHiLoc.resize(numLevels);
  RealVect dxVect = CRDparam::g_domainLength/CRDparam::g_domainBaseSize;
  IntVect gridSize = CRDparam::g_domainBaseSize;
  // Bool to ensure that the dx on the refined region is smaller than the jet
  // length
  // These checks will be redefined in initialization in case mapping exists
  bool jetCheck = false;
  for (int lvl = 0; lvl != numLevels; ++lvl)
    {
      Real dx = dxVect[m_flowDir]/CRDparam::g_refFromBase[lvl];
      IntVect curGridSize = gridSize*CRDparam::g_refFromBase[lvl] -
        IntVect::Unit;
      Box domain(IntVect::Zero, curGridSize);
      // Below jets
      IntVect loBJ(IntVect::Zero);
      IntVect hiBJ(curGridSize);
      // At jets
      IntVect loAtJ(IntVect::Zero);
      IntVect hiAtJ(curGridSize);
      hiAtJ[m_flowDir] = 0;
      // Above jets
      IntVect loAJ(IntVect::Zero);
      IntVect hiAJ(curGridSize);
      // Above viewer
      IntVect loAV(IntVect::Zero);
      IntVect hiAV(curGridSize);
      // Hot spot
      IntVect HSCenterLo(IntVect::Zero);
      IntVect HSCenterHi(IntVect::Zero);
      for (int j = 0; j != curGridSize[m_flowDir]; ++j)
        {
          Real Xj = dx*(j + 0.5);
          Real X = coordOut(Xj);
          if (X < m_jetLoc)
            {
              hiBJ[m_flowDir] = j;
            }
          else if (X > m_jetLoc &&
                   X < m_jetLoc + m_jetLength)
            {
              hiAtJ[m_flowDir] = j;
            }
          else if (X > m_jetLoc + m_jetLength &&
                   X < m_viewerTop)
            {
              hiAJ[m_flowDir] = j;
            }
        }
      loAtJ[m_flowDir] = hiBJ[m_flowDir] + 1;
      loAJ[m_flowDir] = hiAtJ[m_flowDir] + 1;
      loAV[m_flowDir] = hiAJ[m_flowDir] + 1;
      m_belowJets[lvl].define(loBJ, hiBJ);
      if (hiAtJ[m_flowDir] - loAtJ[m_flowDir] >= 2)
        {
          m_jetLoLoc[lvl] = coordOut(dx*(loAtJ[m_flowDir] + 0.5));
          m_jetHiLoc[lvl] = coordOut(dx*(hiAtJ[m_flowDir] + 0.5));
          m_atJets[lvl].define(loAtJ, hiAtJ);
          m_aboveJets[lvl].define(loAJ, hiAJ);
          jetCheck = true;
        }
      else
        {
          m_aboveJets[lvl].define(loAtJ, hiAJ);
        }
      m_slipWalls[lvl].define(loAV, hiAV);
    }
  if (!jetCheck)
    {
      CRD::msg << "Input (BurnerTest IBC): Smallest dx is too large to  "
               << "accomodate 'jet_length'. Either make jets larger or"
               << " increase refinement ratio or levels" << CRD::error;
    }


  m_readInput = true;
  CRD::msg.setTerminateOnError(true);
}

/*--------------------------------------------------------------------*/
//  Fill the parabolic jet profile
/** \param[in]  a_box   Current box
 *  \param[out] a_Wface Face values adjusted based on BC
 *  \param[in]  a_XFab  Face locations
 *  \param[in]  a_inlet Corresponding inlet values
 *  \param[in]  a_loLoc Lower location of the inlet
 *  \param[in]  a_hiLoc Upper location of the inlet
 *//*-----------------------------------------------------------------*/

void
CNSIBCBurnerTest::setJetBCVals(const Box&               a_box,
                               FArrayBox&               a_Wface,
                               const FArrayBox&         a_XFab,
                               const Real&              a_inletT,
                               const Real&              a_inletGamma,
                               const Real&              a_inletVel,
                               const Real&              a_Rgas,
                               const std::vector<Real>& a_inletCn,
                               const Real&              a_loLoc,
                               const Real&              a_hiLoc,
                               const int&               a_tanDir,
                               const int&               a_normDir) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx = velIntv.begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const Interval specIntv = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int specIndx = specIntv.begin();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  Real len = a_hiLoc - a_loLoc;
  Real xm = a_loLoc + len*0.5;
  Real div = a_loLoc - xm;
  Real coeff = -a_inletVel/(div*div);
  Real machSqinv = 1./(sqrt(a_inletGamma*a_inletT*a_Rgas));
  Real alpha = a_inletGamma/(a_inletGamma - 1.);
  MD_ARRAY_RESTRICT(arrX, a_XFab);
  MD_ARRAY_RESTRICT(arrW, a_Wface);
  MD_BOXLOOP(a_box, i)
    {
      Real xloc = arrX[MD_IX(i, a_tanDir)];
      Real pres = arrW[MD_IX(i, presIndx)];
      Real velVal = coeff*(xloc - xm)*(xloc - xm) + a_inletVel;
      Real Msq = std::pow(velVal*machSqinv, 2);
      pres = pres*std::pow(1. + 0.5*(a_inletGamma - 1.)*Msq,alpha);
      arrW[MD_IX(i, presIndx)] = pres;
      arrW[MD_IX(i, tempIndx)] = a_inletT;
      arrW[MD_IX(i, velIndx + a_normDir)] = velVal;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int wComp = sp + specIndx;
          arrW[MD_IX(i, wComp)] = a_inletCn[sp];
        }
      arrW[MD_IX(i, rhoIndx)] = pres/(a_Rgas*a_inletT);
    }
  return;
}  
  
/*--------------------------------------------------------------------*/
//  Return the physical location based on the log stretched mapping
//  FIXME: This type of function should come from the mapping, however,
//         the mapping has not been set yet
/** \param[in]          Computational location
 *  \return             The real location in the flowDir
 *//*-----------------------------------------------------------------*/

Real
CNSIBCBurnerTest::coordOut(const Real a_Xi)
{
  Real realLoc;
  const int d = m_flowDir;
  Real eShift = exp(m_mappedShift[d]) - 1;
  RealVect domL = CRDparam::g_domainLength;
  if (m_mappedShift[d] == 0)
    {
      realLoc = a_Xi;
    }
  else
    {
      realLoc = log(a_Xi/domL[d]*eShift+1.)*domL[d]/m_mappedShift[d];
    }
  return realLoc;
}
