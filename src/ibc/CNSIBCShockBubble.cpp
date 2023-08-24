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
 * \file CNSIBCShockBubble.cpp
 *
 * \brief Member functions for CNSIBCShockBubble
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
 
//----- Internal -----//

#include "CNSIBCShockBubble.H"
#include "CNSIBCFlameF_F.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"
#include "ChordInput.H"

/*******************************************************************************
 *
 * Class CNSIBCShockBubble: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCShockBubble::CNSIBCShockBubble()
  :
  CNSIBCCombustionReference(),
  m_C(-1.),
  m_radius(1.),
  m_shockLoc(-1.),
  m_t0(298.15),
  m_t1(298.15),
  m_p0(101325.0),
  m_p1(101325.0),
  m_U0(RealVect::Zero),
  m_U1(RealVect::Zero),
  m_Center(RealVect::Zero),
  m_O2comp(-1),
  m_H2comp(-1),
  m_N2comp(-1),
  m_loBC(0),
  m_hiBC(0),
  m_bcOrder(1)
{
  readBCInfo();
  BCInfo defBC;
  defBC.m_type = CRDparam::DomainBCTypePeriodic;
  defBC.m_order = m_bcOrder;
  setAllDomainBC(defBC);

//--Define states
//**FIXME Fix input to read a state rather than specific variables

  const Real o2airmix = 0.233;
  const Real n2airmix = 0.767;
  std::string stateName("pre-shock");
  const int cSpecBeg = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  if (CRDState::nameIndex(stateName) == -1) // Define with default values
    {
      CRDState& state = CRDState::get(stateName);
      state.temperature() = m_t0;
      state.pressure()    = m_p0;
      state.velocity()    = m_U0;
      const Real yH2 = 0.5*(1. + std::tanh(m_radius*100./m_C));
      state(m_H2comp + cSpecBeg) = yH2;
      state(m_O2comp + cSpecBeg) = o2airmix*(1. - yH2);
      state(m_N2comp + cSpecBeg) = n2airmix*(1. - yH2);
      state.setExtraThermo();
    }
  stateName = "post-shock";
  if (CRDState::nameIndex(stateName) == -1) // Define with default values
    {
      CRDState& state = CRDState::get(stateName);
      state.temperature() = m_t1;
      state.pressure()    = m_p1;
      state.velocity()    = m_U1;
      state(m_O2comp + cSpecBeg) = o2airmix;
      state(m_N2comp + cSpecBeg) = n2airmix;
      state.setExtraThermo();
    }
  CRDState::writeStateInfo();

  // X-direction, low side
  BoundaryIndex idxB = { 0, 0, Side::Lo };
  if (m_loBC == 0)
    {
      defBC.m_type = CRDparam::DomainBCTypeExtrapolated;
    }
  else if (m_loBC == 1)
    {
      defBC.m_type = CRDparam::DomainBCTypeOutflow;
    }
  else if (m_loBC == 2)
    {
      defBC.m_type = CRDparam::DomainBCTypeDirichlet;
    }
  else
    {
      CRD::msg << "Input (ShockBubble IBC): 'low_bc' must be in range 0:2!"
               << CRD::error;
    }
  defBC.m_idxState = CRDState::nameIndex("pre-shock");
  setDomainBC(idxB, defBC);

  // X-direction, high side
  idxB.m_side = Side::Hi;
  if (m_hiBC == 0)
    {
      defBC.m_type = CRDparam::DomainBCTypeExtrapolated;
    }
  else if (m_hiBC == 1)
    {
      defBC.m_type = CRDparam::DomainBCTypeOutflow;
    }
  else if (m_hiBC == 2)
    {
      defBC.m_type = CRDparam::DomainBCTypeDirichlet;
    }
  else
    {
      CRD::msg << "Input (ShockBubble IBC): 'hi_bc' must be in range 0:2!"
               << CRD::error;
    }
  defBC.m_idxState = CRDState::nameIndex("post-shock");
  setDomainBC(idxB, defBC);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCShockBubble::~CNSIBCShockBubble()
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
CNSIBCShockBubble::IBCName() const
{
  return "Shock bubble combustion problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCShockBubble::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Center location\n" << m_Center << CRD::var;
  if (m_loBC == 0)
    {
      CRD::msg << "Low boundary\nExtrapolated" << CRD::var;
    }
  else if (m_loBC == 1)
    {
      CRD::msg << "Low boundary\nOutflow" << CRD::var;
    }
  else if (m_loBC == 2)
    {
      CRD::msg << "Low boundary\nDirichlet" << CRD::var;
    }
  if (m_hiBC == 0)
    {
      CRD::msg << "High boundary\nExtrapolated" << CRD::var;
    }
  else if (m_hiBC == 1)
    {
      CRD::msg << "High boundary\nOutflow" << CRD::var;
    }
  else if (m_hiBC == 2)
    {
      CRD::msg << "High boundary\nDirichlet" << CRD::var;
    }
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
CNSIBCShockBubble::initialize(LevelData<FArrayBox>&      a_U,
                              LevelGridMetrics&          a_gridMetrics,
                              const LayoutData<FluxBox>& a_unitNormals,
                              const Real                 a_time,
                              const int                  a_level) const
{
  // Set the initial values
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
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
      // Create a FAB of the initial primitive variables
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);
      // Set velocity to zero
      Wc.setVal(0., box2Dom, WvelIndx, SpaceDim);
      FORT_CNSIBCSHOCKBUBBLE(CHF_FRA(Wc),
                             CHF_BOX(box2Dom),
                             CHF_CONST_FRA(XFab),
                             CHF_CONST_INT(rComp),
                             CHF_CONST_INT(pComp),
                             CHF_CONST_INT(tComp),
                             CHF_CONST_INT(WvelIndx),
                             CHF_CONST_INT(wCompStart),
                             CHF_CONST_INT(m_O2comp),
                             CHF_CONST_INT(m_H2comp),
                             CHF_CONST_INT(m_N2comp),
                             CHF_CONST_REALVECT(m_Center),
                             CHF_CONST_REAL(m_radius),
                             CHF_CONST_REALVECT(m_U0),
                             CHF_CONST_REALVECT(m_U1),
                             CHF_CONST_REAL(m_C),
                             CHF_CONST_REAL(m_shockLoc),
                             CHF_CONST_REAL(m_t0),
                             CHF_CONST_REAL(m_t1),
                             CHF_CONST_REAL(m_p0),
                             CHF_CONST_REAL(m_p1));
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
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCShockBubble::haveExactSol() const
{
  return false;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Sets the boundary conditions
/** 
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
 *//*-----------------------------------------------------------------*/

void
CNSIBCShockBubble::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
  const Real o2p = 0.24;
  const Real n2p = 0.76;
  const Real RO2 = CRDparam::g_CRDPhysics->speciesGasConstant(m_O2comp);
  const Real RN2 = CRDparam::g_CRDPhysics->speciesGasConstant(m_N2comp);
  const int numSpecies = CRDparam::g_numSpecies;
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  a_Wface.setVal(0., a_boundaryFaceBox, wCompStart, numSpecies);
  a_Wface.setVal(0., a_boundaryFaceBox, WvelIndx, SpaceDim);
  a_Wface.setVal(o2p, a_boundaryFaceBox, wCompStart+m_O2comp,1);
  a_Wface.setVal(n2p, a_boundaryFaceBox, wCompStart+m_N2comp,1);
  if (a_side == Side::Lo)
    {
      const Real rho = m_p0/((o2p*RO2+n2p*RN2)*m_t0);
      a_Wface.setVal(rho, a_boundaryFaceBox, rComp, 1);
      a_Wface.setVal(m_p0, a_boundaryFaceBox, pComp, 1);
      a_Wface.setVal(m_t0, a_boundaryFaceBox, tComp, 1);
      for (int i = 0; i != SpaceDim; ++i)
        {
          a_Wface.setVal(m_U0[i], a_boundaryFaceBox, WvelIndx + i,1);
        }
    }
  else if (a_side == Side::Hi)
    {
      const Real rho = m_p1/((o2p*RO2+n2p*RN2)*m_t1);
      a_Wface.setVal(rho, a_boundaryFaceBox, rComp, 1);
      a_Wface.setVal(m_p1, a_boundaryFaceBox, pComp, 1);
      a_Wface.setVal(m_t1, a_boundaryFaceBox, tComp, 1);
      for (int i = 0; i != SpaceDim; ++i)
        {
          a_Wface.setVal(m_U1[i], a_boundaryFaceBox, WvelIndx + i,1);
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
CNSIBCShockBubble::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  std::vector<Real> U1(SpaceDim);
  ppIBC.getarr("vel_1", U1, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_U1.dataPtr(),
                                           &U1.front());
  std::vector<Real> U0(SpaceDim);
  ppIBC.getarr("vel_0", U0, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_U0.dataPtr(),
                                           &U0.front());
  std::vector<Real> center(SpaceDim);
  ppIBC.getarr("center", center, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_Center.dataPtr(),
                                           &center.front());
  ppIBC.get("init_const",m_C);
  ppIBC.get("temp_0",m_t0);
  ppIBC.get("temp_1",m_t1);
  ppIBC.get("pres_0",m_p0);
  ppIBC.get("pres_1",m_p1);
  ppIBC.get("shock_loc",m_shockLoc);
  ppIBC.get("radius",m_radius);
  // Retrieve specific components
  for (int i = 0; i != numSpecies; ++i)
    {
      if (CRDparam::g_speciesNames[i] == "H2")
        {
          m_H2comp = i;
        }
      else if (CRDparam::g_speciesNames[i] == "O2")
        {
          m_O2comp = i;
        }
      else if (CRDparam::g_speciesNames[i] == "N2")
        {
          m_N2comp = i;
        }
    }
  if (m_O2comp < 0 || m_H2comp < 0 || m_N2comp < 0)
    {
      CRD::msg << "Input (ShockBubble IBC): Must have O2, H2, and N2!"
               << CRD::error;
    }
  ppIBC.query("low_bc",m_loBC);
  ppIBC.query("hi_bc",m_hiBC);
  ppIBC.query("bc_order",m_bcOrder);
  if (m_bcOrder != 1 && m_bcOrder != 4)
    {
      CRD::msg << "Input (ShockBubble IBC): 'bc_order' must be 1 or 4!"
               << CRD::error;
    }
  m_readInput = true;
}
