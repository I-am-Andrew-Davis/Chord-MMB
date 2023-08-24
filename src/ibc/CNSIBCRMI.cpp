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
 * \file CNSIBCRMI.cpp
 *
 * \brief Member functions for CNSIBCRMI
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "CONSTANTS.H"
 
//----- Internal -----//

#include "CNSIBCRMI.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"
#include "ChordInput.H"

/*******************************************************************************
 *
 * Class CNSIBCRMI: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCRMI::CNSIBCRMI()
  :
  CNSIBCCombustionReference(),
  m_xs(-1.),
  m_lambda(-1.),
  m_rhoo2(-1.),
  m_rhoh2(-1.),
  m_Patm(101325.),
  m_presRatio(1.2),
  m_rhoRatio(1.2),
  m_U0(RealVect::Zero),
  m_O2comp(-1),
  m_H2comp(-1),
  m_loBC(0),
  m_bcOrder(4),
  m_h0(1.),
  m_Delta(1.)
{
  readBCInfo();
  // set default boundary
  BCInfo defBC;
  defBC.m_type = CRDparam::DomainBCTypePeriodic;
  defBC.m_order = m_bcOrder;
  setAllDomainBC(defBC);

  const int cSpecBeg = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();

  std::string stateName("upstream-extrapolated");
  if (CRDState::nameIndex(stateName) == -1)
    {
      CRDState& state = CRDState::get(stateName);
      const Real shockedPres = m_Patm*m_presRatio;
      const Real shockedRho = m_rhoh2*m_rhoRatio;
      state.pressure() = shockedPres;
      state.density() = shockedRho;
      state.velocity() = m_U0;
      state(m_H2comp + cSpecBeg) = 1.;
      state.setExtraThermo();
    }
  stateName = "downstream-slipwall";
  if (CRDState::nameIndex(stateName) == -1)
    {
      CRDState& state = CRDState::get(stateName);
      state.pressure() = m_Patm;
      state.density() = m_rhoo2;
      state.velocity() = RealVect::Zero;
      state(m_O2comp + cSpecBeg) = 1.;
      state.setExtraThermo();
    }

  // sets left side wall
  BoundaryIndex bcIdx;
  bcIdx.m_block = 0;
  bcIdx.m_side = Side::Lo;
  bcIdx.m_dir = 0;
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
      CRD::msg << "Input (RMI IBC): 'low_bc' must be 0, 1, or 2!"
               << CRD::error;
    }
  defBC.m_idxState = CRDState::nameIndex("upstream-extrapolated");
  setDomainBC(bcIdx, defBC);
  // sets right side wall
  bcIdx.m_side = Side::Hi;
  defBC.m_type = CRDparam::DomainBCTypeSlipWall;
  defBC.m_idxState = CRDState::nameIndex("downstream-slipwall");
  setDomainBC(bcIdx, defBC);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCRMI::~CNSIBCRMI()
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
CNSIBCRMI::IBCName() const
{
  return "Richtmyer-Meshkov instability problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCRMI::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Shock location\n" << m_xs << CRD::var;
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
CNSIBCRMI::initialize(LevelData<FArrayBox>&      a_U,
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
  const Real shockedPres = m_Patm*m_presRatio;
  const Real shockedRho = m_rhoh2*m_rhoRatio;
  // FIXME: The actual equation for WVal is 2*abs(erfinv(1 - 2*1.0E-5))
  const Real WVal = 6.03154664028055;
  const Real YH2o = 1.;
  const Real YO2o = 1.;
  // Stoichiometrics mixture of Y_fuel/Y_oxid = s = 8 for H2 and O2
  const Real s = 8.;
  const int h2comp = m_H2comp + wCompStart;
  const int o2comp = m_O2comp + wCompStart;
  D_TERM(
    const int velIndx1 = WvelIndx;,
    const int velIndx2 = WvelIndx + 1;,
    const int velIndx3 = WvelIndx + 2;);
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
      MD_ARRAY_RESTRICT(arrW, Wc);
      MD_ARRAY_RESTRICT(arrX, XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          D_TERM(
            Real x = arrX[MD_IX(i, 0)];,
            Real y = arrX[MD_IX(i, 1)] - m_lambda/2.;,
             (void)0);
          D_SELECT(
            const Real xi = m_xs + m_Delta/2.;,
            const Real xi = m_xs + m_Delta/2. +
            m_h0*(1. - cos(2.*Pi*y/m_lambda));,
            const Real xi = m_xs + m_Delta/2. +
            m_h0*(1. - cos(2.*Pi*y/m_lambda)););
          if (x < m_xs)
            {
              arrW[MD_IX(i, rComp)] = shockedRho;
              arrW[MD_IX(i, pComp)] = shockedPres;
              arrW[MD_IX(i, tComp)] = -1.;
              arrW[MD_IX(i, h2comp)] = 1.;
              D_TERM(
                arrW[MD_IX(i, velIndx1)] = m_U0[0];,
                arrW[MD_IX(i, velIndx2)] = m_U0[1];,
                arrW[MD_IX(i, velIndx3)] = m_U0[2];);
            }
          else
            {
              const Real Z = 0.5*(1. - erf((x - xi)*WVal/m_Delta));
              const Real YO2 = -(Z*(s*YH2o + YO2o) - YO2o - s)/(s + 1);
              const Real YH2 = 1. - YO2;
              const Real rhoMix = 1./(YH2/m_rhoh2 + YO2/m_rhoo2);
              arrW[MD_IX(i, rComp)] = rhoMix;
              arrW[MD_IX(i, pComp)] = m_Patm;
              arrW[MD_IX(i, h2comp)] = YH2;
              arrW[MD_IX(i, o2comp)] = YO2;
              arrW[MD_IX(i, tComp)] = -1.;
            }
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
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCRMI::haveExactSol() const
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
CNSIBCRMI::setImposedBCprimState(
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
  const int numSpecies = CRDparam::g_numSpecies;
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  a_Wface.setVal(0., a_boundaryFaceBox, wCompStart, numSpecies);
  a_Wface.setVal(0., a_boundaryFaceBox, WvelIndx, SpaceDim);
  if (a_side == Side::Lo)
    {
      const Real shockedPres = m_Patm*m_presRatio;
      const Real shockedRho = m_rhoh2*m_rhoRatio;
      a_Wface.setVal(shockedRho, a_boundaryFaceBox, rComp, 1);
      a_Wface.setVal(shockedPres, a_boundaryFaceBox, pComp, 1);
      a_Wface.setVal(1., a_boundaryFaceBox, wCompStart+m_H2comp, 1);
      for (int i = 0; i != SpaceDim; ++i)
        {
          a_Wface.setVal(m_U0[i], a_boundaryFaceBox, WvelIndx + i,1);
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
CNSIBCRMI::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  ppIBC.get("shock_loc",m_xs);
  if (m_xs < CRDparam::g_domainOrigin[0] || m_xs >
     CRDparam::g_domainOrigin[0] + CRDparam::g_domainLength[0])
    {
      CRD::msg << "Input (RMI IBC): 'shock_loc' must be within domain!"
               << CRD::error;
    }
  ppIBC.get("lambda",m_lambda);
  ppIBC.get("rho_o2",m_rhoo2);
  ppIBC.get("rho_h2",m_rhoh2);
  ppIBC.query("low_bc",m_loBC);
  // FIXME: instead of requiring pres ratio, rho ratio, and U0
  // should only need Mach number and calculate the ratios using
  // Rankine Hugoniot relations
  ppIBC.get("pres_ratio",m_presRatio);
  ppIBC.get("rho_ratio",m_rhoRatio);
  std::vector<Real> U0(SpaceDim);
  ppIBC.getarr("shock_vel", U0, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_U0.dataPtr(),
                                           &U0.front());
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
    }
  if (m_O2comp < 0 || m_H2comp < 0)
    {
      CRD::msg << "Input (RMI IBC): Must have O2 and H2!" << CRD::error;
    }
  ppIBC.query("bc_order",m_bcOrder);
  if (m_bcOrder != 1 && m_bcOrder != 4)
    {
      CRD::msg << "Input (RMI IBC): 'bc_order' must be 1 or 4!"
               << CRD::error;
    }
  const Real KVal = 2.*Pi/m_lambda;
  m_h0 = 0.2/KVal;
  m_Delta = 2.*m_h0;
  ppIBC.query("init_amplitude",m_h0);
  ppIBC.query("interface_thickness",m_Delta);
  m_readInput = true;
}
