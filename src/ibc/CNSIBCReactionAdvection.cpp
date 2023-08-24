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
 * \file CNSIBCReactionAdvection.cpp
 *
 * \brief Member functions for CNSIBCReactionAdvection
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCReactionAdvection.H"
#include "CNSIBCCombustionTestF_F.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "ChordInput.H"
#include "CRDState.H"
#include "AMRIO.H"


/*******************************************************************************
 *
 * Class CNSIBCReactionAdvection: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCReactionAdvection::CNSIBCReactionAdvection()
  :
  CNSIBCCombustionReference(),
  m_C(-1.),
  m_x0(0.),
  m_L(-1.),
  m_o2T(298.15),
  m_h2T(298.15),
  m_pressure(101325.),
  m_initVel(RealVect::Zero),
  m_O2comp(-1),
  m_H2comp(-1)
{
  readBCInfo();

  // Define states
  const int cSpecBeg = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  std::string stateName("H2-region");
  if(CRDState::nameIndex(stateName) == -1)
    {
      CRDState& state = CRDState::get(stateName);
      state.temperature() = m_h2T;
      state.pressure() = m_pressure;
      state.velocity() = m_initVel;
      const Real ch2 = 1.;
      const Real co2 = 0.;
      state(m_H2comp + cSpecBeg) = ch2;
      state(m_O2comp + cSpecBeg) = co2;
      state.setExtraThermo();
    }
  stateName = "O2-region";
  if(CRDState::nameIndex(stateName) == -1)
    {
      CRDState& state = CRDState::get(stateName);
      state.temperature() = m_o2T;
      state.pressure() = m_pressure;
      state.velocity() = m_initVel;
      const Real ch2 = 0.;
      const Real co2 = 1.;
      state(m_H2comp + cSpecBeg) = ch2;
      state(m_O2comp + cSpecBeg) = co2;
      state.setExtraThermo();
    }
  CRDState::writeStateInfo();

  // Default for BC
  BCInfo defBC;
  defBC.m_order = 1;
  defBC.m_type = CRDparam::DomainBCTypePeriodic;
  setAllDomainBC(defBC);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCReactionAdvection::~CNSIBCReactionAdvection()
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
CNSIBCReactionAdvection::IBCName() const
{
  return "Advection of H2-O2 flame problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCReactionAdvection::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Center location\n" << m_x0 << CRD::var;
  CRD::msg << "Plateau width\n" << m_L << CRD::var;
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
CNSIBCReactionAdvection::initialize(LevelData<FArrayBox>&      a_U,
                                    LevelGridMetrics&          a_gridMetrics,
                                    const LayoutData<FluxBox>& a_unitNormals,
                                    const Real                 a_time,
                                    const int                  a_level) const
{
  // Set the initial values
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
  const Real RO2 = CRDparam::g_CRDPhysics->speciesGasConstant(m_O2comp);
  const Real RH2 = CRDparam::g_CRDPhysics->speciesGasConstant(m_H2comp);
  const Real rhoO2 = m_pressure/(RO2*m_o2T);
  const Real rhoH2 = m_pressure/(RH2*m_h2T);
  const int o2comp = wCompStart + m_O2comp;
  const int h2comp = wCompStart + m_H2comp;
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
          Real x = arrX[MD_IX(i, 0)];
          const Real ch2 =
            0.5*(1. + std::tanh(m_C*(m_L/2. - std::abs(x - m_x0))*100.));
          const Real co2 = 1. - ch2;
          const Real rhomix = 1./(ch2/rhoH2 + co2/rhoO2);
          arrW[MD_IX(i, rComp)] = rhomix;
          arrW[MD_IX(i, pComp)] = m_pressure;
          arrW[MD_IX(i, tComp)] = -1.;
          D_TERM(
            arrW[MD_IX(i, WvelIndx)] = m_initVel[0];,
            arrW[MD_IX(i, WvelIndx+1)] = m_initVel[1];,
            arrW[MD_IX(i, WvelIndx+2)] = m_initVel[2];);
          arrW[MD_IX(i, o2comp)] = co2;
          arrW[MD_IX(i, h2comp)] = ch2;
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
CNSIBCReactionAdvection::haveExactSol() const
{
  return false;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

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
CNSIBCReactionAdvection::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  std::vector<Real> initVel(SpaceDim);
  ppIBC.getarr("init_velocity", initVel, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_initVel.dataPtr(),
                                           &initVel.front());
  ppIBC.get("temp_o2",m_o2T);
  ppIBC.get("temp_h2",m_h2T);
  ppIBC.get("pressure",m_pressure);
  ppIBC.get("init_const",m_C);
  ppIBC.get("center_loc",m_x0);
  ppIBC.get("plat_width",m_L);
  // Retrieve specific components
  for (int i = 0; i != numSpecies; ++i)
    {
      if (CRDparam::g_speciesNames[i] == "H2" || CRDparam::g_speciesNames[i] == "N2")
        {
          m_H2comp = i;
        }
      else if (CRDparam::g_speciesNames[i] == "O2")
        {
          m_O2comp = i;
        }
    }
  CH_assert(m_O2comp >= 0);
  CH_assert(m_H2comp >= 0);
  m_readInput = true;
}
