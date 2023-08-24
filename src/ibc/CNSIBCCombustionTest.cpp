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
 * \file CNSIBCCombustionTest.cpp
 *
 * \brief Member functions for CNSIBCCombustionTest
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCCombustionTest.H"
#include "CNSIBCCombustionTestF_F.H"
#include "CNSIBCF_F.H"
#include "LGintegrator.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "ViscousTensor4thOrderOp.H"
#include "TagLevel.H"
#include "TagLevelFactory.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCCombustionTest: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCCombustionTest::CNSIBCCombustionTest()
  :
  CNSIBCCombustionReference(),
  m_upRho(2.),
  m_loRho(1.),
  m_grav(-9.8),
  m_pressure(200000.)
{
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }
  readBCInfo();
  // default is all periodic
  //setAllDomainBC(CRDparam::DomainBCTypePeriodic);
  BoundaryIndex bcIdx;
  bcIdx.m_block = 0;
  bcIdx.m_dir = 1;
  BCInfo bc;
  bc.m_order = 1;
  bc.m_type = CRDparam::DomainBCTypeSlipWall;
  // set lo side
  bcIdx.m_side = Side::Lo;
  setDomainBC(bcIdx, bc);
  // set hi side
  bcIdx.m_side = Side::Hi;
  setDomainBC(bcIdx, bc);
  // setDomainBC(1, 0, CRDparam::DomainBCTypeSlipWall);
  // setDomainBC(1, 1, CRDparam::DomainBCTypeSlipWall);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCCombustionTest::~CNSIBCCombustionTest()
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
CNSIBCCombustionTest::IBCName() const
{
  return "Rayleigh-Taylor instability test";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCCombustionTest::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Upper mass fractions\n(";
  CRD::msg << m_upMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_upMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Lower mass fractions\n(";
  CRD::msg << m_loMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_loMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Gravity force\n" << m_grav << CRD::var;
  CRD::msg << "Upper density\n" << m_upRho << CRD::var;
  CRD::msg << "Lower density\n" << m_loRho << CRD::var;
  CRD::msg << "Initial pressure\n" << m_pressure << CRD::var;
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
CNSIBCCombustionTest::initialize(LevelData<FArrayBox>&      a_U,
                                 LevelGridMetrics&          a_gridMetrics,
                                 const LayoutData<FluxBox>& a_unitNormals,
                                 const Real                 a_time,
                                 const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
  // Set the initial values
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int numSpecies = CRDparam::g_numSpecies;
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const Real Ly = CRDparam::g_domainLength[1] -
    CRDparam::g_domainOrigin[1];
  const Real Lx = (CRDparam::g_domainLength[0] -
                   CRDparam::g_domainOrigin[0])/2.;
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(box));
      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);// Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim); // Physical coordinates
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0., box2Dom, WvelIndx, SpaceDim);
      Wc.setVal(m_pressure, pComp);
      Wc.setVal(-1., tComp);
      this->CNSIBC::getCellCoordinates(box2Dom, 
                                       XiFab,
                                       XFab,
                                       blockCoordSys);
      FORT_CNSIBCRTINIT(CHF_FRA(Wc),
                        CHF_BOX(box2Dom),
                        CHF_CONST_FRA(XFab),
                        CHF_CONST_REAL(m_loRho),
                        CHF_CONST_REAL(m_upRho),
                        CHF_CONST_REAL(Ly),
                        CHF_CONST_REAL(Lx),
                        CHF_CONST_INT(rComp),
                        CHF_CONST_INT(wCompStart),
                        CHF_CONST_INT(numSpecies),
                        CHF_CONST_VR(m_loMassFraction),
                        CHF_CONST_VR(m_upMassFraction));
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
    }
#endif 
}

/*--------------------------------------------------------------------*/
//  Add gravity body force
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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCCombustionTest::addSourceTerm(FArrayBox&           a_sourceFab,
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
  const int momcomp = CRDparam::g_CRDPhysics->vectorFluxInterval().begin() + 1;
  const int rcomp = CRDparam::g_CRDPhysics->densityIndex();
  const int engcomp = CRDparam::g_CRDPhysics->energyFluxIndex();
  
  FORT_GRAVITYFORCE(CHF_FRA(a_sourceFab),
                    CHF_BOX(a_solveBox),
                    CHF_CONST_FRA(a_Wcell),
                    CHF_CONST_REAL(m_grav),
                    CHF_CONST_INT(rcomp),
                    CHF_CONST_INT(momcomp),
                    CHF_CONST_INT(engcomp));
  return;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCCombustionTest::haveExactSol() const
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
CNSIBCCombustionTest::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  m_loMassFraction.resize(numSpecies);
  m_loMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  int massCheck = assignMassFractions(m_loMassFraction,
                                      "lower_specs",
                                      "lower_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (CombustionTest IBC): 'lower_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }
  m_upMassFraction.resize(numSpecies);
  m_upMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  massCheck = assignMassFractions(m_upMassFraction,
                                      "upper_specs",
                                      "upper_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (CombustionTest IBC): 'upper_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }
  ppIBC.query("grav_force", m_grav);
  ppIBC.query("upper_rho", m_upRho);
  ppIBC.query("lower_rho", m_loRho);
  if (m_upRho < 0. || m_loRho < 0.)
    {
      CRD::msg << "Input (CombustionTest IBC): Initial upper or lower rho must "
        "greater than 0." << CRD::error;
    }
  ppIBC.query("init_pressure", m_pressure);
  m_readInput = true;
}
