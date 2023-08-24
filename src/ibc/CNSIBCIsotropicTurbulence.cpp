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
 * \file CNSIBCIsotropicTurbulence.cpp
 *
 * \brief Member functions for CNSIBCIsotropicTurbulence
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

//----- Internal -----//

#include "CNSIBCIsotropicTurbulence.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "PatchMappedFunc.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodVMS.H"

/*******************************************************************************
 *
 * Class CNSIBCIsotropicTurbulence: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCIsotropicTurbulence::CNSIBCIsotropicTurbulence()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_P0(101325.),
  m_rho0(-1.),
  m_T0(300.),
  m_Einit(-1.)
{
  CNSIBCIsotropicTurbulence::readBCInfo();
  if(SpaceDim < 3)
    {
      CRD::msg << "Problem not defined for less than 3D!" << CRD::abort;
    }

  // BC are all periodic by default
  // setAllDomainBC(CRDparam::DomainBCTypePeriodic);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCIsotropicTurbulence::~CNSIBCIsotropicTurbulence()
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
CNSIBCIsotropicTurbulence::IBCName() const
{
  return "Homogeneous Isotropic Turbulence (HIT) case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCIsotropicTurbulence::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Domain-mean density\n"        << m_rho0  << CRD::var;
  CRD::msg << "Domain-mean pressure\n"       << m_P0      << CRD::var;
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\>.
 *  \param[inout] a_U   State on the level to be initialized in this
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
CNSIBCIsotropicTurbulence::initialize(LevelData<FArrayBox>&      a_U,
                                      LevelGridMetrics&          a_gridMetrics,
                                      const LayoutData<FluxBox>& a_unitNormals,
                                      const Real                 a_time,
                                      const int                  a_level) const
{
  CRD::msg << CRD::fv2 << "CNSIBCIsotropicTurbulence::initialize" << CRD::end;

  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int preIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int ncells = CRDparam::g_domainBaseSize.product();

  // -----------------------------------------------------------------
  // Calculate the unscaled kinetic energy of the random initial condition
  Real localMeanKE = 0.;
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_U.disjointBoxLayout()[dit];
      FArrayBox& UFab = a_U[dit];
      MD_ARRAY_RESTRICT(arrayU, UFab);
      MD_BOXLOOP(box, i)
        {
          Real& u0 = arrayU[MD_IX(i, velIndx)];
          Real& u1 = arrayU[MD_IX(i, velIndx + 1)];
          Real& u2 = arrayU[MD_IX(i, velIndx + 2)];
          localMeanKE += u0*u0 + u1*u1 + u2*u2;
        }
    }

  Real globalMeanKE = localMeanKE;
#ifdef CH_MPI
  MPI_Allreduce(&localMeanKE, &globalMeanKE, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  globalMeanKE /= ncells;

  // -----------------------------------------------------------------
  // Apply the random IC with the correct kinetic energy
  const Real scale = std::sqrt(m_Einit / (0.5 * m_rho0 * globalMeanKE));
  CRD::msg << CRD::fv2 << "scale factor \n" << scale << CRD::var;

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box& box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2dom = grow(box, 2);
      box2dom &= blockDomain;
      Box box1dom = grow(box, 1);
      box1dom &= blockDomain;

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      CH_assert(UFab.box().contains(box2dom));

      FABSTACKTEMP(WFab, box2dom, numWVar);
      WFab.setVal(0.);

      // Fill WFab with data
      MD_ARRAY_RESTRICT(arrayW, WFab);
      MD_ARRAY_RESTRICT(arrayU, UFab);
      MD_BOXLOOP(box2dom, i)
        {
          arrayW[MD_IX(i, rhoIndx)] = m_rho0;
          arrayW[MD_IX(i, velIndx)] = scale * arrayU[MD_IX(i, velIndx)];
          arrayW[MD_IX(i, velIndx+1)] = scale * arrayU[MD_IX(i, velIndx+1)];
          arrayW[MD_IX(i, velIndx+2)] = scale * arrayU[MD_IX(i, velIndx+2)];
          arrayW[MD_IX(i, preIndx)] = m_P0;
        }

      // overwrite U using W and then convert to cell average
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         WFab,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2dom);
      fourthOrderAverageCell(UFab, blockDomain, box1dom);
    }
}

/*-------------------------------------------------------------------------*//**
 *  \brief no-op source term for now
 *
 *  \param[out] a_sourceFab
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
 *//*-------------------------------------------------------------------------*/
void
CNSIBCIsotropicTurbulence::addSourceTerm(
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
  return;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCIsotropicTurbulence::haveExactSol() const
{
  return false;
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
CNSIBCIsotropicTurbulence::readBCInfo()
{
  ParmParse ppIBC("ibc");
  const Real R = CRDparam::g_R;
  ppIBC.query("init_pressure", m_P0);
  int readT = ppIBC.query("init_temperature", m_T0);
  ppIBC.query("init_density", m_rho0);
  if (m_rho0 > 0.0) {
    if (readT == 1) {
      m_P0 = m_rho0*R*m_T0;
    }
    else {
      m_T0 = m_P0/(R*m_rho0);
    }
  }
  else {
    m_rho0 = m_P0/(R*m_T0);
  }
  ppIBC.query("init_kinetic_energy", m_Einit);
  if (m_Einit < 0.0) {
    CRD::msg << "initial kinetic energy not defined!" << CRD::abort;
  }
  m_readInput = true;
}
