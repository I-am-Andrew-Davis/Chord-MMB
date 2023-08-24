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
 * \file CNSIBCVortex.cpp
 *
 * \brief Member functions for CNSIBCVortex
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCVortex.H"
#include "CNSIBCVortexF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"


/*******************************************************************************
 *
 * Class CNSIBCVortex: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCVortex::CNSIBCVortex()
  :
  CNSIBCReferenceCubeBC(),
  m_beta(0.),
  m_U(RealVect::Zero),
  m_Pinf(101325.),
  m_Tinf(300.)
{
  readBCInfo();
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }
  BCInfo bc;
  bc.m_order = 1;
  bc.m_type = CRDparam::DomainBCTypeCNSCBCOutflow;
  setAllDomainBC(bc);

  BoundaryIndex bcIdx;
  bcIdx.m_block = 0;
  bcIdx.m_dir = 0;
  bcIdx.m_side = Side::Lo;
  bc.m_type = CRDparam::DomainBCTypeCNSCBCInflow;
  setDomainBC(bcIdx, bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCVortex::~CNSIBCVortex()
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
CNSIBCVortex::IBCName() const
{
  return "Vortex case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCVortex::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Vortex strength\n" << m_beta << CRD::var;
  CRD::msg << "Freestream velocity\n" << m_U << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg << "Outflow transverse smoothing beta\n" << m_CBCbeta << CRD::var;
  CRD::msg << "Outflow transverse smoothing sigma\n" << m_CBCsigma << CRD::var;
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
CNSIBCVortex::initialize(LevelData<FArrayBox>&      a_U,
                         LevelGridMetrics&          a_gridMetrics,
                         const LayoutData<FluxBox>& a_unitNormals,
                         const Real                 a_time,
                         const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
  RealVect center(RealVect::Unit);
  for (int i = 0; i != SpaceDim; ++i)
    {
      center[i] = CRDparam::g_domainLength[i]/2. + CRDparam::g_domainOrigin[i];
    }
  const int numSpecies = CRDparam::g_numSpecies;
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

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;
      
      // Get coordinate system for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
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
      Real gasR = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          gasR += m_massFractions[sp]*
            CRDparam::g_CRDPhysics->speciesGasConstant(sp);
        }
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      FORT_CNSIBCVORTEXINITCN(CHF_FRA(Wc),
                              CHF_BOX(box2Dom),
                              CHF_CONST_INT(rComp),
                              CHF_CONST_INT(pComp),
                              CHF_CONST_INT(tComp),
                              CHF_CONST_INT(WvelIndx),
                              CHF_CONST_INT(wCompStart),
                              CHF_CONST_INT(numSpecies),
                              CHF_CONST_FRA(XFab),
                              CHF_CONST_REAL(m_Pinf),
                              CHF_CONST_REAL(m_Tinf),
                              CHF_CONST_REAL(m_beta),
                              CHF_CONST_REAL(gasR),
                              CHF_CONST_REALVECT(center),
                              CHF_CONST_REALVECT(m_U),
                              CHF_CONST_VR(m_massFractions));
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
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCVortex::haveExactSol() const
{
  return false;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Fill in the profiles for the CNSCBC boundary condition.
//  Inflow - Must specify velocity, temperature, and cn
//  Outflow - Must specify pressure
/** \param[out] a_BCProfile
 *                      Profile of boundary values
 *  \param[in]  a_boundaryBox
 *                      Box of locations at the boundary
 *  \param[in]  a_Wcell Cell-averaged primitive variables
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces
 *  \param[in]  a_dir   Normal direction to the boundary
 *  \param[in]  a_side  LoHi side of boundary
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current CNSCBC boundary condition
 *//*-----------------------------------------------------------------*/

void
CNSIBCVortex::setCNSCBCProfiles(
  FArrayBox&                    a_BCProfile,
  const Box&                    a_boundaryBox,
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
  if (CRDparam::DomainBCTypeCNSCBCInflow & a_domT)
    {
      CH_assert(a_BCProfile.nComp() == CRDparam::g_CRDPhysics->numPrimitive());
      const int velComp = CRDparam::g_CRDPhysics->velocityInterval().begin();
      const int tempComp = CRDparam::g_CRDPhysics->temperatureIndex();
      const int numSpecies = CRDparam::g_numSpecies;
      const int wCompStart =
        CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
      const int vComp = velComp;
      a_BCProfile.setVal(0., a_boundaryBox, velComp, SpaceDim);
      a_BCProfile.setVal(m_U[0], a_boundaryBox, vComp, 1);
      a_BCProfile.setVal(m_Tinf, a_boundaryBox, tempComp, 1);
      for (int i = 0; i != numSpecies; ++i)
        {
          a_BCProfile.setVal(m_massFractions[i], a_boundaryBox,
                             wCompStart + i, 1);
        }
    }
  else if (CRDparam::DomainBCTypeCNSCBCOutflow & a_domT)
    {
      const int presComp = CRDparam::g_CRDPhysics->pressureIndex();
      Real outletPressure = m_Pinf;
      a_BCProfile.setVal(outletPressure, a_boundaryBox, presComp, 1);
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
CNSIBCVortex::readBCInfo()
{
  CH_assert(!m_readInput);
  ParmParse ppIBC("ibc");
  m_beta = 0.;
  ppIBC.query("vortex_strength", m_beta);
  m_U = RealVect::Zero;
  std::vector<Real> tU(SpaceDim);
  ppIBC.queryarr("freestream_velocity", tU, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_U.dataPtr(),
                                           &tU.front());
  ppIBC.query("stag_p",m_Pinf);
  ppIBC.query("stag_T",m_Tinf);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      const int numSpecies = CRDparam::g_numSpecies;
      m_massFractions.resize(numSpecies);
      m_massFractions.assign(numSpecies, 0.);
      int massTest = assignMassFractions(m_massFractions,
                                         "init_specs",
                                         "init_mfs");
      if (massTest == 1)
        {
          CRD::msg << "Input (Vortex IBC): 'init_mfs'"
                   << " must be equal to 1!" << CRD::error;
        }
    }
  ppIBC.query("outflow_smoothing_beta", m_CBCbeta);
  ppIBC.query("outflow_smoothing_sigma", m_CBCsigma);
  m_readInput = true;
}
