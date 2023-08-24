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
 * \file CNSIBCGeneralized.cpp
 *
 * \brief Member functions for CNSIBCGeneralized
 *
 *//*+*************************************************************************/

#include <algorithm>

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
 
//----- Internal -----//

#include "CNSIBCGeneralized.H"
#include "CNSIBCCombustionTestF_F.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"


/*******************************************************************************
 *
 * Class CNSIBCGeneralized: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCGeneralized::CNSIBCGeneralized(const bool a_readInputs)
  :
  CNSIBCCombustionReference(),
  m_idxStateInit(-1), // invalid index
  m_bodyForce(RealVect(D_DECL(0., -9.8, 0.))) // Gravity body force
{
  // Default for all domain BC is 4th order periodic, as set in CNSIBC
  
  // read in all other info
  if (a_readInputs)
    {
      readBCInfo(); 
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCGeneralized::~CNSIBCGeneralized()
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
CNSIBCGeneralized::IBCName() const
{
  return "Generalized case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralized::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  CNSIBCCombustionReference::writeTagInfo();
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
CNSIBCGeneralized::initialize(LevelData<FArrayBox>&      a_U,
                              LevelGridMetrics&          a_gridMetrics,
                              const LayoutData<FluxBox>& a_unitNormals,
                              const Real                 a_time,
                              const int                  a_level) const
{
  // Set the initial values
  const int numWVar      = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx     = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx     = CRDparam::g_CRDPhysics->temperatureIndex();
  int speciesComp = 0;
  // if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
  //   {
  //     speciesComp = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  //     if (m_initT*m_initRho*m_initP > 0.)
  //       {
  //         CRD::msg << "Initial values of density, pressure, or temperature"
  //                  << " are not properly specified!" << CRD::error;
  //       }
  //   }
  const int numSpecies = CRDparam::g_numSpecies;
  const CRDState& state = CRDState::get(m_idxStateInit);
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box4Dom = grow(box, 4);
      box4Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      // Pointwise values of W (required on valid +4, for turbInitialize, except
      // at physical boundaries)
      FABSTACKTEMP(Wc, box4Dom, numWVar);
      Wc.setVal(state.density(), rhoIndx);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Wc.setVal(state.velocity()[dir], box4Dom, WvelIndx + dir, 1);
        }
      Wc.setVal(state.pressure(), presIndx);
      Wc.setVal(state.temperature(), tempIndx);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          CH_assert(false);  //**FIXME
          for (int i = 0; i != numSpecies; ++i)
            {
              int comp = speciesComp + i;
              Wc.setVal(state(comp), comp);
            }
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralized::addSourceTerm(FArrayBox&           a_sourceFab,
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
  const int momIndx = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();
  const int engIndx = CRDparam::g_CRDPhysics->energyFluxIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  MD_ARRAY_RESTRICT(arrSource, a_sourceFab);
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_BOXLOOP(a_solveBox, i)
    {
      Real rho = arrWcell[MD_IX(i, rhoIndx)];
      RealVect vel(D_DECL(arrWcell[MD_IX(i, velIndx)],
                          arrWcell[MD_IX(i, velIndx + 1)],
                          arrWcell[MD_IX(i, velIndx + 2)]));
      // Dot product of f.vel
      Real fDotVel = 0.;
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          fDotVel += m_bodyForce[dir]*vel[dir];
          arrSource[MD_IX(i, momIndx + dir)] = rho*m_bodyForce[dir];
        }
      arrSource[MD_IX(i, engIndx)] = rho*fDotVel;
    }
  return;
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
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralized::setImposedBCprimState(
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
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int cVel   = velIntv.begin();
  const int cRho   = CRDparam::g_CRDPhysics->densityIndex();
  const int cPres  = CRDparam::g_CRDPhysics->pressureIndex();
  const int cTemp   = CRDparam::g_CRDPhysics->temperatureIndex();
  const int numPrim = CRDparam::g_CRDPhysics->numPrimitive();

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    case CRDparam::DomainBCTypeFarfield:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      for (int comp = 0; comp != numPrim; ++comp)
        {
          a_Wface.setVal(state(comp), a_boundaryFaceBox, comp);
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      const Real alpha = a_bcInfo.m_relaxCBCStateParam; // 1/0 = non/reflecting
      for (int c = 0; c != numPrim; ++c)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              const Real interior = a_Wcell[MD_IX(i, c)];
              const Real exterior = state(c);
              a_Wface[MD_IX(i, c)] = alpha*interior + (1. - alpha)*exterior;
            }
        }
      break;
    }
    case  CRDparam::DomainBCTypeOutflow:
    {
      // Transform velocity into normal and tangent components
      FABSTACKTEMP(nttVel, a_boundaryFaceBox, velIntv.size());
      nttVel.copy(a_Wface, a_boundaryFaceBox, cVel,
                  a_boundaryFaceBox, 0, velIntv.size());
      PatchMappedFunc::forwardTransform(nttVel,
                                        a_unitNormalBasisFab,
                                        a_boundaryFaceBox);
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          //**FIXME multispecies
          const Real gamma = CRDparam::g_CRDPhysics->gamma();
          //**FIXME Probably need to replace use of cp with h
          const Real cp    = CRDparam::g_CRDPhysics->cp();
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas();
          // Interior estimate of velocity
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, cVel + idxVel)];
            }
          const RealVect intVel1(intVel - state.m_frameVelocity);
          const Real intVel1Sq = stc::dot(intVel1, intVel1);
          Real intT;
          // Temperature state based on characteristics
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_bcIdx.m_side) < 0.)  // Inflow
            {
              intT = state.temperature() - intVel1Sq/(2*cp);
            }
          else                                                      // Outflow
            {
              intT = a_Wface[MD_IX(i, cTemp)];
            }
          const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
          // Whether inflow or outflow, this is the assigned pressure.  Here,
          // assume a frozen mixture with gamma from the interior.
          const Real p = state.pressure()*std::pow(intc1, -gamma/(gamma - 1.0));
          a_Wface[MD_IX(i, cPres)] = p;
          // Temperature and density assigned only if inflow
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_bcIdx.m_side) < 0.)  // Inflow
            {
              const Real T = state.temperature()/intc1;
              a_Wface[MD_IX(i, cRho)]  = p/(Rgas*T);
              a_Wface[MD_IX(i, cTemp)] = T;
            }
        }
      break;
    }
    case  CRDparam::DomainBCTypeRelaxedCBCIn:
    case  CRDparam::DomainBCTypeInflow:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          //**FIXME multispecies
          const Real gamma = CRDparam::g_CRDPhysics->gamma();
          //**FIXME Probably need to replace use of cp with h
          const Real cp    = CRDparam::g_CRDPhysics->cp();
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas();
          const RealVect extVel = state.velocity();
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, cVel + idxVel)];
              a_Wface[MD_IX(i, cVel + idxVel)] = extVel[idxVel];
            }
          if (CRDparam::DomainBCTypeRelaxedCBCIn & a_bcInfo.m_type)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_Wface[MD_IX(i, comp)] = state(comp);
                }
            }
          else
            {
              // Pressure is determined from the interior.  First, find a
              // stagnation pressure for the given reference frame.  Assume
              // gas is frozen at upstream conditions.
              const RealVect intVel1(intVel - state.m_frameVelocity);
              const Real intVel1Sq = stc::dot(intVel1, intVel1);
              const Real intT = a_Wface[MD_IX(i, cTemp)];
              const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
              const Real p0 = a_Wface[MD_IX(i, cPres)]*
                std::pow(intc1, gamma/(gamma - 1.));
              // Find imposed pressure from stagnation pressure
              const RealVect extVel1(extVel - state.m_frameVelocity);
              const Real extVel1Sq = stc::dot(extVel1, extVel1);
              const Real extc1 = 1 + extVel1Sq/(2*cp*state.temperature());
              a_Wface[MD_IX(i, cPres)] =
                p0*std::pow(extc1, -gamma/(gamma - 1.0));
              // Inflow cannot become outflow, so we always set density and
              // temperature
              a_Wface[MD_IX(i, cRho)] = state.density();
              a_Wface[MD_IX(i, cTemp)] = a_Wface[MD_IX(i, cPres)]/
                (state.density()*Rgas);
            }
        }
      break;
    }
    default:
      CH_assert(false);
    }
}

/*--------------------------------------------------------------------*/
//  Fill in the profiles for the CNSCBC boundary condition.
//  Inflow - Must specify velocity, temperature, and c_n
//  Outflow - Must specify pressure
/** \param[out] a_BCProfile
 *                      Profile of boundary values
 *  \param[in]  a_boundaryBox
 *                      Box of locations at the boundary
 *  \param[in]  a_Wcell Cell-averaged primitive variables
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Normal direction to the boundary
 *  \param[in]  a_side  LoHi side of boundary
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current CNSCBC boundary condition
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralized::setCNSCBCProfiles(
  FArrayBox&                    a_BCProfile,
  const Box&                    a_boundaryBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_domT) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const CRDState& bcState = CRDState::get(a_domT.m_idxState);

  if (CRDparam::DomainBCTypeCNSCBCInflow & a_domT.m_type)
    {
      CH_assert(a_BCProfile.nComp() == CRDparam::g_CRDPhysics->numPrimitive());
      const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          a_BCProfile.setVal(bcState(velIndx+dir), a_boundaryBox, velIndx+dir, 1);
        }
      const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
      a_BCProfile.setVal(bcState(tempIndx), a_boundaryBox, tempIndx, 1);
      for (int i = 0; i != numSpecies; ++i)
        {
          const int speciesComp =
            CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
          a_BCProfile.setVal(bcState(speciesComp+i), a_boundaryBox,
                             speciesComp+i, 1);
        }
    }
  else if (CRDparam::DomainBCTypeCNSCBCOutflow & a_domT.m_type)
    {
      const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
      a_BCProfile.setVal(bcState(presIndx), a_boundaryBox, presIndx, 1);
    }
}

/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In general, one should not need to override this routine except
 *  in special circumstances such as for mixed boundaries.
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state corrected for presence of wall
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralized::setWallBCprimState(
  FArrayBox&           a_Wface,
  const Box&           a_boundaryFaceBox,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_boundarySlipVelocity,
  const FArrayBox&     a_unitNormalBasisFab,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  const BCInfo&        a_domT) const
{
  const int lohiSign = sign(a_bcIdx.m_side);
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  const CRDState& bcState = CRDState::get(a_domT.m_idxState);
  RealVect wallVel = bcState.velocity();

  int viscousSlip = 0;
  if (CRDparam::DomainBCTypeSlipWall & a_domT.m_type)
    {
      viscousSlip = 1;
      wallVel = RealVect::Zero;

    }
  else if (CRDparam::DomainBCTypeIsothermalWall & a_domT.m_type)
    {
      Real Temp = bcState.temperature();
      CH_assert(Temp > 0.);
      a_Wface.setVal(Temp, a_boundaryFaceBox, tempIndx);
    }
  
  // Gamma values in cells
  FABSTACKTEMP(gammaCell, a_boundaryFaceBox, 1);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaCell,
                                    a_Wcell);
  computeWallPrimState(a_Wface,
                       a_boundaryFaceBox,
                       a_Wcell,
                       a_boundarySlipVelocity,
                       gammaCell,
                       a_unitNormalBasisFab,
                       wallVel,
                       viscousSlip,
                       a_bcIdx.m_dir,
                       lohiSign);
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
CNSIBCGeneralized::readBCInfo()
{
  CH_assert(!m_readInput);

//--Read the initial state

  ParmParse ppIBC("ibc");

  // State for intitialization
  const char* initStateNameHdr = "initial_state";
  if (ppIBC.contains(initStateNameHdr))
    {
      std::string initStateName;
      ppIBC.get(initStateNameHdr, initStateName);
      m_idxStateInit = CRDState::nameIndex(initStateName);
    }
  else
    {
      CRD::msg << "Input (Generalized IBC): " << initStateNameHdr
               << " must be given" << CRD::error;
    }

  // FIXME Re-design this in some way
  // if ((m_initT*m_initP*m_initRho > 0.) || (m_initT+m_initP+m_initRho < 0))
  //   {
  //     CRD::msg << "Input (Generalized IBC): Make sure only 2 of the 3 "
  //              << "are specified for 'init_pressure', 'init_density', and "
  //              << "'init_temperature'!"
  //              << " They are currently set as \n"
  //              << "\tpressure = " << m_initP
  //              << "\tdensity = " << m_initRho
  //              << "\ttemperature = " << m_initT << CRD::error;
  //   }

  // Check the solution multispecies physics initialization
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      // Call function to assign mass fraction values, this check total mass
      // FIXME Re-enable this
      // const CRDState& initState = CRDState::get(m_idxStateInit);
      // const Interval speciesIntv = CRDparam::g_CRDPhysics->speciesPrimInterval();
      // int massCheck = assignMassFractions(initState(speciesIntv),
      //                                     "init_specs",
      //                                     "init_mfs");
      // if (massCheck == 1)
      //   {
      //     CRD::msg << "Input (Generalized IBC): " << initStateNameHdr
      //              << "must have mass fractions equal to 1." << CRD::error;
      //   }
    }
  
  // Setup the solution uses body force source terms
  if (CRDparam::g_physicsModels & CRDparam::PhysicsSource)
    {
      std::vector<Real> inputBodyForce(SpaceDim, 0.);
      ppIBC.queryarr("body_force", inputBodyForce, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_bodyForce.dataPtr(),
                                               &inputBodyForce.front());
    }

  std::vector<int> periodicity(SpaceDim, 1);
  ppIBC.queryarr("periodicity", periodicity, 0, SpaceDim);

  // --- Read each boundary condition ---
  const std::string typeStr = "bc_type-";
  const std::string stateStr = "bc_state-";
  const std::string orderStr = "bc_order-";
  const std::string stateRelaxStr = "bc_state_relax-";
  const std::string waveRelaxStr = "bc_wave_relax-";
  const std::array<std::string, 3> dirVecStr {"x", "y", "z"};
  const std::array<std::string, 2> sideVecStr {"lo", "hi"};
  // set up default values
  BCInfo bcDefault;
  std::string defStr = "default";
  std::string defT = typeStr + defStr;
  std::string defaultTypeStr;
  ppIBC.query(defT.c_str(), defaultTypeStr);
  std::string defO = orderStr + defStr;
  ppIBC.query(defO.c_str(), bcDefault.m_order);
  std::string defS = stateStr + defStr;
  std::string bcStateNameDef;
  ppIBC.query(defS.c_str(), bcStateNameDef);
  std::string defWR = waveRelaxStr + defStr;
  std::string defSR = stateRelaxStr + defStr;
  ppIBC.query(defWR.c_str(), bcDefault.m_relaxCBCWaveParam);
  ppIBC.query(defSR.c_str(), bcDefault.m_relaxCBCStateParam);
  // set all defaults. This ensures they can be looped over
  setAllDomainBC(bcDefault);
  for (const auto& bcStatePair : m_blockBCInfo)
    {
      // read the state
      const auto& bcIdx = bcStatePair.first;
      BCInfo bc = bcStatePair.second;
      //BlockInfo bInfo = ;
      std::string idxStr = CRDparam::g_coordSys->blockInfo(bcIdx.m_block).m_name
        + "_" + dirVecStr.at(bcIdx.m_dir)
        + "_" + sideVecStr.at(bcIdx.m_side);
      // Takes blockinfo based on Block number *The user should never
      // Have to know this information

    //  std::string idxStr = "block_" + std::to_string(bcIdx.m_block)
    //    + "_" + dirVecStr.at(bcIdx.m_dir)
    //    + "_" + sideVecStr.at(bcIdx.m_side);
      std::string strType = typeStr + idxStr;
      std::string bcTypeStr = defaultTypeStr;
      ppIBC.query(strType.c_str(), bcTypeStr);
      std::string strOrder = orderStr + idxStr;
      ppIBC.query(strOrder.c_str(), bc.m_order);
      std::string strState = stateStr + idxStr;
      std::string bcStateName = bcStateNameDef;
      ppIBC.query(strState.c_str(), bcStateName);
      //bc.m_stateIdx = stateNameIdx.at(bcStateName);
      bc.m_idxState = CRDState::nameIndex(bcStateName);
      std::string strStateRelax = stateRelaxStr + idxStr;
      ppIBC.query(strStateRelax.c_str(), bc.m_relaxCBCStateParam);
      std::string strWaveRelax = waveRelaxStr + idxStr;
      ppIBC.query(strWaveRelax.c_str(), bc.m_relaxCBCWaveParam);

      bool bcSet = false;
      if (iequals(bcTypeStr, "IsothermalWall"))
        {
          bc.m_type = CRDparam::DomainBCTypeIsothermalWall;
          bcSet = true;
        }
      else if (iequals(bcTypeStr, "Wall") || iequals(bcTypeStr, "NoSlipWall"))
        {
          bc.m_type = CRDparam::DomainBCTypeWall;
          bcSet = true;
        }
      else if (iequals(bcTypeStr, "SlipWall"))
        {
          bc.m_type = CRDparam::DomainBCTypeSlipWall;
          bcSet = true;
        }
      else if (iequals(bcTypeStr, "AdiabaticWall"))
        {
          bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
          bcSet = true;
        }
      else if (iequals(bcTypeStr, "Extrapolated"))
        {
          bc.m_type = CRDparam::DomainBCTypeExtrapolated;
          bcSet = true;
        }
      else
        {
          if (iequals(bcTypeStr, "CNSCBCInflow"))
            {
              bc.m_type = CRDparam::DomainBCTypeCNSCBCInflow;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "RelaxedCBCIn"))
            {
              bc.m_type = CRDparam::DomainBCTypeRelaxedCBCIn;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "RelaxedCBCOut"))
            {
              bc.m_type = CRDparam::DomainBCTypeRelaxedCBCOut;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "RelaxedCBCFar"))
            {
              bc.m_type = CRDparam::DomainBCTypeRelaxedCBCFar;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "Dirichlet"))
            {
              bc.m_type = CRDparam::DomainBCTypeDirichlet;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "Inflow"))
            {
              bc.m_type = CRDparam::DomainBCTypeInflow;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "CNSCBCOutflow"))
            {
              bc.m_type = CRDparam::DomainBCTypeCNSCBCOutflow;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "Outflow"))
            {
              bc.m_type = CRDparam::DomainBCTypeOutflow;
              bcSet = true;
            }
          else if (iequals(bcTypeStr, "Farfield"))
            {
              bc.m_type = CRDparam::DomainBCTypeFarfield;
              bcSet = true;
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
            {
              // const auto& bcStateMS = getDomainBCstate(bc);
              // bcStateMS.m_MF.resize(numSpecies);
              // bcStateMS.m_MF.assign(numSpecies,0.);
              // // Call function to assign mass fraction values
              // std::string specStr = lohiVecStr[bcIdx.m_side]
              //   + "_bc_specs_"+ dirVecStr[bcIdx.m_dir];
              // std::string mfStr = lohiVecStr[bcIdx.m_side]
              //   + "_bc_mfs_" + dirVecStr[bcIdx.m_dir];
              // int massCheck = assignMassFractions(bcStateMS.m_MF,
              //                                     specStr,
              //                                     mfStr);
              // if (massCheck == 1)
              //   {
              //     CRD::msg << "Input (Generalized IBC): 'lo_bc_mfs' "
              //              << "must be equal to 1." << CRD::error;
              //   }
            }
        }
      if (!bcSet)
        {
          CRD::msg << "Input (Generalized IBC): Unrecognized bc "
                   << bcTypeStr << CRD::error;
        }
      if (iequals(bcTypeStr, "Farfield"))
        {
          bc.m_type = CRDparam::DomainBCTypeFarfield;
          //const auto& bcStateFF = getDomainBCstate(bc);
          const CRDState& bcStateFF = CRDState::get(bc.m_idxState);
          if (bcStateFF.temperature() < 0. || bcStateFF.pressure() < 0.)
            {
              CRD::msg << "Input (Generalized IBC): must specify "
                       << "temperature and pressure for farfield "
                       << "condition!" << CRD::error;
            } 
        }
      // create the boundary condition
      setDomainBC(bcIdx, bc);
    }
  m_readInput = true;
}
