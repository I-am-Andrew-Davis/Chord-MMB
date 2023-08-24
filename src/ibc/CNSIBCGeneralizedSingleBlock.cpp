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
 * \file CNSIBCGeneralizedSingleBlock.cpp
 *
 * \brief Member functions for CNSIBCGeneralizedSingleBlock
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

#include "CNSIBCGeneralizedSingleBlock.H"
#include "CNSIBCCombustionTestF_F.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"


/*******************************************************************************
 *
 * Class CNSIBCGeneralizedSingleBlock: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCGeneralizedSingleBlock::CNSIBCGeneralizedSingleBlock()
  :
  CNSIBCGeneralized(false) // do no read the default inputs
{
  // read in all info, with different interface than the true
  // generalized case
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCGeneralizedSingleBlock::~CNSIBCGeneralizedSingleBlock()
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
CNSIBCGeneralizedSingleBlock::IBCName() const
{
  return "GeneralizedSingleBlock case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCGeneralizedSingleBlock::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeIBCInfo();
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
CNSIBCGeneralizedSingleBlock::initialize(
  LevelData<FArrayBox>&      a_U,
  LevelGridMetrics&          a_gridMetrics,
  const LayoutData<FluxBox>& a_unitNormals,
  const Real                 a_time,
  const int                  a_level) const
{
  CNSIBCGeneralized::initialize(a_U,
                                a_gridMetrics,
                                a_unitNormals,
                                a_time,
                                a_level);
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
CNSIBCGeneralizedSingleBlock::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  // If the solution uses multispecies physics
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      // Specify the initial mass fractions in the flow field
      m_initMassFraction.resize(numSpecies);
      m_initMassFraction.assign(numSpecies,0.);
      // Call function to assign mass fraction values
      int massCheck = assignMassFractions(m_initMassFraction,
                                          "init_specs",
                                          "init_mfs");
      if (massCheck == 1)
        {
          CRD::msg << "Input (GeneralizedSingleBlock IBC): 'init_mfs' "
                   << "must be equal to 1." << CRD::error;
        }
    }
  // If the solution uses body force source terms
  if (CRDparam::g_physicsModels & CRDparam::PhysicsSource)
    {
      std::vector<Real> inputBodyForce(SpaceDim, 0.);
      ppIBC.queryarr("body_force", inputBodyForce, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_bodyForce.dataPtr(),
                                               &inputBodyForce.front());
    }
  std::vector<Real> inputInitVel(SpaceDim, 0.);
  ppIBC.queryarr("init_velocity", inputInitVel, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_initVel.dataPtr(),
                                           &inputInitVel.front());
  ppIBC.query("init_temperature", m_initT);
  ppIBC.query("init_pressure", m_initP);
  ppIBC.query("init_density", m_initRho);
  if ((m_initT*m_initP*m_initRho > 0.) || (m_initT+m_initP+m_initRho < 0))
    {
      CRD::msg << "Input (GeneralizedSingleBlock IBC): Make sure only 2 of the 3 "
               << "are specified for 'init_pressure', 'init_density', and "
               << "'init_temperature'!"
               << " They are currently set as \n"
               << "\tpressure = " << m_initP
               << "\tdensity = " << m_initRho
               << "\ttemperature = " << m_initT << CRD::error;
    }
  std::vector<int> periodicity(SpaceDim, 1);
  ppIBC.queryarr("periodicity", periodicity, 0, SpaceDim);
  // Loop over each direction, labeled with x, y, z
  const std::array<std::string, 3> dirVecStr {"x", "y", "z"};
  const std::array<std::string, 2> sideVecStr {"lo", "hi"};
  
  // single block only
  BoundaryIndex bcIdx;
  bcIdx.m_block = 0;
  // boundary state index
  //int bcStateIdx = 0;
  // for (const auto dir : EachDir)
  //   {
  //     bcIdx.m_dir = dir;
  //     if (periodicity[dir] == 1)
  //       {
  //         continue;
  //       }
  //     for (const auto side : EachSide)
  //       {
  //         bcIdx.m_side = side;
  //         // define the boundary state
  //         BCstate bcState;
  //         // Read in bc velocity, temperature, or pressure
  //         std::vector<Real> inputVel(SpaceDim, 0.);
  //         std::string velBCstr = lohiVecStr[side] + "_bc_vel_" + dirVecStr[dir];
  //         ppIBC.queryarr(velBCstr.c_str(), inputVel, 0, SpaceDim);
  //         SpaceDimArray<Real, Real>::loadFromArray(bcState.m_vel.dataPtr(),
  //                                                  &inputVel.front());
  //         std::string bcTemp = lohiVecStr[side] + "_bc_temp_" + dirVecStr[dir];
  //         ppIBC.query(bcTemp.c_str(), bcState.m_temp);
  //         std::string bcPres = lohiVecStr[side] + "_bc_pres_" + dirVecStr[dir];
  //         ppIBC.query(bcPres.c_str(), bcState.m_pres);
          
  //         // create the boundary condition
  //         BCInfo bcType;
  //         //bcType.m_stateIdx = bcStateIdx;
  //         // FIXME use CRDstate
  //         //bcType.m_idxState = CRDState::nameIndex(bcStateName);
  //         bcType.m_order = 4;
  //         std::string BCstr;
  //         std::string lookBC = lohiVecStr[side] + "_bc_" + dirVecStr[dir];
  //         ppIBC.get(lookBC.c_str(), BCstr);
  //         std::string lookBCOrder = lookBC + "_order";
  //         ppIBC.query(lookBCOrder.c_str(), bcType.m_order);
  //         bool bcSet = false;
  //         if (iequals(BCstr, "IsothermalWall"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeIsothermalWall;
  //             bcSet = true;
  //           }
  //         else if (iequals(BCstr, "Wall") || iequals(BCstr, "NoSlipWall"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeWall;
  //             bcSet = true;
  //           }
  //         else if (iequals(BCstr, "SlipWall"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeSlipWall;
  //             bcSet = true;
  //           }
  //         else if (iequals(BCstr, "AdiabaticWall"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  //             bcSet = true;
  //           }
  //         else if (iequals(BCstr, "Extrapolated"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeExtrapolated;
  //             bcSet = true;
  //           }
  //         else
  //           {
  //             if (iequals(BCstr, "CNSCBCInflow"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeCNSCBCInflow;
  //                 bcSet = true;
  //               }
  //             else if (iequals(BCstr, "Dirichlet"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeDirichlet;
  //                 bcSet = true;
  //               }
  //             else if (iequals(BCstr, "Inflow"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeInflow;
  //                 bcSet = true;
  //               }
  //             else if (iequals(BCstr, "CNSCBCOutflow"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeCNSCBCOutflow;
  //                 bcSet = true;
  //               }
  //             else if (iequals(BCstr, "Outflow"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeOutflow;
  //                 bcSet = true;
  //               }
  //             else if (iequals(BCstr, "Farfield"))
  //               {
  //                 bcType.m_type = CRDparam::DomainBCTypeFarfield;
  //                 bcSet = true;
  //               }
  //             if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
  //               {
  //                 bcState.m_MF.resize(numSpecies);
  //                 bcState.m_MF.assign(numSpecies,0.);
  //                 // Call function to assign mass fraction values
  //                 std::string specStr = lohiVecStr[side]
  //                   + "_bc_specs_"+ dirVecStr[dir];
  //                 std::string mfStr = lohiVecStr[side]
  //                   + "_bc_mfs_" + dirVecStr[dir];
  //                 int massCheck = assignMassFractions(bcState.m_MF,
  //                                                     specStr,
  //                                                     mfStr);
  //                 if (massCheck == 1)
  //                   {
  //                     CRD::msg << "Input (GeneralizedSingleBlock IBC): 'lo_bc_mfs' "
  //                              << "must be equal to 1." << CRD::error;
  //                   }
  //               }
  //           }
  //         if (!bcSet)
  //           {
  //             CRD::msg << "Input (GeneralizedSingleBlock IBC): Unrecognized bc "
  //                      << BCstr << CRD::error;
  //           }
  //         if (iequals(BCstr, "Farfield"))
  //           {
  //             bcType.m_type = CRDparam::DomainBCTypeFarfield;
  //             if (bcState.m_temp < 0. || bcState.m_pres < 0.)
  //               {
  //                 CRD::msg << "Input (GeneralizedSingleBlock IBC): must specify "
  //                          << "temperature and pressure for farfield "
  //                          << "condition!" << CRD::error;
  //               }
  //           }
  //         // create the state
  //         // FIXME for CRDstate
  //         //setDomainBCstate(bcStateIdx, bcState);
  //         // increment for a new unique state
  //         ++bcStateIdx;
  //         // create the boundary condition
  //         setDomainBC(bcIdx, bcType);
  //       }
  //   }
  m_readInput = true;
}
