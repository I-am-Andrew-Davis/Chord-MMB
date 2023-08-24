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
 * \file CNSIBCWedgeDetonation.cpp
 *
 * \brief Member functions for CNSIBCWedgeDetonation
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

#include "CNSIBCWedgeDetonation.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"

/*******************************************************************************
 *
 * Class CNSIBCWedgeDetonation: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** All required parameters are either taken from CRD parameters or
 *  read from input via readBCInfo();
 *//*-----------------------------------------------------------------*/

CNSIBCWedgeDetonation::CNSIBCWedgeDetonation()
:
  CNSIBCGeneralized(),
  //**FIXME These should be read from input once IBC can construct CS
  m_alpha(30.),
  m_xLead(0.12),
  m_xRamp(0.4),
  m_wallTagCells(2),
  m_refAlongWall(true)
{
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }

  //--Read any BC info
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCWedgeDetonation::~CNSIBCWedgeDetonation()
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
CNSIBCWedgeDetonation::IBCName() const
{
  return "Wedge detonation";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCWedgeDetonation::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Angle of ramp (deg)\n" << 180.*m_alpha/PI << CRD::var;
  CRD::msg << "Lead length before ramp\n" << m_xLead << CRD::var;
  CRD::msg << "Ramp length (projected along x-axis)\n" << m_xRamp << CRD::var;
  CRD::msg << "x-location of corner\n" << 0.0 << CRD::var;
  CRD::msg << "Refine along entire wall\n" << m_refAlongWall << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\> on valid +1 except at physical boundaries.
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
CNSIBCWedgeDetonation::initialize(LevelData<FArrayBox>&      a_U,
                                  LevelGridMetrics&          a_gridMetrics,
                                  const LayoutData<FluxBox>& a_unitNormals,
                                  const Real                 a_time,
                                  const int                  a_level) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();

  // Initial conditions
  const CRDState& state = CRDState::get("initial");

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block;
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Data
      FArrayBox& UFab = a_U[dit];

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);
      Wc.setVal(state.temperature(), box2Dom, tempIndx);
      Wc.setVal(state.pressure(), box2Dom, presIndx);
      Wc.setVal(-1, box2Dom, rhoIndx);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Wc.setVal(state.velocity()[dir], box2Dom, velIndx + dir);
        }
      for (int spec = 0; spec != numSpecies; ++spec)
        {
          const int wComp = wCompStart + spec;
          Wc.setVal(state(wComp), box2Dom, wComp);
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box1Dom));
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
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
CNSIBCWedgeDetonation::setImposedBCprimState(
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
  const int cTemp      = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    case CRDparam::DomainBCTypeFarfield:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      for (int comp = 0, comp_end = CRDparam::g_CRDPhysics->numPrimitive();
           comp != comp_end; ++comp)
        {
          a_Wface.setVal(state(comp), a_boundaryFaceBox, comp);
        }
      break;
    }
    case  CRDparam::DomainBCTypeOutflow:
    {
      // We don't need to do anything because the default is 1st order
      // extrapolation from the interior
      break;
    }
    case  CRDparam::DomainBCTypeInflow:
    {
      // We set the supersonic inflow characteristics, that is all states are
      // set at the inflow
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          for (const int idxVel : EachDir)
            {
              a_Wface[MD_IX(i, cVel + idxVel)] = state.velocity()[idxVel];
            }
           a_Wface[MD_IX(i, cPres)] = state.pressure();
           a_Wface[MD_IX(i, cRho)] = state.density();
          a_Wface[MD_IX(i, cTemp)] = state.temperature();
          for (int j = 0; j != numSpecies; ++j)
            {
              a_Wface[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
            }
        }
      break;
    }
    default:
      CH_assert(false);
    }
}


/*--------------------------------------------------------------------*/
//  Read any information related to the IBC from input
/** \note
 *  <ul>
 *    <li> No output should be printed aside from errors
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBCWedgeDetonation::readBCInfo()
{
  // CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error

//--Geometry

  {
    ParmParse ppCOORDSYS("coordsys");
    ppCOORDSYS.get("ramp_angle", m_alpha);
    // Convert to radians
    m_alpha *= PI/180.;
  }

  ParmParse ppIBC("ibc");
//--Ramp values

  // ppIBC.get("ramp_angle", m_alpha);
  // ppIBC.get("lead_length", m_xLead);
  // ppIBC.get("ramp_length", m_xRamp);
  ppIBC.query("wall_tag_cells", m_wallTagCells);
  ppIBC.query("ref_along_wall", m_refAlongWall);
  // Convert to radians
  m_alpha *= PI/180.;
  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
