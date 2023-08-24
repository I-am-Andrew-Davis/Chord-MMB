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
 * \file CNSIBCShockTube.cpp
 *
 * \brief Member functions for CNSIBCShockTube
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCShockTube.H"
#include "CNSIBCEulerShockBoxF_F.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCShockTube: member definitions
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

CNSIBCShockTube::CNSIBCShockTube()
:
  CNSIBC(),
  // All invalid values which must be corrected
  m_rhoL(-1.),
  m_pressureL(-1.),
  m_rhoR(-1.),
  m_pressureR(-1.),
  m_threshold(-1.),
  m_tubeDir(-1)
{

//--Read any BC info

  readBCInfo();

//--Define states
//**FIXME Fix input to read a state rather than specific variables

  std::string stateName = "zero";
  if (CRDState::nameIndex(stateName) == -1) // Define with default values
    {
      CRDState& state = CRDState::get(stateName);
      (void)state;      
    }
  CRDState::writeStateInfo();

//--Set BC Type

// Normally you would set everything to wall, but to test corner conditions for
// Euler solvers, it is useful to set directions orthogonal to the tube
// direction as periodic (should make no difference)
// #if 0
//   setAllDomainBC(CRDparam::DomainBCTypeSlipWall);
// #else
//   setDomainBC(m_tubeDir, 0, CRDparam::DomainBCTypeWall);
//   setDomainBC(m_tubeDir, 1, CRDparam::DomainBCTypeWall);
// #endif

  BoundaryIndex bcIdx;
  for (int block = 0; block != CRDparam::g_coordSys->numBlocks(); block++)
    {
      bcIdx.m_block = block;

      BCInfo bc;
      bc.m_order = 4;

      for (const int dir : EachDir)
        {
          bcIdx.m_dir = dir;
          for (const auto side : EachSide)
            {
              bcIdx.m_side = side;
              bc.m_type = CRDparam::DomainBCTypeSlipWall;
              bc.m_idxState = CRDState::nameIndex("zero");
              bc.m_order = 4;
              setDomainBC(bcIdx, bc);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCShockTube::~CNSIBCShockTube()
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
CNSIBCShockTube::IBCName() const
{
  return "Sod's shock tube";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCShockTube::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CRD::msg << "Density at left\n" << m_rhoL << CRD::var;
  CRD::msg << "Pressure at left\n" << m_pressureL << CRD::var;
  CRD::msg << "Density at right\n" << m_rhoR << CRD::var;
  CRD::msg << "Pressure at right\n" << m_pressureR << CRD::var;
  CRD::msg << "Tagging threshold\n" << m_threshold << CRD::var;
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
CNSIBCShockTube::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  // Tag density gradient
  tagLevel->appendTagMethod(
    new TagMethodGradient(CRDparam::g_CRDPhysics->densityIndex(),
                          m_threshold,
                          true));
  // Add in tag buffer
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  // Return the factory
  return new TagLevelFactory(tagLevel);
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
CNSIBCShockTube::initialize(LevelData<FArrayBox>&      a_U,
                            LevelGridMetrics&          a_gridMetrics,
                            const LayoutData<FluxBox>& a_unitNormals,
                            const Real                 a_time,
                            const int                  a_level) const
{
  const Real gamma      = CRDparam::g_gamma;
  const int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg);

  // This box defines the left side.  It grown by quite a bit outside the
  // domain.
  RealVect bxLo = CRDparam::g_domainOrigin - CRDparam::g_domainLength;
  RealVect bxHi = CRDparam::g_domainOrigin + 2.*CRDparam::g_domainLength;
  bxHi[m_tubeDir] = CRDparam::g_domainOrigin[m_tubeDir] +
    0.5*CRDparam::g_domainLength[m_tubeDir];

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Data
      FArrayBox& UFab = a_U[dit];

      // Working set boxes
      Box box2Dom = grow(box, numGhosts);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, numGhosts-1);
      box1Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FORT_CNSIBCEULERSHOCKBOXINIT(CHF_FRA(UFab),
                                   CHF_BOX(box2Dom),
                                   CHF_CONST_FRA(XFab),
                                   CHF_CONST_REAL(m_rhoR),
                                   CHF_CONST_REAL(m_pressureR),
                                   CHF_CONST_REAL(gamma),
                                   CHF_CONST_REAL(m_rhoL),
                                   CHF_CONST_REAL(m_pressureL),
                                   CHF_CONST_REALVECT(bxLo),
                                   CHF_CONST_REALVECT(bxHi));

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box1Dom));
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
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
CNSIBCShockTube::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error
  ParmParse ppIBC("ibc");

//--Tube direction

  m_tubeDir = 0;
  ppIBC.query("tube_dir", m_tubeDir);
  if (m_tubeDir < 0 || m_tubeDir >= SpaceDim)
    {
      CRD::msg << "Input (ShockTube IBC): 'tube_dir' must be >= 0 and < "
               << SpaceDim << '!' << CRD::error;
    }

//--Left state

  m_rhoL = 1.0;
  m_pressureL = 1.0;
  ppIBC.query("density_left", m_rhoL);
  if (m_rhoL <= 0.)
    {
      CRD::msg << "Input (ShockTube IBC): 'density_left' must be > 0!"
               << CRD::error;
    }
  ppIBC.query("pressure_left", m_pressureL);
  if (m_pressureL <= 0.)
    {
      CRD::msg << "Input (ShockTube IBC): 'pressure_left' must be > 0!"
               << CRD::error;
    }

//--Right state

  m_rhoR = 0.125;
  m_pressureR = 0.1;
  ppIBC.query("density_right", m_rhoR);
  if (m_rhoR <= 0.)
    {
      CRD::msg << "Input (ShockTube IBC): 'density_right' must be > 0!"
               << CRD::error;
    }
  ppIBC.query("pressure_right", m_pressureR);
  if (m_pressureR <= 0.)
    {
      CRD::msg << "Input (ShockTube IBC): 'pressure_right' must be > 0!"
               << CRD::error;
    }

//--Threshold of relative density gradient for refinement

  m_threshold = 0.15;
  ppIBC.query("tag_threshold", m_threshold);
  if (m_threshold < 0.)
    {
      CRD::msg << "Input (ShockTube IBC): 'tag_threshold' must be > 0!"
               << CRD::error;
    }

  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
