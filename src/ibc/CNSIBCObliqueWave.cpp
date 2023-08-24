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
 * \file CNSIBCMachReflection.cpp
 *
 * \brief Member functions for CNSIBCMachReflection
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"

//----- Internal -----//

#include "CNSIBCObliqueWave.H"
#include "CRDState.H"


/*******************************************************************************
 *
 * Class CNSIBCObliqueWave: member definitions
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

CNSIBCObliqueWave::CNSIBCObliqueWave()
:
  CNSIBC(),
  m_alpha(PI*30./180.),
  m_wallOrder(4)
{
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }

//--Read any BC info

  readBCInfo();

//--Wall state (defined as all zero if not yet set)

  {
    CRDState& wallState = CRDState::get("wall");
    (void)wallState;
  }

//--Set BC Types

  // The state "initial" must exist
  int idxStateInit = CRDState::nameIndex("initial");
  if (idxStateInit == -1)
    {
      CRD::msg << "Oblique wave problem requires state named \"initial\"!"
               << CRD::error;
    }

  // Assumes single block
  if (CRDparam::g_coordSys->numBlocks() != 1)
    {
      CRD::msg << "Oblique wave problem requires single block!" << CRD::error;
    }

  BCInfo bc;

  // Inflow on left side (supersonic so Dirichlet)
  bc.m_type = CRDparam::DomainBCTypeDirichlet;
  bc.m_idxState = idxStateInit;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 0, Side::Lo), bc);
  // Farfield on right side
  bc.m_type = CRDparam::DomainBCTypeFarfield;
  bc.m_idxState = idxStateInit;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 0, Side::Hi), bc);
  // Wall on bottom
  bc.m_type = CRDparam::DomainBCTypeWall;
  bc.m_idxState = CRDState::nameIndex("wall");
  bc.m_order = m_wallOrder;
  setDomainBC(BoundaryIndex(0, 1, Side::Lo), bc);
  // Farfield on top
  bc.m_type = CRDparam::DomainBCTypeFarfield;
  bc.m_idxState = idxStateInit;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 1, Side::Hi), bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCObliqueWave::~CNSIBCObliqueWave()
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
CNSIBCObliqueWave::IBCName() const
{
  return "Oblique wave";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCObliqueWave::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  CRD::msg << "Ramp angle (deg)\n" << m_alpha*180./PI << CRD::var;
  CRD::msg << "x-location of corner\n" << 0.0 << CRD::var;
  CRD::msg << "Order of accuracy at wall BC\n" << m_wallOrder << CRD::var;
  CRD::msg.newline();
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
CNSIBCObliqueWave::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error

//--Geometry

  {
    ParmParse ppCOORDSYS("coordsys");
    ppCOORDSYS.get("ramp_angle", m_alpha);
    // Convert to radians
    m_alpha *= PI/180.;
  }

//--Additional tests on previous input

  {
    const Real dx0 = CRDparam::g_domainLength[0]/CRDparam::g_domainBaseSize[0];
    for (int dir = 1; dir != SpaceDim; ++dir)
      {
        const Real dx =
          CRDparam::g_domainLength[dir]/CRDparam::g_domainBaseSize[dir];
        if (Misc::compare(dx0, dx, std::numeric_limits<Real>::digits10 - 2))
          {
            CRD::msg << "Input (ObliqueWave IBC): 'domain_length' in "
              "direction " << dir << " must be set for same mesh spacing in "
              "all directions.\ndx[0] = " << dx0 << "\ndx[" << dir << "] = "
                     << dx << '!' << CRD::error;
          }
      }
  }
  if (CRDparam::g_domainOrigin != RealVect::Zero)
    {
      CRD::msg << "Input (ObliqueWave IBC): 'domain_origin' must be set to "
        "zero for Mach reflection problem (do not define for valid default)!"
               << CRD::error;
    }

//--IBC

  ParmParse ppIBC("ibc");

//--Order of extrapolation to use at the wall
  
  ppIBC.query("wall_order", m_wallOrder);
  if (m_wallOrder < 0 || m_wallOrder > 4)
    {
      CRD::msg << "Input (ObliqueWave IBC): 'wall_order' must be > 0 and "
               << "<= 4!" << CRD::error;
    }
}
