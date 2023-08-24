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
 * \file CNSIBC.cpp
 *
 * \brief Member functions for CNSIBC
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"

//----- Internal -----//

#include "CNSIBCEulerShockBox.H"
#include "CRDparam.H"
#include "CNSPhysics.H"
#include "ThermPhysics.H"
#include "CRDState.H"
#include "DataTemp.H"

/*******************************************************************************
 *
 * Class CNSIBCEulerShockBox: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCEulerShockBox::CNSIBCEulerShockBox()
:
  CNSIBC(),
  m_bxLo(RealVect::Zero),
  m_bxHi(CRDparam::g_domainLength/2.)
{
  readBCInfo();
  // Default for all BCs
  BCInfo defBC;
  defBC.m_order = 4;
  defBC.m_type = CRDparam::DomainBCTypeSlipWall;
  defBC.m_idxState = CRDState::nameIndex("wall");
  setAllDomainBC(defBC);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCEulerShockBox::~CNSIBCEulerShockBox()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Allocations of new physics states.  Customize for each derivative
//  class.
/**
 *//*-----------------------------------------------------------------*/

std::vector<CRDPhysics*>
CNSIBCEulerShockBox::allocatePhysics()
{
  std::vector<CRDPhysics*> physics(1, nullptr);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      physics[0] = new ThermPhysics;
    }
  else
    {
      physics[0] = new CNSPhysics;
    }
  return physics;
}

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/

const char *const
CNSIBCEulerShockBox::IBCName() const
{
  return "Euler shock box";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCEulerShockBox::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  CRD::msg << "Inner physical box\n" << m_bxLo << ':' << m_bxHi
           << CRD::var;
  CRD::msg.newline();
}

#if 0
/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** The default implementation just sets the buffer and otherwise does
 *  no tagging.
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
CNSIBCEulerShockBox::setTagMethod(const int a_tagBufferSize)
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
#endif

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
CNSIBCEulerShockBox::initialize(LevelData<FArrayBox>&      a_U,
                                LevelGridMetrics&          a_gridMetrics,
                                const LayoutData<FluxBox>& a_unitNormals,
                                const Real                 a_time,
                                const int                  a_level) const
{
  VECSTACKTEMPSIZE(inner, Real, CRDparam::g_CRDPhysics->numPrimitive());
  VECSTACKTEMPSIZE(outer, Real, CRDparam::g_CRDPhysics->numPrimitive());
  const CRDState& innerState = CRDState::get("inner");
  const CRDState& outerState = CRDState::get("outer");

//--Set the primitive state in inner and outer

  inner[CRDparam::g_CRDPhysics->densityIndex()] = innerState.density();
  for (const int idxVel : EachDir)
    {
      inner[CRDparam::g_CRDPhysics->velocityInterval().begin() + idxVel] =
        innerState.velocity()[idxVel];
    }
  inner[CRDparam::g_CRDPhysics->pressureIndex()] = innerState.pressure();
  // Force calculation of T
  inner[CRDparam::g_CRDPhysics->temperatureIndex()] = -1.0;
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      for (int c = CRDparam::g_CRDPhysics->speciesPrimInterval().begin(),
             c_end = CRDparam::g_CRDPhysics->speciesPrimInterval().end() + 1;
           c != c_end; ++c)
        {
          inner[c] = innerState(c);
        }
    }

  outer[CRDparam::g_CRDPhysics->densityIndex()] = outerState.density();
  for (const int idxVel : EachDir)
    {
      outer[CRDparam::g_CRDPhysics->velocityInterval().begin() + idxVel] =
        outerState.velocity()[idxVel];
    }
  outer[CRDparam::g_CRDPhysics->pressureIndex()] = outerState.pressure();
  // Force calculation of T
  outer[CRDparam::g_CRDPhysics->temperatureIndex()] = -1.0;
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      for (int c = CRDparam::g_CRDPhysics->speciesPrimInterval().begin(),
             c_end = CRDparam::g_CRDPhysics->speciesPrimInterval().end() + 1;
           c != c_end; ++c)
        {
          outer[c] = outerState(c);
        }
    }

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
      FArrayBox& Ufab = a_U[dit];

      // Working set boxes
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box1Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box1Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box1Dom, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(Ufab.box().contains(box1Dom));
      // Average values of W
      FABSTACKTEMP(Wfab, box1Dom, CRDparam::g_CRDPhysics->numPrimitive());
      for (int c = 0, c_end = CRDparam::g_CRDPhysics->numPrimitive();
           c != c_end; ++c)
        {
          MD_BOXLOOP(box1Dom, i)
            {
              RealVect X(D_DECL(XFab[MD_IX(i, 0)],
                                XFab[MD_IX(i, 1)],
                                XFab[MD_IX(i, 2)]));
              if ((X - m_bxLo)*(X - m_bxHi) <= RealVect_zero)
                {
                  Wfab[MD_IX(i, c)] = inner[c];
                }
              else
                {
                  Wfab[MD_IX(i, c)] = outer[c];
                }
            }
        }
      CRDparam::g_CRDPhysics->initialize(Ufab,
                                         Wfab,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box1Dom);
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
CNSIBCEulerShockBox::readBCInfo()
{
  CH_assert(!m_readInput);

//--Check the state information

  if (CRDState::nameIndex("inner") == -1)
    {
      CRD::msg << "Input (ShockBox IBC): requires state named 'inner'"
               << CRD::error;
    }

  if (CRDState::nameIndex("outer") == -1)
    {
      CRD::msg << "Input (ShockBox IBC): requires state named 'outer'"
               << CRD::error;
    }

  if (CRDState::nameIndex("wall") == -1)
    {
      CRD::msg << "Input (ShockBox IBC): requires state named 'wall'"
               << CRD::error;
    }

//--Read information particular to the IBC

  ParmParse ppIBC("ibc");

//--Lower corner of inner box

  std::vector<Real> rvec;
  m_bxLo = RealVect_zero;
  if (ppIBC.contains("lower"))
    {
      ppIBC.getarr("lower", rvec, 0, SpaceDim);
      for (const int dir : EachDir)
        {
          m_bxLo = rvec[dir];
        }
      if (m_bxLo < CRDparam::g_domainOrigin)
        {
          CRD::msg << "Input (ShockBox IBC): 'lower' must be >= "
                   << CRDparam::g_domainOrigin << '!' << CRD::error;
        }
      if (m_bxLo > CRDparam::g_domainOrigin + CRDparam::g_domainLength)
        {
          CRD::msg << "Input (ShockBox IBC): 'lower' must be <= "
                   << RealVect(CRDparam::g_domainOrigin +
                               CRDparam::g_domainLength) << '!' << CRD::error;
        }
    }

//--Upper corner of inner box

  m_bxHi = CRDparam::g_domainLength/2.;
  if (ppIBC.contains("upper"))
    {
      ppIBC.getarr("upper", rvec, 0, SpaceDim);
      for (const int dir : EachDir)
        {
          m_bxHi = rvec[dir];
        }
      if (m_bxHi < CRDparam::g_domainOrigin)
        {
          CRD::msg << "Input (ShockBox IBC): 'upper' must be >= "
                   << CRDparam::g_domainOrigin << '!' << CRD::error;
        }
      if (m_bxHi > CRDparam::g_domainOrigin + CRDparam::g_domainLength)
        {
          CRD::msg << "Input (ShockBox IBC): 'upper' must be <= "
                   << RealVect(CRDparam::g_domainOrigin +
                               CRDparam::g_domainLength) << '!' << CRD::error;
        }
      if (m_bxHi < m_bxLo)
        {
          CRD::msg << "Input (ShockBox IBC): 'upper' must be >= 'lower'!"
                   << CRD::error;
        }
    }

#if 0
//--Pressure ratio

  m_presRatio = 4.;
  ppIBC.query("pres_ratio", m_presRatio);
  if (m_presRatio == 0.)
    {
      CRD::msg << "Input (ShockBox IBC): 'pres_ratio' must not = 0!"
               << CRD::error;
    }

//--Tag threshold

  m_threshold = 0.2;
  ppIBC.query("tag_threshold", m_threshold);
  if (m_threshold < 0.)
    {
      CRD::msg << "Input (ShockBox IBC): 'tag_threshold' must be > 0!"
               << CRD::error;
    }
#endif

  m_readInput = true;
}
