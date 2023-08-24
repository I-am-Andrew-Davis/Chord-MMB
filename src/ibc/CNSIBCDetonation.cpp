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
 * \file CNSIBCDetonation.cpp
 *
 * \brief Member functions for CNSIBCDetonation
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCDetonation.H"
#include "CNSIBCCombustionTestF_F.H"
#include "CNSIBCFlameF_F.H"
#include "CRDPhysics.H"
#include "ChordInput.H"


/*******************************************************************************
 *
 * Class CNSIBCDetonation: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCDetonation::CNSIBCDetonation()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_highT(1000.),
  m_lowT(300.),
  m_gradLength(0.)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCDetonation::~CNSIBCDetonation()
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
CNSIBCDetonation::IBCName() const
{
  return "Detonation problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCDetonation::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "High temperature\n" << m_highT << CRD::var;
  CRD::msg << "Low temperature\n" << m_lowT << CRD::var;
  CRD::msg << "Initial mass fractions\n(" << m_initMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_initMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
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
CNSIBCDetonation::initialize(LevelData<FArrayBox>&      a_U,
                             LevelGridMetrics&          a_gridMetrics,
                             const LayoutData<FluxBox>& a_unitNormals,
                             const Real                 a_time,
                             const int                  a_level) const
{
  // Set the initial values
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

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;
      Box initBox(box2Dom);

      // Get physical coordinates
      FABSTACKTEMP(XiFab, initBox, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, initBox, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(initBox, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(initBox));
      // Create a FAB of the initial primitive variables
      FABSTACKTEMP(Wc, initBox, numWVar);
      Wc.setVal(0.);
      // Set velocity to zero
      Wc.setVal(0., initBox, WvelIndx, SpaceDim);
      MD_ARRAY_RESTRICT(arrW, Wc);
      MD_ARRAY_RESTRICT(arrX, XFab);
      MD_BOXLOOP(initBox, i)
        {
          Real xloc = arrX[MD_IX(i,0)];
          if (xloc < m_gradLength)
            {
              Real xratio = xloc/m_gradLength;
              arrW[MD_IX(i, tComp)] = m_highT - (m_highT - m_lowT)*xratio;
            }
          else
            {
              arrW[MD_IX(i, tComp)] = m_lowT;
            }
          arrW[MD_IX(i, rComp)] = -1.;
          arrW[MD_IX(i, pComp)] = m_initP;
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              const int wComp = wCompStart + spec;
              arrW[MD_IX(i, wComp)] = m_initMassFraction[spec];
            }
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         initBox);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
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
CNSIBCDetonation::readBCInfo()
{
  // Must be thermally perfect physics
  CH_assert(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf);
  ParmParse ppIBC("ibc");
  // Temperature gradient info
  ppIBC.query("gradient_length", m_gradLength);
  // Highest temperature in range
  ppIBC.query("high_temperature", m_highT);
  m_lowT = m_initT;
  ppIBC.query("low_temperature", m_lowT);
}
