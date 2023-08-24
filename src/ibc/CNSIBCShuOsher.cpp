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
 * \file CNSIBCShuOsher.cpp
 *
 * \brief Member functions for CNSIBCShuOsher
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCShuOsher.H"
#include "CRDPhysics.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class CNSIBCShuOsher: member definitions
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

CNSIBCShuOsher::CNSIBCShuOsher()
  :
  CNSIBCGeneralized(),
  // All invalid values which must be corrected
  m_shockLoc(0.1*CRDparam::g_domainLength[0]),
  m_rhoR(-1.),
  m_presR(-1.),
  m_tempR(-1.),
  m_velR(0.),
  m_rhoAmp(0.2),
  m_rhoFreq(5.)
{

//--Read any BC info

  readBCInfo();

//--Set BC Type

  setAllDomainBC(CRDparam::DomainBCTypePeriodic);
  setDomainBC(0, Side::Lo, CRDparam::DomainBCTypeExtrapolated);
  setDomainBC(0, Side::Hi, CRDparam::DomainBCTypeExtrapolated);
  setAllDomainBCOrder(1);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCShuOsher::~CNSIBCShuOsher()
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
CNSIBCShuOsher::IBCName() const
{
  return "Shock tube with species transport";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCShuOsher::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Density at left\n" << m_initRho << CRD::var;
  CRD::msg << "Pressure at left\n" << m_initP << CRD::var;
  CRD::msg << "Density at right\n" << m_rhoR << CRD::var;
  CRD::msg << "Pressure at right\n" << m_presR << CRD::var;
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
CNSIBCShuOsher::initialize(LevelData<FArrayBox>&      a_U,
                           LevelGridMetrics&          a_gridMetrics,
                           const LayoutData<FluxBox>& a_unitNormals,
                           const Real                 a_time,
                           const int                  a_level) const
{
  Real pL = m_initP;
  Real rhoL = m_initRho;
  Real TL = m_initT;
  Real uL = m_initVel[0];

  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int pIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Get coordinate system for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));

      // Data
      FArrayBox& UFab = a_U[dit];

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

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FABSTACKTEMP(Wc, initBox, numWVar);
      Wc.setVal(0., initBox, WvelIndx, SpaceDim);
      MD_ARRAY_RESTRICT(arrW, Wc);
      MD_ARRAY_RESTRICT(arrX, XFab);
      MD_BOXLOOP(initBox, i)
        {
          Real x = arrX[MD_IX(i, 0)];
          Real pres = pL;
          Real rho = rhoL;
          Real temp = TL;
          arrW[MD_IX(i, WvelIndx)] = uL;
          if (x > m_shockLoc)
            {
              rho = m_rhoR*(1. + m_rhoAmp*std::sin(m_rhoFreq*x));
              temp = m_tempR;
              pres = m_presR;
              arrW[MD_IX(i, WvelIndx)] = m_velR;
            }    
          arrW[MD_IX(i, rIndx)] = rho;
          arrW[MD_IX(i, pIndx)] = pres;
          arrW[MD_IX(i, tIndx)] = temp;
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         initBox);
      //fourthOrderAverageCell(UFab, blockDomain, box1Dom);
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
CNSIBCShuOsher::readBCInfo()
{
  ParmParse ppIBC("ibc");

  ppIBC.query("shock_loc", m_shockLoc);

//--Right state is defined as the initial values

//--Left state
  ppIBC.query("density_right", m_rhoR);
  ppIBC.query("pressure_right", m_presR);
  ppIBC.query("temperature_right", m_tempR);
  if (m_presR*m_tempR*m_rhoR > 0. || (m_presR + m_tempR + m_rhoR) == -3.)
    {
      CRD::msg << "Input (ShuOsher IBC): Must specify 2 of the 3 for"
               << " density, pressure, and temperature_left!"
               << CRD::error;
    }
  ppIBC.query("velocity_right", m_velR);
  ppIBC.query("rho_amp", m_rhoAmp);
  ppIBC.query("rho_freq", m_rhoFreq);
  m_readInput = true;
}
