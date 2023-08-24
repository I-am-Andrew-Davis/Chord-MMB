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
 * \file CNSIBCRiemannCube.cpp
 *
 * \brief Member functions for CNSIBCRiemannCube
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

#include "CNSIBCRiemannCube.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCRiemannCube: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCRiemannCube::CNSIBCRiemannCube()
  :
  CNSIBCGeneralized(),
  m_density(1,0.),
  m_temperature(1,0.),
  m_pressure(1,0.),
  m_velocity(1,0.),
  m_case(-1)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCRiemannCube::~CNSIBCRiemannCube()
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
CNSIBCRiemannCube::IBCName() const
{
  return "3D Riemann problem case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCRiemannCube::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  // CRD::msg << "Density    \n" << m_density     << CRD::var;
  // CRD::msg << "Temperature\n" << m_temperature << CRD::var;
  // CRD::msg << "Pressure   \n" << m_pressure    << CRD::var;
  // CRD::msg << "Velocity   \n" << m_velocity    << CRD::var;
  CRD::msg << "Case number\n" << m_case        << CRD::var;
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
CNSIBCRiemannCube::initialize(LevelData<FArrayBox>&      a_U,
                              LevelGridMetrics&          a_gridMetrics,
                              const LayoutData<FluxBox>& a_unitNormals,
                              const Real                 a_time,
                              const int                  a_level) const
{
  // Set the initial values
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int scalSize = std::pow(2,SpaceDim);
  const RealVect halfPt(CRDparam::g_domainOrigin
                        + 0.5*CRDparam::g_domainLength);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      CRD::msg << "Riemann cube does not currently support"
               << " thermally perfect physics!" << CRD::abort;
    }
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

      // Get coordinate system and domain for the block
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
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);

      D_TERM(const unsigned char x_flag = 1 << 0;,
             const unsigned char y_flag = 1 << 1;,
             const unsigned char z_flag = 1 << 2;);
      MD_ARRAY_RESTRICT(arrW,Wc);
      MD_ARRAY_RESTRICT(arrX,XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          // identify what quadrant we're in
          unsigned char quadrant = 0;
          D_TERM(
            if (arrX[MD_IX(i,0)] >= halfPt[0])
              {
                quadrant |= x_flag;
              },
            if (arrX[MD_IX(i,1)] >= halfPt[1])
              {
                quadrant |= y_flag;
              },
            if (arrX[MD_IX(i,2)] >= halfPt[2])
              {
                quadrant |= z_flag;
              });
          const int quad = static_cast<int>(quadrant);

          arrW[MD_IX(i,rhoIndx)]  = m_density[quad];
          arrW[MD_IX(i,presIndx)] = m_pressure[quad];
          arrW[MD_IX(i,tempIndx)] = m_temperature[quad];
          // assign velocities -- a little more thought required
          D_TERM(arrW[MD_IX(i,WvelIndx)]   = m_velocity[quad];,
                 arrW[MD_IX(i,WvelIndx+1)] = m_velocity[quad+scalSize];,
                 arrW[MD_IX(i,WvelIndx+2)] = m_velocity[quad+(2*scalSize)];);
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      //fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCRiemannCube::haveExactSol() const
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
CNSIBCRiemannCube::readBCInfo()
{
  ParmParse ppIBC("ibc");
  const int vecSize = std::pow(2,SpaceDim)*SpaceDim;
  const int strd = std::pow(2,SpaceDim); // stride for number of quadrants
  m_density.resize(strd);
  m_density.assign(strd,0.);
  m_temperature.resize(strd);
  m_temperature.assign(strd, 0.);
  m_pressure.resize(strd);
  m_pressure.assign(strd, 0.);
  m_velocity.resize(vecSize);
  m_velocity.assign(vecSize,0.);
  for (int idx = 0; idx != strd; ++idx)
    {
      ppIBC.query("density", m_density[idx], idx);
      ppIBC.query("temperature", m_temperature[idx], idx);
      ppIBC.query("pressure", m_pressure[idx], idx);
      std::ostringstream vel;
      vel << "velocity_Region" << idx;
      std::string velRegion(vel.str());
      for (int jdx = 0; jdx != SpaceDim; ++jdx)
        {
          const int kdx = strd*jdx + idx;
          ppIBC.query(velRegion.c_str(), m_velocity[kdx], jdx);
        }
    }
  ppIBC.query("case", m_case);
  switch(m_case)
    {
    case -1:
      break;
    case 1:
      D_TERM(m_density[0] = 0.1072;
             m_density[1] = 0.2579;,
             m_density[2] = 0.5197;
             m_density[3] = 1.;,
             m_density[4] = 0.1072;
             m_density[5] = 0.2579;
             m_density[6] = 0.5197;
             m_density[7] = 1.;);
      D_TERM(m_pressure[0] = 0.0439;
             m_pressure[1] = 0.15;,
             m_pressure[2] = 0.4;
             m_pressure[3] = 1.;,
             m_pressure[4] = 0.0439;
             m_pressure[5] = 0.15;
             m_pressure[6] = 0.4;
             m_pressure[7] = 1.;);
      D_TERM(
        D_TERM(m_velocity[0]          = -0.7259;,
               m_velocity[0+strd]     = -1.4045;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = -1.4045;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = -0.7259;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = -0.7259;,
               m_velocity[4+strd]     = -1.4045;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = -1.4045;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = -0.7259;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 2:
      D_TERM(m_density[0] = 1.;
             m_density[1] = 0.5197;,
             m_density[2] = 0.5197;
             m_density[3] = 1.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 1.;
             m_pressure[1] = 0.4;,
             m_pressure[2] = 0.4;
             m_pressure[3] = 1.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = -0.7259;,
               m_velocity[0+strd]     = -0.7259;,
               m_velocity[0+(2*strd)] = 0;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = -0.7259;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = -0.7259;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = -0.7259;,
               m_velocity[4+strd]     = -0.7259;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = -0.7259;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = -0.7259;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 3:
      D_TERM(m_density[0] = 0.138;
             m_density[1] = 0.5323;,
             m_density[2] = 0.5323;
             m_density[3] = 1.5;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.029;
             m_pressure[1] = 0.3;,
             m_pressure[2] = 0.3;
             m_pressure[3] = 1.5;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 1.206;,
               m_velocity[0+strd]     = 1.206;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 1.206;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 1.206;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 1.206;,
               m_velocity[4+strd]     = 1.206;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 1.206;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 1.206;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 4:
      D_TERM(m_density[0] = 1.1;
             m_density[1] = 0.5065;,
             m_density[2] = 0.5065;
             m_density[3] = 1.1;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 1.1;
             m_pressure[1] = 0.35;,
             m_pressure[2] = 0.35;
             m_pressure[3] = 1.1;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 5:
      D_TERM(m_density[0] = 1.;
             m_density[1] = 2.;,
             m_density[2] = 3.;
             m_density[3] = 1.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 1.;
             m_pressure[1] = 1.;,
             m_pressure[2] = 1.;
             m_pressure[3] = 1.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 6:
      D_TERM(m_density[0] = 1.;
             m_density[1] = 3.;,
             m_density[2] = 2.;
             m_density[3] = 1.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 1.;
             m_pressure[1] = 1.;,
             m_pressure[2] = 1.;
             m_pressure[3] = 1.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 7:
      D_TERM(m_density[0] = 0.8;
             m_density[1] = 0.5197;,
             m_density[2] = 0.5197;
             m_density[3] = 1.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.4;
             m_pressure[1] = 0.4;,
             m_pressure[2] = 0.4;
             m_pressure[3] = 1.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 8:
      D_TERM(m_density[0] = 0.8;
             m_density[1] = 1.;,
             m_density[2] = 1.;
             m_density[3] = 0.5197;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 1.;
             m_pressure[1] = 1.;,
             m_pressure[2] = 1.;
             m_pressure[3] = 0.4;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 9:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 10:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 11:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 12:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 13:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 14:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 15:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 16:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 17:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 18:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    case 19:
      D_TERM(m_density[0] = 0.;
             m_density[1] = 0.;,
             m_density[2] = 0.;
             m_density[3] = 0.;,
             m_density[4] = 0.;
             m_density[5] = 0.;
             m_density[6] = 0.;
             m_density[7] = 0.;);
      D_TERM(m_pressure[0] = 0.;
             m_pressure[1] = 0.;,
             m_pressure[2] = 0.;
             m_pressure[3] = 0.;,
             m_pressure[4] = 0.;
             m_pressure[5] = 0.;
             m_pressure[6] = 0.;
             m_pressure[7] = 0.;);
      D_TERM(
        D_TERM(m_velocity[0]          = 0.;,
               m_velocity[0+strd]     = 0.;,
               m_velocity[0+(2*strd)] = 0.;);
        D_TERM(m_velocity[1]          = 0.;,
               m_velocity[1+strd]     = 0.;,
               m_velocity[1+(2*strd)] = 0.;);,
        D_TERM(m_velocity[2]          = 0.;,
               m_velocity[2+strd]     = 0.;,
               m_velocity[2+(2*strd)] = 0.;);
        D_TERM(m_velocity[3]          = 0.;,
               m_velocity[3+strd]     = 0.;,
               m_velocity[3+(2*strd)] = 0.;);,
        D_TERM(m_velocity[4]          = 0.;,
               m_velocity[4+strd]     = 0.;,
               m_velocity[4+(2*strd)] = 0.;);
        D_TERM(m_velocity[5]          = 0.;,
               m_velocity[5+strd]     = 0.;,
               m_velocity[5+(2*strd)] = 0.;);
        D_TERM(m_velocity[6]          = 0.;,
               m_velocity[6+strd]     = 0.;,
               m_velocity[6+(2*strd)] = 0.;);
        D_TERM(m_velocity[7]          = 0.;,
               m_velocity[7+strd]     = 0.;,
               m_velocity[7+(2*strd)] = 0.;););
      break;
    }
  m_readInput = true;
}
