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
 * \file CNSIBCSpecShock.cpp
 *
 * \brief Member functions for CNSIBCSpecShock
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCSpecShock.H"
#include "CRDPhysics.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class CNSIBCSpecShock: member definitions
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

CNSIBCSpecShock::CNSIBCSpecShock()
  :
  CNSIBCCombustionReference(),
  // All invalid values which must be corrected
  m_rhoL(-1.),
  m_presL(-1.),
  m_tempL(-1.),
  m_rhoR(-1.),
  m_presR(-1.),
  m_tempR(-1.),
  m_velL(0.),
  m_velR(0.),
  m_Cr(1.E8),
  m_tubeDir(-1)
{

//--Read any BC info

  readBCInfo();

//--Set BC Type

  setAllDomainBC(CRDparam::DomainBCTypeSlipWall);
  setAllDomainBCOrder(1);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCSpecShock::~CNSIBCSpecShock()
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
CNSIBCSpecShock::IBCName() const
{
  return "Shock tube with species transport";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpecShock::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Density at left\n" << m_rhoL << CRD::var;
  CRD::msg << "Pressure at left\n" << m_presL << CRD::var;
  CRD::msg << "Density at right\n" << m_rhoR << CRD::var;
  CRD::msg << "Pressure at right\n" << m_presR << CRD::var;
  CRD::msg << "Initial smearing parameter\n" << m_Cr << CRD::var;
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
CNSIBCSpecShock::initialize(LevelData<FArrayBox>&      a_U,
                            LevelGridMetrics&          a_gridMetrics,
                            const LayoutData<FluxBox>& a_unitNormals,
                            const Real                 a_time,
                            const int                  a_level) const
{
  Real len = CRDparam::g_domainLength[0];
  Real center = len/2.;
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  Real sameCnTest = 10.;
  bool sameCnBool = true;
  for (int spec = 0; spec != numSpecies; ++spec)
    {
      sameCnTest = std::min(
        std::abs(m_leftMassFraction[spec] - m_rightMassFraction[spec]),
        sameCnTest);
    }
  if (sameCnTest == 0.)
    {
      sameCnBool = false;
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
          Real gap = x - center;
          Real pres = m_presL;
          Real rho = m_rhoL;
          Real temp = m_tempL;
          if (sameCnBool)
            {
              Real perc = 0.5*(1. + std::tanh(m_Cr*gap/len));
              Real summf = 0.;
              Real sumside = 0.;
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  const int wComp = wCompStart + spec;
                  const Real cnL = m_leftMassFraction[spec];
                  const Real cnR = m_rightMassFraction[spec];
                  Real cnval = cnL + perc*(cnR - cnL);
                  arrW[MD_IX(i, wComp)] = cnval;
                  summf += (cnval - cnR)/(cnL - cnR);
                  sumside += 1.;
                }
              summf /= sumside;
              pres = 1./(summf/m_presL + (1. - summf)/m_presR);
              rho = 1./(summf/m_rhoL + (1. - summf)/m_rhoR);
              temp = 1./(summf/m_tempL + (1. - summf)/m_tempR);
            }
          else
            {
              if (x > center)
                {
                  pres = m_presR;
                  rho = m_rhoR;
                  temp = m_tempR;
                  for (int spec = 0; spec != numSpecies; ++spec)
                    {
                      const int wComp = wCompStart + spec;
                      arrW[MD_IX(i, wComp)] = m_rightMassFraction[spec];
                    }
                }
              else
                {
                  for (int spec = 0; spec != numSpecies; ++spec)
                    {
                      const int wComp = wCompStart + spec;
                      arrW[MD_IX(i, wComp)] = m_leftMassFraction[spec];
                    }
                }
            }
          arrW[MD_IX(i, rComp)] = rho;
          arrW[MD_IX(i, pComp)] = pres;
          arrW[MD_IX(i, tComp)] = temp;
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
CNSIBCSpecShock::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(true);  // Disable terminate on error
  ParmParse ppIBC("ibc");
  const int numSpecies = CRDparam::g_numSpecies;
//--Tube direction

  m_tubeDir = 0;
  ppIBC.query("tube_dir", m_tubeDir);
  if (m_tubeDir < 0 || m_tubeDir >= SpaceDim)
    {
      CRD::msg << "Input (SpecShock IBC): 'tube_dir' must be >= 0 and < "
               << SpaceDim << '!' << CRD::error;
    }

//--Left state
  ppIBC.query("density_left", m_rhoL);
  ppIBC.query("pressure_left", m_presL);
  ppIBC.query("temperature_left", m_tempL);
  if (m_presL*m_tempL*m_rhoL > 0. || (m_presL + m_tempL + m_rhoL) == -3.)
    {
      CRD::msg << "Input (SpecShock IBC): Must specify 2 of the 3 for"
               << " P_L, T_L, or rho_L!" << CRD::error;
    }
  ppIBC.query("velocity_left", m_velL);

//--Right state
  ppIBC.query("density_right", m_rhoR);
  ppIBC.query("pressure_right", m_presR);
  ppIBC.query("temperature_right", m_tempR);
  if (m_presR*m_tempR*m_rhoR > 0. || (m_presR + m_tempR + m_rhoR) == -3.)
    {
      CRD::msg << "Input (SpecShock IBC): Must specify 2 of the 3 for"
               << " P_R, T_R, or rho_R!" << CRD::error;
    }
  ppIBC.query("velocity_right", m_velR);
  m_leftMassFraction.resize(numSpecies);
  m_leftMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  int massCheck = assignMassFractions(m_leftMassFraction,
                                      "left_specs",
                                      "left_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecShock IBC): 'left_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }
  m_rightMassFraction.resize(numSpecies);
  m_rightMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  massCheck = assignMassFractions(m_rightMassFraction,
                                  "right_specs",
                                  "right_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecShock IBC): 'right_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }

  ppIBC.query("init_smearing_parameter", m_Cr);
  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  // Figure out if the temperature behind the contact surface is less than
  // 280
  if ((ppIBC.contains("temperature_left") && m_tempL < 280.) ||
      (ppIBC.contains("temperature_right") && m_tempR < 280.))
    {
      CRD::msg << "Temperature values must be > 280!" << CRD::error;
    }
  m_readInput = true;
}
