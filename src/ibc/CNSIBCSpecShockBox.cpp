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
 * \file CNSIBCSpecShockBox.cpp
 *
 * \brief Member functions for CNSIBCSpecShockBox
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCSpecShockBox.H"
#include "CRDPhysics.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class CNSIBCSpecShockBox: member definitions
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

CNSIBCSpecShockBox::CNSIBCSpecShockBox()
  :
  CNSIBCCombustionReference(),
  // All invalid values which must be corrected
  m_rhoL(-1.),
  m_presL(-1.),
  m_tL(-1.),
  m_rhoU(-1.),
  m_presU(-1.),
  m_tU(-1.)
{
//--Read any BC info

  readBCInfo();

//--Set BC Type

// Normally you would set everything to wall
  setAllDomainBC(CRDparam::DomainBCTypeWall);
  setAllDomainBCOrder(1);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCSpecShockBox::~CNSIBCSpecShockBox()
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
CNSIBCSpecShockBox::IBCName() const
{
  return "Shock box with species transport";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpecShockBox::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Lower density\n" << m_rhoL << CRD::var;
  CRD::msg << "Lower pressure\n" << m_presL << CRD::var;
  CRD::msg << "Upper density\n" << m_rhoU << CRD::var;
  CRD::msg << "Upper pressure\n" << m_presU << CRD::var;
  CRD::msg << "Initial upper mass fractions\n(" << m_upperMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_upperMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Initial lower mass fractions\n(" << m_lowerMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_lowerMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
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
CNSIBCSpecShockBox::initialize(LevelData<FArrayBox>&      a_U,
                               LevelGridMetrics&          a_gridMetrics,
                               const LayoutData<FluxBox>& a_unitNormals,
                               const Real                 a_time,
                               const int                  a_level) const
{
  // This box defines the lower side.  It grown by quite a bit outside the
  // domain.
  IntVect middle(IntVect::Zero);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      middle[dir] = CRDparam::g_domainBaseSize[dir]/2-1;
    }
  Box lowerBox(-5*IntVect::Unit, middle);
  lowerBox.refine(CRDparam::g_refFromBase[a_level]);
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

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Data
      FArrayBox& UFab = a_U[dit];

      // Working set boxes
      Box box3Dom = grow(box, 3);
      box3Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;
      Box initBox(box3Dom);

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(initBox));
      FABSTACKTEMP(Wc, initBox, numWVar);
      Wc.setVal(m_rhoU, rComp);
      Wc.setVal(0., initBox, WvelIndx, SpaceDim);
      Wc.setVal(m_presU, pComp);
      Wc.setVal(m_tU, tComp);
      for (int i = 0; i != numSpecies; ++i)
        {
          int comp = wCompStart + i;
          Wc.setVal(m_upperMassFraction[i], comp);
        }
      Box cBoxL(initBox);
      cBoxL&=lowerBox;
      if (!(cBoxL.isEmpty()))
        {
          Wc.setVal(m_rhoL, cBoxL, rComp);
          Wc.setVal(m_presL, cBoxL, pComp);
          Wc.setVal(m_tL, cBoxL, tComp);
          for (int i = 0; i != numSpecies; ++i)
            {
              int comp = wCompStart + i;
              Wc.setVal(m_lowerMassFraction[i], cBoxL, comp);
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
CNSIBCSpecShockBox::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(true);  // Disable terminate on error
  ParmParse ppIBC("ibc");
  const int numSpecies = CRDparam::g_numSpecies;

//--Lower state

  ppIBC.query("density_lower", m_rhoL);
  ppIBC.query("pressure_lower", m_presL);
  ppIBC.query("temperature_lower", m_tL);
  if (m_presL*m_tL*m_rhoL > 0.)
    {
      CRD::msg << "Input (SpecShockBox IBC): Must specify only 2 of the 3 for"
               << " P, T, or rho for lower region!" << CRD::error;
    }

//--Upper state

  ppIBC.query("density_upper", m_rhoU);
  ppIBC.query("pressure_upper", m_presU);
  ppIBC.query("temperature_upper", m_tU);
  if (m_presU*m_tU*m_rhoU > 0.)
    {
      CRD::msg << "Input (SpecShockBox IBC): Must specify only 2 of the 3 for"
               << " P, T, or rho for upper region!" << CRD::error;
    }
  m_lowerMassFraction.resize(numSpecies);
  m_lowerMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  int massCheck = assignMassFractions(m_lowerMassFraction,
                                      "lower_specs",
                                      "lower_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecShock IBC): 'lower_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }
  
  m_upperMassFraction.resize(numSpecies);
  m_upperMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  massCheck = assignMassFractions(m_upperMassFraction,
                                  "upper_specs",
                                  "upper_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecShock IBC): 'upper_mass_fractions' must be "
        "equal to 1." << CRD::error;
    }

  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
