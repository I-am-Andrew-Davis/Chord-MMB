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
 * \file CNSIBCGaussPulse.cpp
 *
 * \brief Member functions for CNSIBCGaussPulse
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "CONSTANTS.H"
 
//----- Internal -----//

#include "CNSIBCGaussPulse.H"
#include "CNSIBCVortexF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"


/*******************************************************************************
 *
 * Class CNSIBCGaussPulse: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCGaussPulse::CNSIBCGaussPulse()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_deltaRho(0.),
  m_center(-1.*RealVect::Unit)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCGaussPulse::~CNSIBCGaussPulse()
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
CNSIBCGaussPulse::IBCName() const
{
  return "Gaussian pulse case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCGaussPulse::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Delta density\n" << m_deltaRho << CRD::var;
  CRD::msg << "Center of Gaussian (domain  Dim.)\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << m_center[dir];
    }
  CRD::msg << ')' << CRD::var;
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
CNSIBCGaussPulse::initialize(LevelData<FArrayBox>&      a_U,
                             LevelGridMetrics&          a_gridMetrics,
                             const LayoutData<FluxBox>& a_unitNormals,
                             const Real                 a_time,
                             const int                  a_level) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg);
  RealVect timeCenter(D_DECL(m_center[0] + m_initVel[0]*a_time,
                             m_center[1] + m_initVel[1]*a_time,
                             m_center[2] + m_initVel[2]*a_time));
  D_TERM(timeCenter[0] -= std::floor(timeCenter[0]);,
         timeCenter[1] -= std::floor(timeCenter[1]);,
         timeCenter[2] -= std::floor(timeCenter[2]););
  std::vector<Real> cn(numSpecies);
  if (m_initRho < 0.)
    {
      CRD::msg << "Must specify 'initial_density'!" << CRD::error;
    }
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, numGhosts);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, numGhosts - 1);
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
      MD_ARRAY_RESTRICT(arrW, Wc);
      MD_ARRAY_RESTRICT(arrX, XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          RealVect xloc(D_DECL(arrX[MD_IX(i,0)],
                               arrX[MD_IX(i,1)],
                               arrX[MD_IX(i,2)]));
          RealVect pt = xloc - timeCenter;
          Real radsq = D_TERM(pt[0]*pt[0],+ pt[1]*pt[1], + pt[2]*pt[2]);
          Real radius = std::sqrt(radsq);
          Real smoo = 0.;
          if (std::abs(radius) <= 0.5)
            {
              smoo = m_deltaRho*std::exp(-16.*radsq)*
                std::pow(std::cos(Pi*radius),6);
            }
          Real rho = m_initRho + smoo;
          arrW[MD_IX(i,rhoIndx)] = rho;
          Real rgas = CRDparam::g_R;
          if (numSpecies > 0)
            {
              rgas = 0.;
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  int wComp = wCompStart + spec;
                  cn[spec] = m_initMassFraction[spec];
                  arrW[MD_IX(i,wComp)] = cn[spec];
                  rgas += cn[spec]*
                    CRDparam::g_CRDPhysics->speciesGasConstant(spec);
                }
            }
          Real T = std::max(m_initT, m_initP/(rgas*m_initRho));
          if (T < 300.)
            {
              CRD::msg << "T must be greater than 300!" << CRD::error;
            }
          Real gamma = CRDparam::g_CRDPhysics->gamma(T, cn.data());
          Real pres = std::max(m_initP, m_initRho*T*rgas);
          arrW[MD_IX(i,presIndx)] = pres*std::pow(rho/m_initRho, gamma);
          arrW[MD_IX(i,tempIndx)] = -1.;
          D_TERM(arrW[MD_IX(i, velIndx)] = m_initVel[0];,
                 arrW[MD_IX(i, velIndx + 1)] = m_initVel[1];,
                 arrW[MD_IX(i, velIndx + 2)] = m_initVel[2];);
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

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCGaussPulse::haveExactSol() const
{
  return false;
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
CNSIBCGaussPulse::readBCInfo()
{
  ParmParse ppIBC("ibc");

  m_deltaRho = 0.14;
  ppIBC.query("delta_density", m_deltaRho);

//--Center of the Gaussian at time t=0, on a domain from 0 to 1
//--(default (0.5,0.5,0.5))

  m_center = 0.5*RealVect::Unit;
  if (ppIBC.contains("center"))
    {
      std::vector<Real> IBCcenter(SpaceDim);
      ppIBC.getarr("center", IBCcenter, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_center.dataPtr(),
                                               &IBCcenter.front());
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (m_center[dir] < 0. || m_center[dir] > 1.)
            {
              CRD::msg << "Input (GaussPulse IBC): 'center' must be >= 0.0 "
                "and < 1.0!" << CRD::error;
            }
        }
      m_center = m_center*CRDparam::g_domainLength;
    }

  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
