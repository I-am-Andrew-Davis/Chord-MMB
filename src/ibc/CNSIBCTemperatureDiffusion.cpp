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
 * \file CNSIBCTemperatureDiffusion.cpp
 *
 * \brief Member functions for CNSIBCTemperatureDiffusion
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "CRDparam.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "ViscousTensor4thOrderOp.H"
#include "ViscousTensor4thOrderOpF_F.H"
 
//----- Internal -----//

#include "CNSIBCTemperatureDiffusion.H"
#include "CNSIBCTemperatureDiffusionF_F.H"
#include "CNSIBCF_F.H"
#include "CRDState.H"
#include "CRDPhysics.H"
#include "LGintegrator.H"
#include "ChordInput.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"

/*******************************************************************************
 *
 * Class CNSIBCTemperatureDiffusion: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCTemperatureDiffusion::CNSIBCTemperatureDiffusion()
  :
  CNSIBCGeneralized()
{
  CNSIBCTemperatureDiffusion::readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTemperatureDiffusion::~CNSIBCTemperatureDiffusion()
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
CNSIBCTemperatureDiffusion::IBCName() const
{
  return "Temperature diffusion test";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCTemperatureDiffusion::writeIBCInfo() const
{
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
CNSIBCTemperatureDiffusion::initialize(LevelData<FArrayBox>&      a_U,
                                       LevelGridMetrics&          a_gridMetrics,
                                       const LayoutData<FluxBox>& a_unitNormals,
                                       const Real                 a_time,
                                       const int                  a_level) const
{
  // Global indices and variables
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cRho  = CRDparam::g_CRDPhysics->densityIndex();
  const int cPres = CRDparam::g_CRDPhysics->pressureIndex();
  // Obtain the mean values and the wave-amplitude values
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real rho = state.density();
  const RealVect vel = state.velocity();
  // const Real pres = state.pressure();
  
  // Define R
  Real R = CRDparam::g_R;

// Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box3Dom = grow(box, 3);
      box3Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      // Get physical coordinates
      FABSTACKTEMP(XiFab, box3Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box3Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box3Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box2Dom));
      // Pointwise values of W
      FABSTACKTEMP(Wc, box3Dom, numWVar);
      Wc.setVal(-1.);

      MD_ARRAY_RESTRICT(arrX, XFab);
      // Fill Wc with data
      
      MD_BOXLOOP(box3Dom, i)
        {
          RealVect xloc(D_DECL(arrX[MD_IX(i,0)],
                               arrX[MD_IX(i,1)],
                               arrX[MD_IX(i,2)]));
          RealVect pt = xloc - m_sourceLoc;
          Real radsq = D_TERM(pt[0]*pt[0],+ pt[1]*pt[1], + pt[2]*pt[2]);
          Real radius = std::sqrt(radsq);
          if (std::abs(radius) <= m_initRadius)
            {
              Wc[MD_IX(i, cPres)] = rho * R * m_tHi;
            }
          else
            {
              Wc[MD_IX(i, cPres)] = rho * R * m_tLo;
            }
          
          Wc[MD_IX(i, cRho)] = rho;
          D_TERM(
            Wc[MD_IX(i, cVel)] = vel[0];,
            Wc[MD_IX(i, cVel+1)] = vel[1];,
            Wc[MD_IX(i, cVel+2)] = vel[2];);
          //Wc[MD_IX(i, cPres)] = pres;

        }
      // Initialize the values in UFab
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
//  Add body force
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg
 *                      Cell-averaged conservative variables
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Current level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_didx  Index of current box in layout
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCTemperatureDiffusion::addSourceTerm(
  FArrayBox&           a_sourceFab,
  FArrayBox&           a_invDtFab,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_UcellAvg,
  const FluxBox&       a_WfaceAvgFxb,
  const ProblemDomain& a_domain,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const Real           a_stageWeight,
  const int            a_level,
  const Box&           a_disjointBox,
  const Box&           a_solveBox,
  const DataIndex&     a_didx,
  const Real           a_globalKE,
  const Real           a_globalHelicity) const
{
  // Add the source terms to a_sourceFab
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_solveBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_solveBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(a_solveBox, XiFab, XFab, blockCoordSys);

  // Global indices and variables
  // const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int cPres    = CRDparam::g_CRDPhysics->pressureIndex();
  const Real R       = CRDparam::g_R;
  const Real gamma   = CRDparam::g_gamma;
  const Real cv      = R / (gamma-1);  // Only works when R is constant (?)

  // Grab the state
  const CRDState& state = CRDState::get(m_idxStateInit);
  const Real rho        = state.density();
  const RealVect vel    = state.velocity();
  // const Real pres       = state.pressure();
  
  const Real sourceCoeff = m_sourceAlpha / (m_sourceSigma * 2.506628); // 2.5 is sqrt(2*pi)
  const Real multDenom   = 1/(2*m_sourceSigma*m_sourceSigma);
  
  MD_ARRAY_RESTRICT(arrSource, a_sourceFab);
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_BOXLOOP(a_solveBox, i)
    {
      if (a_time <= m_timeSourceEnd)
        {
          RealVect xloc(D_DECL(XFab[MD_IX(i,0)],
                               XFab[MD_IX(i,1)],
                               XFab[MD_IX(i,2)]));
          RealVect pt = xloc - m_sourceLoc;
          Real absPtSq(D_TERM( pt[0]*pt[0],
                               +pt[1]*pt[1],
                               +pt[2]*pt[2]));
          Real tempAdd = sourceCoeff * std::exp( - absPtSq * multDenom);
          // Must convert to pressure, the primitive var in arrSource
          // the tempAdd term above is the temperature to add
          arrSource[MD_IX(i, cPres)] = rho* cv* tempAdd;
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTemperatureDiffusion::haveExactSol() const
{
  return false;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/** 
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update
 *  \param[in]  a_problemDomain
 *                      The problem domain on the level
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx  Current DataIndex on the current disjoint box
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Grid level
 *  \return             0 - Successfully computed exact solution
 *                      1 - Exact solution is not known
 *//*-----------------------------------------------------------------*/

int
CNSIBCTemperatureDiffusion::exactSol(FArrayBox&              a_Ux,
                                     const Box&              a_box,
                                     const Box&              a_disjointBox,
                                     const LevelGridMetrics& a_gridMetrics,
                                     const FluxBox&          a_unitNormals,
                                     const DataIndex&        a_didx,
                                     const Real              a_time,
                                     const int               a_level) const
{
  return 0;
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
CNSIBCTemperatureDiffusion::readBCInfo()
{
  //CH_assert(!m_readInput); :'(
  ParmParse ppIBC("ibc");

  std::vector<Real> trueSourceLoc(SpaceDim);
  ppIBC.getarr("true_source_location", trueSourceLoc, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_trueSourceLoc.dataPtr(),
                                           &trueSourceLoc.front());
  m_trueSourceLoc = m_trueSourceLoc*CRDparam::g_domainLength;
  
  std::vector<Real> sourceLoc(SpaceDim);
  ppIBC.getarr("source_location", sourceLoc, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sourceLoc.dataPtr(),
                                           &sourceLoc.front());
  m_sourceLoc = m_sourceLoc*CRDparam::g_domainLength;

  m_sourceSigma = 0.1;
  ppIBC.query("source_sigma", m_sourceSigma);

  m_sourceAlpha = 2.;
  ppIBC.query("source_alpha", m_sourceAlpha);

  m_timeSourceEnd = 0.05;
  ppIBC.query("source_time_end", m_timeSourceEnd);

  m_tHi = 1600;
  ppIBC.query("initial_temp_hi", m_tHi);

  m_tLo = 298;
  ppIBC.query("initial_temp_lo", m_tLo);

  m_initRadius = 0.1;
  ppIBC.query("initial_radius_size", m_initRadius);
  
  m_readInput = true;
}
