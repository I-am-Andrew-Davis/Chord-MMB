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
 * \file CNSIBCTransientCouette.cpp
 *
 * \brief Member functions for CNSIBCTransientCouette
 *
 *//*+*************************************************************************/

//----- Standard -----//
//
#include <random>
#include <fstream>

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCTransientCouette.H"
#include "CNSIBCTransientCouetteF_F.H"
#include "CNSIBCF_F.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodBaseBox.H"
#include "CRDState.H"


/*******************************************************************************
 *
 * Class CNSIBCTransientCouette: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCTransientCouette::CNSIBCTransientCouette()
  :
  CNSIBCGeneralized(),
  m_wallNormalDir(1),
  m_movWall(0),
  m_tagPct(0.),
  m_threshold(0.0)
{
  readBCInfo();
  const int wallPct = std::ceil(m_tagPct*
                                CRDparam::g_domainBaseSize[m_wallNormalDir]);
  {
    // Set lower box for tagging, grow by 5 in periodic directions
    IntVect lo(IntVect::Zero);
    IntVect hi(CRDparam::g_domainBaseSize + 5*IntVect::Unit);
    hi[m_wallNormalDir] = wallPct;
    m_loTagBox.define(lo, hi);
  }
  {
    // Set upper box for tagging, grow by 5 in periodic directions
    IntVect lo(IntVect::Zero);
    IntVect hi(CRDparam::g_domainBaseSize + 5*IntVect::Unit);
    hi[m_wallNormalDir] = CRDparam::g_domainBaseSize[m_wallNormalDir] - 1;
    lo[m_wallNormalDir] = CRDparam::g_domainBaseSize[m_wallNormalDir] - 1
      - wallPct;
    m_hiTagBox.define(lo, hi);
  }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTransientCouette::~CNSIBCTransientCouette()
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
CNSIBCTransientCouette::IBCName() const
{
  return "Navier-Stokes transient Couette";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCTransientCouette::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Wall velocity\n(";
  CRD::msg << ')' << CRD::var;
  CRD::msg << "Wall normal direction\n" << m_wallNormalDir << CRD::var;
  CRD::msg << "Moving wall\n";
  if (m_movWall == 0)
    {
      CRD::msg << "low" << CRD::var;
    }
  else
    {
      CRD::msg << "high" << CRD::var;
    }
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
CNSIBCTransientCouette::initialize(LevelData<FArrayBox>&      a_U,
                                   LevelGridMetrics&          a_gridMetrics,
                                   const LayoutData<FluxBox>& a_unitNormals,
                                   const Real                 a_time,
                                   const int                  a_level) const
{
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

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      this->exactSol(UFab,
                     box1Dom,
                     box,
                     a_gridMetrics,
                     a_unitNormals[dit],
                     dit(),
                     a_time,
                     a_level);
    }
}

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
CNSIBCTransientCouette::setTagMethod(const int a_tagBufferSize)
{
  if(m_steadyState)
    {
      return setTagMethodSteadyState(a_tagBufferSize);
    }

  TagLevel* tagLevel = new TagLevel;
  CNSIBCCombustionReference::setTagMethodLevel(a_tagBufferSize, tagLevel);
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  // Return the factory
  return new TagLevelFactory(tagLevel);
}

TagLevelFactory*
CNSIBCTransientCouette::setTagMethodSteadyState(const int a_tagBufferSize)
{
  const int maxAMRLevel = CRDparam::maxAMRLevel();
  TagLevel* tagLevel;
  std::vector<TagLevel*> tagLevelVec;
  std::vector<int> levelMapVec;
  const RealVect dxVect = CRDparam::g_domainLength/CRDparam::g_domainBaseSize;

  switch (maxAMRLevel)
    {
    case 0:
    case 1:
      {
        // One TagLevel
        tagLevelVec.resize(1);
        levelMapVec.resize(1);
        // Apply this TagLevel to all levels
        levelMapVec[0] = 0;
        // Tag a box grown by the max displacement of the BL
        tagLevel = new TagLevel;
        tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
        // Refine region near top wall
        {
          IntVect lo(IntVect::Zero);
          IntVect hi = CRDparam::g_domainBaseSize;
          hi[m_wallNormalDir] = std::ceil(CRDparam::g_domainLength[m_wallNormalDir]/
                                        dxVect[m_wallNormalDir]*0.25);
          Box tagBox(lo,hi);
          tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        }
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[0] = tagLevel;
        break;
      }
    case 2:
      {
        // One TagLevel
        tagLevelVec.resize(2);
        levelMapVec.resize(2);
        // Apply this TagLevel to all levels
        levelMapVec[0] = 0;
        levelMapVec[1] = 1;

        // Tag a box grown by the max displacement of the BL
        tagLevel = new TagLevel;
        tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
        IntVect lo(IntVect::Zero);
        IntVect hi = CRDparam::g_domainBaseSize;
        hi[m_wallNormalDir] = std::ceil((CRDparam::g_domainLength[m_wallNormalDir]/
                                       dxVect[m_wallNormalDir])*0.25);
        Box tagBox(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[0] = tagLevel;

        // Second level
        tagLevel = new TagLevel;
        hi[m_wallNormalDir] = std::ceil((CRDparam::g_domainLength[m_wallNormalDir]/
                                       dxVect[m_wallNormalDir])*0.15);
        tagBox.define(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[1] = tagLevel;
        break;
      }
    case 3:
      {
        // One TagLevel
        tagLevelVec.resize(2);
        levelMapVec.resize(2);
        // Apply this TagLevel to all levels
        levelMapVec[0] = 0;
        levelMapVec[1] = 1;

        // Tag a box grown by the max displacement of the BL
        tagLevel = new TagLevel;
        tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
        IntVect lo(IntVect::Zero);
        IntVect hi = CRDparam::g_domainBaseSize;
        hi[m_wallNormalDir] = std::ceil((CRDparam::g_domainLength[m_wallNormalDir]/
                                       dxVect[m_wallNormalDir])*0.25);
        Box tagBox(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[0] = tagLevel;

        // Second level
        tagLevel = new TagLevel;
        hi[m_wallNormalDir] = std::ceil((CRDparam::g_domainLength[m_wallNormalDir]/
                                       dxVect[m_wallNormalDir])*0.15);
        tagBox.define(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[1] = tagLevel;

        //  Third level
        tagLevel = new TagLevel;
        hi[m_wallNormalDir] = std::ceil((CRDparam::g_domainLength[m_wallNormalDir]/
                                       dxVect[m_wallNormalDir])*0.05);
        tagBox.define(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
        // Add in tag buffer
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[1] = tagLevel;
        break;
      }
    }

  // Return the factory
  return new TagLevelFactory(tagLevelVec, levelMapVec);
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTransientCouette::haveExactSol() const
{
  // no exact solution exists for the turbulent case
  // return ~m_turbulent;
  return !m_turbulent;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/** The exact solution is for \f$Re = 200\f$, air at 300K
 *  \f$\nu = 1.56e-5 m^2/s\f$ \f$h = 0.1m\f$, solution time = 500s,
 *  the number n in the series is 500.  The exact solution form
 *  (moving wall is on the low side of y direction) is
 *  \f[
 *     U(y,t) = U_0\left(1 - \frac{y}{h}\right) -
 *       2 U_0\sum\limits_{n=1}^{N}\left(\frac{\sin(n\pi \frac{y}{h})}
 *       {n\pi}e^{-\frac{\displaystyle n^2\pi^2\nu t}
 *       {\displaystyle h^2}}\right)
 *  \f]
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
CNSIBCTransientCouette::exactSol(FArrayBox&              a_Ux,
                                 const Box&              a_box,
                                 const Box&              a_disjointBox,
                                 const LevelGridMetrics& a_gridMetrics,
                                 const FluxBox&          a_unitNormals,
                                 const DataIndex&        a_didx,
                                 const Real              a_time,
                                 const int               a_level) const
{
  // Initial state
  const CRDState& state = CRDState::get(m_idxStateInit);
  //
  if(m_steadyState)
    {
      this->steadyState(a_Ux,
                        a_box,
                        a_disjointBox,
                        a_gridMetrics,
                        a_time,
                        a_level);
      return 0;
    }

  const int  nTerms = 5000;  // 40K at y = 0.0001*h

  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  // Working set boxes
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;
  Box initBox(box1Dom);

  // Get physical coordinates
  FABSTACKTEMP(XiFab, box1Dom, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, box1Dom, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(box1Dom, XiFab, XFab, blockCoordSys);
  
  // Set all the components to state with zero velocity before computing the
  // exact solution.
  Real rho = state.density();
  Real pres = state.pressure();

  a_Ux.setVal(rho, URHO);
  a_Ux.setVal(0., a_Ux.box(), UMOMX, SpaceDim);
  const Real gamma = CRDparam::g_gamma;
  const Real nu = CRDparam::g_mu/rho;
  const Real height = CRDparam::g_physicalLength[m_wallNormalDir];
  a_Ux.setVal(pres/(CRDparam::g_gamma - 1.), UENG);
  const RealVect wallSpeed = state.velocity();
  const Real tstar = nu*a_time/(height*height);
  // At time 0, the velocity everywhere is 0
  if (tstar <= 0.)
    {
      fourthOrderAverageCell(a_Ux, blockDomain, a_box);
      return 0;
    }

  // Between 0 and 0.01, we have not validated the exact solution
  if (tstar < 0.01)
    {
      Real min_t = 0.01*height*height/nu;
      CRD::msg << "Unvalidated time for analytic Couette flow solution"
               << "t < " << min_t << CRD::warn;
    }

  // comp is the component to solve for in a_Ux
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  int velDir = m_wallNormalDir - 1;
  if (velDir < 0)
    {
      velDir = SpaceDim - 1;
    }
  if (wallSpeed[velDir] == 0.)
    {
      return 0;
    }
  int velComp = velDir + velIndx;

  const int engComp = UENG;
  FORT_CNSIBCTRANSIENTCOUETTEEXACTSOL(
    CHF_FRA(a_Ux),
    CHF_BOX(box1Dom),
    CHF_CONST_FRA1(XFab,m_wallNormalDir),
    CHF_CONST_REAL(a_time),
    CHF_CONST_REAL(height),
    CHF_CONST_REAL(nu),
    CHF_CONST_REAL(wallSpeed[velDir]),
    CHF_CONST_REAL(rho),
    CHF_CONST_REAL(pres),
    CHF_CONST_REAL(gamma),
    CHF_CONST_INT(m_movWall),
    CHF_CONST_INT(nTerms),
    CHF_CONST_INT(velComp),
    CHF_CONST_INT(engComp));

  fourthOrderAverageCell(a_Ux, blockDomain, a_box);
  return 0;
}

/*--------------------------------------------------------------------*/
//  Compute the exact steady state solution \<U\> in the cells
/** The exact solution is for \f$Re = 200\f$, air at 300K
 *  \f$\nu = 1.56e-5 m^2/s\f$ \f$h = 0.1m\f$, solution time = 500s,
 *  the number n in the series is 500.  The exact solution form
 *  (moving wall is on the low side of y direction) for point centered
 *  values is
 *  \f[
 *     U(y) = U_0\frac{y}{h}
 *  \f]
 *  which is then integrated to get the fourth order cell average value.
 *
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update
 *  \param[in]  a_problemDomain
 *                      The problem domain on the level
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Grid level
 *  \return             0 - Successfully computed exact solution
 *                      1 - Exact solution is not known
 *//*-----------------------------------------------------------------*/

int
CNSIBCTransientCouette::steadyState(FArrayBox&              a_Ux,
                                    const Box&              a_box,
                                    const Box&              a_disjointBox,
                                    const LevelGridMetrics& a_gridMetrics,
                                    const Real              a_time,
                                    const int               a_level) const
{
  // Initial state
  const CRDState& state = CRDState::get(m_idxStateInit);

  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box);

  // Working set boxes
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;
  Box initBox(box1Dom);

  // Get physical coordinates
  FABSTACKTEMP(XiFab, box1Dom, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, box1Dom, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(box1Dom, XiFab, XFab, blockCoordSys);

  // Set all the components to state with zero velocity before computing the
  // exact solution.
  Real rho = state.density();
  Real pres = state.pressure();

  a_Ux.setVal(rho, URHO);
  a_Ux.setVal(0., a_Ux.box(), UMOMX, SpaceDim);
  const Real gamma = CRDparam::g_gamma;
  const Real height = CRDparam::g_domainLength[m_wallNormalDir];
  a_Ux.setVal(pres/(CRDparam::g_gamma - 1.), UENG);
  const RealVect wallSpeed = state.velocity();

  // comp is the component to solve for in a_Ux
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  int velDir = m_wallNormalDir - 1;
  if(velDir < 0)
    {
      velDir = SpaceDim - 1;
    }
  if (wallSpeed[velDir] == 0.)
    {
      return 0;
    }
  int velComp = velDir + velIndx;

  const int engComp = UENG;

  // At this point we only need to set momentum parallel to the wall
  // and energy
  MD_ARRAY_RESTRICT(Uvals, a_Ux);
  MD_ARRAY_RESTRICT(Xvals, XFab);
  MD_BOXLOOP(box1Dom, cell)
    {
      Real y = Xvals[MD_IX(cell, m_wallNormalDir)];
      if(m_movWall == 0)
        {
          y = height - y;
        }
      Real vel = wallSpeed[velDir]*(y/height);
      Uvals[MD_IX(cell, velComp)] = rho*vel;
      Uvals[MD_IX(cell, engComp)] = pres/(gamma - 1.) + rho*vel*vel;
    }

  fourthOrderAverageCell(a_Ux, blockDomain, a_box);
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
CNSIBCTransientCouette::readBCInfo()
{
  ParmParse ppIBC("ibc");
  // Determine if steady state or transient
  m_steadyState = false;
  ppIBC.query("steady_state", m_steadyState);
  // set the wall normal direction
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      m_wallNormalDir = 1;
    }
  if (m_wallNormalDir < 0 || m_wallNormalDir >= SpaceDim)
    {
      CRD::msg << "Input (TransientCouette IBC): 'wall_normal_dir' must be "
        ">= 0 and < " << SpaceDim << '!' << CRD::error;
    }
  // Determine if low wall is moving or high wall
  ppIBC.query("moving_wall", m_movWall);

  ppIBC.query("tag_threshold", m_threshold);

  m_readInput = true;
}
