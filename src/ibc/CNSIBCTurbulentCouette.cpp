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
 * \file CNSIBCTurbulentCouette.cpp
 *
 * \brief Member functions for CNSIBCTurbulentCouette
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "NumericIntegral.H"
#include "CONSTANTS.H"
 
//----- Internal -----//

#include "CNSIBCTurbulentCouette.H"
#include "CNSIBCTransientCouetteF_F.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "ChordInput.H"

/*******************************************************************************
 *
 * Class CNSIBCTurbulentCouette: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCTurbulentCouette::CNSIBCTurbulentCouette()
  :
  CNSIBCGeneralized(),
  m_tagPct(0.)
{
  CH_assert(false);
  //**FIXME This needs to be update to the new BC format (address other fixmes
  //**      in this file as well)
#if 0
  readBCInfo();
  int wallNormDir = 0;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      if (m_domainBC[dir][0].type & CRDparam::DomainBCTypeAllWall)
        {
          wallNormDir = dir;
          break;
        }
    }
  const int wallPct = std::ceil(m_tagPct*
                                CRDparam::g_domainBaseSize[wallNormDir]);
  {
    // Set lower box for tagging, grow by 5 in periodic directions
    IntVect lo(IntVect::Zero);
    IntVect hi(CRDparam::g_domainBaseSize + 5*IntVect::Unit);
    hi[wallNormDir] = wallPct;
    m_loTagBox.define(lo, hi);
  }
  {
    // Set upper box for tagging, grow by 5 in periodic directions
    IntVect lo(IntVect::Zero);
    IntVect hi(CRDparam::g_domainBaseSize + 5*IntVect::Unit);
    hi[wallNormDir] = CRDparam::g_domainBaseSize[wallNormDir] - 1;
    lo[wallNormDir] = CRDparam::g_domainBaseSize[wallNormDir] - 1 - wallPct;
    m_hiTagBox.define(lo, hi);
  }
#endif
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTurbulentCouette::~CNSIBCTurbulentCouette()
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
CNSIBCTurbulentCouette::IBCName() const
{
  return "Turbulent Couette flow";
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/**
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
CNSIBCTurbulentCouette::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  CNSIBCCombustionReference::setTagMethodLevel(a_tagBufferSize, tagLevel);
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  tagLevel->appendTagMethod(new TagMethodBaseBox(m_loTagBox));
  tagLevel->appendTagMethod(new TagMethodBaseBox(m_hiTagBox));
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  return new TagLevelFactory(tagLevel);
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
CNSIBCTurbulentCouette::initialize(LevelData<FArrayBox>&      a_U,
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

      // Get domain for the block
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
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTurbulentCouette::haveExactSol() const
{
  return false;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/** No exact solution exists for turbulent couette flow, but an
 *  analytic one does exist which gives a reasonable attempt at a
 *  steady state solution. It is not perfect however and will take
 *  some time to reach true steady state.
 *
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
CNSIBCTurbulentCouette::exactSol(FArrayBox&              a_Ux,
                                 const Box&              a_box,
                                 const Box&              a_disjointBox,
                                 const LevelGridMetrics& a_gridMetrics,
                                 const FluxBox&          a_unitNormals,
                                 const DataIndex&        a_didx,
                                 const Real              a_time,
                                 const int               a_level) const
{
  const int nTerms = 500;
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));

  // Get domain for the block
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
  
  // Set all the components to start with zero velocity before computing the
  // exact solution.
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  FABSTACKTEMP(Wc, initBox, numWVar);
  Real pres = m_initP;
  Real temp = m_initT;
  Real rho = m_initRho;
  Real Rval = CRDparam::g_R;
  if (pres < 0. && rho > 0)
    {
      pres = rho*Rval*temp;
    }
  else if (rho < 0. && pres > 0.)
    {
      rho = pres/(Rval*temp);
    }
  //**FIXME: this setup only works with walls in one direction: if there are
  //         more than two directions with walls, only one of them is used
  //         and it may not be used correctly
  int wallNormalDir = 0;
//**FIXME (new BC setup)
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     if (m_domainBC[dir][0].type & CRDparam::DomainBCTypeAllWall)
  //       {
  //         wallNormalDir = dir;
  //       }
  //   }
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      CRD::msg << "Cannot use thermally perfect physics with this problem"
               << CRD::error;
    }
  
  const Real height = CRDparam::g_domainLength[wallNormalDir];
  Wc.setVal(m_initRho, rhoIndx);
  Wc.setVal(0., box1Dom, velIndx, SpaceDim);
  Wc.setVal(m_initP, presIndx);
  Wc.setVal(m_initT, tempIndx);

  Real wallSpeedRel = 0.;
//  Real wallSpeedRelTemp = 0.;
  int  wallDir = 0;
  Real wallSpeedLo = 0;
  // Real wallSpeedLo, wallSpeedHi, wallSpeedLoTemp, wallSpeedHiTemp;
  // Set to 0 if lower wall is moving and 1 if upper wall is moving
 //  for (int dir = 0; dir != SpaceDim; ++dir)
//     {
//       if (dir != wallNormalDir)
//         {
// //**FIXME (new BC setup)
//           // wallSpeedLoTemp  = m_loBCVel[wallNormalDir][dir];
//           // wallSpeedHiTemp  = m_hiBCVel[wallNormalDir][dir];
   
//           wallSpeedRelTemp = wallSpeedHiTemp - wallSpeedLoTemp;
//           if (wallSpeedRelTemp != 0.)
//             {
//               wallSpeedRel = wallSpeedRelTemp;
//               wallSpeedLo = wallSpeedLoTemp;
//               wallSpeedHi = wallSpeedHiTemp;
//               wallDir = dir;
//             }
//         }
//     }
  // comp is the component to solve for in a_Ux
  int velComp = velIndx + wallDir;

  // solve for the analytic solution of the turbulent Couette flow
  {
    Real T = 17.5; // a empirical solution constant
    Real f = 2.3041e-04; // for Re = 100000
    Real tau = (8.0/f)*std::pow(T/CRDparam::g_Re, 2.0);
    Real y, yp, phi, cn, c4, c, a, ta2;
    // first term in solution is phi = normalized height
    Wc.copy(XFab, wallNormalDir, velComp);
    Wc.mult(1./height, velComp);
    MD_ARRAY_RESTRICT(arrW, Wc);
    MD_ARRAY_RESTRICT(arrY, XFab);
    for (int n = 1; n <= nTerms; ++n)
      {
        const anlySolTermC4& eq_c4 = anlySolTermC4(f, n);
        // integrate from 0 to 1, except the function approaches -Inf at
        // 0 so use the 0+tolerance instead
        c4 = NumericIntegral::adaptiveSimpsons(
          eq_c4,
          1.0e-6, // lo
          1.,     // hi
          NumericIntegral::NITr<Real>::tolerance(),
          10000);
        a = n*Pi;
        cn = 2.0*(std::pow(-1, n)/a + c4);
        ta2 = tau*std::pow(a, 2.0);
        c = (1.0-std::exp(-ta2))*(cn/ta2);
        MD_BOXLOOP(initBox, i)
          {
            y = arrY[MD_IX(i, wallNormalDir)];
            yp = std::min(y, height - y)/height;
            phi = std::sin(a*yp)*c;
            // mirror solution
            if ((y/height) > 0.5)
            {
              phi = -phi;
            }
            arrW[MD_IX(i, velComp)] += phi;
          }
      }
    // scale up
    Wc.mult(wallSpeedRel, velComp);
    Wc.plus(wallSpeedLo, velComp);
  }
  
  // convert to cell averaged conservative solution
  CRDparam::g_CRDPhysics->initialize(a_Ux,
                                     Wc,
                                     a_gridMetrics,
                                     a_unitNormals,
                                     a_didx,
                                     a_disjointBox,
                                     initBox);
  fourthOrderAverageCell(a_Ux, blockDomain, a_box);
  return 0;
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
CNSIBCTurbulentCouette::readBCInfo()
{
  ParmParse ppIBC("ibc");
  // Percent into the domain to tag the region near the walls
  m_tagPct = 0.1;
  ppIBC.query("tag_pct", m_tagPct);
}
