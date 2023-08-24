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
 * \file CNSIBCTransientPoiseuille.cpp
 *
 * \brief Member functions for CNSIBCTransientPoiseuille
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "CRDparam.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
//**FIXME
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCTransientPoiseuille.H"
#include "CNSIBCTransientPoiseuilleF_F.H"
#include "CNSIBCF_F.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"

/*******************************************************************************
 *
 * Class CNSIBCTransientPoiseuille: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_wallNormalDir
 *                      Direction normal to the wall.  Default is
 *                      SpaceDim-1
 *//*-----------------------------------------------------------------*/

CNSIBCTransientPoiseuille::CNSIBCTransientPoiseuille(
  const int a_wallNormalDir)
:
  CNSIBCReferenceCubeBC(),
  m_wallNormalDir(a_wallNormalDir),
  m_forceDir(-1),
  m_forceVal(1.),
  m_refPerc(0.5)
{
  readBCInfo();
  const Real reynolds = CRDparam::g_Re;
  if (reynolds > 2000.)
    {
      CRD::msg << "Input (TransientPoiseuille IBC): 'Re' must be "
        "> 0 and < 2000" << '!' << CRD::error;
    }
  const Real nu = CRDparam::g_mu/CRDparam::g_rho;
  Real height = CRDparam::g_domainLength[m_wallNormalDir];
  m_forceVal = 12.*reynolds*nu*nu/(pow(height,3));
  m_U = m_forceVal*height*height/(12.*nu);
  // Default BC type is periodic 
  //setAllDomainBC(CRDparam::DomainBCTypePeriodic);
  BoundaryIndex bcIdx;
  bcIdx.m_block = 0; // single block only
  bcIdx.m_dir = m_wallNormalDir;
  
  BCInfo bc;
  bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  bc.m_order = 1;
  
  // Low side is stationary
  bcIdx.m_side = Side::Lo;
  setDomainBC(bcIdx, bc);
  // Hi side is stationary
  bcIdx.m_side = Side::Hi;
  setDomainBC(bcIdx, bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTransientPoiseuille::~CNSIBCTransientPoiseuille()
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
CNSIBCTransientPoiseuille::IBCName() const
{
  return "Navier-Stokes transient Poiseuille";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCTransientPoiseuille::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Wall normal direction\n" << m_wallNormalDir << CRD::var;
  CRD::msg << "Forcing direction\n" << m_forceDir << CRD::var;
  CRD::msg << "Forcing value\n" << m_forceVal << CRD::var;
  CRD::msg << "Bulk velocity\n" << m_U << CRD::var;
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
CNSIBCTransientPoiseuille::initialize(LevelData<FArrayBox>&      a_U,
                                      LevelGridMetrics&          a_gridMetrics,
                                      const LayoutData<FluxBox>& a_unitNormals,
                                      const Real                 a_time,
                                      const int                  a_level) const
{
  const Real rho      = m_refPrimStateBC(m_wallNormalDir, 1, WRHO);
  const Real pressure = rho*CRDparam::g_R*CRDparam::g_T;
  const Real gamma    = CRDparam::g_gamma;

  const RealVect vel(D_DECL(0.,0.,0.));
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
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      UFab.setVal(rho, URHO);
      UFab.setVal(pressure/(gamma - 1.), UENG);
      CH_assert(UFab.box().contains(box2Dom));
      if (a_time <= 0.)
        {                  
          FORT_CNSIBCINIT(CHF_FRA(UFab),
                          CHF_BOX(box2Dom),
                          CHF_CONST_REAL(rho),
                          CHF_CONST_REAL(pressure),
                          CHF_CONST_REALVECT(vel),
                          CHF_CONST_REAL(gamma));
          // Average values of U (required on valid +1 except at physical
          // boundaries)
          CH_assert(UFab.box().contains(box1Dom));
          fourthOrderAverageCell(UFab, blockDomain, box1Dom);
        }
      else
        {
          this->exactSol(UFab,
                         box2Dom,
                         box,
                         a_gridMetrics,
                         a_unitNormals[dit],
                         dit(),
                         a_time,
                         a_level);
        }
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
CNSIBCTransientPoiseuille::setTagMethod(const int a_tagBufferSize)
{
  const int maxAMRLevel = CRDparam::maxAMRLevel();
  TagLevel* tagLevel;
  std::vector<TagLevel*> tagLevelVec;
  std::vector<int> levelMapVec;

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
      IntVect lo(IntVect::Zero);
      IntVect hi(CRDparam::g_domainBaseSize);
      hi[m_wallNormalDir] = CRDparam::g_domainBaseSize[m_wallNormalDir]
        *m_refPerc;
      {  
        Box tagBox(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      }
      lo = IntVect::Zero;
      lo[m_wallNormalDir] = (CRDparam::g_domainBaseSize[m_wallNormalDir])
        *(1. - m_refPerc);
      hi = CRDparam::g_domainBaseSize;
      {
        Box tagBox(lo,hi);
        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      }
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[0] = tagLevel;
      break;
    }
    }

    // Return the factory
    return new TagLevelFactory(tagLevelVec, levelMapVec);
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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCTransientPoiseuille::addSourceTerm(
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
  const DataIndex&     a_dataIndx,
  const Real           a_globalKE,
  const Real           a_globalHelicity) const
{
  // Add the source terms to a_sourceFab
  int momComp = m_forceDir + 1;
  const int engIndx = CRDparam::g_CRDPhysics->energyFluxIndex();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  FORT_ADDPOISEUILLEFORCE(CHF_FRA(a_sourceFab),
                          CHF_BOX(a_solveBox),
                          CHF_CONST_FRA(a_Wcell),
                          CHF_CONST_REAL(m_forceVal),
                          CHF_CONST_INT(rhoIndx),
                          CHF_CONST_INT(engIndx),
                          CHF_CONST_INT(momComp));
  return;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTransientPoiseuille::haveExactSol() const
{
  return true;
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
CNSIBCTransientPoiseuille::exactSol(FArrayBox&              a_Ux,
                                    const Box&              a_box,
                                    const Box&              a_disjointBox,
                                    const LevelGridMetrics& a_gridMetrics,
                                    const FluxBox&          a_unitNormals,
                                    const DataIndex&        a_didx,
                                    const Real              a_time,
                                    const int               a_level) const
{
  const int maxIter = 2000;

  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  // Working set boxes
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;

  // Get physical coordinates
  FABSTACKTEMP(XiFab, box1Dom, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, box1Dom, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(box1Dom, XiFab, XFab, blockCoordSys);
  
  // Set all the components to state with zero velocity before computing the
  // exact solution.
  a_Ux.setVal(CRDparam::g_rho, URHO);
  a_Ux.setVal(0., a_Ux.box(), UMOMX, SpaceDim);
  const Real rho = CRDparam::g_rho;
  const Real nu = CRDparam::g_mu/CRDparam::g_rho;
  Real height = CRDparam::g_domainLength[m_wallNormalDir];
  // comp is the component to solve for in a_Ux
  int comp = -1;
  D_TERM(
  if (m_forceDir == 0)
    {
      comp = UMOMX;
    },
  else if (m_forceDir == 1)
    {
      comp = UMOMY;
    },
  else if (m_forceDir == 2)
    {
      comp = UMOMZ;
    });
  const Real presVal = rho*CRDparam::g_R*CRDparam::g_T;
  const int eComp = CRDparam::g_CRDPhysics->energyFluxIndex();
  FORT_CNSIBCPOISEUILLEEXACTSOL(
    CHF_FRA(a_Ux),
    CHF_BOX(box1Dom),
    CHF_CONST_FRA1(XFab, m_wallNormalDir),
    CHF_CONST_INT(comp),
    CHF_CONST_INT(eComp),
    CHF_CONST_REAL(a_time),
    CHF_CONST_REAL(height),
    CHF_CONST_REAL(nu),
    CHF_CONST_REAL(rho),
    CHF_CONST_REAL(presVal),
    CHF_CONST_REAL(m_forceVal),
    CHF_CONST_INT(maxIter));

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
CNSIBCTransientPoiseuille::readBCInfo()
{
  CH_assert(!m_readInput);
  ParmParse ppIBC("ibc");
  ppIBC.query("wall_normal_dir", m_wallNormalDir);
  if (m_wallNormalDir < 0 || m_wallNormalDir >= SpaceDim)
    {
      CRD::msg << "Input (TransientPoiseuille IBC): 'wall_normal_dir' must be "
        ">= 0 and < " << SpaceDim << '!' << CRD::error;
    }
  ppIBC.query("force_dir", m_forceDir);
  if (m_forceDir < 0 || m_forceDir >= SpaceDim)
    {
      CRD::msg << "Input (TransientPoiseuille IBC): 'force_dir' must be "
        ">= 0 and < " << SpaceDim << '!' << CRD::error;
    }
  m_refPerc = 0.5;
  ppIBC.query("refine_perc", m_refPerc);
  if (m_refPerc > 0.5)
    {
      CRD::msg << "Input (TransientPoiseuille IBC): 'refine_perc' must be "
        "< 0.5'!" << CRD::error;
    }
  m_readInput = true;
}
