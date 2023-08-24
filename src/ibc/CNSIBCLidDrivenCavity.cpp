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
 * \file CNSIBCLidDrivenCavity.cpp
 *
 * \brief Member functions for CNSIBCLidDrivenCavity
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCLidDrivenCavity.H"
#include "CNSIBCCombustionTestF_F.H"
#include "CRDPhysics.H"
#include "CNSIBCF_F.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"


/*******************************************************************************
 *
 * Class CNSIBCLidDrivenCavity: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCLidDrivenCavity::CNSIBCLidDrivenCavity()
  :
  CNSIBCReferenceCubeBC(),
  m_wallDir(-1),
  m_wallNormDir(-1),
  m_wallVel(-1.),
  m_wallOrder(-1)
{
  readBCInfo();
  BCInfo bcType;
  bcType.m_order = m_wallOrder;
  bcType.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  setAllDomainBC(bcType);
  if (SpaceDim == 3)
    {
      BoundaryIndex bcIdx;
      bcIdx.m_block = 0;
      bcIdx.m_dir = 2;
      bcIdx.m_side = Side::Lo;
      BCInfo bc = getDomainBC(bcIdx);
      bc.m_type = CRDparam::DomainBCTypeWall;
      setDomainBC(bcIdx, bc);
      bcIdx.m_side = Side::Hi;
      bc = getDomainBC(bcIdx);
      bc.m_type = CRDparam::DomainBCTypeWall;
      setDomainBC(bcIdx, bc);
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCLidDrivenCavity::~CNSIBCLidDrivenCavity()
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
CNSIBCLidDrivenCavity::IBCName() const
{
  return "Mixing lid-driven cavity flow test";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCLidDrivenCavity::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Wall direction\n" << m_wallDir << CRD::var;
  CRD::msg << "Wall normal direction\n" << m_wallNormDir << CRD::var;
  CRD::msg << "Wall velocity\n" << m_wallVel << CRD::var;
  CRD::msg << "Initial mass fractions\n(";
  CRD::msg << m_leftMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_leftMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Right side mass fractions\n(";
  CRD::msg << m_rightMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_rightMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Wall extrapolation order\n" << m_wallOrder << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** For the shock tube, tags are set based on gradients of density.
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
CNSIBCLidDrivenCavity::setTagMethod(const int a_tagBufferSize)
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
        lo[m_wallNormDir] = std::ceil(CRDparam::g_domainLength[m_wallNormDir]/
                                      dxVect[m_wallNormDir]*0.85);
        IntVect hi = CRDparam::g_domainBaseSize;
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
      lo[m_wallNormDir] = std::ceil((CRDparam::g_domainLength[m_wallNormDir]/
                                     dxVect[m_wallNormDir])*0.85);
      IntVect hi = CRDparam::g_domainBaseSize;
      Box tagBox(lo,hi);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[0] = tagLevel;

      // Second level
      tagLevel = new TagLevel;
      lo[m_wallNormDir] = std::ceil(0.9*CRDparam::g_domainLength[m_wallNormDir]/
                                    dxVect[m_wallNormDir]);
      hi[m_wallDir] = std::floor(CRDparam::g_domainLength[m_wallDir]/
                                 dxVect[m_wallDir]*1/8);
      tagBox.define(lo,hi);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      lo[m_wallDir] = std::ceil(CRDparam::g_domainLength[m_wallDir]/
                               dxVect[m_wallDir]*7/8)-1;
      hi = CRDparam::g_domainBaseSize;
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
CNSIBCLidDrivenCavity::initialize(LevelData<FArrayBox>&      a_U,
                                  LevelGridMetrics&          a_gridMetrics,
                                  const LayoutData<FluxBox>& a_unitNormals,
                                  const Real                 a_time,
                                  const int                  a_level) const
{
  // Define the right side box
  IntVect lo(-5*IntVect::Unit);
  lo[m_wallDir] = CRDparam::g_domainBaseSize[m_wallDir]*0.5;
  IntVect hi(CRDparam::g_domainBaseSize);
  hi += IntVect::Unit*5;
  Box rightBox(lo, hi);
  rightBox.refine(CRDparam::g_refFromBase[a_level]);
  // Define initial values
  const Real rho = CRDparam::g_rho;
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

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(rho, rComp);
      Wc.setVal(0., box2Dom, WvelIndx, SpaceDim);
      Wc.setVal(101325., pComp);
      Wc.setVal(-1., tComp);
      for (int i = 0; i != numSpecies; ++i)
        {
          int comp = wCompStart + i;
          Wc.setVal(m_leftMassFraction[i], comp);
        }
      Box mixingBox(box2Dom);
      mixingBox&=rightBox;
      if (!mixingBox.isEmpty() && numSpecies > 1)
        {
          for (int i = 0; i != numSpecies; ++i)
            {
              int comp = wCompStart + i;
              Wc.setVal(m_rightMassFraction[i], mixingBox,  comp);
            }
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      
      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box1Dom));
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
CNSIBCLidDrivenCavity::haveExactSol() const
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
CNSIBCLidDrivenCavity::exactSol(FArrayBox&              a_Ux,
                                const Box&              a_box,
                                const Box&              a_disjointBox,
                                const LevelGridMetrics& a_gridMetrics,
                                const FluxBox&          a_unitNormals,
                                const DataIndex&        a_didx,
                                const Real              a_time,
                                const int               a_level) const
{
  a_Ux.setVal(0.5);
  return 0;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In this override, the wallVelocity is not zero on low-side walls
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCLidDrivenCavity::setWallBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_boundarySlipVelocity,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_domT) const
{
  const int lohiSign = sign(a_bcIdx.m_side);
  RealVect wallVelocity = RealVect::Zero;
  if (a_bcIdx.m_side == Side::Hi && a_bcIdx.m_dir == m_wallNormDir)
    {
      wallVelocity[m_wallDir] = m_wallVel;
    }
  int viscousSlip = 0;
  if (CRDparam::DomainBCTypeSlipWall & a_domT.m_type)
    {
      viscousSlip = 1;
    }
  // Gamma values in cells
  FABSTACKTEMP(gammaCell, a_boundaryFaceBox, 1);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaCell,
                                    a_Wcell);
  computeWallPrimState(a_Wface,
                       a_boundaryFaceBox,
                       a_Wcell,
                       a_boundarySlipVelocity,
                       gammaCell,
                       a_unitNormalBasisFab,
                       wallVelocity,
                       viscousSlip,
                       a_bcIdx.m_dir,
                       lohiSign);
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
CNSIBCLidDrivenCavity::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  m_wallDir = 0;
  ppIBC.query("wall_dir",m_wallDir);
  m_wallNormDir = m_wallDir + 1;
  if (m_wallNormDir > SpaceDim-1)
    {
      m_wallNormDir = 0;
    }
  // We assume the whole domain is the first species listed in the input file
  // unless otherwise specified with left or right_mass_fractions
  m_leftMassFraction.resize(numSpecies);
  m_leftMassFraction.assign(numSpecies,0.);
  m_leftMassFraction[0] = 1.;
  if (ppIBC.contains("left_mfs"))
    {
      int massTest = assignMassFractions(m_leftMassFraction,
                                         "left_specs",
                                         "left_mfs");
      if (massTest == 1)
        {
          CRD::msg << "Input (LidDrivenCavity IBC): 'left_mfs'"
                   << " must be equal to 1!" << CRD::error;
        }
    }
  m_rightMassFraction.resize(numSpecies);
  m_rightMassFraction.assign(numSpecies,0.);
  m_rightMassFraction[0] = 1.;
  if (ppIBC.contains("right_mfs"))
    {
      int massTest = assignMassFractions(m_rightMassFraction,
                                         "right_specs",
                                         "right_mfs");
      if (massTest == 1)
        {
          CRD::msg << "Input (LidDrivenCavity IBC): 'right_mfs'"
                   << " must be equal to 1!" << CRD::error;
        }
    }
  m_wallVel = 1.;
  ppIBC.query("wall_velocity", m_wallVel);
//--Order of extrapolation to use at the wall
  m_wallOrder = 4;
  ppIBC.query("wall_order", m_wallOrder);
  if (m_wallOrder != 1 && m_wallOrder != 4)
    {
      CRD::msg << "Input (LidDrivenCavity IBC): 'wall_order' must be 1 or 4 "
               << CRD::error;
    }
  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
