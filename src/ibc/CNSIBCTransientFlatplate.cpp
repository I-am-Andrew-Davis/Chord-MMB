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
 * \file CNSIBCTransientFlatplate.cpp
 *
 * \brief Member functions for CNSIBCTransientFlatplate
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
//**FIXME
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "CanonicalCasesF_F.H"
#include "LoHiCenter.H"

//----- Internal -----//

#include "CNSIBCTransientFlatplate.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "LGintegrator.H"
#include "DataTemp.H"
#include "ChordInput.H"
#include "CNSIBCF_F.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "ViscousTensor4thOrderOp.H"
#include "ViscousTensor4thOrderOpF_F.H"


/*******************************************************************************
 *
 * Class CNSIBCTransientFlatplate: member definitions
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

CNSIBCTransientFlatplate::CNSIBCTransientFlatplate()
  :
  CNSIBCReferenceCubeBC(),
  // All invalid values which must be corrected
  m_farfieldVelocity(-1.*RealVect::Unit),
  m_plateLength(-1.),
  m_maxBLdisplacement(-1.),
  m_wallNormalDir(-1),
  m_flowDir(-1)
{
//--Read any BC info

  readBCInfo();

//--Set BC Type

  // All boundaries periodic by default
  //setAllDomainBC(CRDparam::DomainBCTypePeriodic);

  BoundaryIndex bcIdx;
  bcIdx.m_block = 0; // single block only
  
  
  BCInfo bc;
  bc.m_order = 1;

  // Wall normal dir
  bcIdx.m_dir = m_wallNormalDir;
  bcIdx.m_side = Side::Lo;
  bc.m_type = CRDparam::DomainBCTypeMixed;
  setDomainBC(bcIdx, bc);
  // setDomainBC(m_wallNormalDir,
  //             0,
  //             CRDparam::DomainBCTypeMixed |
  //             CRDparam::DomainBCTypeAdiabaticWall |
  //             CRDparam::DomainBCTypeSlipWall);
  bcIdx.m_side = Side::Hi;
  bc.m_type = CRDparam::DomainBCTypeFarfield;
  setDomainBC(bcIdx, bc);
  // setDomainBC(m_wallNormalDir,
  //             1,
  //             CRDparam::DomainBCTypeFarfield);

  // Flow dir
  bcIdx.m_dir = m_flowDir;
  bcIdx.m_side = Side::Lo;
  bc.m_type = CRDparam::DomainBCTypeInflow;
  setDomainBC(bcIdx, bc);
  // setDomainBC(m_flowDir,
  //             0,
  //             CRDparam::DomainBCTypeInflow);
  bcIdx.m_side = Side::Hi;
  bc.m_type = CRDparam::DomainBCTypeOutflow;
  setDomainBC(bcIdx, bc);
  // setDomainBC(m_flowDir,
  //             1,
  //             CRDparam::DomainBCTypeOutflow);

  // // Order         Direction        Si Order
  // setDomainBCOrder(m_flowDir,       0, 1);  // Left
  // setDomainBCOrder(m_flowDir,       1, 1);  // Right
  // //setDomainBCOrder(m_wallNormalDir, 0, 4);  // Wall
  // setDomainBCOrder(m_wallNormalDir, 0, 1);
  // setDomainBCOrder(m_wallNormalDir, 1, 1);  // Top

//--Define the plate box

  {
    IntVect lo(IntVect::Zero);
    lo[m_flowDir] = CRDparam::g_domainBaseSize[m_flowDir] - m_cellPlateLength;
    IntVect hi(CRDparam::g_domainBaseSize - IntVect::Unit);
    hi[m_wallNormalDir] = 0;
    m_plateBox.define(lo, hi);

    // Adjust the plate box for periodic directions (by a safe 5 cells)
    bcIdx.m_side = Side::Lo;
    for (const auto dir : EachDir)
      {
        bcIdx.m_dir = dir;
        if (getDomainBC(bcIdx).m_type & CRDparam::DomainBCTypePeriodic)
          {
            m_plateBox.grow(dir, 5);
          }
      }
  }

//--Define the state at BC

  const Real rho = CRDparam::g_rho;
  const Real pres = CRDparam::g_rho*CRDparam::g_R*CRDparam::g_T;
  // Inflow
  setReferenceBCState(m_flowDir, 0, WRHO, rho);
  D_TERM(
    setReferenceBCState(m_flowDir, 0, WVELX, m_farfieldVelocity[0]);,
    setReferenceBCState(m_flowDir, 0, WVELY, m_farfieldVelocity[1]);,
    setReferenceBCState(m_flowDir, 0, WVELZ, m_farfieldVelocity[2]);)
  setReferenceBCState(m_flowDir, 0, WPRES, pres);
  // Outflow
  setReferenceBCState(m_flowDir, 1, WRHO, rho);
  D_TERM(
    setReferenceBCState(m_flowDir, 1, WVELX, m_farfieldVelocity[0]);,
    setReferenceBCState(m_flowDir, 1, WVELY, m_farfieldVelocity[1]);,
    setReferenceBCState(m_flowDir, 1, WVELZ, m_farfieldVelocity[2]);)
  setReferenceBCState(m_flowDir, 1, WPRES, pres);
  // Top (farfield)
  setReferenceBCState(m_wallNormalDir, 1, WRHO, rho);
  D_TERM(
    setReferenceBCState(m_wallNormalDir, 1, WVELX, m_farfieldVelocity[0]);,
    setReferenceBCState(m_wallNormalDir, 1, WVELY, m_farfieldVelocity[1]);,
    setReferenceBCState(m_wallNormalDir, 1, WVELZ, m_farfieldVelocity[2]);)
  setReferenceBCState(m_wallNormalDir, 1, WPRES, pres);

//--Set the maximum height of the BL (assumes Blasius)
  const Real ReN = CRDparam::g_Re;
  if (ReN < 5.E5)
    {
      m_maxBLdisplacement = 5./sqrt(CRDparam::g_Re)*m_plateLength;
    }
  else
    {
      m_maxBLdisplacement = 0.37*m_plateLength*pow(CRDparam::g_Re, -1./5.);
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTransientFlatplate::~CNSIBCTransientFlatplate()
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
CNSIBCTransientFlatplate::IBCName() const
{
  return "Navier-Stokes transient flat plate";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCTransientFlatplate::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Flow direction\n" << m_flowDir << CRD::var;
  CRD::msg << "Wall normal direction\n" << m_wallNormalDir << CRD::var;
  CRD::msg << "Farfield velocity\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << m_farfieldVelocity[dir];
    }
  CRD::msg << ')' << CRD::var;
  CRD::msg << "Box defining no-slip plate\n" << m_plateBox << CRD::var;
  CRD::msg << "Discretized no-slip plate length\n" << m_plateLength << CRD::var;
  CRD::msg << "Max height of BL (Blasius)\n" << m_maxBLdisplacement << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
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
CNSIBCTransientFlatplate::setTagMethod(const int a_tagBufferSize)
{
  const int maxAMRLevel = CRDparam::maxAMRLevel();
  // Cell size in wallNormalDir
  const Real dxWallNormal = CRDparam::g_domainLength[m_wallNormalDir]/
    CRDparam::g_domainBaseSize[m_wallNormalDir];
  const int cellMaxBLdisplacement = std::ceil(m_maxBLdisplacement/dxWallNormal);

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
      Box tagBox = grow(m_plateBox, cellMaxBLdisplacement);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[0] = tagLevel;
      break;
    }
    case 2:
    {
      // Two TagLevel
      tagLevelVec.resize(2);
      levelMapVec.resize(2);
      // Apply the first to level 0 and the second to level 1
      levelMapVec[0] = 0;
      levelMapVec[1] = 1;

      // First TagLevel is a box grown by 4 times the max displacement of the BL
      // This is supposed to do a better resolution of the inviscid flow
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(0);
      Box tagBox = grow(m_plateBox, 4*cellMaxBLdisplacement);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // No tag buffer
      tagLevelVec[0] = tagLevel;

      // Second TagLevel is a box grown by the max displacement of the BL
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
      tagBox = grow(m_plateBox, cellMaxBLdisplacement);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[1] = tagLevel;
      break;
    }
    default:
    {
      // Three TagLevel
      tagLevelVec.resize(3);
      levelMapVec.resize(3);
      // Apply the first to level 0, the second to level 1:max-2, and the third
      // to level max-1
      levelMapVec[0] = 0;
      levelMapVec[1] = 1;
      levelMapVec[2] = maxAMRLevel-1;

      // First TagLevel is a box grown by 4 times the max displacement of the BL
      // This is supposed to do a better resolution of the inviscid flow
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
      Box tagBox = grow(m_plateBox, 4*cellMaxBLdisplacement);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // No tag buffer
      tagLevelVec[0] = tagLevel;

      // Middle TagLevel is a box grown by the max displacement of the BL
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
      tagBox = grow(m_plateBox, cellMaxBLdisplacement);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[1] = tagLevel;

      // Last TagLevel attempts better resolution at the first sixth of the
      // plate
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(3);
      tagBox = m_plateBox;
      // tagBox.setBig(m_flowDir,
      //               tagBox.smallEnd(m_flowDir) +
      //               (m_plateBox.bigEnd(m_flowDir) -
      //                m_plateBox.smallEnd(m_flowDir))/6);
      // By making it face-centered, it will only capture the first layer of
      // cells along the boundary on the level of refinement where tagging is
      // performed
      tagBox.shiftHalf(m_wallNormalDir, -1);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(3));
      tagLevelVec[2] = tagLevel;
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
CNSIBCTransientFlatplate::initialize(LevelData<FArrayBox>&      a_U,
                                     LevelGridMetrics&          a_gridMetrics,
                                     const LayoutData<FluxBox>& a_unitNormals,
                                     const Real                 a_time,
                                     const int                  a_level) const
{
  const Real rho      = m_refPrimStateBC(m_wallNormalDir, 1, WRHO);
  const Real pressure = m_refPrimStateBC(m_wallNormalDir, 1, WPRES);
  const Real gamma    = CRDparam::g_gamma;
  const RealVect vel(D_DECL(m_refPrimStateBC(m_wallNormalDir, 1, WVELX),
                            m_refPrimStateBC(m_wallNormalDir, 1, WVELY),
                            m_refPrimStateBC(m_wallNormalDir, 1, WVELZ)));

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
      Box box2 = grow(box, 2);
      box2 &= blockDomain;
      Box box1 = grow(box, 1);
      box1 &= blockDomain;

      // Example if state depends on location
      // FABSTACKTEMP(XiFab, dataBox, SpaceDim);  // Cartesian coordinates
      // FABSTACKTEMP(XFab, dataBox, SpaceDim);   // Physical coordinates
      // getCellCoordinates(dataBox, XiFab, XFab, blockCoordSys);

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2));
      FORT_CNSIBCINIT(CHF_FRA(UFab),
                      CHF_BOX(box2),
                      CHF_CONST_REAL(rho),
                      CHF_CONST_REAL(pressure),
                      CHF_CONST_REALVECT(vel),
                      CHF_CONST_REAL(gamma));

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box1));
      fourthOrderAverageCell(UFab, blockDomain, box1);
    }
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTransientFlatplate::haveExactSol() const
{
  return false;  //**FIXME turn back on when working
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state in the cells (defaults to a
//  no-op)
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
 *
 *  \note
 *  <ul>
 *    <li> Only velocity components are set.  Mass and energy are set
 *         to zero
 *  </ul>
 *//*-----------------------------------------------------------------*/

int
CNSIBCTransientFlatplate::exactSol(FArrayBox&              a_Ux,
                                   const Box&              a_box,
                                   const Box&              a_disjointBox,
                                   const LevelGridMetrics& a_gridMetrics,
                                   const FluxBox&          a_unitNormals,
                                   const DataIndex&        a_didx,
                                   const Real              a_time,
                                   const int               a_level) const
{
  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  // Working set boxes
  Box box2Dom = grow(a_box, 2);
  box2Dom &= blockDomain;
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;

  // Data
  FArrayBox YFab(a_disjointBox, 3);
  YFab.setVal(0.0);

  // Get physical coordinates
  FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

  const RealVect a_dxVect = a_gridMetrics.dxVect();
  const Real dx = a_dxVect[m_wallNormalDir];
  const Real uInfty = m_farfieldVelocity[m_flowDir];
  
  const Real nu = CRDparam::g_mu/CRDparam::g_rho;
  FORT_BLASIUS(CHF_FRA(YFab),
               CHF_BOX(a_disjointBox),
               CHF_CONST_FRA(XFab),
               CHF_CONST_REAL(dx),
               CHF_CONST_REAL(uInfty),
               CHF_CONST_REAL(nu),
               CHF_CONST_INT(m_wallNormalDir));

  return 0;
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In this override, the wallVelocity is not zero on low-side walls
 *  \param[in]  a_W     Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_W     Primitive state corrected for presence of wall
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCTransientFlatplate::setWallBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_boundarySlipVelocity,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
  const int lohiSign = sign(a_side);
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  // Gamma values in cells
  FABSTACKTEMP(gammaCell, a_boundaryFaceBox, 1);
  // Gamma values on faces
  FABSTACKTEMP(gammaFace, a_boundaryFaceBox, 1);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaCell,
                                    a_Wcell);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaFace,
                                    a_Wface);
  if (a_dir == m_wallNormalDir && a_side == Side::Lo)
    // Possibly apply mixed BC here
    {
      CH_assert(m_plateBox.size(a_dir) == 1);

      // Match input box in terms of centering and refinement
      Box boundaryBoxPlate(m_plateBox);
      boundaryBoxPlate.shiftHalf(a_dir, -1);
      boundaryBoxPlate.refine(CRDparam::g_refFromBase[a_level]);

      // This is the part before the no-slip plate
      Box boundaryBoxPrePlate;
      if (boundaryBoxPlate.smallEnd(m_flowDir) >
          a_boundaryFaceBox.smallEnd(m_flowDir))
        // a_boundaryFaceBox is partially or fully before the no-slip plate
        {
          boundaryBoxPrePlate = a_boundaryFaceBox;
          boundaryBoxPrePlate.setBig(
            m_flowDir,
            std::min(a_boundaryFaceBox.bigEnd(m_flowDir),
                     boundaryBoxPlate.smallEnd(m_flowDir) - 1));
        }

      // The part overlapping the plate is simply determined through
      // intersection
      boundaryBoxPlate &= a_boundaryFaceBox;

      RealVect wallVelocity = RealVect::Zero;

      if (!boundaryBoxPrePlate.isEmpty())
        {
          int viscousSlip = 1;
          FORT_CNSIBCPRIMSTATEWALL(CHF_FRA(a_Wface),
                                   CHF_BOX(boundaryBoxPrePlate),
                                   CHF_CONST_FRA(a_Wcell),
                                   CHF_CONST_FRA1(gammaFace,0),
                                   CHF_CONST_FRA1(gammaCell,0),
                                   CHF_CONST_REALVECT(wallVelocity),
                                   CHF_CONST_INT(rhoIndx),
                                   CHF_CONST_INT(presIndx),
                                   CHF_CONST_INT(velIndx),
                                   CHF_CONST_INT(viscousSlip),
                                   CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(lohiSign));
        }
      if (!boundaryBoxPlate.isEmpty())
        {
          int viscousSlip = 0;
          FORT_CNSIBCPRIMSTATEWALL(CHF_FRA(a_Wface),
                                   CHF_BOX(boundaryBoxPlate),
                                   CHF_CONST_FRA(a_Wcell),
                                   CHF_CONST_FRA1(gammaFace,0),
                                   CHF_CONST_FRA1(gammaCell,0),
                                   CHF_CONST_REALVECT(wallVelocity),
                                   CHF_CONST_INT(rhoIndx),
                                   CHF_CONST_INT(presIndx),
                                   CHF_CONST_INT(velIndx),
                                   CHF_CONST_INT(viscousSlip),
                                   CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(lohiSign));
        }
    }
  else
    // Regular wall BC
    {
      int viscousSlip = 0;
      BoundaryIndex bcIdx;
      bcIdx.m_block = 0;
      bcIdx.m_dir = a_dir;
      bcIdx.m_side = a_side;
      if (CRDparam::DomainBCTypeSlipWall & getDomainBC(bcIdx).m_type)
        {
          viscousSlip = 1;
        }

      RealVect wallVelocity = RealVect::Zero;

      FORT_CNSIBCPRIMSTATEWALL(CHF_FRA(a_Wface),
                               CHF_BOX(a_boundaryFaceBox),
                               CHF_CONST_FRA(a_Wcell),
                               CHF_CONST_FRA1(gammaFace,0),
                               CHF_CONST_FRA1(gammaCell,0),
                               CHF_CONST_REALVECT(wallVelocity),
                               CHF_CONST_INT(rhoIndx),
                               CHF_CONST_INT(presIndx),
                               CHF_CONST_INT(velIndx),
                               CHF_CONST_INT(viscousSlip),
                               CHF_CONST_INT(a_dir),
                               CHF_CONST_INT(lohiSign));
    }
}

/*--------------------------------------------------------------------*/
//  Set the boxes that are inflow or outflow if mixed BC is being used
/** \param[in]  a_boundaryFaceBox
 *                      Box of boundary being operated on
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_boxVect
 *                      Vector of boxes
 *//*-----------------------------------------------------------------*/

int
CNSIBCTransientFlatplate::setMixedBC(
  const Box&            a_boundaryFaceBox,
  const int             a_dir,
  const Side::LoHiSide& a_side,
  const Box&            a_disjointBox,
  LevelGridMetrics&     a_gridMetrics,
  const Real            a_time,
  const int             a_level,
  Vector<Box>&          a_boxVect,
  Vector<BCInfo>&       a_domainType) const
{
  // NOTE: this is not the intended use for mixed boundaries, just a fix
  a_boxVect.resize(1);
  a_domainType.resize(1);
  a_boxVect[0].define(a_boundaryFaceBox);
  BCInfo bc;
  bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  bc.m_order = 1;
  a_domainType[0] = bc;
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
CNSIBCTransientFlatplate::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error
  ParmParse ppIBC("ibc");

//--Wall normal direction

  m_wallNormalDir = SpaceDim - 1;
  ppIBC.query("wall_normal_dir", m_wallNormalDir);
  if (m_wallNormalDir < 0 || m_wallNormalDir >= SpaceDim)
    {
      CRD::msg << "Input (TransientFlatplate IBC): 'wall_normal_dir' must be "
        ">= 0 and < " << SpaceDim << '!' << CRD::error;
    }

//--Flow direction

  m_flowDir = 0;

//--Plate length

  m_plateLength = CRDparam::g_domainLength[0];
  ppIBC.query("plate_length", m_plateLength);
  if (m_plateLength <= 0. || m_plateLength > CRDparam::g_domainLength[0])
    {
      CRD::msg << "Input (TransientFlatplate IBC): 'plate_length' must be "
        ">= 0 and < " << CRDparam::g_domainLength[0] << '!' << CRD::error;
    }
  // Discretize the plate length in terms of base cells
  const Real flowDx = CRDparam::g_domainLength[m_flowDir]/
    CRDparam::g_domainBaseSize[m_flowDir];
  m_cellPlateLength = std::min((int)(m_plateLength/flowDx),
                               CRDparam::g_domainBaseSize[m_flowDir]);
  m_plateLength = m_cellPlateLength*flowDx;

//--Farfield velocity (and Re)

  if (ppIBC.contains("farfield_velocity"))
    {
      std::vector<Real> IBCfarfieldVelocity(SpaceDim);
      ppIBC.getarr("farfield_velocity", IBCfarfieldVelocity, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_farfieldVelocity.dataPtr(),
                                               &IBCfarfieldVelocity.front());
      // If the velocity is specified, modify the Re
      CRDparam::CRDP.adjustFluidSpeed(m_farfieldVelocity.vectorLength(),
                                     m_plateLength);
    }
  else
    {
      m_farfieldVelocity =
        RealVect(D_DECL(CRDparam::g_Re*CRDparam::g_mu/
                        (CRDparam::g_rho*m_plateLength),
                        0.,
                        0.));
      // If Re is specified, modify the speed
      CRDparam::CRDP.adjustFluidRe(CRDparam::g_Re, m_plateLength);
    }

  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}




