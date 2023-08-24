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
 * \file CNSIBCSpecMachReflection.cpp
 *
 * \brief Member functions for CNSIBCSpecMachReflection
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

#include "CNSIBCSpecMachReflection.H"
#include "CNSIBCMachReflectionF_F.H"
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
 * Class CNSIBCSpecMachReflection: member definitions
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

CNSIBCSpecMachReflection::CNSIBCSpecMachReflection()
:
  CNSIBCReferenceCubeBC(),
  m_lowXBound(0),
  m_hiYBound(0),
  //**FIXME These should be read from input once IBC can construct CS
  m_alpha(PI*30./180.),
  m_xLead(0.12),
  m_xRamp(0.4),
  m_xBase(-0.05),
  m_xloclow(-0.05),
  // All invalid values which must be corrected
  m_p0(-1.),
  m_r0(-1.),
  m_Mach(-1.),
  m_threshold(-1.),
  m_wallOrder(4),
  m_tagComp(-1),
  m_tagWall(false)
{
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }

//--Read any BC info

  readBCInfo();

//--Set BC Type

  setAllDomainBC(CRDparam::DomainBCTypePeriodic);
  // Set lower x boundary
  if (m_lowXBound == 0)
    {
      setDomainBC(0, 0, CRDparam::DomainBCTypeDirichlet);
    }
  else if (m_lowXBound == 1)
    {
      setDomainBC(0, 0, CRDparam::DomainBCTypeExtrapolated);
    }
  // Set upper y boundary
  if (m_hiYBound == 0)
    {
      setDomainBC(1, 1, CRDparam::DomainBCTypeDirichlet);
    }
  else if (m_hiYBound == 1)
    {
      setDomainBC(1, 1, CRDparam::DomainBCTypeExtrapolated);
    }
  else if (m_hiYBound == 2)
    {
      setDomainBC(1, 1, CRDparam::DomainBCTypeFarfield);
    }
  setDomainBC(0, 1, CRDparam::DomainBCTypeDirichlet);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
    {
      setDomainBC(1, 0, CRDparam::DomainBCTypeAdiabaticWall);
    }
  else
    {
      setDomainBC(1, 0, CRDparam::DomainBCTypeWall);
    }

  // Order         Di Si Order
  setDomainBCOrder(0, 0, 1);            // Left
  setDomainBCOrder(0, 1, 1);            // Right
  setDomainBCOrder(1, 0, m_wallOrder);  // Wall
  setDomainBCOrder(1, 1, 1);            // Top
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCSpecMachReflection::~CNSIBCSpecMachReflection()
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
CNSIBCSpecMachReflection::IBCName() const
{
  return "Multispecies Mach reflection";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpecMachReflection::writeIBCInfo() const
{
  if (m_lowXBound == 0)
    {
      CRD::msg << "Lower X boundary\n" << "Dirichlet" << CRD::var;
    }
  else if (m_lowXBound == 1)
    {
      CRD::msg << "Lower X boundary\n" << "Extrapolated" << CRD::var;
    }
  if (m_hiYBound == 0)
    {
      CRD::msg << "Upper Y boundary\n" << "Dirichlet" << CRD::var;
    }
  else if (m_hiYBound == 1)
    {
      CRD::msg << "Upper Y boundary\n" << "Extrapolated" << CRD::var;
    }
  else if (m_hiYBound == 2)
    {
      CRD::msg << "Upper Y boundary\n" << "Farfield" << CRD::var;
    }
  CRD::msg << "Angle of ramp (deg)\n" << 180.*m_alpha/PI << CRD::var;
  CRD::msg << "Lead length before ramp\n" << m_xLead << CRD::var;
  CRD::msg << "Ramp length (projected along x-axis)\n" << m_xRamp << CRD::var;
  CRD::msg << "x-location of corner\n" << 0.0 << CRD::var;
  CRD::msg << "Starting x-location of shock\n" << m_xBase << CRD::var;
  CRD::msg << "End of slip wall portion\n" << m_xloclow << CRD::var;
  CRD::msg << "Shock Mach number\n" << m_Mach << CRD::var;
  CRD::msg.setPrecFloatSN(5);
  CRD::msg << "Conditions in front of shock:" << CRD::body;
  CRD::msg << "Pressure\n" << m_p0 << CRD::var;
  CRD::msg << "Density\n" << m_r0 << CRD::var;
  CRD::msg << "Sound speed\n" << m_c0 << CRD::var;
  CRD::msg << "Conditions behind the shock:" << CRD::body;
  CRD::msg << "Pressure\n" << m_p1 << CRD::var;
  CRD::msg << "Density\n" << m_r1 << CRD::var;
  CRD::msg << "Velocity\n" << m_u1 << CRD::var;
  CRD::msg << "Order of accuracy at wall BC\n" << m_wallOrder << CRD::var;
  CRD::msg << "Tagging threshold\n" << m_threshold << CRD::var;
  CRD::msg << "Tagging component\n" << m_tagComp << CRD::var;
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
CNSIBCSpecMachReflection::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  // Tag the component gradient
  tagLevel->appendTagMethod(new TagMethodGradient(m_tagComp,
                                                  m_threshold,
                                                  true));
  // Add in tag buffer
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  // Tag the bottom boundary
  if (m_tagWall)
    {
      Box tagBox(IntVect(D_DECL(4, 0, 0)),
                 IntVect(D_DECL(CRDparam::g_domainBaseSize[0], 0,
                                CRDparam::g_domainBaseSize[2])));
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
    }
  // Return the factory
  return new TagLevelFactory(tagLevel);
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
CNSIBCSpecMachReflection::initialize(LevelData<FArrayBox>&      a_U,
                                     LevelGridMetrics&          a_gridMetrics,
                                     const LayoutData<FluxBox>& a_unitNormals,
                                     const Real                 a_time,
                                     const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
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

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);
      // Pointwise values of W (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FORT_CNSIBCSPECMACHREFLECTIONINIT(CHF_FRA(Wc),
                                        CHF_BOX(box2Dom),
                                        CHF_CONST_FRA(XFab),
                                        CHF_CONST_INT(rComp),
                                        CHF_CONST_INT(pComp),
                                        CHF_CONST_INT(tComp),
                                        CHF_CONST_INT(WvelIndx),
                                        CHF_CONST_INT(wCompStart),
                                        CHF_CONST_INT(numSpecies),
                                        CHF_CONST_REAL(m_xBase),
                                        CHF_CONST_REAL(m_p0),
                                        CHF_CONST_REAL(m_r0),
                                        CHF_CONST_REAL(m_p1),
                                        CHF_CONST_REAL(m_r1),
                                        CHF_CONST_REAL(m_u1),
                                        CHF_CONST_VR(m_downstreamMF),
                                        CHF_CONST_VR(m_upstreamMF));
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
#endif
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the exterior (or farfield) primitive state at flow BC
/** This override is used to impose conditions before and after the
 *  shock in the "farfield" boundary at ymax.  Otherwise the base
 *  scheme is used to set constant values.
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
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
CNSIBCSpecMachReflection::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
#if (CH_SPACEDIM != 1)
  // Set all the component variables and intervals
  const int numSpecies = CRDparam::g_numSpecies;
  Interval velInt = CRDparam::g_CRDPhysics->velocityInterval();
  const int velComp = velInt.begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  Interval specInt = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int speciesComp = specInt.begin();
  const int lohiSign = sign(a_side);

  if (a_dir == 1 && a_side == Side::Hi)
    {
      // Get a disjoint interior box
      Box interiorBox = adjCellBox(a_boundaryFaceBox,
                                   a_dir,
                                   Side::flip(a_side),
                                   1);
      //**FIXME - this should be single-block so any box should work.  Ideally
      //**        the parent block domain should be an annotation in the box
      //**FIXME - this should work as long as min box size is >= 8
      // Grow by -5 in transverse directions
      const int numGhostW =
        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry);
      interiorBox.grow(a_dir, numGhostW);
      interiorBox.grow(-numGhostW);
      CH_assert(!interiorBox.isEmpty());

      // Get coordinate system for the block
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(a_disjointBox));
      
      // Get physical coordinates
      FABSTACKTEMP(XiFab, a_boundaryFaceBox, SpaceDim);  // Cart. coordinates
      FABSTACKTEMP(XFab, a_boundaryFaceBox, SpaceDim);   // Phys. coordinates
      this->CNSIBC::getFaceCoordinates(a_boundaryFaceBox,
                                       XiFab,
                                       XFab,
                                       a_dir,
                                       blockCoordSys);

      // Set the state on the boundary
      int viscousSlip = 0;  // Only required for lo side so value doesn't matter
      const Real xsloc = m_xBase + m_Mach * m_c0 * a_time;

      FORT_CNSIBCMACHREFLECTIONBCY(CHF_FRA(a_Wface),
                                   CHF_BOX(a_boundaryFaceBox),
                                   CHF_CONST_FRA(a_Wcell),
                                   CHF_CONST_FRA(XFab),
                                   CHF_CONST_INT(viscousSlip),
                                   CHF_CONST_INT(lohiSign),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_p0),
                                   CHF_CONST_REAL(m_r0),
                                   CHF_CONST_REAL(m_p1),
                                   CHF_CONST_REAL(m_r1),
                                   CHF_CONST_REAL(m_u1),
                                   CHF_CONST_REAL(xsloc));
      // FIXME: This should change based on where we are in the domain
      for (int i = 0; i != numSpecies; ++i)
        {
          const int scomp = speciesComp + i;
          a_Wface.setVal(m_upstreamMF[i], a_boundaryFaceBox, scomp);
        }
    }
  else if (a_dir == 0 && a_side == Side::Lo)
    {
      a_Wface.setVal(m_r1, a_boundaryFaceBox, rComp);
      a_Wface.setVal(m_p1, a_boundaryFaceBox, pComp);
      a_Wface.setVal(0., a_boundaryFaceBox, velComp, SpaceDim);
      a_Wface.setVal(m_u1, a_boundaryFaceBox, velComp);
      for (int i = 0; i != numSpecies; ++i)
        {
          const int scomp = speciesComp + i;
          a_Wface.setVal(m_upstreamMF[i], a_boundaryFaceBox, scomp);
        }
    }
  else if (a_dir == 0 && a_side == Side::Hi)
    {
      a_Wface.setVal(m_r0, a_boundaryFaceBox, rComp);
      a_Wface.setVal(m_p0, a_boundaryFaceBox, pComp);
      a_Wface.setVal(0., a_boundaryFaceBox, velComp, SpaceDim);
      for (int i = 0; i != numSpecies; ++i)
        {
          const int scomp = speciesComp + i;
          a_Wface.setVal(m_downstreamMF[i], a_boundaryFaceBox, scomp);
        }
    }
  CRDparam::g_CRDPhysics->temperature(a_Wface, a_boundaryFaceBox);
#endif
}

//**FIXME Default implementation should be able to get wall normal from metrics
/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In this override, the velocity normal to the wall is determined
 *  from the known geoemtry
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state corrected for presence of wall
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
CNSIBCSpecMachReflection::setWallBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_boundarySlipVelocity,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
#if (CH_SPACEDIM != 1)
  CH_assert(a_dir == 1 && a_side == Side::Lo);

  const int lohiSign = sign(a_side);

  // Get a disjoint interior box
  Box interiorBox = adjCellBox(a_boundaryFaceBox,
                               a_dir,
                               Side::flip(a_side),
                               1);
  //**FIXME - this should be single-block so any box should work.  Ideally
  //**        the parent block domain should be an annotation in the box
  //**FIXME - this should work as long as min box size is >= 8
  // Grow by -5 in transverse directions
  const int numGhostW =
    CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry);
  interiorBox.grow(a_dir, numGhostW);
  interiorBox.grow(-numGhostW);
  CH_assert(!interiorBox.isEmpty());

  // Get coordinate system for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
      
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_boundaryFaceBox, SpaceDim);  // Cart. coordinates
  FABSTACKTEMP(XFab, a_boundaryFaceBox, SpaceDim);   // Phys. coordinates
  this->CNSIBC::getFaceCoordinates(a_boundaryFaceBox,
                                   XiFab,
                                   XFab,
                                   a_dir,
                                   blockCoordSys);
  
  const Real xsloc = 0.;  // Only required for lo side value doesn't matter
  FORT_CNSIBCCOMBMACHREFLECTIONBCY(CHF_FRA(a_Wface),
                                   CHF_BOX(a_boundaryFaceBox),
                                   CHF_CONST_FRA(a_Wcell),
                                   CHF_CONST_FRA(XFab),
                                   CHF_CONST_INT(lohiSign),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_p0),
                                   CHF_CONST_REAL(m_r0),
                                   CHF_CONST_REAL(m_p1),
                                   CHF_CONST_REAL(m_r1),
                                   CHF_CONST_REAL(m_u1),
                                   CHF_CONST_REAL(xsloc),
                                   CHF_CONST_REAL(m_xloclow));
  CRDparam::g_CRDPhysics->temperature(a_Wface, a_boundaryFaceBox);
#endif
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
CNSIBCSpecMachReflection::readBCInfo()
{
  CH_assert(!m_readInput);
  const int numSpecies = CRDparam::g_numSpecies;
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error
  ParmParse ppIBC("ibc");
//--Ramp values

  ppIBC.get("ramp_angle", m_alpha);
  ppIBC.get("lead_length", m_xLead);
  ppIBC.get("ramp_length", m_xRamp);
  // Convert to radians
  m_alpha *= PI/180.;
  ppIBC.query("base_loc", m_xBase);

//--Determine the boundary conditions
  ppIBC.query("lower_x_bound", m_lowXBound);
  ppIBC.query("upper_y_bound", m_hiYBound);
  if (m_lowXBound < 0 || m_lowXBound > 1 || m_hiYBound < 0 || m_hiYBound > 2)
    {
      CRD::msg << "Input (SpecMachReflection IBC): lower and upper bound "
               << " are not specified correctly!" << CRD::error;
    }

//--Additional tests on previous input

  if (CRDparam::g_domainLength[0] != 1.0)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'domain_length' "
               << "must be set to 1.0 in x-direction for Mach reflection"
               << " problem!" << CRD::error;
    }
  {
    const Real dx0 = CRDparam::g_domainLength[0]/CRDparam::g_domainBaseSize[0];
    for (int dir = 1; dir != SpaceDim; ++dir)
    {
      const Real dx =
        CRDparam::g_domainLength[dir]/CRDparam::g_domainBaseSize[dir];
      if (Misc::compare(dx0, dx, std::numeric_limits<Real>::digits10 - 2))
        {
          CRD::msg << "Input (SpecMachReflection IBC): 'domain_length' in "
            "direction " << dir << " must be set for same mesh spacing in "
            "all directions.\ndx[0] = " << dx0 << "\ndx[" << dir << "] = "
                   << dx << '!' << CRD::error;
        }
    }
  }
  if (CRDparam::g_domainOrigin != RealVect::Zero)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'domain_origin' must be "
               << "set to zero for Mach reflection problem "
               << "(do not define for valid default)!"
               << CRD::error;
    }

//--Quiescent conditions in front of the shock

  m_p0 = 1.0;
  m_r0 = CRDparam::g_gamma;
  ppIBC.query("pressure", m_p0);
  if (m_p0 <= 0.)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'pressure' must be > 0!"
               << CRD::error;
    }
  ppIBC.query("density", m_r0);
  if (m_r0 <= 0.)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'density' must be > 0!"
               << CRD::error;
    }
  Real gamma = CRDparam::g_gamma;
  m_c0 = std::sqrt(gamma*m_p0/m_r0);

//--Conditions behind of the shock

  m_p1 = -1.;
  m_r1 = -1.;
  m_u1 = -1.;
  ppIBC.get("pressure_1", m_p1);
  ppIBC.get("density_1", m_r1);
  ppIBC.get("u_1", m_u1);

//--Shock Mach number

  m_Mach = 10.0;
  ppIBC.query("shock_mach", m_Mach);
  if (m_Mach <= 1.)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'shock_mach' must be > 1!"
               << CRD::error;
    }

//--Threshold of relative density gradient for refinement

  m_threshold = 0.2;
  ppIBC.query("tag_threshold", m_threshold);
  if (m_threshold < 0.)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'tag_threshold' must be > 0!"
               << CRD::error;
    }

//--Distance between shock and start of ramp
  
  ppIBC.query("base_size", m_xBase);

//--Distance between the start of ramp and the slip wall portion end
  m_xloclow = m_xBase;
  ppIBC.query("slip_wall_end", m_xloclow);

//--Order of extrapolation to use at the wall
  
  ppIBC.query("wall_order", m_wallOrder);
  if (m_wallOrder < 0 || m_wallOrder > 4)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'wall_order' must be > 0 "
               << "and <= 4!" << CRD::error;
    }

//--Set the component to do the tag for AMR

  m_tagComp = 0;
  ppIBC.query("tag_comp", m_tagComp);
  if (m_tagComp < 0 || m_tagComp > WNUM)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'tag_comp' must be > 0 and "
               << "< the number of primitive variables" << CRD::error;
    }

//--Determine if to do base tagging near wall
  ppIBC.query("tag_wall", m_tagWall);

  // Mass fractions upstream of shock
  m_upstreamMF.resize(numSpecies);
  m_upstreamMF.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  int massCheck = assignMassFractions(m_upstreamMF,
                                      "upstream_specs",
                                      "upstream_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'upstream_mass_fractions' "
               << "must be equal to 1!" << CRD::error;
    }
  // Mass fractions downstream of shock
  m_downstreamMF.resize(numSpecies);
  m_downstreamMF.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  massCheck = assignMassFractions(m_downstreamMF,
                                  "downstream_specs",
                                  "downstream_mfs");
  if (massCheck == 1)
    {
      CRD::msg << "Input (SpecMachReflection IBC): 'downstream_mass_fractions' "
               << "must be equal to 1!" << CRD::error;
    }
  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
