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
 * \file CNSIBCMachReflection.cpp
 *
 * \brief Member functions for CNSIBCMachReflection
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

#include "CNSIBCMachReflection.H"
#include "CNSIBCMachReflectionF_F.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCMachReflection: member definitions
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

CNSIBCMachReflection::CNSIBCMachReflection()
:
  CNSIBC(),
  m_alpha(PI*30./180.),
  m_xBase(-0.05),
  // All invalid values which must be corrected
  m_Mach(-1.),
  m_wallOrder(4),
  m_idxStatePreShock(-1),
  m_idxStatePostShock(-1)
{
  if (SpaceDim < 2)
    {
      CRD::msg << "Problem not defined for 1D!" << CRD::error;
    }

//--Read any BC info

  readBCInfo();

//--Wall state (defined as all zero if not yet set)

  {
    CRDState& wallState = CRDState::get("wall");
    (void)wallState;
  }

//--State in front of shock

  if (CRDState::nameIndex("pre-shock") == -1) // Define with default values
    {
      CRDState& state0 = CRDState::get("pre-shock");
      state0.density()  = CRDparam::g_gamma;
      state0.pressure() = 1.0;
      state0.velocity() = RealVect_zero;
      state0.setExtraThermo();
    }
  m_idxStatePreShock = CRDState::nameIndex("pre-shock");
  Real r0;  // Saved since addition of new states makes refs to state0 invalid
  Real p0;
  Real gamma;
  {
    const CRDState& state0 = CRDState::get(m_idxStatePreShock);
    r0 = state0.density();
    p0 = state0.pressure();

    // Compute additional state info in front of the shock
    const Real* cn = nullptr;
    if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
      // Note however that this code assumes a frozen gas for computing the
      // shock
      {
        cn = &state0(CRDparam::g_CRDPhysics->speciesPrimInterval().begin());
      }
    gamma = CRDparam::g_CRDPhysics->gamma(state0.temperature(), cn);
    m_c0 = std::sqrt(gamma*p0/r0);
  }

//--Compute state behind the shock

  // If state behind the shock is not set from input file, solve for the values
  // using the following equations
  if (CRDState::nameIndex("post-shock") == -1) // Define with default values
    {
      CRDState& state1 = CRDState::get("post-shock");
      const Real p1 = ((std::pow(m_Mach, 2) - 1.0)*(2*gamma)/
                       (gamma + 1.0) + 1.0)*p0;
      const Real prat = p1/p0;
      const Real grat = (gamma + 1.0)/(gamma - 1.0);
      const Real r1 = (1.0 + grat*prat)/(grat + prat)*r0;
      const Real u1 = (m_c0/gamma)*(prat - 1.0)
        *std::sqrt((2*gamma/(gamma + 1.0))/(prat + 1.0/grat));
      state1.density()  = r1;
      state1.pressure() = p1;
      state1.velocity() = u1*RealVect_basis(0);
      state1.setExtraThermo();
    }
  m_idxStatePostShock = CRDState::nameIndex("post-shock");
  const CRDState& state1 = CRDState::get(m_idxStatePostShock);

//--Set common reference state

  CNSIBC::defineCommonRef(r0, state1.velocity()[0]);

//--Set BC Types

  // Assumes single block
  if (CRDparam::g_coordSys->numBlocks() != 1)
    {
      CRD::msg << "Mach reflection problem requires single block!"
               << CRD::error;
    }

  BCInfo bc;

  // Inflow on left side (supersonic so Dirichlet)
  bc.m_type = CRDparam::DomainBCTypeDirichlet;
  bc.m_idxState = m_idxStatePostShock;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 0, Side::Lo), bc);
  // Quiescent on right side (intersection not supported so Dirichlet)
  bc.m_type = CRDparam::DomainBCTypeDirichlet;
  bc.m_idxState = m_idxStatePreShock;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 0, Side::Hi), bc);
  // Wall on bottom
  bc.m_type = CRDparam::DomainBCTypeWall;
  bc.m_idxState = CRDState::nameIndex("wall");
  bc.m_order = m_wallOrder;
  setDomainBC(BoundaryIndex(0, 1, Side::Lo), bc);
  // Farfield on top (but we will customize the implementation)
  bc.m_type = CRDparam::DomainBCTypeFarfield;
  bc.m_idxState = m_idxStatePostShock;
  bc.m_order = 1;
  setDomainBC(BoundaryIndex(0, 1, Side::Hi), bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCMachReflection::~CNSIBCMachReflection()
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
CNSIBCMachReflection::IBCName() const
{
  return "Mach reflection";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCMachReflection::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  CRD::msg << "x-location of corner\n" << 0.0 << CRD::var;
  CRD::msg << "Starting x-location of shock\n" << m_xBase << CRD::var;
  CRD::msg << "Shock Mach number\n" << m_Mach << CRD::var;
  CRD::msg << "Order of accuracy at wall BC\n" << m_wallOrder << CRD::var;
  CRD::msg.newline();
}

#if 0
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
CNSIBCMachReflection::setTagMethod(const int a_tagBufferSize)
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
#endif

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
CNSIBCMachReflection::initialize(LevelData<FArrayBox>&      a_U,
                                 LevelGridMetrics&          a_gridMetrics,
                                 const LayoutData<FluxBox>& a_unitNormals,
                                 const Real                 a_time,
                                 const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
  const Real gamma      = CRDparam::g_gamma;

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

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CRDState& state0 = CRDState::get(m_idxStatePreShock);
      const Real r0 = state0.density();
      const Real p0 = state0.pressure();
      CRDState& state1 = CRDState::get(m_idxStatePostShock);
      const Real r1 = state1.density();
      const Real p1 = state1.pressure();
      const Real u1 = state1.velocity()[0];
      CH_assert(UFab.box().contains(box2Dom));
      FORT_CNSIBCMACHREFLECTIONINIT(CHF_FRA(UFab),
                                    CHF_BOX(box2Dom),
                                    CHF_CONST_FRA(XFab),
                                    CHF_CONST_REAL(gamma),
                                    CHF_CONST_REAL(m_xBase),
                                    CHF_CONST_REAL(p0),
                                    CHF_CONST_REAL(r0),
                                    CHF_CONST_REAL(p1),
                                    CHF_CONST_REAL(r1),
                                    CHF_CONST_REAL(u1));

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
 *  \param[in]  a_domT  Unused in this class
 *//*-----------------------------------------------------------------*/

void
CNSIBCMachReflection::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_domT) const
{
#if (CH_SPACEDIM != 1)
  const int lohiSign = sign(a_bcIdx.m_side);

  if (a_bcIdx.m_dir == 1 && a_bcIdx.m_side == Side::Hi)
    {
      // Get a disjoint interior box
      Box interiorBox = adjCellBox(a_boundaryFaceBox,
                                   a_bcIdx.m_dir,
                                   Side::flip(a_bcIdx.m_side),
                                   1);
      //**FIXME - this should be single-block so any box should work.  Ideally
      //**        the parent block domain should be an annotation in the box
      //**FIXME - this should work as long as min box size is >= 8
      // Grow by -5 in transverse directions
      const int numGhostW =
        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry);
      interiorBox.grow(a_bcIdx.m_dir, numGhostW);
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
                                       a_bcIdx.m_dir,
                                       blockCoordSys);

      // Set the state on the boundary
      int viscousSlip = 0;  // Only required for lo side so value doesn't matter
      const Real xsloc = m_xBase + m_Mach * m_c0 * a_time;

      CRDState& state0 = CRDState::get(m_idxStatePreShock);
      const Real r0 = state0.density();
      const Real p0 = state0.pressure();
      CRDState& state1 = CRDState::get(m_idxStatePostShock);
      const Real r1 = state1.density();
      const Real p1 = state1.pressure();
      const Real u1 = state1.velocity()[0];
      FORT_CNSIBCMACHREFLECTIONBCY(CHF_FRA(a_Wface),
                                   CHF_BOX(a_boundaryFaceBox),
                                   CHF_CONST_FRA(a_Wcell),
                                   CHF_CONST_FRA(XFab),
                                   CHF_CONST_INT(viscousSlip),
                                   CHF_CONST_INT(lohiSign),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(p0),
                                   CHF_CONST_REAL(r0),
                                   CHF_CONST_REAL(p1),
                                   CHF_CONST_REAL(r1),
                                   CHF_CONST_REAL(u1),
                                   CHF_CONST_REAL(xsloc));
    }
  else // Call base implementation that sets constant values
    {
      CNSIBC::setImposedBCprimState(a_Wface,
                                    a_boundaryFaceBox,
                                    a_Wcell,
                                    a_unitNormalBasisFab,
                                    a_bcIdx,
                                    a_disjointBox,
                                    a_gridMetrics,
                                    a_time,
                                    a_level,
                                    a_domT);
    }
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
CNSIBCMachReflection::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error

//--Geometry

  {
    ParmParse ppCOORDSYS("coordsys");
    ppCOORDSYS.get("ramp_angle", m_alpha);
    // Convert to radians
    m_alpha *= PI/180.;
  }

//--IBC

  ParmParse ppIBC("ibc");
  ppIBC.query("base_loc", m_xBase);

//--Additional tests on previous input

  {
    const Real dx0 = CRDparam::g_domainLength[0]/CRDparam::g_domainBaseSize[0];
    for (int dir = 1; dir != SpaceDim; ++dir)
      {
        const Real dx =
          CRDparam::g_domainLength[dir]/CRDparam::g_domainBaseSize[dir];
        if (Misc::compare(dx0, dx, std::numeric_limits<Real>::digits10 - 2))
          {
            CRD::msg << "Input (MachReflection IBC): 'domain_length' in "
              "direction " << dir << " must be set for same mesh spacing in "
              "all directions.\ndx[0] = " << dx0 << "\ndx[" << dir << "] = "
                     << dx << '!' << CRD::error;
          }
      }
  }
  if (CRDparam::g_domainOrigin != RealVect::Zero)
    {
      CRD::msg << "Input (MachReflection IBC): 'domain_origin' must be set to "
        "zero for Mach reflection problem (do not define for valid default)!"
               << CRD::error;
    }

//--Shock Mach number

  m_Mach = 10.0;
  ppIBC.query("shock_mach", m_Mach);
  if (m_Mach <= 1.)
    {
      CRD::msg << "Input (MachReflection IBC): 'shock_mach' must be > 1!"
               << CRD::error;
    }

//--Distance between shock and start of ramp
  
  ppIBC.query("base_size", m_xBase);

//--Order of extrapolation to use at the wall
  
  ppIBC.query("wall_order", m_wallOrder);
  if (m_wallOrder < 0 || m_wallOrder > 4)
    {
      CRD::msg << "Input (MachReflection IBC): 'wall_order' must be > 0 and "
               << "<= 4!" << CRD::error;
    }
}
