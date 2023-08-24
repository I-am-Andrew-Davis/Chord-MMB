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
 * \file CNSIBC.cpp
 *
 * \brief Member functions for CNSIBC
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
#include "PolytropicPhysicsF_F.H"
#include "LGintegrator.H"
#include "GodunovUtilitiesF_F.H"

//----- Internal -----//

#include "CNSIBC.H"
#include "CNSIBCF_F.H"
#include "CRDparam.H"
#include "CNSPhysics.H"
#include "CRDState.H"
#include "DataTemp.H"
#include "PatchMappedFunc.H"
#include "ViscousTensor4thOrderOp.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "BdryCharacteristics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "SetCentersF_F.H"

/*--------------------------------------------------------------------*
 * Definition of statics
 *--------------------------------------------------------------------*/

bool CNSIBC::s_firstRKStage = false;
bool CNSIBC::s_lastRKStage = false;

/*--------------------------------------------------------------------*/
//  Output of a boundary index
/** \param[in]  a_os    Output stream
 *  \param[in]  a_idxB  Boundary index to output
 *//*-----------------------------------------------------------------*/

std::ostream& operator<<(std::ostream&       a_os,
                         const BoundaryIndex a_idxB)
{
  static constexpr const char* dirName[] =
    { "x", "y", "z" };
  static constexpr const char* sideName[] =
    { "min", "max" };
  a_os << "block[" << a_idxB.m_block
       << "] " << dirName[a_idxB.m_dir]
       << '-' << sideName[a_idxB.m_side];
  return a_os;
}

/*******************************************************************************
 *
 * Class CNSIBC: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor reads state, sets reference, and default BC
/** The reference value are used to construct small values.
 *//*-----------------------------------------------------------------*/

CNSIBC::CNSIBC()
:
  m_readInput(false),
  m_idxStateRef(-1)
{

//--Read the state information

  CRDState::readStateInfo();
  if (CRDState::nameIndex("reference") == -1) // Define with normalized default
                                              // values
    {
      CRDState& stateRef = CRDState::get("reference");
      stateRef.density()  = CRDparam::g_gamma;
      stateRef.pressure() = 1.0;
      stateRef.velocity() = RealVect_basis(0);
      stateRef.setExtraThermo();
    }
  m_idxStateRef = CRDState::nameIndex("reference");
  {
    const CRDState& stateRef = CRDState::get(m_idxStateRef);
    const Real speedRef = stc::mag(stateRef.velocity());
    if (speedRef <= 0.)
      {
        CRD::msg << "CNSIBC::CNSIBC: magnitude of reference velocity must be "
          "greater than zero.  Value given was " << speedRef << CRD::error;
      }
    defineCommonRef(stateRef.density(), speedRef);
  }

  // Default for all domain BC is 4th order periodic
  BCInfo defBC;
  defBC.m_order = 4;
  defBC.m_type = CRDparam::DomainBCTypePeriodic;
  setAllDomainBC(defBC);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBC::~CNSIBC()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Allocations of new physics states.  Customize for each derivative
//  class.
/**
 *//*-----------------------------------------------------------------*/

std::vector<CRDPhysics*>
CNSIBC::allocatePhysics()
{
  std::vector<CRDPhysics*> physics(1, nullptr);
  physics[0] = new CNSPhysics;
  return physics;
}

/*--------------------------------------------------------------------*/
//  Define common reference state
/** WARNING: This should only be called in a constructor as the
 *  values set are assumed to be constant for a run
 *//*-----------------------------------------------------------------*/

void
CNSIBC::defineCommonRef(const Real a_rhoRef, const Real a_uRef)
{
  Real smallr, smallp;
  FORT_CNSIBCSETCOMMONREF(CHF_REAL(smallr),
                          CHF_REAL(smallp),
                          CHF_CONST_REAL(CRDparam::g_gamma),
                          CHF_CONST_REAL(a_rhoRef),
                          CHF_CONST_REAL(a_uRef));
  CRDparam::CRDP.defineSmall(smallr, smallp);
}

/*--------------------------------------------------------------------*/
//  Set domainBC on all boundaries of all blocks
/** \param[in]  a_domainBC
 *                      Domain BC type
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setAllDomainBC(const BCInfo& a_domainBC)
{
  BoundaryIndex bdIdx;
  for (int block = 0; block != CRDparam::g_coordSys->numBlocks(); ++block)
    {
      bdIdx.m_block = block;
      for (const int dir : EachDir)
        {
          bdIdx.m_dir = dir;
          for (const auto side : EachSide)
            {
              bdIdx.m_side = side;
              setDomainBC(bdIdx, a_domainBC);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set domainBC info on one boundary of
/** \param[in]  a_bcIndex 
 *                      The multi-block boundary
 *  \param[in]  a_domainBC
 *                      Domain BC info on this side of block domain
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setDomainBC(const BoundaryIndex& a_bcIndex,
                    const BCInfo&        a_domainBC)
{
  a_bcIndex.assertValid();
  auto boundary = CRDparam::g_coordSys->boundary(a_bcIndex.m_block,
                                                 a_bcIndex.m_dir,
                                                 a_bcIndex.m_side);
  if (boundary.isDomainBoundary())
    {
      m_blockBCInfo[a_bcIndex] = a_domainBC;
    }
  else if (CRDparam::g_verbosity > 2)
    {
      pout() << "SetDomainBC not valid for " << a_bcIndex
             << " because it is a multi-block interface.\n";
    }
}

/*--------------------------------------------------------------------*/
//  Set domain state on one boundary
/** \param[in]  a_stateIndex 
 *                      The index for the boundary state (can be anything)
 *                      where each state is indexed by a unique number
 *  \param[in]  a_stateVal
 *                      Boundary values for this state
 *//*-----------------------------------------------------------------*/

// void
// CNSIBC::setDomainBCstate(const int a_stateIndex,
//                          const BCstate& a_stateVal)
// {
//   m_blockBCstate[a_stateIndex] = a_stateVal;
// }

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBC::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  //**FIXME working to avoid use of g_rho, etc, in preference of reference state
  // CRD::msg.setPrecFloatSN(5);
  // CRD::msg << "Uniform density\n" <<  CRDparam::g_rho << CRD::var;
  // CRD::msg << "Uniform temperature\n" << CRDparam::g_T << CRD::var;
  // CRD::msg << "Uniform pressure\n"
  //          << CRDparam::g_rho*CRDparam::g_R*CRDparam::g_T << CRD::var;
  // CRD::msg << "Velocity\n(";
  // for (int dir = 0; dir != SpaceDim; ++dir)
  //   {
  //     if (dir != 0)
  //       {
  //         CRD::msg << ',';
  //       }
  //     CRD::msg << 0.;
  //   }
  // CRD::msg << ')' << CRD::var;
  // CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Write each specified boundary state and boundary condition
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBC::writeBoundaryConditions() const
{
  CRD::msg.setPrecFloatSN(5);
  //CRD::msg << "Defined boundary states:" << CRD::body;
  // for (const auto &bcStatePair : m_blockBCstate) {
  //   CRD::msg << "State Number\n" << bcStatePair.first << CRD::var;
  //   CRD::msg << "  Temperature\n" << bcStatePair.second.m_temp << CRD::var;
  //   CRD::msg << "  Pressure\n" << bcStatePair.second.m_pres << CRD::var;
  //   CRD::msg << "  Velocity\n(";
  //   for (const auto dir : EachDir)
  //     {
  //       if (dir != 0)
  //         {
  //           CRD::msg << ',';
  //         }
  //       CRD::msg << bcStatePair.second.m_vel[dir];
  //     }
  //   CRD::msg << ')' << CRD::var;
  //   bool first = true;
  //   CRD::msg << "  Mass Fraction\n(";
  //   for (const auto &mf : bcStatePair.second.m_MF)
  //     {
  //       if (!first)
  //         {
  //           CRD::msg << ',';
  //         }
  //       first = false;
  //       CRD::msg << mf;
  //     }
  //   CRD::msg << ')' << CRD::var;
  // }
  CRD::msg << "Defined boundary types:" << CRD::body;
  for (int idxBlk = 0, idxBlk_end = CRDparam::g_coordSys->numBlocks();
       idxBlk != idxBlk_end; ++idxBlk)
    {
      for (const auto &bcInfoPair : m_blockBCInfo)
        {
          if (bcInfoPair.first.m_block == idxBlk)
            {
              CRD::msg << "Boundary index\n" << bcInfoPair.first << CRD::var;
              CRD::msg << "  Type\n"
                       << CRDparam::domainBCname(bcInfoPair.second.m_type)
                       << CRD::var;
              CRD::msg << "  Order\n" << bcInfoPair.second.m_order << CRD::var;
              // CRD::msg << "  State\n" << bcInfoPair.second.m_stateIdx
              //          << CRD::var;
              CRD::msg << "  State\n" << bcInfoPair.second.m_idxState
                       << CRD::var;
              CRD::msg << "  CBC state relaxation parameter\n"
                       << bcInfoPair.second.m_relaxCBCStateParam << CRD::var;
              CRD::msg << "  CBC wave relaxation parameter\n"
                       << bcInfoPair.second.m_relaxCBCWaveParam << CRD::var;
            }
        }
    }
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
 *    <li> Allocate new methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory*
CNSIBC::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  // Return the factory
  return new TagLevelFactory(tagLevel);
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\>.  This default implementation initializes to a state named
 *  initial.  If it is not found, the reference state is used.
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
CNSIBC::initialize(LevelData<FArrayBox>&      a_U,
                   LevelGridMetrics&          a_gridMetrics,
                   const LayoutData<FluxBox>& a_unitNormals,
                   const Real                 a_time,
                   const int                  a_level) const
{
  CH_TIME("CNSIBC::initialize");
  int idxStateInit = CRDState::nameIndex("initial");
  if (idxStateInit == -1)
    {
      CH_assert(m_idxStateRef != -1);
      idxStateInit = m_idxStateRef;
    }
  const CRDState& stateInit = CRDState::get(idxStateInit);

  const Real density = stateInit.density();
  const RealVect momentum = density*stateInit.velocity();
  const Real energy = stateInit.pressure()/(CRDparam::g_gamma - 1.) +
    0.5*density*stc::dot(stateInit.velocity(), stateInit.velocity());

  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const Interval cMomIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int cEng = CRDparam::g_CRDPhysics->energyFluxIndex();

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {

// Do not delete commented out items.  This is a a guide for derived classes.

      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block (using multiblock coordinate system)
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Data
      FArrayBox& UFab = a_U[dit];

      // Working set boxes
      // Box box2Dom = grow(box, 2);
      // box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // // Get physical coordinates
      // FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      // FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      // this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      // CH_assert(UFab.box().contains(box2Dom));

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box1Dom));
      // fourthOrderAverageCell(UFab, blockDomain, box1Dom);

      // Here, the average values are set directly
      UFab.setVal(density, box1Dom, cRho);
      for (const int dir : EachDir)
        {
          UFab.setVal(momentum[dir], box1Dom, cMomIntv.begin() + dir);
        }
      UFab.setVal(energy, box1Dom, cEng);
      for (int c = cEng+1, c_end = CRDparam::g_CRDPhysics->numConservative();
           c < c_end; ++c)
        {
          UFab.setVal(0., box1Dom, c);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Initialize wall-model state after restart
/** After restarting without the wall-model, we need a way to initialize
 *  the wall-model eta_0 value to something sensible
 *  \param[out] a_JU    State on the level to be initialized in this
 *                      routine
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *//*-----------------------------------------------------------------*/

void
CNSIBC::restartWallModel(LevelData<FArrayBox>&      a_JU,
                         LevelGridMetrics&          a_gridMetrics,
                         const LayoutData<FluxBox>& a_unitNormals) const
{
  CH_TIME("CNSIBC::restartWallModel");

  // Iterate over all boxes on this level
  for (DataIterator dit = a_JU.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_JU.disjointBoxLayout()[dit];

      // Data
      FArrayBox& JUFab = a_JU[dit];

      // Call restartWallModel
      CRDparam::g_CRDPhysics->initWallModelAfterRestart(
        JUFab, a_unitNormals[dit], a_gridMetrics, dit(), box);
    }
}

/*--------------------------------------------------------------------*/
//  Initialize the inlet region data structures
/** \param[in] a_U     Conservative state
 *  \param[in] a_level Current level
 *  \param[in]  a_disjointBoxLayout
 *                      Disjoint box layout
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *//*-----------------------------------------------------------------*/

void
CNSIBC::initializeInletDataStructures(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics,
  const bool                  a_hasFinerGrid) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Initialize the inlet region data structures
/** \param[in]  a_U     Conservative state
 *  \param[in]  a_level Current level
 *  \param[in]  a_t     Solution time
 *  \param[in]  a_dt    Current time-step size
 *  \param[in]  a_stage Current stage of the time-marching method
 *  \param[in]  a_disjointBoxLayout
 *                      Disjoint box layout
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *//*-----------------------------------------------------------------*/

void
CNSIBC::copyInteriorToInletAndRescale(
  const LevelData<FArrayBox>& a_U,
  const int                   a_level,
  const Real                  a_t,
  const Real                  a_dt,
  const int                   a_stage,
  const DisjointBoxLayout&    a_disjointBoxLayout,
  const LevelGridMetrics&     a_gridMetrics) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Set boundary slopes
/** The boundary slopes in a_dW are already set to one sided
 *  difference approximations.  If this function doesn't change them
 *  they will be used for the slopes at the boundaries.
 *  \param[out] a_dW
 *  \param[in]  a_W
 *  \param[in]  a_dir
 *  \param[in]  a_time
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setBdrySlopes(FArrayBox&       a_dW,
                      const FArrayBox& a_W,
                      const int&       a_dir,
                      const Real&      a_time)
{
  CRD::msg << "CNSIBC::setBdrySlopes is not supported." << CRD::error;
}

/*--------------------------------------------------------------------*/
//  Adjust boundary fluxes to account for artificial viscosity
/** Use at walls with the Euler equations
 *  \param[in]  a_NtFdir
 *                      All components at indices in a_loFaceBox and
 *                      a_hiFaceBox must be set to zero
 *  \param[out] a_NtFdir
 *                      Fluxes due to artificial viscosity added to
 *                      'a_dir' faces along the domain boundary
 *  \param[in]  a_Nctg  N, with components stored contiguously, on the
 *                      'a_dir' faces
 *  \param[in]  a_U     Conservative state \<U\> in the cells
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_divVel
 *                      Physical space divergence of the velocity on
 *                      the 'a_dir' faces
 *  \param[in]  a_csq   Min of c^2/gamma (over +1 layer of cells in
 *                      all transverse directions), stored on cells
 *  \param[in]  a_dxFace
 *                      Distance between cell centers in direction
 *                      'a_dir' stored on faces
 *  \param[in]  a_momIntv
 *                      Interval of momentum components
 *  \param[in]  a_alpha Artificial viscosity coefficient
 *  \param[in]  a_beta  Artificial viscosity coefficient
 *  \param[in]  a_loFaceBox
 *                      Domain faces on the low side
 *  \param[in]  a_hasLo Low-side box exists
 *  \param[in]  a_hiFaceBox
 *                      Domain faces on the high side
 *  \param[in]  a_hasHi High-side box exists
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void CNSIBC::artViscBC(
  FArrayBox&                                   a_NtFdir,
  const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
  const FArrayBox&                             a_U,
  const FArrayBox&                             a_unitNormalBasis,
  const FArrayBox&                             a_divVel,
  const FArrayBox&                             a_csq,
  const FArrayBox&                             a_dxFace,
  const Interval&                              a_momIntv,
  const Real                                   a_alpha,
  const Real                                   a_beta,
  const Box&                                   a_loFaceBox,
  const int                                    a_hasLo,
  const Box&                                   a_hiFaceBox,
  const int                                    a_hasHi,
  const int                                    a_dir,
  const Box&                                   a_disjointBox,
  LevelGridMetrics&                            a_gridMetrics,
  const Real                                   a_time,
  const int                                    a_level) const
{
  CH_TIME("CNSIBC::artViscBC");
  // Const cast just for shifting
  FArrayBox &csq = const_cast<FArrayBox&>(a_csq);
  FArrayBox &U = const_cast<FArrayBox&>(a_U);
  Vector<BCInfo> bcTypes;
  Vector<Box> bcBoxes;
  // The boundary index
  BoundaryIndex bcIdx;
  bcIdx.m_block = a_gridMetrics.getCoordSys().whichBlock(a_disjointBox);
  bcIdx.m_dir = a_dir;
  // Size of the stencil used adjacent to the wall
  const int stenSize = 1;

  // Create lambda for call to setting either side
  auto operateOnSide =
    [&] (const BoundaryIndex& a_bcIdx, const Box& a_faceBox)
    { 
      const BCInfo& domBC = getDomainBC(a_bcIdx);
      const CRDState& domBCstate = CRDState::get(domBC.m_idxState);
      const int lohiSign = sign(a_bcIdx.m_side);
      // All cell centered quantities need to be shifted to the faces
      csq.shiftHalf(a_dir, lohiSign);
      U.shiftHalf(a_dir, lohiSign);
      if (CRDparam::DomainBCTypeMixed & domBC.m_type)
        {
          int setBC = setMixedBC(a_faceBox,
                                 bcIdx,
                                 a_disjointBox,
                                 a_gridMetrics,
                                 a_time,
                                 a_level,
                                 bcBoxes,
                                 bcTypes);
          CH_assert(setBC == 0);
        }
      else
        {
          bcBoxes.resize(1);
          bcTypes.resize(1);
          bcBoxes[0].define(a_faceBox);
          bcTypes[0] = domBC;
        }
      CH_assert(a_NtFdir.contains(a_faceBox));
      CH_assert(U.contains(a_faceBox));
      CH_assert(a_unitNormalBasis.contains(a_faceBox));
      CH_assert(a_divVel.contains(a_faceBox));
      CH_assert(csq.contains(a_faceBox));
      CH_assert(a_dxFace.contains(a_faceBox));

      // Alias to the momentum components
      FArrayBox NtFdirMom(a_momIntv, a_NtFdir);  // Alias
      // Make box to include more than just boundary adjacent values
      Box nearBoundaryBox(a_faceBox);
      nearBoundaryBox.growDir(a_dir, flip(a_bcIdx.m_side), stenSize);
      // For U, copy the momentum components so it can be altered
      FArrayBox UMom(nearBoundaryBox, SpaceDim);
      UMom.copy(a_U, a_momIntv.begin(), 0, SpaceDim);
      for (int BCNum = 0; BCNum != bcTypes.size(); ++BCNum)
        {
          Box modBox(bcBoxes[BCNum]);
          // Do not do the operations if the box is empty
          if (modBox.isEmpty()) continue;
          auto domT = bcTypes[BCNum];
          if (domT.m_type & CRDparam::DomainBCTypeAllWall)
            {
              Box stenBox(modBox);
              stenBox.growDir(a_dir, flip(a_bcIdx.m_side), stenSize);
              RealVect wallVel(domBCstate.velocity());
              // Modified from MOLPhysicsMappedArtViscF.ChF
              FORT_MAPPEDARTVISCINVISCIDHO(
                CHF_FRA(NtFdirMom),
                CHF_FRA(UMom),
                CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, a_Nctg),
                CHF_CONST_FRA(a_unitNormalBasis),
                CHF_CONST_FRA1(a_divVel, 0),
                CHF_CONST_FRA1(csq, 0),
                CHF_CONST_FRA1(a_dxFace, 0),
                CHF_CONST_REAL(a_alpha),
                CHF_CONST_REAL(a_beta),
                CHF_BOX(modBox),
                CHF_BOX(stenBox),
                CHF_CONST_INT(lohiSign),
                CHF_CONST_INT(a_dir));
            }
        }

      // Shift back csq
      csq.shiftHalf(a_dir, -lohiSign);
      U.shiftHalf(a_dir, -lohiSign);
    }; // lambda close

  // Now call the function on each respective side
  if (a_hasLo)
    {
      // Get the block boundary
      bcIdx.m_side = Side::Lo;
      bcIdx.assertValid();
      // Call the function with the corresponding box
      operateOnSide(bcIdx, a_loFaceBox);
    }
  if (a_hasHi)
    {
      // Get the block boundary
      bcIdx.m_side = Side::Hi;
      bcIdx.assertValid();
      // Call the function with the corresponding box
      operateOnSide(bcIdx, a_hiFaceBox);
    }
}

/*--------------------------------------------------------------------*/
//  This function is called to get boundary faces of a_box.
/** The boundary faces are cropped at all physical domain boundaries
 *  \param[out] a_boundaryBox
 *                      Face-centered box of boundary faces to fill
 *                      in.  Sub-box of a_dataFaceBox
 *  \param[in]  a_box   Box to find faces of.  This can be either
 *                      face-centered or cell-centered.
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_dir   Normal direction of faces
 *  \param[in]  a_side  Side of the domain where boundary faces are
 *                      needed
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getBoundaryFaces(Box&                  a_boundaryBox,
                         const Box&            a_box,
                         const BlockDomain&    a_domain,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side) const
{
  a_boundaryBox = Box();  // Return the empty box by default.
  if (!a_domain.isPeriodic(a_dir))
    {
      Box box = a_box;
      box &= a_domain;
      // Determine which side and thus shifting directions
      int lohisign = sign(a_side);
      box.shift(a_dir, lohisign);
      // Is there a domain boundary next to this grid?
      if (!a_domain.contains(box))
        {
          box &= a_domain;            // After cropping external cells
          CH_assert(!box.isEmpty());  // must still have interior cells
          if (a_side == Side::Lo)
            {
              a_boundaryBox = bdryLo(box, a_dir);
            }
          else  // (a_side == Side::Hi)
            {
              a_boundaryBox = bdryHi(box, a_dir);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  This function returns the cell centered box that intersects 2 physical
//  boundaries.
/**
 *  \param[out] a_boundaryBox
 *                      Cell-centered box of boundary adjacent cells
 *  \param[in]  a_box   Disjoint box
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_dir1  Normal direction of faces in first direction
 *  \param[in]  a_side1 Side of the domain for the first direction
 *  \param[in]  a_dir2  Normal direction of faces in second direction
 *  \param[in]  a_side2 Side of the domain for the second direction
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getBoundaryEdgeBox(Box&                  a_boundaryBox,
                           const Box&            a_box,
                           const BlockDomain&    a_domain,
                           const int&            a_dir1,
                           const Side::LoHiSide& a_side1,
                           const int&            a_dir2,
                           const Side::LoHiSide& a_side2) const
{
  a_boundaryBox = Box();  // Return the empty box by default.
  if (!a_domain.isPeriodic(a_dir1) && !a_domain.isPeriodic(a_dir2))
    {
      Box lDomain = grow(a_domain.domainBox(), 10);
      lDomain &= a_domain;
      Box box1 = adjCellBox(lDomain, a_dir1, a_side1, -1);
      Box box2 = adjCellBox(lDomain, a_dir2, a_side2, -1);
      a_boundaryBox = box1&box2;
      a_boundaryBox &= a_box;
// #if CH_SPACEDIM==3
//       // Get the third direction
//       int dir3 = SpaceDim - a_dir1 - a_dir2;
//       // Remove cells where box intersects 3 boundaries
//       if (!a_domain.isPeriodic(dir3))
//         {
//           lDomain.grow(dir3, -1);
//           a_boundaryBox &= lDomain;
//         }
// #endif
    }
}

/*--------------------------------------------------------------------*/
//  This function is called to get interior boundary faces of a_box.
/** The boundary faces are shrunk by one at all physical domain boundaries
 *  \param[out] a_boundaryBox
 *                      Face-centered box of boundary faces to fill in
 *  \param[in]  a_box   Box to find faces of.  This can be either
 *                      face-centered or cell-centered.
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_dir   Normal direction of faces
 *  \param[in]  a_side  Side of the domain where boundary faces are
 *                      needed
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getInteriorBoundaryFaces(Box&                  a_boundaryBox,
                                 const Box&            a_box,
                                 const BlockDomain&    a_domain,
                                 const int&            a_dir,
                                 const Side::LoHiSide& a_side) const
{
  a_boundaryBox = Box();  // Return the empty box by default.
  if (!a_domain.isPeriodic(a_dir))
    {
      Box box = a_box;
      box &= a_domain;
      // Determine which side and thus shifting directions
      int lohisign = sign(a_side);
      box.shift(a_dir, lohisign);
      // Is there a domain boundary next to this grid?
      if (!a_domain.contains(box))
        {
          Box modDomainBox = a_domain.domainBox();

          for (int tandir = 0; tandir != SpaceDim; ++tandir)
            {
              if (!a_domain.isPeriodic(tandir) && tandir != a_dir)
                {
                  modDomainBox.grow(tandir, -1);
                }
              else if (a_domain.isPeriodic(tandir))
                {
                  modDomainBox.grow(tandir, 8);
                }
            }
          box &= modDomainBox;        // After cropping external cells
          CH_assert(!box.isEmpty());  // must still have interior cells
          if (a_side == Side::Lo)
            {
              a_boundaryBox = bdryLo(box, a_dir);
            }
          else  // (a_side == Side::Hi)
            {
              a_boundaryBox = bdryHi(box, a_dir);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  This fills a box for boundary adjacent cells for specific boundary types
/** \param[out] a_boundAdjCells
 *                      Vector of boxes of boundary adjacent cells
 *  \param[in]  a_box   Input box
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_dir   Direction of interest normal to the boundary
 *  \param[in]  a_side  Current side
 *  \param[in]  a_level Current level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_numCells
 *                      Number of cells adjacent to the boundary to
 *                      return
 *  \param[in]  a_gridMetrics
 *                      The level grid metrics
 *  \param[in]  a_boundaryType
 *                      Boundary type of interest
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getBoundaryAdjCells(Vector<Box>&            a_boundAdjCells,
                            const Box&              a_box,
                            const BlockDomain&      a_domain,
                            const int&              a_dir,
                            const Side::LoHiSide&   a_side,
                            const int&              a_level,
                            const Real&             a_time,
                            const int               a_numCells,
                            LevelGridMetrics&       a_gridMetrics,
                            CRDparam::DomainBCType& a_boundaryType) const
{
  const int lohiSign = Side::sign(a_side);
  const MultiBlockCoordSys& mbcs = a_gridMetrics.getCoordSys();
  const int idxBlk = mbcs.whichBlock(a_domain.domainBox());
  // Check that this is a multi-block boundary
  const BlockBoundary& blkBdry = mbcs.boundary(idxBlk, a_dir, a_side);
  if (!blkBdry.isDomainBoundary())
    {
      return;
    }
  // The boundary index and type
  BoundaryIndex bcIdx;
  bcIdx.define(idxBlk, a_dir, a_side);
  auto domBC = getDomainBC(bcIdx);
  // Get the boundary box
  Box boundaryBox;
  getBoundaryFaces(boundaryBox, a_box, a_domain, a_dir, a_side);
  bool boundCheck = false;
  if (!boundaryBox.isEmpty())
    {
      if (CRDparam::DomainBCTypeMixed & domBC.m_type)
        {
          Vector<BCInfo> bcInfo;
          setMixedBC(boundaryBox, bcIdx,
                     a_domain.domainBox(), a_gridMetrics, a_time,
                     a_level, a_boundAdjCells, bcInfo);
          // Check if any boundaries are of the specified type boundaries
          for (int BCNum = 0; BCNum != bcInfo.size(); ++BCNum)
            {
              // If types match and box is not empty
              if ((a_boundaryType & bcInfo[BCNum].m_type)
                  && !a_boundAdjCells[BCNum].isEmpty())
                {
                  boundCheck = true;
                }
              else
                {
                  Box emptyBox;
                  a_boundAdjCells[BCNum].define(emptyBox);
                }
            }
        }
      else if (a_boundaryType & domBC.m_type)
        {
          a_boundAdjCells.resize(1);
          a_boundAdjCells[0] = boundaryBox;
          boundCheck = true;
        }
    }
  if (boundCheck)
    {
      // Opposite side for growDir
      Side::LoHiSide oppSide = Side::flip(a_side);
      for (int BCNum = 0; BCNum != a_boundAdjCells.size(); ++BCNum)
        {
          if (!a_boundAdjCells[BCNum].isEmpty())
            {
              a_boundAdjCells[BCNum].shiftHalf(a_dir, -lohiSign);
              a_boundAdjCells[BCNum].growDir(a_dir, oppSide, a_numCells-1);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  This function is called to get a list of wall boundary boxes
/** 
 *  \param[out] a_boundaryBoxes
 *                      Vector of face-centered boxes of the boundary
 *  \param[out] a_faceDirs
 *                      Vector of face normal directions, corresponding
 *                      to elements in a_boundaryBox
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_time  The current time
 *  \param[in]  a_level Current level
 *  \param[in]  a_gridMetrics
 *                      The level grid metrics
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getAllWallBoundaryFaces(Vector<Box>&         a_boundaryBoxes,
                                Vector<int>&         a_faceDirs,
                                const Real           a_time,
                                const Real           a_level,
                                LevelGridMetrics&    a_gridMetrics) const
{
  // Make sure vectors are empty
  a_boundaryBoxes.clear();
  a_faceDirs.clear();

  // The boundary index
  BoundaryIndex bcIdx;
  // Loop over all blocks
  for (bcIdx.m_block = 0;
       bcIdx.m_block != a_gridMetrics.getCoordSys().numBlocks();
       ++bcIdx.m_block)
    {
      const BlockDomain& domain =
        a_gridMetrics.getCoordSys().problemDomain(bcIdx.m_block);
      
      // Loop over each direction and sideto check boundaries
      for (bcIdx.m_dir = 0; bcIdx.m_dir != SpaceDim; ++bcIdx.m_dir)
        {
          Box dbox(domain.domainBox());
          Box bndFaceBox;
          for (auto side : EachSide)
            {
              bcIdx.m_side = side;
              auto domBC = getDomainBC(bcIdx);
              getBoundaryFaces(bndFaceBox, dbox, domain, bcIdx.m_dir, side);
              if (!bndFaceBox.isEmpty())
                {
                  // Check for a mixed boundary
                  if (CRDparam::DomainBCTypeMixed & domBC.m_type)
                    {
                      // Get the state for a mixed boundary
                      Vector<BCInfo> bcInfo;
                      Vector<Box> bcBoxes;
                      setMixedBC(bndFaceBox,
                                 bcIdx,
                                 dbox,
                                 a_gridMetrics,
                                 a_time,
                                 a_level,
                                 bcBoxes,
                                 bcInfo);
                      // Check each part of mixBC for walls
                      for (int bIdx = 0; bIdx != bcBoxes.size(); ++bIdx)
                        {
                          if (CRDparam::DomainBCTypeSWall
                              & bcInfo[bIdx].m_type)
                            {
                              // Append face box
                              a_boundaryBoxes.push_back(bcBoxes[bIdx]);
                              a_faceDirs.push_back(bcIdx.m_dir);
                            }
                        }
                    }
                  // Any other kind of wall
                  else if (CRDparam::DomainBCTypeSWall & domBC.m_type)
                    {
                      // Append face box
                      a_boundaryBoxes.push_back(bndFaceBox);
                      a_faceDirs.push_back(bcIdx.m_dir);
                    }
                }
            }  // Loop over sides
        }  // Loop over directions
    }  // Loop over blocks
}

/*--------------------------------------------------------------------*/
//  Set the primitive state at a domain boundary
/** In general, one should not need to override this routine.  To
 *  customize the BC imposed on the boundary face, use routine
 *  setImposedBCprimState.  This will not necessarily be the final
 *  state set on the face of flow BC (except for supersonic inflow),
 *  because a Riemann solution is used to resolve the difference
 *  between the imposed state and that determined by the interior
 *  scheme.
 *
 *  \param[in]  a_WavgFace
 *                      Should have an interpolation of the average
 *                      primitive state on the face from the interior
 *                      scheme
 *  \param[out] a_WavgFace
 *                      Average primitive state on the face computed
 *                      either through wall considerations or as a
 *                      Riemann solution between the interpolated
 *                      state and an imposed state.
 *  \param[in]  a_WavgCell
 *                      Average primitive state in core (interior)
 *                      cells
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_bndryCellFab
 *                      Previous-time data for cells immediately adjacent
 *                      to boundaries
 *  \param[in]  a_boundaryFaceBox
 *                      A box of the boundary faces to be adjusted
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box boundary faces to set adjacent ghost cells
 *                      If external ghost cells are not used, it is the
 *                      same as a_boundaryFaceBox
 *  \param[in]  a_disjointBox
 *                      Disjoint box of cells adjacent to boundary
 *  \param[in]  a_unitNormalBasisFab
 *                      A unit normal basis for this direction
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of the domain
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_prevDt
 *                      Time-step size for past time-step
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setBCprimState(FArrayBox&            a_WfaceAvgBndryDirFab,
                       FArrayBox&            a_WfaceAvgDirFab,
                       FArrayBox&            a_WcellAvgFab,
                       const FArrayBox&      a_bndryFaceAvgSlipVelDirFab,
                       FArrayBox&            a_bndryCellFab,
                       const Box&            a_boundaryFaceBox,
                       const Box&            a_boundaryFaceGhostBox,
                       const Box&            a_disjointBox,
                       const FArrayBox&      a_unitNormalBasisFab,
                       BCInfo&               a_bcInfo,
                       Box&                  a_bcBox,
                       const int             a_dir,
                       const Side::LoHiSide& a_side,
                       LevelGridMetrics&     a_gridMetrics,
                       const Real            a_time,
                       const Real            a_prevDt,
                       const int             a_level) const
{
  CH_TIME("CNSIBC::setBCprimState");
  //**FIXME: We are temporarily discarding the idea of mixed boundaries
  //         To bring this back, we should consider
  //         a) allowing multiple boundary types per boundary and a method to
  //            retrieve these using getDomainBC(bcIdx, bcNum) for example
  //         b) segmenting boundaries into arbitrary sub-boundaries, setting
  //            boundary types and states for each of these, and subsequently
  //            looping over every sub-boundary (similar to above except that
  //            this may be more complicated)

  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir,
               a_side);

  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int lohiSign = Side::sign(a_side);

  auto domBC = getDomainBC(bcIdx);

  // Set totalFaceBox over which the face-averaged primitive state will be set
  Box totalFaceBox = a_boundaryFaceBox;
  // If using viscous physics, we need to set one more tangential ghost face
  if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
    {
      totalFaceBox = a_boundaryFaceGhostBox;
    }

  a_bcBox.define(totalFaceBox);
  a_bcInfo = domBC;

  //**FIXME: check on these and set a lot more CH_asserts
  // CH_assert(a_WavgFace.contains(a_boundaryFaceBox));
  // CH_assert(a_boundarySlipVelocity.contains(a_boundaryFaceBox));

  // First, fill a_WfaceAvgBndryDirFab with the first layer of interior cells
  // Shift cells to align first interior layer with boundary faces
  a_WcellAvgFab.shiftHalf(a_dir, lohiSign);
  // Copy the first interior cell values into the complete boundary face box
  a_WfaceAvgBndryDirFab.copy(a_WcellAvgFab, totalFaceBox);

  // Next, fill a_WfaceAvgBndryDirFab with the limited fourth-order face-values
  // If this boundary is first order, leave it as interior cell values
  if (a_bcInfo.m_order != 1)
    {
      a_WfaceAvgBndryDirFab.copy(a_WfaceAvgDirFab, a_boundaryFaceBox);
    }

  // We need to set all boundary primitive states here
  if (CRDparam::DomainBCTypeAllWall & domBC.m_type) // any wall boundary
    {
      setWallBCprimState(a_WfaceAvgBndryDirFab,
                         totalFaceBox,
                         a_WcellAvgFab,
                         a_bndryFaceAvgSlipVelDirFab,
                         a_unitNormalBasisFab,
                         bcIdx,
                         a_disjointBox,
                         a_gridMetrics,
                         a_time,
                         a_level,
                         domBC);
    }
  else if (CRDparam::DomainBCTypeDirichlet & domBC.m_type)
    {
      // All comps specified and fixed, no Riemann solve should be necessary
      setImposedBCprimState(a_WfaceAvgBndryDirFab,
                            totalFaceBox,
                            a_WcellAvgFab,
                            a_unitNormalBasisFab,
                            bcIdx,
                            a_disjointBox,
                            a_gridMetrics,
                            a_time,
                            a_level,
                            domBC);
    }
  else if ((CRDparam::DomainBCTypeInflow & domBC.m_type) ||
           (CRDparam::DomainBCTypeOutflow & domBC.m_type) ||
           (CRDparam::DomainBCTypeFarfield & domBC.m_type) ||
           (CRDparam::DomainBCTypeRelaxedCBC & domBC.m_type))
    {
      // Temporarily store the values in a_WfaceAvgBndryDirFab for Riemann solve
      FABSTACKTEMP(WfaceAvgInterior, totalFaceBox, numWcomp);
      WfaceAvgInterior.copy(a_WfaceAvgBndryDirFab);
      // Temporarily create exterior state fab
      FABSTACKTEMP(WfaceAvgExterior, totalFaceBox, numWcomp);
      WfaceAvgExterior.copy(a_WfaceAvgBndryDirFab);
      // Set the exterior state
      setImposedBCprimState(WfaceAvgExterior,
                            totalFaceBox,
                            a_WcellAvgFab,
                            a_unitNormalBasisFab,
                            bcIdx,
                            a_disjointBox,
                            a_gridMetrics,
                            a_time,
                            a_level,
                            domBC);
      // Modify the exterior state if this is relaxed characteristic BCs
      if (CRDparam::DomainBCTypeRelaxedCBC & domBC.m_type)
        {
          // DEBUG ONLY
          if (!a_prevDt && CNSIBC::s_firstRKStage)
            {
              a_bndryCellFab.setVal(0.);
            }
          // Shift cells back to interior
          a_WcellAvgFab.shiftHalf(a_dir, -lohiSign);
          // Shift faces to exterior cell
          WfaceAvgExterior.shiftHalf(a_dir, lohiSign);
          setRelaxedCBCPrimState(
            WfaceAvgExterior, // update the exterior face-averaged state
            a_bndryCellFab,   // provides interior/exterior temporal-gradients
            a_WcellAvgFab,    // provides spatial-gradients for testing
            a_WfaceAvgDirFab, // could help modify BC using first interior-face
            a_unitNormalBasisFab,
            bcIdx,
            domBC,
            a_disjointBox,
            totalFaceBox,     // need to update on faces of totalFaceBox
            a_gridMetrics,
            a_dir,
            a_side,           // we need to determine inside/outside cells
            a_time,
            a_prevDt,
            a_level);
          // Shift cells to align first interior layer with boundary faces
          a_WcellAvgFab.shiftHalf(a_dir, lohiSign);
          // Shift faces back to boundary-face location
          WfaceAvgExterior.shiftHalf(a_dir, -lohiSign);
        }
      // Then solve a Riemann problem
      // Figure out if we can simplify this
      FArrayBox* WleftPtr;
      FArrayBox* WrightPtr;
      if (a_side == Side::Lo)
        {
          // Left is exterior
          WleftPtr = &WfaceAvgExterior;
          // Right is interior
          WrightPtr = &WfaceAvgInterior;
        }
      else  // High side
        {
          // Left is interior
          WleftPtr = &WfaceAvgInterior;
          // Right is exterior
          WrightPtr = &WfaceAvgExterior;
        }
      FArrayBox& Wleft  = *WleftPtr;
      FArrayBox& Wright = *WrightPtr;
      CRDparam::g_CRDPhysics->riemannBC(a_WfaceAvgBndryDirFab,
                                        Wleft,
                                        Wright,
                                        a_unitNormalBasisFab,
                                        a_dir,
                                        a_side,
                                        totalFaceBox);
      if (!(CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
        {
          // Compute extra primitive state
          CRDparam::g_CRDPhysics->extraPrimitiveState(a_WfaceAvgBndryDirFab,
                                                      totalFaceBox);
        }
    }
  else if (CRDparam::DomainBCTypeExtrapolated & domBC.m_type)
    {
      // Use the already extrapolated state
    }
  else if (CRDparam::DomainBCTypeMixed & domBC.m_type)
    {
      CRD::msg << "CNSIBC::setBCprimState: DomainBCTypeMixed is not supported!"
               << CRD::error;
    }
  else if (CRDparam::DomainBCTypeJet & domBC.m_type)
    {
      CRD::msg << "CNSIBC::setBCprimState: DomainBCTypeJet is not supported!"
               << CRD::error;
    }
  else if (CRDparam::DomainBCTypeCNSCBCInflow & domBC.m_type)
    {
      CRD::msg << "CNSIBC::setBCprimState: DomainBCTypeCNSCBCInflow"
               << " is not supported!" << CRD::error;
    }
  else if (CRDparam::DomainBCTypeCNSCBCOutflow & domBC.m_type)
    {
      CRD::msg << "CNSIBC::setBCprimState: DomainBCTypeCNSCBCOutflow"
               << " is not supported!" << CRD::error;
    }
  // Set any turbulent variables as necessary
  if (numTurbVar > 0)
    {
      CRDparam::g_CRDPhysics->setTurbulentBC(totalFaceBox,
                                             a_WfaceAvgBndryDirFab,
                                             a_WcellAvgFab,
                                             a_unitNormalBasisFab,
                                             a_dir,
                                             a_side,
                                             a_gridMetrics,
                                             a_time,
                                             a_level,
                                             domBC.m_type);
    }
  // Copy a_WfaceAvgBndryDirFab to a_WfaceAvgDirFab
  a_WfaceAvgDirFab.copy(a_WfaceAvgBndryDirFab);
  // Shift cells back to interior
  a_WcellAvgFab.shiftHalf(a_dir, -lohiSign);
}


/*--------------------------------------------------------------------*/
//  Set exterior ghost cells and final boundary values
/** \param[in]  a_WfaceBdryFab
 *                      Temporary storage of final boundary values. Needed
 *                      Outside this function for turbulence models
 *  \param[in]  a_WavgFace
 *                      Face-averaged primitive variables
 *  \param[out] a_WavgFace
 *                      Face values corrected on boundary faces where necessary
 *  \param[out] a_WavgCell
 *                      Cell-averaged values with external ghost cells filled
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_boundaryFaceBox
 *                      Face-centered box of boundary
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box boundary faces to set adjacent ghost cells
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_unitNormalBasisFab
 *                      Fab of unit normal
 *  \param[in]  a_dir   Directions normal to boundary
 *  \param[in]  a_side  Low or high side
 *  \param[in]  a_bcInfo
 *                      Boundary types corresponding to a_bcBoxes
 *  \param[in]  a_bcBoxes
 *                      Vector of boundary boxes corresponding to a_bcInfo
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBC::extrapolateDomainGhostCells(
  const FArrayBox&      a_WfaceAvgBndryDirFab,
  const FArrayBox&      a_WfaceAvgDirFab,
  FArrayBox&            a_WcellAvgFab,
  const Box&            a_boundaryFaceBox,
  const Box&            a_boundaryFaceGhostBox,
  const Box&            a_disjointBox,
  const FArrayBox&      a_unitNormalBasisFab,
  const int             a_dir,
  const Side::LoHiSide& a_side,
  const BCInfo&         a_bcInfo,
  const Box&            a_bcBox,
  LevelGridMetrics&     a_gridMetrics,
  const Real            a_time,
  const int             a_level) const
{
  CH_TIME("CNSIBC::extrapolateDomainGhostCells");
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir, a_side);
  const int lohiSign = sign(a_side);
  const Real dx = a_gridMetrics.dxVect()[a_dir];
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();

  //**FIXME: Check these and add a lot more CH_asserts
  // CH_assert(a_WavgFace.contains(a_boundaryFaceBox));
  // CH_assert(a_boundarySlipVelocity.contains(a_boundaryFaceBox));

  // Box of faces adjacent to ghost cells that must be filled
  Box ghostFaceBox(a_bcBox);
  ghostFaceBox &= a_boundaryFaceGhostBox;
  // Box of faces of boundary values that are adjusted where necessary
  Box faceAvgBox(a_bcBox);
  faceAvgBox &= a_boundaryFaceBox;

  // Vector of components to be set with Neumann conditions
  Vector<int> NeumannComp(1, -1);
  // Vector of components to extrapolate ghost cells directly from interior
  Vector<int> interiorExtrap(1, -1);
  // Vector of components to extrapolate using values set on faces
  Vector<int> faceValExtrap(1, -1);

  // Set the boundary condition info
  auto domT = a_bcInfo;

  // FArrayBox for the face gradients needed in Neumann conditions
  FABSTACKTEMP(WgradFace, ghostFaceBox, numWcomp);
  // Fill any imposed boundary conditions
  if (CRDparam::DomainBCTypeCNSCBC & domT.m_type)
    {
      //**FIXME: Originally, this method determined the gradient of W on the
      //         faces and then used this in conjunction with interior cell
      //         data to set the exterior ghost cells. If this is to be
      //         repeated in the future, one of the following should be done
      //         a) set the boundary derivatives here and then use this to
      //            determine the domain ghost cells
      //         b) set the boundary derivatives using another function in
      //            Inertial4thOrderOp, provide them here, and then specify
      //            the domain ghost cells here
      //         As it is, this BC type never truly worked anyway.
      CRD::msg << "CNSIBC::extrapolateDomainGhostCells: DomainBCTypeCNSCBC"
               << " is not supported!" << CRD::error;
    }

  int solveVar = setExtraBCs(NeumannComp, interiorExtrap, faceValExtrap,
                             WgradFace, bcIdx,
                             a_disjointBox, a_gridMetrics,
                             a_time, a_level, domT);
  // Align first external ghost cells with boundary faces
  a_WcellAvgFab.shiftHalf(a_dir, -lohiSign);

  //**FIXME: what's going on here? Isn't this redundant?
  // Box for the boundary face values for a_WfaceBdryFab,
  //  the first external ghost layer for a_WavgCell
  Box ghostBox1(ghostFaceBox);
  // Box for the boundary face values for a_WavgFace
  Box faceAvgGhostBox1(faceAvgBox);

  //**FIXME: Are the Neumann conditions imposed correctly? Why would we try to
  //         overwrite the face values? Maybe the limited interpolated state
  //         computed in Inertial4thOrderOp should take boundary conditions into
  //         account when computing the extrapolation of boundary values? Maybe
  //         it should even take boundary conditions into account when computing
  //         the interpolated face values on the interior? In any case, the face
  //         values should be completely set before this and never touched again
  // Extrapolate any Neumann conditions
  if (NeumannComp[0] != -1)
    {
      //**FIXME: Rather than allocate more memory, we should decide whether to
      //         make a_WfaceAvgBndryDirFab non-const as it's brought into this
      //         function (extrapolateDomainGhostCells) or if it should be const
      //         as it goes into extrapGhostNeumann.
      FABSTACKTEMP(WfaceTemp, ghostFaceBox, numWcomp);
      WfaceTemp.copy(a_WfaceAvgBndryDirFab);
      extrapGhostNeumann(ghostBox1, a_WcellAvgFab, WfaceTemp,
                         WgradFace, a_side, a_dir, dx, NeumannComp);
    }
  // Extrapolate any conditions from interior
  if (interiorExtrap[0] != -1)
    {
      extrapGhostFromInterior(ghostBox1, a_WcellAvgFab, a_side, a_dir,
                              interiorExtrap);
    }
  // Extrapolate any ghost cells from face values
  if (faceValExtrap[0] != -1)
    {
      extrapGhostFromFaceVals(ghostBox1, a_WcellAvgFab, a_WfaceAvgBndryDirFab,
                              a_side, a_dir, faceValExtrap);
    }
  // Both layers of ghost cells
  Box ghostCellBox(ghostBox1);
  ghostCellBox.shift(a_dir, lohiSign);
  ghostCellBox.minBox(ghostBox1);
  // Solve for density, pressure, or temperature depending on setExtraVarBCs
  // FIXME: Solving for the 3rd variable means the BC is consistent
  //        but no longer stable. Must find solution to this
  // FIXME: Additionally, a_WfaceAvgBndryDirFab and a_WfaceAvgDirFab are now
  //        const inputs.
  solveVar = -1;
  switch(solveVar)
    {
    case 1:
      // CRDparam::g_CRDPhysics->density(a_WfaceAvgBndryDirFab, ghostBox1);
      // CRDparam::g_CRDPhysics->density(a_WfaceAvgDirFab, faceAvgGhostBox1);
      CRDparam::g_CRDPhysics->density(a_WcellAvgFab, ghostCellBox);
      break;
    case 2:
      // CRDparam::g_CRDPhysics->pressure(a_WfaceAvgBndryDirFab, ghostBox1);
      // CRDparam::g_CRDPhysics->pressure(a_WfaceAvgDirFab, faceAvgGhostBox1);
      CRDparam::g_CRDPhysics->pressure(a_WcellAvgFab, ghostCellBox);
      break;
    case 0:
      // CRDparam::g_CRDPhysics->temperature(a_WfaceAvgBndryDirFab, ghostBox1);
      // CRDparam::g_CRDPhysics->temperature(a_WfaceAvgDirFab, faceAvgGhostBox1);
      CRDparam::g_CRDPhysics->temperature(a_WcellAvgFab, ghostCellBox);
      break;
    default:
      break;
    }
  // Shift cells back to exterior
  a_WcellAvgFab.shiftHalf(a_dir, lohiSign);
}

/*--------------------------------------------------------------------*/
//  Set the primitive state at a domain boundary -- DEPRECATED
/** In general, one should not need to override this routine.  To
 *  customize the BC imposed on the boundary face, use routine
 *  setImposedBCprimState.  This will not necessarily be the final
 *  state set on the face of flow BC (except for supersonic inflow),
 *  because a Riemann solution is used to resolve the difference
 *  between the imposed state and that determined by the interior
 *  scheme.
 *
 *  \param[in]  a_WavgFace
 *                      Should have an interpolation of the average
 *                      primitive state on the face from the interior
 *                      scheme
 *  \param[out] a_WavgFace
 *                      Average primitive state on the face computed
 *                      either through wall considerations or as a
 *                      Riemann solution between the interpolated
 *                      state and an imposed state.
 *  \param[in]  a_WavgCell
 *                      Average primitive state in core (interior)
 *                      cells
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_boundaryFaceBox
 *                      A box of the boundary faces to be adjusted
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box boundary faces to set adjacent ghost cells
 *                      If external ghost cells are not used, it is the
 *                      same as a_boundaryFaceBox
 *  \param[in]  a_disjointBox
 *                      Disjoint box of cells adjacent to boundary
 *  \param[in]  a_unitNormalBasisFab
 *                      A unit normal basis for this direction
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of the domain
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

//**DEPRECATED

void
CNSIBC::setBCprimState(
  FArrayBox&                      a_WavgFace,
  const FArrayBox&                a_WavgCell,
  const FArrayBox&                a_boundarySlipVelocity,
  const Box&                      a_boundaryFaceBox,
  const Box&                      a_boundaryFaceGhostBox,
  const Box&                      a_disjointBox,
  const FArrayBox&                a_unitNormalBasisFab,
  Vector<BCInfo>&                 a_bcInfo,
  Vector<Box>&                    a_bcBoxes,
  const int                       a_dir,
  const Side::LoHiSide&           a_side,
  LevelGridMetrics&               a_gridMetrics,
  const Real                      a_time,
  const int                       a_level) const
{
  CH_TIME("CNSIBC::setBCprimState");
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir,
               a_side);

  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int lohiSign = Side::sign(a_side);
  
  auto domBC = getDomainBC(bcIdx);

//--Shift cell array to overlap boundary faces.  For all setWallBC* and
//--setImposedBC*, the first interior cells is indexed the same as the face.

  FArrayBox& Wcell = (FArrayBox&)a_WavgCell;
  Wcell.shiftHalf(a_dir, lohiSign);   // Shift cell to overlap face

//--Set the state on the face (often, a derived virtual function will be invoked
//--here)

  // Count the number of BC's on a single boundary
  Box totalFaceBox = a_boundaryFaceBox;
  if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
    {
      totalFaceBox = a_boundaryFaceGhostBox;
    }
  if (CRDparam::DomainBCTypeMixed & domBC.m_type)
    {
      // Initialize the order to what is set from before
      int setBC = setMixedBC(totalFaceBox,
                             bcIdx,
                             a_disjointBox,
                             a_gridMetrics,
                             a_time,
                             a_level,
                             a_bcBoxes,
                             a_bcInfo);
      CH_assert(setBC == 0);
    }
  else
    {
      a_bcBoxes.resize(1);
      a_bcInfo.resize(1);
      a_bcBoxes[0].define(totalFaceBox);
      a_bcInfo[0] = domBC;
    }

  for (int BCNum = 0; BCNum != a_bcInfo.size(); ++BCNum)
    {
      Box modBox(a_bcBoxes[BCNum]);
      // Make sure the box types are the same, even if they are empty
      CH_assert(modBox.ixType() == a_boundaryFaceBox.ixType());
      modBox &= a_boundaryFaceBox;
      // Do not do the operations if the box is empty
      if (modBox.isEmpty()) continue;
      auto domT = a_bcInfo[BCNum];

//--Reduce the order of the solution if required.  By default, a_WavgFace has
//--a limited fourth-order interpolation.  If first-order is requested, just use
//--piece-wise constant data from the cell.

      if (a_bcInfo[BCNum].m_order == 1)
        {
          a_WavgFace.copy(Wcell, modBox, 0, modBox, 0, numWcomp);
        }
      // FIXME: This is only necessary for jet, inflow, outflow should change
      // since we want to set the ghost cells to be the BC values before
      // the Riemann solution
      FABSTACKTEMP(Wexterior, modBox, numWcomp);
      Wexterior.copy(a_WavgFace, 0, 0, numWcomp);
      // Ensure the mixed type is not being used as an actual BC type
      CH_assert(!(CRDparam::DomainBCTypeMixed & domT.m_type));
      if (CRDparam::DomainBCTypeCNSCBC & domT.m_type)
        {
          setCNSCBCProfiles(a_WavgFace, modBox, Wcell, a_unitNormalBasisFab,
                            bcIdx, a_disjointBox, a_gridMetrics,
                            a_time, a_level, domT);
        }
      else if (CRDparam::DomainBCTypeAllWall & domT.m_type)
        // A wall BC
        {
          setWallBCprimState(a_WavgFace,
                             modBox,
                             Wcell,
                             a_boundarySlipVelocity,
                             a_unitNormalBasisFab,
                             bcIdx,
                             a_disjointBox,
                             a_gridMetrics,
                             a_time,
                             a_level,
                             domT);
        }
      else if (CRDparam::DomainBCTypeDirichlet & domT.m_type)
        // Fixed BC
        {
          setImposedBCprimState(a_WavgFace,
                                modBox,
                                Wcell,
                                a_unitNormalBasisFab,
                                bcIdx,
                                a_disjointBox,
                                a_gridMetrics,
                                a_time,
                                a_level,
                                domT);
        }
      else if (CRDparam::DomainBCTypeExtrapolated & domT.m_type)
        {
          // Use the already extrapolated state
        }
      else
        // Solve a Riemann problem to select between interior and exterior state
        // This applies for inlet, outlet, and jet BCs
        {
          setImposedBCprimState(Wexterior,
                                modBox,
                                Wcell,
                                a_unitNormalBasisFab,
                                bcIdx,
                                a_disjointBox,
                                a_gridMetrics,
                                a_time,
                                a_level,
                                domT);
          // Average the exterior BC if requested
          // if (CRDparam::g_averageBC < 1.0)
          //   {
          //     averageBC
          //   }
          // Interior is a copy of a_WavgFace
          FABSTACKTEMP(Winterior, modBox, numWcomp);
          Winterior.copy(a_WavgFace);

          FArrayBox* WleftPtr;
          FArrayBox* WrightPtr;
          if (a_side == Side::Lo)
            {
              // Left is exterior
              WleftPtr = &Wexterior;
              // Right is interior
              WrightPtr = &Winterior;
            }
          else  // High side
            {
              // Left is interior
              WleftPtr = &Winterior;
              // Right is exterior
              WrightPtr = &Wexterior;
            }
          FArrayBox& Wleft  = *WleftPtr;
          FArrayBox& Wright = *WrightPtr;
          CRDparam::g_CRDPhysics->riemannBC(a_WavgFace,
                                            Wleft,
                                            Wright,
                                            a_unitNormalBasisFab,
                                            a_dir,
                                            a_side,
                                            modBox);
          // If inflow BC, the velocity is fixed
          //**FIXME - this isn't correct.  Need a wave to fix BC.
          // if (CRDparam::DomainBCTypeInflow & m_domainBC[a_dir][lohiIdx].m_type)
          //   {
          //     FORT_REVERSETRANSFORMF(CHF_FRA(Wexterior),
          //                            CHF_CONST_FRA(a_unitNormalBasisFab),
          //                            CHF_BOX(a_boundaryFaceBox));
          //     a_WavgFace.copy(Wexterior, WVELX, WVELX, SpaceDim);
          //   }
        }
      if (numTurbVar > 0)
        {
          CRDparam::g_CRDPhysics->setTurbulentBC(modBox,
                                                 a_WavgFace,
                                                 Wcell,
                                                 a_unitNormalBasisFab,
                                                 a_dir,
                                                 a_side,
                                                 a_gridMetrics,
                                                 a_time,
                                                 a_level,
                                                 domT.m_type);
        }
    }
  Wcell.shiftHalf(a_dir, -lohiSign);
}

/*--------------------------------------------------------------------*/
//  Set exterior ghost cells and final boundary values -- DEPRECATED
/** \param[in]  a_WfaceBdryFab
 *                      Temporary storage of final boundary values. Needed
 *                      Outside this function for turbulence models
 *  \param[in]  a_WavgFace
 *                      Face-averaged primitive variables
 *  \param[out] a_WavgFace
 *                      Face values corrected on boundary faces where necessary
 *  \param[out] a_WavgCell
 *                      Cell-averaged values with external ghost cells filled
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_boundaryFaceBox
 *                      Face-centered box of boundary
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box boundary faces to set adjacent ghost cells
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_unitNormalBasisFab
 *                      Fab of unit normal
 *  \param[in]  a_dir   Directions normal to boundary
 *  \param[in]  a_side  Low or high side
 *  \param[in]  a_bcInfo
 *                      Boundary types corresponding to a_bcBoxes
 *  \param[in]  a_bcBoxes
 *                      Vector of boundary boxes corresponding to a_bcInfo
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

//**DEPRECATED

void
CNSIBC::applyBCStates(
  FArrayBox&                            a_WfaceBdryFab,
  FArrayBox&                            a_WavgFace,
  FArrayBox&                            a_WavgCell,
  const FArrayBox&                      a_boundarySlipVelocity,
  const Box&                            a_boundaryFaceBox,
  const Box&                            a_boundaryFaceGhostBox,
  const Box&                            a_disjointBox,
  const FArrayBox&                      a_unitNormalBasisFab,
  const int                             a_dir,
  const Side::LoHiSide&                 a_side,
  const Vector<BCInfo>&                 a_bcInfo,
  const Vector<Box>&                    a_bcBoxes,
  LevelGridMetrics&                     a_gridMetrics,
  const Real                            a_time,
  const int                             a_level) const
{
  CH_TIME("CNSIBC::applyBCStates");
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir, a_side);
  const int lohiSign = sign(a_side);
  const Real dx = a_gridMetrics.dxVect()[a_dir];
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  CH_assert(a_WavgFace.contains(a_boundaryFaceBox));
  CH_assert(a_boundarySlipVelocity.contains(a_boundaryFaceBox));
  // Shift cells to align first interior layer with boundary faces
  a_WavgCell.shiftHalf(a_dir, lohiSign);
  // Copy the first interior cell values into the boundary faces
  a_WfaceBdryFab.copy(a_WavgCell, a_boundaryFaceGhostBox);

  for (int bcNum = 0; bcNum != a_bcInfo.size(); ++bcNum)
    {
      // Box of faces adjacent to ghost cells that must be filled
      Box ghostFaceBox(a_bcBoxes[bcNum]);
      ghostFaceBox &= a_boundaryFaceGhostBox;
      // Box of faces of boundary values that are adjusted where necessary 
      Box faceAvgBox(a_bcBoxes[bcNum]);
      faceAvgBox &= a_boundaryFaceBox;
      if (ghostFaceBox.isEmpty()) continue;
      // Vector of components to be set with Neumann conditions
      Vector<int> NeumannComp(1, -1);
      // Vector of components to extrapolate ghost cells directly from interior
      Vector<int> interiorExtrap(1, -1);
      // Vector of components to extrapolate using values set on faces
      Vector<int> faceValExtrap(1, -1);
      auto domT = a_bcInfo[bcNum];
      // FArrayBox for the face gradients needed in Neumann conditions
      FABSTACKTEMP(WgradFace, ghostFaceBox, numWcomp);
      // Fill any imposed boundary conditions
      if (CRDparam::DomainBCTypeCNSCBC & domT.m_type)
        {
          setCNSCBCProfiles(a_WfaceBdryFab, ghostFaceBox, a_WavgCell,
                            a_unitNormalBasisFab, bcIdx,
                            a_disjointBox, a_gridMetrics,
                            a_time, a_level, domT);
          a_WavgCell.shiftHalf(a_dir, -lohiSign);
          applyCNSCBC(a_WfaceBdryFab, ghostFaceBox, a_WavgCell,
                      a_unitNormalBasisFab, bcIdx,
                      a_disjointBox, a_gridMetrics,
                      a_time, a_level, false, domT);
          a_WavgCell.shiftHalf(a_dir, lohiSign);
          a_WavgFace.copy(a_WfaceBdryFab, faceAvgBox, 0,
                          faceAvgBox, 0, numWcomp);
          continue;
        }
      if (CRDparam::DomainBCTypeAllWall & domT.m_type)
        // A wall BBC
        {
          setWallBCprimState(a_WfaceBdryFab, ghostFaceBox, a_WavgCell,
                             a_boundarySlipVelocity, a_unitNormalBasisFab,
                             bcIdx, a_disjointBox, a_gridMetrics,
                             a_time, a_level, domT);
        }
      else if (CRDparam::DomainBCTypeDirichlet & domT.m_type)
        // Fixed BC
        {
          setImposedBCprimState(a_WfaceBdryFab, ghostFaceBox, a_WavgCell,
                                a_unitNormalBasisFab, bcIdx,
                                a_disjointBox, a_gridMetrics,
                                a_time, a_level, domT);
        }
      int solveVar = setExtraBCs(NeumannComp, interiorExtrap, faceValExtrap,
                                 WgradFace, bcIdx,
                                 a_disjointBox, a_gridMetrics,
                                 a_time, a_level, domT);
      // Align first external ghost cells with boundary faces
      a_WavgCell.shift(a_dir, -lohiSign);
      // Box for the boundary face values for a_WfaceBdryFab,
      //  the first external ghost layer for a_WavgCell
      Box ghostBox1(ghostFaceBox);
      // Box for the boundary face values for a_WavgFace
      Box faceAvgGhostBox1(faceAvgBox);

      // Extrapolate any Neumann conditions and fill face values
      if (NeumannComp[0] != -1)
        {
          extrapGhostNeumann(ghostBox1, a_WavgCell, a_WfaceBdryFab, WgradFace,
                             a_side, a_dir, dx, NeumannComp);
          // FIXME: This leads to instabilities at walls
          // Copy corrected face values into shifted a_WavgFace
          // for (int neuComp = 0; neuComp != NeumannComp.size(); ++neuComp)
          //   {
          //     int comp = NeumannComp[neuComp];
          //     a_WavgFace.copy(a_WfaceBdryFab, faceAvgGhostBox1, comp,
          //                     faceAvgGhostBox1, comp, 1);
          //   }
        }
      // Extrapolate any conditions from interior
      if (interiorExtrap[0] != -1)
        {
          extrapGhostFromInterior(ghostBox1, a_WavgCell, a_side, a_dir,
                                  interiorExtrap);
        }
      // Extrapolate any ghost cells from face values
      if (faceValExtrap[0] != -1)
        {
          extrapGhostFromFaceVals(ghostBox1, a_WavgCell, a_WfaceBdryFab,
                                  a_side, a_dir, faceValExtrap);
        }
      // Both layers of ghost cells
      Box ghostCellBox(ghostBox1);
      ghostCellBox.shift(a_dir, lohiSign);
      ghostCellBox.minBox(ghostBox1);
      // Solve for density, pressure, or temperature depending on setExtraVarBCs
      // FIXME: Solving for the 3rd variable means the BC is consistent
      //        but no longer stable. Must find solution to this
      solveVar = -1;
      switch(solveVar)
        {
        case 1:
          CRDparam::g_CRDPhysics->density(a_WfaceBdryFab, ghostBox1);
          CRDparam::g_CRDPhysics->density(a_WavgFace, faceAvgGhostBox1);
          CRDparam::g_CRDPhysics->density(a_WavgCell, ghostCellBox);
          break;
        case 2:
          CRDparam::g_CRDPhysics->pressure(a_WfaceBdryFab, ghostBox1);
          CRDparam::g_CRDPhysics->pressure(a_WavgFace, faceAvgGhostBox1);
          CRDparam::g_CRDPhysics->pressure(a_WavgCell, ghostCellBox);
          break;
        case 0:
          CRDparam::g_CRDPhysics->temperature(a_WfaceBdryFab, ghostBox1);
          CRDparam::g_CRDPhysics->temperature(a_WavgFace, faceAvgGhostBox1);
          CRDparam::g_CRDPhysics->temperature(a_WavgCell, ghostCellBox);
          break;
        default:
          break;
        }
      // Shift back to align first interior layer with boundary faces
      a_WavgCell.shift(a_dir, lohiSign);
    }
  // Shift cell back
  a_WavgCell.shiftHalf(a_dir, -lohiSign);
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBC::haveExactSol() const
{
  return false;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state in the cells (defaults to a
//  no-op)
/** \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update with \<U\>.  One
 *                      assumes that this is larger than the disjoint
 *                      box so that \<JU\> can be computed.
 *  \param[in]  a_disjointBox
 *                      The disjoint box
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
 *    <li> The average conservative physical state \<U\> is returned
 *         in a_Ux.
 *  </ul>
 *//*-----------------------------------------------------------------*/

int 
CNSIBC::exactSol(FArrayBox&              a_Ux,
                 const Box&              a_box,
                 const Box&              a_disjointBox,
                 const LevelGridMetrics& a_gridMetrics,
                 const FluxBox&          a_unitNormals,
                 const DataIndex&        a_didx,
                 const Real              a_time,
                 const int               a_level) const
{
  return 1;
}

/*--------------------------------------------------------------------*/
//  Add source term (defaults to a no-op)
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
 *                      Block domain
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Grid level
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
CNSIBC::addSourceTerm(FArrayBox&           a_sourceFab,
                      FArrayBox&           a_invDtFab,
                      const FArrayBox&     a_Wcell,
                      const FArrayBox&     a_UcellAvg,
                      const FluxBox&       a_WfaceAvgFxb,
                      const BlockDomain&   a_domain,
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
  return;
}

/*--------------------------------------------------------------------*/
//  Solve the transverse derivatives in the exterior ghost cells
/** \param[in]  a_PhicellAvgFab
 *                      Cell-averaged phi known at interior cells and
 *                      ghost cells
 *  \param[out] a_GradPhicellAvgFab
 *                      The gradient of phi FAB
 *  \param[in]  a_box   The canonical box that we are working on.  It
 *                      must be adjacent to a domain extent where BC
 *                      need to be applied
 *  \param[in]  a_domain
 *                      Block domain
 *  \param[in]  a_dxVect
 *                      The computational grid spacing
 *  \param[in]  a_solvePhi
 *                      true - For solving single variable
 *                      false - For solving tensor variables
 *  \param[in]  a_level Grid level
 *  \note
 *  <ul>
 *    <li> This routine should only be called for boxes adjacent to a
 *         domain extent where BC need to be applied
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBC::cellAvgGradBC(const FArrayBox&     a_PhicellAvgFab,
                      FArrayBox&           a_GradPhicellAvgFab,
                      const Box&           a_box,
                      const BlockDomain&   a_domain,
                      const RealVect&      a_dxVect,
                      const bool           a_solvePhi,
                      const int            a_level) const
{
  CH_TIME("CNSIBC::cellAvgGradBC");
  CRD::msg << CRD::fv4 <<"Solve cell-averaged gradients at boundaries "
           << CRD::end;

  // Set up all boxes needed for corner computations
  // Vectors are arranged by component as
  // (boundary direction, tan direction)
  // For example: if we have a 3D box that is non-periodic in any direction
  // X bdrys: boxLo1[0]: Box that is adjacent to the lower Y boundary
  //          that includes 2 ghost ghost cells in the X direction
  //          boxLo1[1]: Box that is adjacent to the lower Z boundary
  //          that includes 2 ghost ghost cells in the X direction
  //
  // Y bdrys: boxLo1[2]: Box that is adjacent to the lower X boundary
  //          that includes 2 ghost ghost cells in the Y direction
  //          boxLo1[3]: Box that is adjacent to the lower Z boundary
  //          that includes 2 ghost ghost cells in the Y direction
  //
  // Z bdrys: boxLo1[4]: Box that is adjacent to the lower X boundary
  //          that includes 2 ghost ghost cells in the Z direction
  //          boxLo1[5]: Box that is adjacent to the lower Y boundary
  //          that includes 2 ghost ghost cells in the Z direction
  // If the Y boundary were periodic, components 0,2,3 and 5 would be empty
  stc::Vector<int, SpaceDim*(SpaceDim-1)> tanGradDirs;
  stc::Vector<Box, SpaceDim*(SpaceDim-1)> boxLo1;
  stc::Vector<Box, SpaceDim*(SpaceDim-1)> boxLo2;
  stc::Vector<Box, SpaceDim*(SpaceDim-1)> boxHi1;
  stc::Vector<Box, SpaceDim*(SpaceDim-1)> boxHi2;

  // Test to see if the box is adjacent to two physical domain boundaries
  // The artifacts are just used for cropping later
  int cornerTest = 0;
  for (const int dir : EachDir)
    {
      if ((!a_domain.isConnected(dir, Side::Lo)) ||
          (!a_domain.isConnected(dir, Side::Hi)))
        {
          Box dBox(a_domain.domainBox());
          dBox.grow(dir, -2);
          if (!dBox.contains(a_box))
            {
              dBox.grow(2);       // Domain box +2 in all directions
              dBox &= a_domain;   // ... except at all domain boundaries
              dBox.grow(dir, 2);  // ... and in direction of interest
              int k = 0;
              for (const int tanDir : EachDir)
                {
                  if (tanDir != dir)
                    {
                      Box dBoxTm1(dBox);
                      dBoxTm1.grow(tanDir, -1);
                      Box dBoxTm2(dBox);
                      dBoxTm2.grow(tanDir, -2);
                      if (!dBoxTm1.contains(a_box))
                        {
                          const int comp = dir*(SpaceDim-1) + k;
                          for (const auto tanSide : EachSide)
                            {
                              if (!a_domain.isConnected(tanDir, tanSide))
                                {
                                  // Have two physical domain boundaries at
                                  // the edge of dir and tanDir+tanSide
                                  ++cornerTest;
                                  tanGradDirs[comp] = tanDir;
                                  if (tanSide == Side::Lo)
                                    {
                                      boxLo1[comp] =
                                        adjCellLo(dBoxTm1, tanDir, 1);
                                      boxLo2[comp] =
                                        adjCellLo(dBoxTm2, tanDir, 1);
                                    }
                                  else
                                    {
                                      boxHi1[comp] =
                                        adjCellHi(dBoxTm1, tanDir, 1);
                                      boxHi2[comp] =
                                        adjCellHi(dBoxTm2, tanDir, 1);
                                    }
                                }
                            }
                        }  // Box adjacent to tangential boundary
                      ++k;
                    }
                }  // Loop over tangential directions
            }
        }
    }

  int phiComp = SpaceDim;
  int numGhostTan = 2;
  // This checks to see if we are solving the velocity gradient or 
  // phi (single variable) gradient
  if (a_solvePhi)
    {
      numGhostTan = 2;
      phiComp = 1;
    }

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if ((!a_domain.isConnected(dir, Side::Lo)) ||
          (!a_domain.isConnected(dir, Side::Hi)))
        {
          const IntVect growVect = IntVect_unit - IntVect_basis(dir);
          Box loBox;
          Box hiBox;
          Box centerBox;
          Box entireBox;
          int hasLo;
          int hasHi;
          Box bdryBox;
          Box combinedBox;
          Box inBox = a_box;
          Side::LoHiSide side;
          inBox.grow(dir,1);
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, inBox,
                     a_domain, dir);
          if (hasLo)
            {
              combinedBox = adjCellLo(entireBox,dir,2);
              combinedBox.grow(numGhostTan*growVect);
              BlockDomain limBox = a_domain;
              limBox.growLo(dir,2);
              combinedBox &= limBox;
              side = Side::Lo;
              // Solve for cell-averaged gradients in ghost cells
              computeCellAvgGradBdry(a_GradPhicellAvgFab,
                                     a_PhicellAvgFab,
                                     combinedBox,
                                     dir,
                                     side,
                                     phiComp,
                                     a_dxVect);
              // Check and solve for corners
              if (cornerTest > 0)
                {
                  Side::LoHiSide tanSide;
                  for (int k = 0; k != SpaceDim - 1; ++k)
                    {
                      int comp = dir*(SpaceDim-1) + k;
                      int gradDir = tanGradDirs[comp];
                      // Solve for low corners
                      Box cornerBox = combinedBox;
                      cornerBox &= boxLo1[comp];
                      if (!cornerBox.isEmpty())
                        {
                          tanSide = Side::Lo;
                          Box biasedBox = combinedBox;
                          biasedBox &= boxLo2[comp];
                          computeCornerGrad(a_GradPhicellAvgFab,
                                            a_PhicellAvgFab,
                                            cornerBox, biasedBox,
                                            gradDir, tanSide,
                                            phiComp, a_dxVect[gradDir]);
                        }
                      // Solve for high corners
                      cornerBox = combinedBox;
                      cornerBox &= boxHi1[comp];
                      if (!cornerBox.isEmpty())
                        {
                          tanSide = Side::Hi;
                          Box biasedBox = combinedBox;
                          biasedBox &= boxHi2[comp];
                          computeCornerGrad(a_GradPhicellAvgFab,
                                            a_PhicellAvgFab,
                                            cornerBox, biasedBox,
                                            gradDir, tanSide,
                                            phiComp, a_dxVect[gradDir]);
                        }
                    }
                }
            }
          if (hasHi)
            {
              combinedBox = adjCellHi(entireBox,dir,2);
              combinedBox.grow(numGhostTan*growVect);
              BlockDomain limBox = a_domain;
              limBox.growHi(dir, 2);
              combinedBox &= limBox;
              side = Side::Hi;
              // Solve for cell-averaged gradients in ghost cells
              computeCellAvgGradBdry(a_GradPhicellAvgFab,
                                     a_PhicellAvgFab,
                                     combinedBox,
                                     dir,
                                     side,
                                     phiComp,
                                     a_dxVect);
              // Check and solve for corners
              if (cornerTest > 0)
                {
                  Side::LoHiSide tanSide;
                  for (int k = 0; k != SpaceDim - 1; ++k)
                    {
                      int comp = dir*(SpaceDim-1) + k;
                      int gradDir = tanGradDirs[comp];
                      // Solve for low corners
                      Box cornerBox = combinedBox;
                      cornerBox &= boxLo1[comp];
                      if (!cornerBox.isEmpty())
                        {
                          tanSide = Side::Lo;
                          Box biasedBox = combinedBox;
                          biasedBox &= boxLo2[comp];
                          computeCornerGrad(a_GradPhicellAvgFab,
                                            a_PhicellAvgFab,
                                            cornerBox, biasedBox,
                                            gradDir, tanSide,
                                            phiComp, a_dxVect[gradDir]);
                        }
                      // Solve for high corners
                      cornerBox = combinedBox;
                      cornerBox &= boxHi1[comp];
                      if (!cornerBox.isEmpty())
                        {
                          tanSide = Side::Hi;
                          Box biasedBox = combinedBox;
                          biasedBox &= boxHi2[comp];
                          computeCornerGrad(a_GradPhicellAvgFab,
                                            a_PhicellAvgFab,
                                            cornerBox, biasedBox,
                                            gradDir, tanSide,
                                            phiComp, a_dxVect[gradDir]);
                        }
                    }
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Get the computational coordinates on cells
/** \param[in]  a_box   Box on which to compute coordinates
 *  \param[out] a_XiFab Cartesian computational space coordinates
 *  \param[in]  a_blockCoordSys
 *                      Coordinate system for the block
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getCellCompCoordinates(const Box&           a_box,
                               FArrayBox&           a_XiFab,
                               const BlockCoordSys& a_blockCoordSys) const
{
  const RealVect& dx = a_blockCoordSys.dx();
  // Computational space
  FORT_SETCELLCENTERSVEC(CHF_FRA(a_XiFab),
                         CHF_CONST_REALVECT(dx),
                         CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Get the computational and physical space coordinates on cells
/** \param[in]  a_box   Box on which to compute coordinates
 *  \param[out] a_XiFab Cartesian computational space coordinates
 *  \param[out] a_XFab  Physical coordinates
 *  \param[in]  a_blockCoordSys
 *                      Coordinate system for the block
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getCellCoordinates(const Box&           a_box,
                           FArrayBox&           a_XiFab,
                           FArrayBox&           a_XFab,
                           const BlockCoordSys& a_blockCoordSys) const
{
  CH_TIME("CNSIBC::getCellCoordinates");
  // Computational space
  getCellCompCoordinates(a_box, a_XiFab, a_blockCoordSys);

  // Physical space
  a_blockCoordSys.realCoord(a_XFab, a_XiFab, a_box);
}

/*--------------------------------------------------------------------*/
//  Get the computational and physical space coordinates on faces
/** \param[in]  a_box   Face box on which to compute coordinates
 *  \param[out] a_XiFab Cartesian computational space coordinates
 *  \param[out] a_XFab  Physical coordinates
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_blockCoordSys
 *                      Coordinate system for the block
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getFaceCompCoordinates(const Box&           a_box,
                               FArrayBox&           a_XiFab,
                               const int            a_dir,
                               const BlockCoordSys& a_blockCoordSys) const
{
  CH_TIME("CNSIBC::getFaceCompCoordinates");
  const RealVect& dx = a_blockCoordSys.dx();

  // Computational space
  FORT_SETFACECENTERSVEC(CHF_FRA(a_XiFab),
                         CHF_CONST_REALVECT(dx),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Get the computational and physical space coordinates on faces
/** \param[in]  a_box   Face box on which to compute coordinates
 *  \param[out] a_XiFab Cartesian computational space coordinates
 *  \param[out] a_XFab  Physical coordinates
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_blockCoordSys
 *                      Coordinate system for the block
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getFaceCoordinates(const Box&           a_box,
                           FArrayBox&           a_XiFab,
                           FArrayBox&           a_XFab,
                           const int            a_dir,
                           const BlockCoordSys& a_blockCoordSys) const
{
  CH_TIME("CNSIBC::getFaceCoordinates");
  const RealVect& dx = a_blockCoordSys.dx();

  // Computational space
  FORT_SETFACECENTERSVEC(CHF_FRA(a_XiFab),
                         CHF_CONST_REALVECT(dx),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_box));

  // Physical space
  a_blockCoordSys.realCoord(a_XFab, a_XiFab, a_box);
}

/*--------------------------------------------------------------------*/
//  Get the computational and physical space coordinates at nodes
/** \param[in]  a_box   Box on which to compute coordinates
 *  \param[out] a_XiFab Cartesian computational space coordinates
 *  \param[out] a_XFab  Physical coordinates
 *  \param[in]  a_blockCoordSys
 *                      Coordinate system for the block
 *//*-----------------------------------------------------------------*/

void
CNSIBC::getNodeCoordinates(const Box&           a_box,
                           FArrayBox&           a_XiFab,
                           FArrayBox&           a_XFab,
                           const BlockCoordSys& a_blockCoordSys) const
{
  CH_TIME("CNSIBC::getNodeCoordinates");
  CH_assert(a_XiFab.box().contains(a_box));
  const RealVect& dx = a_blockCoordSys.dx();

  // Computational space
  FORT_SETCORNERSVEC(CHF_FRA(a_XiFab),
                     CHF_CONST_REALVECT(dx),
                     CHF_BOX(a_box));

  // Physical space
  a_blockCoordSys.realCoord(a_XFab, a_XiFab, a_box);
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the exterior (or farfield) primitive state at flow BC
/** This routine is used to set both Dirichlet and the exterior state
 *  for characteristic boundaries.  This particular implementation
 *  sets the state to that defined in the global parameters.  The
 *  velocity is set to (g_speed, 0, 0).  In general, one should
 *  override this routine to specify the exact conditions you need.
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
 *  \param[in]  a_bcIdx Index to be MB block, direction, and side
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_bcInfo) const
{
  const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
  for (int c = 0, c_end = CRDparam::g_CRDPhysics->numPrimitive(); c != c_end;
       ++c)
    {
      a_Wface.setVal(state(c), a_boundaryFaceBox, c);
    }
}

/*--------------------------------------------------------------------*/
//  Set the CNSCBC boundary values and ghost cells
/** \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.
 *  \param[out] a_Wcell Primitive state with exterior ghost cells filled
 *  \param[in]  a_ghostCellBox1, a_ghostCellBox2
 *                      First and second layers of exterior ghost cells
 *  \param[in]  a_bcIdx Index to be MB block, direction, and side
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBC::applyCNSCBC(FArrayBox&                    a_WavgFace,
                    const Box&                    a_boundaryFaceBox,
                    FArrayBox&                    a_WavgCell,
                    const FArrayBox&              a_unitNormalBasisFab,
                    const BoundaryIndex&          a_bcIdx,
                    const Box&                    a_disjointBox,
                    LevelGridMetrics&             a_gridMetrics,
                    const Real                    a_time,
                    const int                     a_level,
                    const bool                    a_edgeBool,
                    const BCInfo&                 a_domT) const
{
  // Set all the component variables and intervals
  const Side::LoHiSide side = a_bcIdx.m_side;
  const int dir = a_bcIdx.m_dir;
  const int lohiSign = Side::sign(side);
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  // Get the block index and boundary conditions
  a_bcIdx.assertValid();
  // Shift face to match first interior cell
  a_WavgFace.shiftHalf(dir, -lohiSign);
  // Coefficients
  RealVect dxVect = a_gridMetrics.dxVect();
  const Real dxN = dxVect[dir];
  // Get a disjoint interior box
  Box interiorBox = adjCellBox(a_boundaryFaceBox, dir, Side::flip(side), 1);
  Box ghostCellBox1 = interiorBox;
  ghostCellBox1.shift(dir, lohiSign);
  Box ghostCellBox2 = ghostCellBox1;
  ghostCellBox2.shift(dir, lohiSign);
  Box ghostCells(ghostCellBox1);
  ghostCells.minBox(ghostCellBox2);
  const int startIndx = 0;
  const int lastIndx = numWcomp - 1;
  // Normal derivatives
  FABSTACKTEMP(d1FabN, interiorBox, numWcomp);
  MD_ARRAY_RESTRICT(arrD1, d1FabN);
  MD_ARRAY_RESTRICT(arrWCell, a_WavgCell);
  const int MD_ID(o, dir);
  for (int comp = startIndx; comp != lastIndx + 1; ++comp)
    {
      MD_BOXLOOP(interiorBox, i)
        {
          Real Wjp3 = arrWCell[MD_OFFSETIX(i,-,2*lohiSign*o,comp)];
          Real Wjp2 = arrWCell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
          Real Wjp1 = arrWCell[MD_IX(i,comp)];
          Real dphidxn = lohiSign*(2.*Wjp1 - 3.*Wjp2 + Wjp3)/dxN;
          arrD1[MD_IX(i, comp)] = dphidxn;
        }
    }
  FABSTACKTEMP(d1FabT, interiorBox, (SpaceDim - 1)*numWcomp);
  d1FabT.setVal(0.);
  // d1FabT components are organized as comp contiguous,
  // For example, if we are on an y boundary in 3D the comps are
  // 0: d u_0/dy
  // 1: d u_1/dy
  // 2: d u_2/dy
  // ...
  // numWcomp+1: d u_0/dz
  // numWcomp+2: d u_1/dz etc
  //int stride = 0;
  // FIXME: Permanently refuses to use transverse operations for CNSCBC
  // for (int tandir = 0; tandir != SpaceDim; ++tandir)
  //   {
  //     if (tandir == dir)
  //       {
  //         continue;
  //       }
  //     Box solveBox(interiorBox);
  //     const Real dxT = dxVect[tandir];
  //     const int fST = stride*numWcomp;
  //     FORT_CNSCBCTRANSGRAD(CHF_FRA(d1FabT),
  //                          CHF_CONST_FRA(a_WavgCell),
  //                          CHF_CONST_INT(tandir),
  //                          CHF_CONST_INT(fST),
  //                          CHF_CONST_INT(startIndx),
  //                          CHF_CONST_INT(lastIndx),
  //                          CHF_CONST_REAL(dxT),
  //                          CHF_BOX(interiorBox));
  //     ++stride;
  //   }
  FABSTACKTEMP(gammaFab, interiorBox, 1);
  CRDparam::g_CRDPhysics->calcGamma(interiorBox, gammaFab, a_WavgCell);
  FABSTACKTEMP(dqdx, interiorBox, numWcomp);
  // FIXME: Not using transverse components at inflow boundary
  Real betain = 0.;
  // FIXME: Use boundary specific etaMax, etaCN sigma and beta
  if (CRDparam::DomainBCTypeCNSCBCInflow & a_domT.m_type)
    {
      CNSCBCInflow(dqdx, a_WavgCell, a_WavgFace, d1FabN, d1FabT, gammaFab,
                   interiorBox, m_etaMax, m_etaCN, betain, dir, side);
    }
  else if (CRDparam::DomainBCTypeCNSCBCOutflow & a_domT.m_type)
    {
      CNSCBCOutflow(dqdx, a_WavgCell, a_WavgFace, d1FabN, d1FabT, gammaFab,
                    interiorBox, m_CBCsigma, m_CBCbeta, dir, side);
    }
  Vector<int> compVect(numWcomp);
  for (int comp = 0; comp != numWcomp; ++comp)
    {
      compVect[comp] = comp;
    }
  MD_ARRAY_RESTRICT(arrWFace, a_WavgFace);
  MD_ARRAY_RESTRICT(arrDQ, dqdx);
  for (int comp = 0; comp != compVect.size(); ++comp)
    {
      const int curComp = compVect[comp];
      MD_BOXLOOP(interiorBox, i)
        {
          Real dqdxV = arrDQ[MD_IX(i, curComp)];
          Real Wjp  = arrWCell[MD_OFFSETIX(i,-,lohiSign*o,curComp)];
          Real Wj  = arrWCell[MD_IX(i,curComp)];
          arrWFace[MD_IX(i, curComp)] =
            (-Wjp + 9.*Wj + lohiSign*3.*dxN*dqdxV)/8.;
          Real Wjm = Wj + lohiSign*dqdxV*dxN;
          arrWCell[MD_OFFSETIX(i,+,lohiSign*o,curComp)] = Wjm;
          arrWCell[MD_OFFSETIX(i,+,2*lohiSign*o,curComp)] =
            -Wjp + 27.*(Wj + Wjm) + 24.*dqdxV*dxN;
        }
    }
  a_WavgFace.shiftHalf(dir, lohiSign);
  CRDparam::g_CRDPhysics->temperature(a_WavgCell, ghostCells);
  CRDparam::g_CRDPhysics->temperature(a_WavgFace, a_boundaryFaceBox);
}

/*--------------------------------------------------------------------*/
//  Set wall velocity on entire box
/** This function is not intended to require interior information
 *  \param[out] a_wallVelocity 
 *                      Velocity (face-averaged) of the wall boundary
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_disjointBox
 *                      Disjoint box of the current block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_dir   Normal direction for current wall boundary
 *  \param[in]  a_side  Side of the block the boundary is on
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setBoxWallVelocity(FArrayBox&           a_wallVelocity,
                           const FArrayBox&     a_unitNormalBasisFab,
                           const Box&           a_boundaryFaceBox,
                           const Box&           a_disjointBox,
                           LevelGridMetrics&    a_gridMetrics,
                           const ProblemDomain& a_domain,
                           const Real           a_time,
                           const int            a_level,
                           const int            a_dir,
                           const Side::LoHiSide a_side) const
{
  CH_TIME("CNSIBC::setBoxWallVelocity");
  // The boundary index and type
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir,
               a_side);
  const BCInfo& domBC = getDomainBC(bcIdx);

  // If this is not a wall, don't set the wall velocity
  if (!(CRDparam::DomainBCTypeAllWall & domBC.m_type))
    {
      return;
    }

  setBndryWallVelocity(a_wallVelocity,
                       a_unitNormalBasisFab,
                       a_boundaryFaceBox,
                       a_disjointBox,
                       a_gridMetrics,
                       a_domain,
                       a_time,
                       a_level,
                       a_dir,
                       a_side,
                       domBC);
}

/*--------------------------------------------------------------------*/
//  Set wall velocity on specific boundary -- user can specialize this
/** This function is not intended to require interior information
 *  \param[out] a_wallVelocity 
 *                      Velocity (face-averaged) of the wall boundary
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_disjointBox
 *                      Disjoint box of the current block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Current-block domain
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_dir   Normal direction for current wall boundary
 *  \param[in]  a_side  Side of the block the boundary is on
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setBndryWallVelocity(FArrayBox&           a_wallVelocity,
                             const FArrayBox&     a_unitNormalBasisFab,
                             const Box&           a_boundaryFaceBox,
                             const Box&           a_disjointBox,
                             LevelGridMetrics&    a_gridMetrics,
                             const ProblemDomain& a_domain,
                             const Real           a_time,
                             const int            a_level,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const BCInfo&        a_domT) const
{
  CH_TIME("CNSIBC::setBndryWallVelocity");
  // Wall boundary primitive state
  const CRDState& domBCstate = CRDState::get(a_domT.m_idxState);

  // Set the wall-velocity to the primitive state set on this wall boundary
  RealVect wallVelocity = domBCstate.velocity();
  // Now set the boundary-face values
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          a_wallVelocity[MD_IX(i, comp)] = wallVelocity[comp];
        }
    }
}

/*--------------------------------------------------------------------*/
//  Spatially average data if necessary
/** \param[out] a_avgData
 *                      Spatially-averaged data
 *  \param[in]  a_avgData
 *                      Data to spatially-average
 *  \param[in]  a_stage Current stage of the time-marching method
 *  \param[in]  a_t     Current time
 *//*-----------------------------------------------------------------*/

void
CNSIBC::spatiallyAverageData(LevelData<FArrayBox>& a_avgData,
                             const int             a_stage,
                             const Real            a_t) const
{
  CH_TIME("CNSIBC::spatiallyAverageData");
}

/*--------------------------------------------------------------------*/
//  Output pnt-values and domain-sums specific to simulation
/** \param[in] a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in] a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in] a_WcellPntFab
 *                      Cell-centered primitive state
 *  \param[in] a_unitNormalsFxb
 *                      Unit-normal basis for faces
 *  \param[in] a_gridMetrics
 *                      Grid metrics for current level
 *  \param[in] a_time   Current time on current level
 *  \param[in] a_level  Current level
 *//*-----------------------------------------------------------------*/

void
CNSIBC::inSituSumPntProcessing(
  const LevelData<FluxBox>&   a_faceAvgPlotData,
  const LevelData<FluxBox>&   a_WfaceAvgFxb,
  const LevelData<FArrayBox>& a_WcellAvgFab,
  const LevelData<FArrayBox>& a_WcellPntFab,
  const LayoutData<FluxBox>&  a_unitNormalsFxb,
  const LevelGridMetrics&     a_gridMetrics,
  const Real                  a_time,
  const int                   a_stage,
  const int                   a_level) const
{
  CH_TIME("CNSIBC::inSituSumPntProcessing");
}

/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In general, one should not need to override this routine except
 *  in special circumstances such as for mixed boundaries.
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
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setWallBCprimState(FArrayBox&                    a_Wface,
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
  // Sign to shift boxes out of the domain
  const int lohiSign = Side::sign(a_bcIdx.m_side);

  // Wall boundary primitive state
  const CRDState& domBCstate = CRDState::get(a_domT.m_idxState);

  // Need to specify whether or not this is a slip-wall
  int viscousSlip = 0;
  if (CRDparam::DomainBCTypeSlipWall & a_domT.m_type)
    {
      viscousSlip = 1;
    }
  // Set the wall-velocity to the primitive state set on this wall boundary
  //**FIXME: it would be really nice to wrap this up in computeWallPrimState
  //         rather than have it as another input
  RealVect wallVelocity = domBCstate.velocity();
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
  CRDparam::g_CRDPhysics->temperature(a_Wface, a_boundaryFaceBox);
}

/*--------------------------------------------------------------------*/
//  Fill in the profiles for the CNSCBC boundary condition.
//  Inflow - Must specify velocity, temperature, and cn. a_BCProfile
//           has numWcomp components
//  Outflow - Must specify pressure. a_BCProfile only has 1 component
/** \param[out] a_BCProfile
 *                      Profile of boundary values
 *  \param[in]  a_boundaryBox
 *                      Box of locations at the boundary
 *  \param[in]  a_Wcell Cell-averaged primitive variables
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current CNSCBC boundary condition
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setCNSCBCProfiles(FArrayBox&                    a_BCProfile,
                          const Box&                    a_boundaryBox,
                          const FArrayBox&              a_Wcell,
                          const FArrayBox&              a_unitNormalBasisFab,
                          const BoundaryIndex&          a_bcIdx,
                          const Box&                    a_disjointBox,
                          LevelGridMetrics&             a_gridMetrics,
                          const Real                    a_time,
                          const int                     a_level,
                          const BCInfo&                 a_domT) const
{
  CRD::msg << "CNSIBC::setCNSCBCProfiles is not defined!" << CRD::error;
}

/*--------------------------------------------------------------------*/
//  Set exterior face-state for relaxed characteristic BCs
/** \param[out] a_WfaceAvgExterior
 *                      Exterior face-averaged value for Riemann solve
 *  \param[in]  a_bndryCellFab
 *                      Cell-averaged values exterior/interior to bndry.
 *                      These values are at current and previous
 *                      time-points
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged values (including interior cells)
 *  \param[in]  a_WfaceAvgDirFab
 *                      Face-averaged values (including interior faces)
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_bcInfo
 *                      Current boundary type
 *  \param[in]  a_disjointBox
 *                      Comment
 *  \param[in]  a_totalFaceBox
 *                      Comment
 *  \param[in]  a_gridMetrics
 *                      Grid-metrics on the current level
 *  \param[in]  a_dir   Face-normal direction
 *  \param[in]  a_side  Side of box on which bndry resides (low/high)
 *  \param[in]  a_time  Current time
 *  \param[in]  a_prevDt
 *                      Time-step size at the previous time-step
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBC::setRelaxedCBCPrimState(
  FArrayBox&              a_WfaceAvgExterior,
  FArrayBox&              a_bndryCellFab,
  FArrayBox&              a_WcellAvgFab,
  const FArrayBox&        a_WfaceAvgDirFab,
  const FArrayBox&        a_unitNormalBasisFab,
  const BoundaryIndex&    a_bcIdx,
  const BCInfo&           a_bcInfo,
  const Box&              a_disjointBox,
  const Box&              a_totalFaceBox,
  const LevelGridMetrics& a_gridMetrics,
  const int               a_dir,
  const Side::LoHiSide    a_side,
  const Real              a_time,
  const Real              a_prevDt,
  const int               a_level) const
{
  CH_TIME("CNSIBC::setRelaxedCBCPrimState");
  // Constants
  const int cRho       = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel       = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cPres      = CRDparam::g_CRDPhysics->pressureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const Real beta = a_bcInfo.m_relaxCBCWaveParam; // 1/0 = non/reflecting
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  D_TERM(const int idnorm = a_dir;,
         const int idtan1 = (a_dir + 1) % SpaceDim;,
         const int idtan2 = (a_dir + 2) % SpaceDim;);
  //**NOTE: While a_bndryCellFab doesn't change throughout a time-step,
  //        the update is performed using the user-defined exterior state.
  //        Since this can change through a time-step, it does here as well.
  // Set up some boxes for determining interior and exterior cells
  Box interiorCellBox = a_totalFaceBox;
  Box exteriorCellBox = a_totalFaceBox;
  Box totalCellBox = a_totalFaceBox;
  totalCellBox.grow(a_dir, 1); // Grow both the low and high sides
  totalCellBox.enclosedCells(a_dir); // Cell-box covering in/exterior cells
  if (a_side == Side::Lo)
    {
      interiorCellBox.growHi(a_dir, 1);
      interiorCellBox.enclosedCells(a_dir);
      exteriorCellBox.growLo(a_dir, 1);
      exteriorCellBox.enclosedCells(a_dir);
    }
  else
    {
      interiorCellBox.growLo(a_dir, 1);
      interiorCellBox.enclosedCells(a_dir);
      exteriorCellBox.growHi(a_dir, 1);
      exteriorCellBox.enclosedCells(a_dir);
    }
  int numPrim = CRDparam::g_CRDPhysics->numPrimitive();
  if (CNSIBC::s_firstRKStage)
    {
      // Fill interior bndryCellFab cells with a_WcellAvgFab
      for (int comp = 0; comp != numPrim; ++comp)
        {
          MD_BOXLOOP(interiorCellBox, i)
            {
              a_bndryCellFab[MD_IX(i, comp)] = a_WcellAvgFab[MD_IX(i, comp)];
              CH_assert(a_bndryCellFab[MD_IX(i, comp)] < hiTol);
              CH_assert(a_bndryCellFab[MD_IX(i, comp)] > loTol);
            }
        }
      // Fill exterior bndryCellFab cells with a_WfaceAvgExterior (1st-order)
      for (int comp = 0; comp != numPrim; ++comp)
        {
          MD_BOXLOOP(exteriorCellBox, i)
            {
              a_bndryCellFab[MD_IX(i, comp)] =
                a_WfaceAvgExterior[MD_IX(i, comp)];
              CH_assert(a_bndryCellFab[MD_IX(i, comp)] < hiTol);
              CH_assert(a_bndryCellFab[MD_IX(i, comp)] > loTol);
            }
        }
    }
  else if (CNSIBC::s_lastRKStage)
    {
      // Copy currBndryState to prevBndryState (interior and exterior cells)
      for (int comp = 0; comp != numPrim; ++comp)
        {
          int pastIdx = comp + numPrim; // Comp is the index of current time
          MD_BOXLOOP(totalCellBox, i)
            {
              a_bndryCellFab[MD_IX(i, pastIdx)] =
                a_bndryCellFab[MD_IX(i, comp)];
              CH_assert(a_bndryCellFab[MD_IX(i, pastIdx)] < hiTol);
              CH_assert(a_bndryCellFab[MD_IX(i, pastIdx)] > loTol);
            }
        }
    }

  // Test for first stage of first time-step
  if (!a_prevDt) { return; } // nothing else is required in this function

  // Inverse of the previous time-step size
  CH_assert(a_prevDt > 0.);
  const Real invDt = 1./a_prevDt;

  // Boundary-side indicator
  const int lohiSign = Side::sign(a_side);

  // Alias the velocity in a_bndryCellFab so we can transform it
  Interval velCurrIntv = CRDparam::g_CRDPhysics->velocityInterval();
  Interval velPastIntv(velCurrIntv.begin()+numPrim, velCurrIntv.end()+numPrim);
  FArrayBox velPastFab(velPastIntv, a_bndryCellFab);
  // First we shift velPastFab exterior to the face so unitNormalBasis works
  velPastFab.shiftHalf(a_dir, -lohiSign);
  // Transform the exterior velocity to be normal to the faces
  PatchMappedFunc::forwardTransform(
    velPastFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift it back so interior is on face
  velPastFab.shift(a_dir, lohiSign);
  // Transform the interior velocity to be normal to the faces
  PatchMappedFunc::forwardTransform(
    velPastFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift it back so everything is cell-based again
  velPastFab.shiftHalf(a_dir, -lohiSign);

  // Alias the velocity in a_WfaceAvgExterior
  FArrayBox WfaceAvgExteriorVelFab(velCurrIntv, a_WfaceAvgExterior);
  // Alias the velocity in a_WcellAvgFab
  FArrayBox WcellAvgInteriorVelFab(velCurrIntv, a_WcellAvgFab);
  // Shift WfaceAvgExteriorVelFab to the boundary face
  WfaceAvgExteriorVelFab.shiftHalf(a_dir, -lohiSign);
  // Transform the velocities to be normal to the faces
  PatchMappedFunc::forwardTransform(
    WfaceAvgExteriorVelFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift WfaceAvgExteriorVelFab back to the exterior cell location
  WfaceAvgExteriorVelFab.shiftHalf(a_dir, lohiSign);
  // Shift WcellAvgInteriorVelFab to the boundary face
  WcellAvgInteriorVelFab.shiftHalf(a_dir, lohiSign);
  // Transform the velocity to be normal to the faces
  PatchMappedFunc::forwardTransform(
    WcellAvgInteriorVelFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift WcellAvgInteriorVelFab back to the interior cell
  WcellAvgInteriorVelFab.shiftHalf(a_dir, -lohiSign);

  // (a) Compute the wave-amplitudes for the interior state at the
  //     current time using current-time-data (time-step n) and
  //     past-time-data (time-step n-1)

  // Compute the derivative and assign the reference state
  FABSTACKTEMP(derivState1, interiorCellBox, numPrim);
  FABSTACKTEMP(refW1, interiorCellBox, numPrim);
  for (int comp = 0; comp != numPrim; ++comp)
    {
      int pastIdx = comp + numPrim; // Comp is the index of current time
      MD_BOXLOOP(interiorCellBox, i)
        {
          const Real currVal = a_WcellAvgFab[MD_IX(i, comp)];
          const Real pastVal = a_bndryCellFab[MD_IX(i, pastIdx)];
          derivState1[MD_IX(i, comp)] = invDt*(currVal - pastVal);
          // The reference state is open to experimentation
          refW1[MD_IX(i, comp)] = pastVal; // Upwind values for now

          CH_assert(currVal < hiTol);
          CH_assert(currVal > loTol);
          CH_assert(pastVal < hiTol);
          CH_assert(pastVal > loTol);
        }
    }
  // Compute the wave-amplitude
  FABSTACKTEMP(waveAmpInt1, interiorCellBox, numPrim);
  computeCharWaveAmp(waveAmpInt1, refW1, derivState1, a_dir, interiorCellBox);

  // (b) Compute the wave-amplitudes for the exterior state at the
  //     current time using current-time-data (time-step n) and
  //     past-time-data (time-step n-1)

  // Compute the derivative and assign the reference state
  // The reference state is open to experimentation
  FABSTACKTEMP(derivState2, exteriorCellBox, numPrim);
  FABSTACKTEMP(refW2, exteriorCellBox, numPrim);
  for (int comp = 0; comp != numPrim; ++comp)
    {
      int pastIdx = comp + numPrim; // Comp is the index of current time
      MD_BOXLOOP(exteriorCellBox, i)
        {
          const Real currVal = a_WfaceAvgExterior[MD_IX(i, comp)];
          const Real pastVal = a_bndryCellFab[MD_IX(i, pastIdx)];
          derivState2[MD_IX(i, comp)] = invDt*(currVal - pastVal);

          refW2[MD_IX(i, comp)] = pastVal; // Upwind values for now

          CH_assert(currVal < hiTol);
          CH_assert(currVal > loTol);
          CH_assert(pastVal < hiTol);
          CH_assert(pastVal > loTol);
        }
    }
  // Compute the wave-amplitude
  FABSTACKTEMP(waveAmpExt1, exteriorCellBox, numPrim);
  computeCharWaveAmp(waveAmpExt1, refW2, derivState2, a_dir, exteriorCellBox);

  // Alias the primitive state at the previous time
  Interval pastIntv(numPrim, (2*numPrim) - 1);
  FArrayBox WcellAvgPreviousTime(pastIntv, a_bndryCellFab);
  // Compute reference speed of sound
  FABSTACKTEMP(interiorC, interiorCellBox, 1);
  CRDparam::g_CRDPhysics->soundSpeed(
    interiorC, refW1, interiorCellBox);
  FABSTACKTEMP(exteriorC, exteriorCellBox, 1);
  CRDparam::g_CRDPhysics->soundSpeed(
    exteriorC, refW2, exteriorCellBox);

  // Shift interior cell-averaged state to the exterior
  a_WcellAvgFab.shift(a_dir, lohiSign);
  // Shift interior speed-of-sound to the exterior
  interiorC.shift(a_dir, lohiSign);
  // Shift refW1 to the exterior
  refW1.shift(a_dir, lohiSign);
  // Shift interior wave amplitudes to the exterior
  waveAmpInt1.shift(a_dir, lohiSign);
  // Now update the exterior state using the wave-amplitudes
  MD_BOXLOOP(exteriorCellBox, i)
    {
      // Check the direction of interior velocity
      RealVect intVel(D_DECL(a_WcellAvgFab[MD_IX(i, cVel)],
                             a_WcellAvgFab[MD_IX(i, cVel + 1)],
                             a_WcellAvgFab[MD_IX(i, cVel + 2)]));
      // Set the default ref state to be the previous exterior state
      Real refRho = refW2[MD_IX(i, cRho)];
      RealVect refVel(D_DECL(refW2[MD_IX(i, cVel)],
                             refW2[MD_IX(i, cVel + 1)],
                             refW2[MD_IX(i, cVel + 2)]));
      Real refPress = refW2[MD_IX(i, cPres)];
      std::vector<Real> refSpecies(numSpecies, 0.);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              refSpecies[spec] = refW2[MD_IX(i, wSpecIndex)];
            }
        }
      Real refC = exteriorC[MD_IX(i, 0)];
      // If an outflow set the ref state to be the previous interior state
      if (intVel[idnorm]*lohiSign >= 0.) // Outflow
        {
          refRho = refW1[MD_IX(i, cRho)];
          refVel = RealVect(D_DECL(refW1[MD_IX(i, cVel)],
                                   refW1[MD_IX(i, cVel + 1)],
                                   refW1[MD_IX(i, cVel + 2)]));
          refPress = refW1[MD_IX(i, cPres)];
          if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
            {
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  int wSpecIndx = spec + cSpecies;
                  refSpecies[spec] = refW1[MD_IX(i, wSpecIndx)];
                }
            }
          refC = interiorC[MD_IX(i, 0)];
        }
      // Set the default wave-amplitudes to come from the exterior
      Real LWave = waveAmpExt1[MD_IX(i, cRho)];
      Real RWave = waveAmpExt1[MD_IX(i, cPres)];
      RealVect EntWaves(D_DECL(waveAmpExt1[MD_IX(i, cVel + idnorm)],
                               waveAmpExt1[MD_IX(i, cVel + idtan1)],
                               waveAmpExt1[MD_IX(i, cVel + idtan2)]));
      std::vector<Real> speciesWaves(numSpecies, 0.);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              speciesWaves[spec] = waveAmpExt1[MD_IX(i, wSpecIndex)];
            }
        }
      // Check if the interior bndry-normal velocity is supersonic
      Real intBndryNormalMach =
        (std::abs(intVel[idnorm]))/(interiorC[MD_IX(i, 0)]);
      if (intBndryNormalMach >= 1.)
        {
          // If it is and it's outflow, set everything to interior waves
          if (intVel[idnorm]*lohiSign > 0.) // Outflow
            {
              LWave = waveAmpInt1[MD_IX(i, cRho)];
              RWave = waveAmpInt1[MD_IX(i, cPres)];
              EntWaves = RealVect(D_DECL(waveAmpInt1[MD_IX(i, cVel + idnorm)],
                                         waveAmpInt1[MD_IX(i, cVel + idtan1)],
                                         waveAmpInt1[MD_IX(i, cVel + idtan2)]));
              if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
                {
                  for (int spec = 0; spec != numSpecies; ++spec)
                    {
                      int wSpecIndex = spec + cSpecies;
                      speciesWaves[spec] = waveAmpInt1[MD_IX(i, wSpecIndex)];
                    }
                }
            }
        }
      else // Subsonic
        {
          // If it's the low side, the left moving wave must exit
          if (a_side == Side::Lo)
            {
              LWave = waveAmpInt1[MD_IX(i, cRho)];
              // Also, if it's an outflow, use the interior entropy waves
              if (intVel[idnorm]*lohiSign >= 0.) // Outflow
                {
                  EntWaves = RealVect(
                    D_DECL(waveAmpInt1[MD_IX(i, cVel + idnorm)],
                           waveAmpInt1[MD_IX(i, cVel + idtan1)],
                           waveAmpInt1[MD_IX(i, cVel + idtan2)]));
                  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
                    {
                      for (int spec = 0; spec != numSpecies; ++spec)
                        {
                          int wSpecIndex = spec + cSpecies;
                          speciesWaves[spec] =
                            waveAmpInt1[MD_IX(i, wSpecIndex)];
                        }
                    }
                }
            }
          else // If it's the high side, the right moving wave must exit
            {
              RWave = waveAmpInt1[MD_IX(i, cPres)];
              // Also, if it's an outflow, use the interior entropy waves
              if (intVel[idnorm]*lohiSign >= 0.) // Outflow
                {
                  EntWaves = RealVect(
                    D_DECL(waveAmpInt1[MD_IX(i, cVel + idnorm)],
                           waveAmpInt1[MD_IX(i, cVel + idtan1)],
                           waveAmpInt1[MD_IX(i, cVel + idtan2)]));
                  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
                    {
                      for (int spec = 0; spec != numSpecies; ++spec)
                        {
                          int wSpecIndex = spec + cSpecies;
                          speciesWaves[spec] =
                            waveAmpInt1[MD_IX(i, wSpecIndex)];
                        }
                    }
                }
            }
        }
      // Relax the non-reflecting state towards the imposed state if necessary
      refRho = beta*refRho + (1. - beta)*refW2[MD_IX(i, cRho)];
      D_TERM(
        refVel[0] = beta*refVel[0] + (1. - beta)*refW2[MD_IX(i, cVel)];,
        refVel[1] = beta*refVel[1] + (1. - beta)*refW2[MD_IX(i, cVel + 1)];,
        refVel[2] = beta*refVel[2] + (1. - beta)*refW2[MD_IX(i, cVel + 2)];);
      refPress = beta*refPress + (1. - beta)*refW2[MD_IX(i, cPres)];
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              refSpecies[spec] = beta*refSpecies[spec]
                + (1. - beta)*refW2[MD_IX(i, wSpecIndex)];
            }
        }
      refC = beta*refC + (1. - beta)*exteriorC[MD_IX(i, 0)];
      // Relax the non-reflecting waves towards the imposed waves if necessary
      LWave = beta*LWave + (1. - beta)*waveAmpExt1[MD_IX(i, cRho)];
      RWave = beta*RWave + (1. - beta)*waveAmpExt1[MD_IX(i, cPres)];
      D_TERM(
        EntWaves[0] = beta*EntWaves[0]
        + (1. - beta)*waveAmpExt1[MD_IX(i, cVel + idnorm)];,
        EntWaves[1] = beta*EntWaves[1]
        + (1. - beta)*waveAmpExt1[MD_IX(i, cVel + idtan1)];,
        EntWaves[2] = beta*EntWaves[2]
        + (1. - beta)*waveAmpExt1[MD_IX(i, cVel + idtan2)];);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              speciesWaves[spec] = beta*speciesWaves[spec]
                + (1. - beta)*waveAmpExt1[MD_IX(i, wSpecIndex)];
            }
        }

      // Finally, get to performing the update
      a_WfaceAvgExterior[MD_IX(i, cRho)] = refRho
        - a_prevDt*((1./(refC*refC))*(EntWaves[idnorm] + 0.5*(LWave + RWave)));
      a_WfaceAvgExterior[MD_IX(i, cPres)] = refPress
        - a_prevDt*0.5*(LWave + RWave);
      D_TERM(
        a_WfaceAvgExterior[MD_IX(i, cVel + idnorm)] = refVel[idnorm]
        - a_prevDt*(1./(2.*refRho*refC))*(RWave - LWave);
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idnorm)] < hiTol);
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idnorm)] > loTol);,
        a_WfaceAvgExterior[MD_IX(i, cVel + idtan1)] = refVel[idtan1]
        - a_prevDt*EntWaves[idtan1];
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idtan1)] < hiTol);
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idtan1)] > loTol);,
        a_WfaceAvgExterior[MD_IX(i, cVel + idtan2)] = refVel[idtan2]
        - a_prevDt*EntWaves[idtan2];
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idtan2)] < hiTol);
        CH_assert(a_WfaceAvgExterior[MD_IX(i, cVel + idtan2)] > loTol););
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              a_WfaceAvgExterior[MD_IX(i, wSpecIndex)] =
                refSpecies[spec] - a_prevDt*speciesWaves[spec];
              CH_assert(a_WfaceAvgExterior[MD_IX(i, wSpecIndex)] < hiTol);
              CH_assert(a_WfaceAvgExterior[MD_IX(i, wSpecIndex)] > loTol);
            }
        }
      CH_assert(a_WfaceAvgExterior[MD_IX(i, cRho)] < hiTol);
      CH_assert(a_WfaceAvgExterior[MD_IX(i, cRho)] > loTol);
      CH_assert(a_WfaceAvgExterior[MD_IX(i, cPres)] < hiTol);
      CH_assert(a_WfaceAvgExterior[MD_IX(i, cPres)] > loTol);
    }
  // Now set temperature -- maybe unnecessary, but it probably doesn't hurt
  CRDparam::g_CRDPhysics->temperature(a_WfaceAvgExterior, exteriorCellBox);

  // Shift interior cell-averaged state back to the interior
  a_WcellAvgFab.shift(a_dir, -lohiSign);

  // First we shift velPastFab exterior to the face so unitNormalBasis works
  velPastFab.shiftHalf(a_dir, -lohiSign);
  // Transform the velocities back to the mapped space
  PatchMappedFunc::reverseTransform(
    velPastFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift it back so interior is on face
  velPastFab.shift(a_dir, lohiSign);
  // Transform the interior velocity to be normal to the faces
  PatchMappedFunc::reverseTransform(
    velPastFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift it back so everything is cell-based again
  velPastFab.shiftHalf(a_dir, -lohiSign);
  // Shift WfaceAvgExteriorVelFab to the boundary face
  WfaceAvgExteriorVelFab.shiftHalf(a_dir, -lohiSign);
  // Transform the velocity
  PatchMappedFunc::reverseTransform(
    WfaceAvgExteriorVelFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift WfaceAvgExteriorVelFab back to the exterior
  WfaceAvgExteriorVelFab.shiftHalf(a_dir, lohiSign);
  // Shift WcellAvgInterior to the boundary face
  WcellAvgInteriorVelFab.shiftHalf(a_dir, lohiSign);
  // Transform the velocity
  PatchMappedFunc::reverseTransform(
    WcellAvgInteriorVelFab, a_unitNormalBasisFab, a_totalFaceBox);
  // Shift WcellAvgInteriorVelFab back to the interior
  WcellAvgInteriorVelFab.shiftHalf(a_dir, -lohiSign);
}

/*--------------------------------------------------------------------*/
//  Compute characteristic wave amplitudes
/** \param[out] a_waveAmp
 *                      Characteristic wave amplitudes
 *  \param[in]  a_refState
 *                      Cell-averaged reference state values
 *  \param[in]  a_derivState
 *                      Primitive-state derivatives (space or time)
 *  \param[in]  a_dir   Face-normal direction for boundary faces
 *  \param[in]  a_cellBox
 *                      Cells on which to compute the wave-amplitudes
 *//*-----------------------------------------------------------------*/

void
CNSIBC::computeCharWaveAmp(FArrayBox&       a_waveAmp,
                           const FArrayBox& a_refState,
                           const FArrayBox& a_derivState,
                           const int        a_dir,
                           const Box&       a_cellBox) const
{
  const int cRho       = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel       = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cPres      = CRDparam::g_CRDPhysics->pressureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  D_TERM(const int idnorm = a_dir;,
         const int idtan1 = (a_dir + 1) % SpaceDim;,
         const int idtan2 = (a_dir + 2) % SpaceDim;);
  // Compute reference speed of sound
  FABSTACKTEMP(speedOfSound, a_cellBox, 1);
  CRDparam::g_CRDPhysics->soundSpeed(speedOfSound, a_refState, a_cellBox);
  MD_BOXLOOP(a_cellBox, i)
    {
      // Derivatives
      Real rhoDeriv = a_derivState[MD_IX(i, cRho)];
      RealVect velDeriv(D_DECL(a_derivState[MD_IX(i, cVel)],
                               a_derivState[MD_IX(i, cVel + 1)],
                               a_derivState[MD_IX(i, cVel + 2)]));
      Real pressDeriv = a_derivState[MD_IX(i, cPres)];
      std::vector<Real> specDeriv(numSpecies, 0.);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              specDeriv[spec] = a_derivState[MD_IX(i, wSpecIndex)];
            }
        }

      // Reference state
      Real rho = a_refState[MD_IX(i, cRho)];
      RealVect vel(D_DECL(a_refState[MD_IX(i, cVel)],
                          a_refState[MD_IX(i, cVel + 1)],
                          a_refState[MD_IX(i, cVel + 2)]));
      std::vector<Real> species(numSpecies, 0.);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              species[spec] = a_refState[MD_IX(i, wSpecIndex)];
            }
        }
      // Reference speed of sound
      Real refC = speedOfSound[MD_IX(i, 0)];

      // Wave-amplitudes (Lodato et al. 2008 JCP 227 -- inverse of eqn. A.20)
      // Note: For now, assume there are no bndry-parallel waves in the domain
      //       interior. Assuming otherwise, the system lacks sufficient eqns
      //       to determine the waves in multiple dirs. One logical step is
      //       implicitly linking cells in the bndry-parallel dir and solving a
      //       system of eqns with reasonable constraints. This should improve
      //       non-bndry-normal flow.
      // Note: Another improvement-concept is to use previous fluxes instead of
      //       cell-avgs. This would provide a sense of where information is
      //       flowing. It could eliminate the need for these wave-amplitudes.
      // Note: In fact, considering 1D contributions from only bndry-normal
      //       faces, we isolate bndry-tngnt waves (and can exclude them).
      // Note: Considering just diffusion (parabolic), CBCs and MOC
      //       (method-of-characteristics) are useless. A possible work-around
      //       is to just try to determine how much conserved-quantity diffused
      //       into the exterior ghost-cell. This is hard to do well though.
      // Note: A logical method for improving the BCs is to ditch
      //       characteristics and specify how much conserved-quantity flows
      //       into and out of a face. Then, this along with the interior cell
      //       value (or bndry + interior face values) would provide a final,
      //       desired face value. Basically, it updates the bndry flux to get
      //       as close to the desired incoming/outgoing states as possible
      Real leftMovingWaveAmp = rho*refC*velDeriv[idnorm] - pressDeriv;
      Real rightMovingWaveAmp = -pressDeriv - rho*refC*velDeriv[idnorm];
      RealVect entropyWaveAmps(D_DECL(pressDeriv - refC*refC*rhoDeriv,
                                      -velDeriv[idtan1], -velDeriv[idtan2]));
      std::vector<Real> speciesWaveAmps(numSpecies, 0.);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              speciesWaveAmps[spec] = -specDeriv[spec];
            }
        }
      a_waveAmp[MD_IX(i, cRho)] = leftMovingWaveAmp;
      a_waveAmp[MD_IX(i, cPres)] = rightMovingWaveAmp;
      D_TERM(a_waveAmp[MD_IX(i, cVel + idnorm)] = entropyWaveAmps[idnorm];,
             a_waveAmp[MD_IX(i, cVel + idtan1)] = entropyWaveAmps[idtan1];,
             a_waveAmp[MD_IX(i, cVel + idtan2)] = entropyWaveAmps[idtan2];);
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              int wSpecIndex = spec + cSpecies;
              a_waveAmp[MD_IX(i, wSpecIndex)] = speciesWaveAmps[spec];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set the boundary values and fill the ghost cells
/** \param[out] a_NeumannComp
 *                      Components for applying Neumann conditions
 *  \param[out] a_interiorExtrap
 *                      Components to extrapolate directly from interior
 *  \param[out] a_faceValExtrap
 *                      Components to extrapolate based on boundary values
 *  \param[out] a_WgradFace
 *                      Face boundary gradients for Neumann conditions
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set
 *  \return     solveV  Determines which variable is solved at the boundary
 *                      0  - Solve for temperature from rho and p
 *                      1  - Solve for density from p and T
 *                      2  - Solve for pressure from rho and T
 *                      -1 - Use all specified variables for rho, p, and T
 *//*-----------------------------------------------------------------*/

int
CNSIBC::setExtraBCs(Vector<int>&                  a_NeumannComp,
                    Vector<int>&                  a_interiorExtrap,
                    Vector<int>&                  a_faceValExtrap,
                    FArrayBox&                    a_WgradFace,
                    const BoundaryIndex&          a_bcIdx,
                    const Box&                    a_disjointBox,
                    LevelGridMetrics&             a_gridMetrics,
                    const Real                    a_time,
                    const int                     a_level,
                    const BCInfo&                 a_domT) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const Interval specIntv = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  const int numWpComp = CRDparam::g_CRDPhysics->numNativePrimitive();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const Interval wCompIntv(0, numWcomp - 1);
  const Interval WpCompIntv(0, numWpComp - 1);

  // Default is to solve temperature at boundary and in ghost cells
  // Set to 1 to solve for density instead
  // Set to 2 to solve for pressure instead
  // Set to -1 to not overwrite the value prescribed
  int solveV = 0; // FIXME: Not actually used due to stability issues
  if (CRDparam::DomainBCTypeAllWall & a_domT.m_type)
    {
      // Assume we are not solving for any variables
      solveV = -1;
      a_interiorExtrap[0] = presIndx; // Extrapolate p from interior
      a_interiorExtrap.push_back(rhoIndx); // Extrapolate rho from interior
      if (CRDparam::DomainBCTypeSlipWall & a_domT.m_type)
        {
          a_faceValExtrap[0] = velIntv.begin() + a_bcIdx.m_dir;
          a_NeumannComp[0] = tempIndx; // Assume adiabatic slip wall
          for (int velDir = 0; velDir != SpaceDim; ++velDir)
            {
              if (velDir != a_bcIdx.m_dir)
                {
                  a_interiorExtrap.push_back(velIntv.begin() + velDir);
                }
            }
        }
      else
        {
          // Set temperature at wall for isothermal walls
          if (CRDparam::DomainBCTypeIsothermalWall & a_domT.m_type)
            {
              a_faceValExtrap[0] = tempIndx; // Use set temperature values
            }
          // Otherwise, assume wall is adiabatic
          else
            {
              a_NeumannComp[0] = tempIndx; // Zero grad temperature
            }
          // Use face values of 0
          addIntervaltoVec(a_faceValExtrap, velIntv);
        }
      if (numSpecies > 0)
        {
          addIntervaltoVec(a_NeumannComp, specIntv);
        }
    }
  else if (CRDparam::DomainBCTypeDirichlet & a_domT.m_type)
    {
      addIntervaltoVec(a_faceValExtrap, wCompIntv);
      solveV = -1;
    }
  else if (CRDparam::DomainBCTypeInflow & a_domT.m_type ||
           CRDparam::DomainBCTypeFarfield & a_domT.m_type ||
           CRDparam::DomainBCTypeRelaxedCBCIn & a_domT.m_type ||
           CRDparam::DomainBCTypeRelaxedCBCFar & a_domT.m_type )
    {
      a_NeumannComp[0] = presIndx; // Use zero gradient Neumann
      a_faceValExtrap[0] = tempIndx; // Use set temperature and velocity values
      addIntervaltoVec(a_faceValExtrap, velIntv);
      if (numSpecies > 0)
        {
          addIntervaltoVec(a_faceValExtrap, specIntv);
        }
      solveV = 1;
    }
  else if (CRDparam::DomainBCTypeOutflow & a_domT.m_type ||
           CRDparam::DomainBCTypeRelaxedCBCOut & a_domT.m_type)
    {
      // Use zero gradient Neumann for viscous BC
      a_NeumannComp[0] = presIndx;
      a_NeumannComp.push_back(rhoIndx);
      a_NeumannComp.push_back(tempIndx);
      addIntervaltoVec(a_NeumannComp, velIntv);
      solveV = 0; // Solve for temperature for calorically perfect
      if (numSpecies > 0)
        {
          addIntervaltoVec(a_NeumannComp, specIntv);
          solveV = 1;
        }
    }
  else if (CRDparam::DomainBCTypeExtrapolated & a_domT.m_type)
    {
      addIntervaltoVec(a_interiorExtrap, wCompIntv);
      solveV = 0; // Solve for temperature
      if (numSpecies > 0)
        {
          solveV = 1; // Otherwise, solve for density
        }
    }
  // Default is to set the gradients to 0. This can be overwritten
  a_WgradFace.setVal(0.);
  return solveV;
}

/*--------------------------------------------------------------------*/
//  Compute the cell averaged velocity gradients at the physical 
//  boundaries
/** 
 *  \param[out] a_GradPhicellAvgFab
 *                      Cell-averaged gradient of phi
 *  \param[in]  a_PhicellAvgFab   
 *                      Cell-averaged phi
 *  \param[in]  a_combinedBox
 *                      Box of 2 layers of ghost cells outside the domain
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  Side of the boundary
 *  \param[in]  a_numComp
 *                      Number of components to solve for
 *                      a_numComp = SpaceDim for velocity gradients
 *                      a_numComp = 1 for temperature gradients
 *  \param[in]  a_dxVect
 *                      Vector of computational grid spacing
 *//*-----------------------------------------------------------------*/

void
CNSIBC::computeCellAvgGradBdry(FArrayBox&           a_GradPhicellAvgFab,
                               const FArrayBox&     a_PhicellAvgFab,
                               const Box&           a_combinedBox,
                               const int            a_dir,
                               const Side::LoHiSide a_side,
                               const int            a_numComp,
                               const RealVect&      a_dxVect) const
{
  CH_assert(a_PhicellAvgFab.box().contains(a_combinedBox));
  CH_assert(a_GradPhicellAvgFab.box().contains(a_combinedBox));
  for (int phiComp = 0; phiComp != a_numComp; ++phiComp)
    {
      for (int gradDir = 0; gradDir != SpaceDim; ++gradDir)
        {
          const int gradVelComp = 
            ViscousTensor4thOrderOp::tensorIdxRowOrder(phiComp,gradDir);

          if (gradDir != a_dir)
            {
              // Solve for gradients tangent to the wall in the ghost cells
              FORT_CELLGRADDIR4THO(CHF_FRA1(a_GradPhicellAvgFab,gradVelComp),
                                   CHF_CONST_FRA1(a_PhicellAvgFab,phiComp),
                                   CHF_BOX(a_combinedBox),
                                   CHF_CONST_INT(gradDir),
                                   CHF_CONST_REAL(a_dxVect[gradDir]));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the cell averaged velocity gradients at ghost cells adjacent
//  to the corners
/** 
 *  \param[out] a_GradPhicellAvgFab
 *                      Cell-averaged gradient of phi
 *  \param[out] a_PhicellAvgFab   
 *                      Cell-averaged phi
 *  \param[in]  a_cornerBox
 *                      Box to solve gradients using one-sided stencil
 *  \param[in]  a_biasedBox
 *                      Box to solve gradients using biased stencil
 *  \param[in]  a_gradDir
 *                      Direction to solve the gradient
 *  \param[in]  a_side  Side of the boundary
 *  \param[in]  a_numComp
 *                      Number of components to solve for
 *                      a_numComp = SpaceDim for velocity gradients
 *                      a_numComp = 1 for temperature gradients
 *  \param[in]  a_dxVect
 *                      Vector of computational grid spacing
 *//*-----------------------------------------------------------------*/

void
CNSIBC::computeCornerGrad(FArrayBox&           a_GradPhicellAvgFab,
                          const FArrayBox&     a_PhicellAvgFab,
                          const Box&           a_cornerBox,
                          const Box&           a_biasedBox,
                          const int            a_gradDir,
                          const Side::LoHiSide a_side,
                          const int            a_numComp,
                          const Real           a_dx) const
{
  CH_assert(a_PhicellAvgFab.box().contains(a_cornerBox));
  CH_assert(a_PhicellAvgFab.box().contains(a_biasedBox));
  for (int phiComp = 0; phiComp != a_numComp; ++phiComp)
    {
      const int gradComp = 
        ViscousTensor4thOrderOp::tensorIdxRowOrder(phiComp, a_gradDir);
      const int lohiSign = Side::sign(a_side);
      FORT_ONESIDEDCELLAVGGRAD(
        CHF_FRA1(a_GradPhicellAvgFab,gradComp),
        CHF_CONST_FRA1(a_PhicellAvgFab,phiComp),
        CHF_BOX(a_cornerBox),
        CHF_CONST_INT(a_gradDir),
        CHF_CONST_INT(lohiSign),
        CHF_CONST_REAL(a_dx));

      FORT_BIASEDCELLAVGGRAD(
        CHF_FRA1(a_GradPhicellAvgFab,gradComp),
        CHF_CONST_FRA1(a_PhicellAvgFab,phiComp),
        CHF_BOX(a_biasedBox),
        CHF_CONST_INT(a_gradDir),
        CHF_CONST_INT(lohiSign),
        CHF_CONST_REAL(a_dx));
    }
}

/*--------------------------------------------------------------------*/
//  Compute the primative values on the wall faces
//  This works for mapped walls
/** 
 *  \param[out] a_Wface
 *  \param[in]  a_Wface
 *                      Primitive state to be corrected at wall
 *                      (face-centered average)
 *  \param[in]  a_boundaryFaceBox
 *                      Box of boundary faces to adjust
 *  \param[in]  a_Wcell             
 *                      Average primitive state in the cell
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip-velocity for boundary faces provided by
 *                      turbulence model wall-models
 *  \param[in]  a_gammaFaceFab      
 *                      Gamma values from Wface
 *  \param[in]  a_gammaCellFab      
 *                      Gamma values from Wcell
 *  \param[in]  a_wallVel           
 *                      Velocity of the wall (for no-slip)
 *  \param[in]  a_viscousSlip       
 *                      1 - Use slip conditions, 0 - No-slip
 *  \param[in]  a_dir               
 *                     Direction of the boundary
 *  \param[in]  a_lohiSign          
 *                     Sign indicating side of the boundary
 *                     (-1 = lo; +1 = hi)
 *//*-----------------------------------------------------------------*/

void
CNSIBC::computeWallPrimState(FArrayBox&       a_Wface,
                             const Box        a_boundaryFaceBox,
                             const FArrayBox& a_Wcell,
                             const FArrayBox& a_boundarySlipVelocity,
                             const FArrayBox& a_gammaCellFab,
                             const FArrayBox& a_unitNormalBasisFab,
                             const RealVect&  a_wallVel,
                             const int        a_viscousSlip,
                             const int        a_dir,
                             const int        a_lohiSign) const
{
  CH_assert(a_Wface.box().contains(a_boundaryFaceBox));
  CH_assert(a_Wcell.box().contains(a_boundaryFaceBox));
  CH_assert(a_gammaCellFab.box().contains(a_boundaryFaceBox));
  CH_assert(a_unitNormalBasisFab.box().contains(a_boundaryFaceBox));
  if (a_viscousSlip != 1)
    {
      CH_assert(a_boundarySlipVelocity.box().contains(a_boundaryFaceBox));
    }
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int numSpecies = CRDparam::g_numSpecies;
  const int specIndx = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  // Only use the temperature values if thermally perfect physics is used
  int tempIndx = 0;
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
    }
  const Interval velInterval = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx = velInterval.begin();

  // A flag for if the wall velocity is being set in normal-tangent components
  // or in x-y-z components.
  // Should probably be a input, but this needs a little work still in 3D
  const bool setNormTanVel = (a_viscousSlip == 1) ? true : false;

  // Transform velocity into normal and tangent components
  FArrayBox velFabFace(velInterval, a_Wface); // this is an alias
  PatchMappedFunc::forwardTransform(velFabFace,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  // Copy and transform the velocity instead of modifying cell velocities
  const int tempVelIdx = 0;
  FABSTACKTEMP(WcellVel, a_boundaryFaceBox, SpaceDim);
  PatchMappedFunc::forwardTransform(WcellVel,
                                    tempVelIdx,
                                    a_Wcell,
                                    velIndx,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);

  MD_BOXLOOP(a_boundaryFaceBox, i)
    {
      // Thermodynamics
      Real unormFace = a_Wface[MD_IX(i, velIndx+a_dir)];
      Real unormCell = WcellVel[MD_IX(i, a_dir)];
      Real gammaCell = a_gammaCellFab[MD_IX(i, 0)];

      Real rgas = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int wComp = specIndx + sp;
          rgas += a_Wcell[MD_IX(i, wComp)]*
            CRDparam::g_CRDPhysics->speciesGasConstant(sp);
        }

      // Limiting
      // unormFace is interpolated. Logically, it must be between 0 and
      // unormCell. Limit unormFace so it is in this range.
      if (unormFace*unormCell < 0.)
        {
          unormFace = 0.0;
        }
      else if (std::fabs(unormFace) > std::fabs(unormCell))
        {
          unormFace = unormCell;
        }
      
      // Knowing unormFace and unormCell, the limit of the thermodynamic
      // state at the face is estimated. An acoustic correction is used
      // based on the difference between the normal velocity at the cell
      // center and the face
      Real unormDiff = unormCell - unormFace;
      Real presCell = a_Wcell[MD_IX(i,presIndx)];
      Real rhoCell = a_Wcell[MD_IX(i,rhoIndx)];
      Real TCell = a_Wcell[MD_IX(i,tempIndx)];
      Real rho    = std::max(rhoCell, CRDparam::g_smallr);
      Real pres   = std::max(presCell, CRDparam::g_smallp);
      Real c      = std::sqrt(gammaCell*pres/rho);
      Real T      = TCell;
      Real presLimit = pres + a_lohiSign*rho*unormDiff*c;
      Real rhoLimit = rho*std::pow(presLimit/pres, 1./gammaCell);
      Real TLimit = T*std::pow(presLimit/pres, (gammaCell - 1.)/gammaCell);
      
      // Limit the density and pressure at the face
      if (unormFace*a_lohiSign > 0.) // Compression
        {
          // Pressure and density should increase approaching the wall
          rho  = std::max(rhoLimit, a_Wface[MD_IX(i, rhoIndx)]);
          pres = std::max(presLimit, a_Wface[MD_IX(i, presIndx)]);
          T    = std::max(TLimit, a_Wface[MD_IX(i, tempIndx)]);
        }
      else  // Expansion
        {
          // Pressure and density should decrease approaching the wall
          rho  = std::min(rhoLimit, a_Wface[MD_IX(i, rhoIndx)]);
          pres = std::min(presLimit, a_Wface[MD_IX(i, presIndx)]);
          T    = std::min(TLimit, a_Wface[MD_IX(i, tempIndx)]);
        }
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          pres = std::max(pres, CRDparam::g_smallp);
          rho  = pres/(rgas*T);
          c    = std::sqrt(gammaCell*pres/rho);

          if (CRDparam::g_wallCorrections)
            {
              // Acoustic correction: delta u = delta p/(rho*c)
              Real presFace = pres + a_lohiSign*rho*unormFace*c;
              a_Wface[MD_IX(i, presIndx)] = presFace;
              Real TFace =
                T*std::pow(presFace/pres, (gammaCell - 1.)/gammaCell);
              a_Wface[MD_IX(i, tempIndx)] = TFace;
              a_Wface[MD_IX(i, rhoIndx)] = presFace/(rgas*TFace);
            }
          // First-order extrapolation of mass fractions and tangential velocity
          //**FIXME: There should be an option here to keep the existing
          //         4th-order wall face state. It should only be overwritten if
          //         the boundary is 1st-order.
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              const int wComp = specIndx + sp;
              a_Wface[MD_IX(i, wComp)] = a_Wcell[MD_IX(i, wComp)];
            }
          for (int tandir = 0; tandir != SpaceDim; ++tandir)
            {
              if (tandir != a_dir)
                {
                  a_Wface[MD_IX(i, velIndx+tandir)] =
                    WcellVel[MD_IX(i, tandir)];
                }
            }
        }
      else
        {
          rho  = std::max(rho, CRDparam::g_smallr);
          pres = std::max(pres, CRDparam::g_smallp);
          c    = std::sqrt(gammaCell*pres/rho);

          if (CRDparam::g_wallCorrections)
            {
              // Acoustic correction: delta u = delta p/(rho*c)
              Real presFace = pres + a_lohiSign*rho*unormFace*c;
              presFace = std::max(presFace, CRDparam::g_smallp);
              a_Wface[MD_IX(i, presIndx)] = presFace;
              // Isentropic correction: rho2/rho1 = (p2/p1)^(1/gamma)
              a_Wface[MD_IX(i, rhoIndx)] =
                std::max(rho*std::pow(presFace/pres, 1./gammaCell),
                         CRDparam::g_smallr);
            }
        }

      // Momentum
      // Specified velocity is x-y-z by default, but if we apply the reverse
      // transformation we instead specify normal-tangent velocity
      if (a_viscousSlip == 1)
        {
          // Normal velocity is zero (make sure to apply transform later)
          a_Wface[MD_IX(i, velIndx + a_dir)] = 0.;
        }
      else
        {
          D_TERM(
            a_Wface[MD_IX(i, velIndx)]   =
            a_boundarySlipVelocity[MD_IX(i, 0)];,
            a_Wface[MD_IX(i, velIndx + 1)] =
            a_boundarySlipVelocity[MD_IX(i, 1)];,
            a_Wface[MD_IX(i, velIndx + 2)] =
            a_boundarySlipVelocity[MD_IX(i, 2)];);
        }
    }

  // Do the velocity transform if norm-tan velocity was specified, this
  // gets us back to x-y-z space. Be careful in 3D
  if (setNormTanVel)
    {
      PatchMappedFunc::reverseTransform(velFabFace,
                                        a_unitNormalBasisFab,
                                        a_boundaryFaceBox);
    }
}

/*--------------------------------------------------------------------*/
//  Set the boxes that are inflow or outflow if mixed BC is being used
/** To use this function you must first set the boundary in question 
    in the problem specific CNSIBC file constructors 
    setDomainBC(0, 
    0, 
    CRDparam::DomainBCTypeMixed |
    CRDparam::DomainBCTypeInflow |
    CRDparam::DomainBCTypeWall);
    This creates a lower x boundary that is both an inflow portion and a wall
    portion. Then, within your problem specific CNSIBC file, you must create
    this setMixedBC function. Inside the function, you must set each box in
    a_boxVect to be a box on the boundary of interest and, respectively, 
    set the domain type in a_domainType to match the box. Here is an example
    where we want to set the lower boundary to both a jet defined by member data
    called m_jetBoxLo and a wall defined by the rest of the plate.

    Box boundaryJet;
    if (a_side == Side::Lo)
    {
      boundaryJet = m_jetBoxLo;
      boundaryJet.shiftHalf(a_dir, -1);
      boundaryJet.refine(CRDparam::g_refFromBase[a_level]);
      a_boxVect[0].define(boundaryJet);
      a_boxVect[0]&=a_boundaryFaceBox;
      a_domainType[0] = CRDparam::DomainBCTypeInflow;
      a_orderBC[0] = 1;
    {
    if (boundaryJet.smallEnd(0) >
        a_boundaryFaceBox.smallEnd(0))
    {
      a_boxVect[1].define(a_boundaryFaceBox);
      a_boxVect[1].setBig(
      m_jetWallDir,
      std::min(a_boundaryFaceBox.bigEnd(0),
      boundaryJet.smallEnd(0) - 1));
      a_domainType[1] = CRDparam::DomainBCTypeWall;
      a_orderBC[1] = 4;
    }
    
    We set a_boxVect[0] to be the box representing the jet so we also have to
    set a_domainType[0] to be the inflow boundary. We then set a_boxVect[1] to
    be the wall so we must also set a_domainType[1] to be the wall. This order
    could also be reversed; we could set a_boxVect[1] to be the jet box as long
    as a_domainType[1] is also set to inflow.
 * \param[in]  a_boundaryFaceBox
 *                      Box of boundary being operated on
 *  \param[in]  a_bcIdx Index of the multi-block boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_boxVect
 *                      Vector of boxes (order corresponds to a_domainType
 *                      vector
 *  \param[out] a_domainType
 *                      Vector of the domain types corresponding to the boxVect
 *  \return             0 - Successfully applied multiple BCs
 *                      1 - Boundary only has one BC
 *//*-----------------------------------------------------------------*/

int
CNSIBC::setMixedBC(const Box&           a_boundaryFaceBox,
                   const BoundaryIndex& a_bcIdx,
                   const Box&           a_disjointBox,
                   LevelGridMetrics&    a_gridMetrics,
                   const Real           a_time,
                   const int            a_level,
                   Vector<Box>&         a_boxVect,
                   Vector<BCInfo>&      a_bcInfo) const
{
  return 1;
}

/*--------------------------------------------------------------------*/
//  Extrapolate ghost values based on imposed face values
/** \param[in]  a_ghostBox1
 *                      Box of first layer of exterior ghost cells
 *  \param[in]  a_Wcell Cell primitive state in interior
 *  \param[out] a_Wcell Cell primitive state with filled ghosts
 *  \param[in]  a_Wface Face values shifted to align with first layer
 *                      of external ghost cells
 *  \param[in]  a_side  Lo or hi side
 *  \param[in]  a_dir   Normal direction to boundary
 *  \param[in]  a_compList
 *                      Vector of ints that correspond to components
 *//*-----------------------------------------------------------------*/

void
CNSIBC::extrapGhostFromFaceVals(const Box&            a_ghostBox1,
                                FArrayBox&            a_Wcell,
                                const FArrayBox&      a_Wface,
                                const Side::LoHiSide& a_side,
                                const int             a_dir,
                                const Vector<int>&    a_compList) const
{
  //**FIXME: The order of this extrapolation should be reducible
  const int lohiSign = sign(a_side);
  CH_assert(a_Wcell.contains(a_ghostBox1));
  CH_assert(a_Wface.contains(a_ghostBox1));
  const int MD_ID(o, a_dir);
  for (int ref = 0; ref != a_compList.size(); ++ref)
    {
      const int comp = a_compList[ref];
      MD_BOXLOOP(a_ghostBox1, i)
        {
          if (CRDparam::g_diffusiveDerivativeOrder == 2)
            {
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wjph1 = a_Wface[MD_IX(i,comp)];
              Real Wj = 2.*Wjph1 - 1.*Wjp1;
              Real Wjm1 =  4.*Wjph1 - 3.*Wjp1;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              CH_assert(Wj < 1.E300);
              CH_assert(Wjm1 < 1.E300);
            }
          else
            {
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wjp2 = a_Wcell[MD_OFFSETIX(i,-,2*lohiSign*o,comp)];
              Real Wjp3 = a_Wcell[MD_OFFSETIX(i,-,3*lohiSign*o,comp)];
              Real Wjph1 = a_Wface[MD_IX(i,comp)];
              Real Wj = 4.*Wjph1 - (Wjp3 - 5.*Wjp2 + 13.*Wjp1)/3.;
              Real Wjm1 =  20.*Wjph1 - (8.*Wjp3 - 37.*Wjp2 + 83.*Wjp1)/3. - Wj;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              CH_assert(Wj < 1.E300);
              CH_assert(Wjm1 < 1.E300);
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Extrapolate ghost values and overwrite face values for Neumann
//  boundary conditions
/** \param[in]  a_ghostBox1
 *                      Box of first layer of exterior ghost cells
 *  \param[in]  a_Wcell Cell primitive state in interior
 *  \param[out] a_Wcell Cell primitive state with filled ghosts
 *  \param[out] a_Wface Face values shifted to align with first layer
 *                      of external ghost cells
 *  \param[in]  a_WgradFace
 *                      Gradient values on boundary faces shifted to align
 *                      with first layer of external ghost cells
 *  \param[in]  a_side  Lo or hi side
 *  \param[in]  a_dir   Normal direction to boundary
 *  \param[in]  a_dx    Grid spacing in a_dir
 *  \param[in]  a_compList
 *                      Vector of ints that correspond to components
 *//*-----------------------------------------------------------------*/

void
CNSIBC::extrapGhostNeumann(const Box&            a_ghostBox1,
                           FArrayBox&            a_Wcell,
                           FArrayBox&            a_Wface,
                           const FArrayBox&      a_WgradFace,
                           const Side::LoHiSide& a_side,
                           const int             a_dir,
                           const Real            a_dx,
                           const Vector<int>&    a_compList) const
{
  //**FIXME: The order of this extrapolation should be reducible
  const int lohiSign = sign(a_side);
  const Real ca1 = lohiSign*12.*a_dx;
  const Real ca2 = lohiSign*a_dx;
  const int MD_ID(o, a_dir);
  CH_assert(a_Wcell.contains(a_ghostBox1));
  CH_assert(a_Wface.contains(a_ghostBox1));
  for (int ref = 0; ref != a_compList.size(); ++ref)
    {
      const int comp = a_compList[ref];
      MD_BOXLOOP(a_ghostBox1, i)
        {
          if (CRDparam::g_diffusiveDerivativeOrder == 2)
            {
              Real gradVal = a_WgradFace[MD_IX(i,comp)];
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wj = ca2*gradVal + Wjp1;
              Real Wjm1 = 2.*ca2*gradVal + Wjp1;
              Real Wface = (Wjp1 + Wj)/2.;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              a_Wface[MD_IX(i,comp)] = Wface;
              CH_assert(Wj < 1.E300);
              CH_assert(Wjm1 < 1.E300);
              CH_assert(Wface < 1.E300);
            }
          else
            {
              Real gradVal = a_WgradFace[MD_IX(i,comp)];
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wjp2 = a_Wcell[MD_OFFSETIX(i,-,2*lohiSign*o,comp)];
              Real Wjp3 = a_Wcell[MD_OFFSETIX(i,-,3*lohiSign*o,comp)];
              Real Wj = (ca1*gradVal + 9.*Wjp1 + 3.*Wjp2 - Wjp3)/11.;
              Real Wjm1 = -ca1*gradVal + 15.*Wj - 15.*Wjp1 + Wjp2;
              Real Wface = (7.*(Wjp1 + Wj) - (Wjm1 + Wjp2))/12.;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              a_Wface[MD_IX(i,comp)] = Wface;
              CH_assert(Wj < 1.E300);
              CH_assert(Wjm1 < 1.E300);
              CH_assert(Wface < 1.E300);
            } 
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Extrapolate ghost cells using interior cell-averaged values
/** \param[in]  a_ghostBox1
 *                      Box of first layer of exterior ghost cells
 *  \param[in]  a_Wcell Cell primitive state in interior
 *  \param[out] a_Wcell Cell primitive state with filled ghosts
 *  \param[in]  a_side  Lo or hi side
 *  \param[in]  a_dir   Normal direction to boundary
 *  \param[in]  a_compList
 *                      Vector of ints that correspond to components
 *//*-----------------------------------------------------------------*/

void
CNSIBC::extrapGhostFromInterior(const Box&            a_ghostBox1,
                                FArrayBox&            a_Wcell,
                                const Side::LoHiSide& a_side,
                                const int             a_dir,
                                const Vector<int>&    a_compList) const
{
  //**FIXME: The order of this extrapolation should be reducible
  const int lohiSign = sign(a_side);
  CH_assert(a_Wcell.contains(a_ghostBox1));
  const int MD_ID(o, a_dir);
  for (int ref = 0; ref != a_compList.size(); ++ref)
    {
      const int comp = a_compList[ref];
      MD_BOXLOOP(a_ghostBox1, i)
        {
          if (CRDparam::g_diffusiveDerivativeOrder == 2)
            {
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wjp2 = a_Wcell[MD_OFFSETIX(i,-,2*lohiSign*o,comp)];
              Real Wj = 2.*Wjp1 - 1.*Wjp2;
              Real Wjm1 =  3.*Wjp1 - 2.*Wjp2;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              CH_assert(Wj < 1.E300);
              CH_assert(Wjm1 < 1.E300);
            }
          else
            {
              Real Wjp1 = a_Wcell[MD_OFFSETIX(i,-,lohiSign*o,comp)];
              Real Wjp2 = a_Wcell[MD_OFFSETIX(i,-,2*lohiSign*o,comp)];
              Real Wjp3 = a_Wcell[MD_OFFSETIX(i,-,3*lohiSign*o,comp)];
              Real Wjp4 = a_Wcell[MD_OFFSETIX(i,-,4*lohiSign*o,comp)];
              Real Wj = 4.*Wjp1 - 6.*Wjp2 + 4.*Wjp3 - Wjp4;
              Real Wjm1 = 10.*Wjp1 - 20.*Wjp2 + 15.*Wjp3 - 4.*Wjp4;
              a_Wcell[MD_IX(i,comp)] = Wj;
              a_Wcell[MD_OFFSETIX(i,+,lohiSign*o,comp)] = Wjm1;
              CH_assert(Wjm1 < 1.E300);
              CH_assert(Wj < 1.E300);
            }
        }
    }
  return;
}

/*==============================================================================
 * Deprecated member functions
 * These override some base class functions and force an error.  Do not use
 * them.
 *============================================================================*/

// Initialize a level (DEPRECATED!)
void
CNSIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CRD::msg << "CNSIBC::initialize() is deprecated." << CRD::error;
}

/// Set domainBC on all boundaries to a single type
void CNSIBC::setAllDomainBC(const int a_domainBC)
{
  CRD::msg << "CNSIBC::setAllDomainBC is deprecated, use the updated function arguments" << CRD::error;
}

/// Set domainBC on all boundaries
void CNSIBC::setAllDomainBC(const CRDparam::DomainBCType* a_domainBC)
{
  CRD::msg << "CNSIBC::setAllDomainBC is deprecated, use the updated function arguments" << CRD::error;
}


/// Set domainBC on one boundary
void CNSIBC::setDomainBC(const int  a_dir,
                         const int  a_side,
                         const int  a_domainBC)
{
  CRD::msg << "CNSIBC::setDomainBC is deprecated, use the updated function arguments" << CRD::error;
}

  
/// Set BC order on all boundaries
void CNSIBC::setAllDomainBCOrder(const int a_order)
{
  CRD::msg << "CNSIBC::setAllDomainBCOrder is deprecated, use the updated function arguments" << CRD::error;
}


/// Set BC order on one boundary
void CNSIBC::setDomainBCOrder(const int  a_dir,
                              const int  a_side,
                              const int  a_order)
{
  CRD::msg << "CNSIBC::setDomainBCOrder is deprecated, use the updated function arguments" << CRD::error;
}


/// (DEPRECATED) use updated interface instead
/// Set the imposed (exterior or farfield) primitive state at flow BC
void
CNSIBC::setImposedBCprimState(
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
  CRD::msg << "CNSIBC::setImposedBCprimState function call has been updated!" << CRD::warn;
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir, a_side);
  BCInfo bc;
  bc.m_type = a_domT;
  
  setImposedBCprimState(
    a_Wface,
    a_boundaryFaceBox,
    a_Wcell,
    a_unitNormalBasisFab,
    bcIdx,
    a_disjointBox,
    a_gridMetrics,
    a_time,
    a_level,
    bc);
}

/// (DEPRECATED) use updated interface instead
/// Set the primitive state at wall BC
void
CNSIBC::setWallBCprimState(
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
  CRD::msg << "CNSIBC::setImposedBCprimState function call has been updated!" << CRD::warn;
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir, a_side);
  BCInfo bc;
  bc.m_type = a_domT;
  
  setWallBCprimState(
    a_Wface,
    a_boundaryFaceBox,
    a_Wcell,
    a_boundarySlipVelocity,
    a_unitNormalBasisFab,
    bcIdx,
    a_disjointBox,
    a_gridMetrics,
    a_time,
    a_level,
    bc);
}

/// (DEPRECATED) use updated interface instead
/// Set boxes for inflow or outflow conditions in mixed boundaries
int
CNSIBC::setMixedBC(const Box&                      a_boundaryFaceBox,
                   const int                       a_dir,
                   const Side::LoHiSide&           a_side,
                   const Box&                      a_disjointBox,
                   LevelGridMetrics&               a_gridMetrics,
                   const Real                      a_time,
                   const int                       a_level,
                   Vector<Box>&                    a_boxVect,
                   Vector<int>&                    a_orderBC,
                   Vector<CRDparam::DomainBCType>& a_domainType) const
{
  CRD::msg << "CNSIBC::setImposedBCprimState function call has been updated!" << CRD::warn;
  BoundaryIndex bcIdx;
  bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
               a_dir, a_side);

  Vector<BCInfo> bc(a_orderBC.size());
  for(int i=0; i!=bc.size(); ++i)
    {
      auto& thisBC = bc[i];
      thisBC.m_order = a_orderBC[i];
      thisBC.m_type = a_domainType[i];
    }

  return setMixedBC(a_boundaryFaceBox,
                    bcIdx,
                    a_disjointBox,
                    a_gridMetrics,
                    a_time,
                    a_level,
                    a_boxVect,
                    bc);
}
