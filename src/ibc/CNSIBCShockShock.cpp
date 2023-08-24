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
 * \file CNSIBCShockShock.cpp
 *
 * \brief Member functions for shock-shock intersections IBC
 *
 *//*+*************************************************************************/

#include <algorithm>

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "CNSIBCShockShock.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class CNSIBCShockShock: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCShockShock::CNSIBCShockShock()
  :
  CNSIBC(),
  m_beta(30.0*PI/180.0),
  m_m(std::tan(m_beta)),
  //      y            x intersection
  m_b(-0.02 - m_m*(-0.75))
{

//--Read any BC info

  readBCInfo();

//--Wall state (defined as all zero if not yet set)

  {
    CRDState& wallState = CRDState::get("wall");
    (void)wallState;
  }

//--Initial state (pre-shock state)

  if (CRDState::nameIndex("initial") == -1) // Define with default values
    {
      CRDState& state0 = CRDState::get("initial");
      state0.density()  = CRDparam::g_gamma;
      state0.pressure() = 1.0;
      state0.velocity() = RealVect{10., 0., 0.};
      state0.setExtraThermo();
    }
  m_idxStatePreShock = CRDState::nameIndex("initial");

//--Set common reference state

  Real r0;  // Saved since addition of new states makes refs to state0 invalid
  Real p0;
  Real c0;
  Real u0;
  Real gamma;
  {
    const CRDState& state0 = CRDState::get(m_idxStatePreShock);
    r0 = state0.density();
    p0 = state0.pressure();
    u0 = state0.velocity()[0];
    CNSIBC::defineCommonRef(r0, u0);

    // Compute additional state info in front of the shock
    const Real* cn = nullptr;
    if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
      // Note however that this code assumes a frozen gas for computing the
      // shock
      {
        cn = &state0(CRDparam::g_CRDPhysics->speciesPrimInterval().begin());
      }
    gamma = CRDparam::g_CRDPhysics->gamma(state0.temperature(), cn);
    c0 = std::sqrt(gamma*p0/r0);
  }

//--Set the post-shock state

  if (CRDState::nameIndex("post-shock") == -1) // Define with default values
    {
      CRDState& state1 = CRDState::get("post-shock");
      const Real cosBeta = std::cos(m_beta);
      const Real sinBeta = std::sin(m_beta);
      const Real M0 = u0/c0;
      pout() << "M0: " << M0 << std::endl;
      pout() << "c0: " << c0 << std::endl;
      CH_assert(M0 > 1.);
      const Real Mn0 = std::fabs(M0*sinBeta);
      const Real Mn0sq = Mn0*Mn0;
      const Real theta = 2*cosBeta/sinBeta*(Mn0sq - 1.0)/
        (M0*M0*(gamma + std::cos(2*m_beta)) + 2.0);
      const Real Mn1 = std::sqrt((Mn0sq + (2.0/(gamma - 1.0)))/
                                 (2.0*gamma*Mn0sq/(gamma - 1.0) - 1.0));
      const Real M1 = std::fabs(Mn1/(std::sin(m_beta - theta)));
      pout() << "M1: " << M1 << std::endl;
      CH_assert(M1 > 1.);
      state1.density() = r0*(gamma + 1.0)*Mn0sq/((gamma - 1.0)*Mn0sq + 2.0);
      state1.pressure() = p0*(1.0 + 2.0*gamma/(gamma + 1.0)*(Mn0sq - 1.0));
      const Real c1 = std::sqrt(gamma*state1.pressure()/state1.density());
      const Real Vmag1 = M1*c1;
      pout() << "c1: " << c1 << std::endl;
      state1.velocity() = RealVect{
        Vmag1*std::cos(theta),
        Vmag1*std::sin(theta),
        0.0
      };
      state1.setExtraThermo();
    }
  m_idxStatePostShock = CRDState::nameIndex("post-shock");

//--Set BC Types

  // Assumes single block
  if (CRDparam::g_coordSys->numBlocks() != 1)
    {
      CRD::msg << "Joukowski airfoil problem requires single block!"
               << CRD::error;
    }

  BCInfo bc;
  // Wall at y-min
  bc.m_type = CRDparam::DomainBCTypeAdiabaticWall;
  bc.m_idxState = CRDState::nameIndex("wall");
  bc.m_order = 4;
  setDomainBC(BoundaryIndex(0, 1, Side::Lo), bc);
  // Farfield at y-max
  bc.m_type = CRDparam::DomainBCTypeFarfield;
  bc.m_idxState = m_idxStatePreShock;
  bc.m_order = 4;
  setDomainBC(BoundaryIndex(0, 1, Side::Hi), bc);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCShockShock::~CNSIBCShockShock()
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
CNSIBCShockShock::IBCName() const
{
  return "Shock-shock interaction";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCShockShock::writeIBCInfo() const
{
  CRDState::writeStateInfo();
  CRD::msg.newline();
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
CNSIBCShockShock::initialize(LevelData<FArrayBox>&      a_U,
                             LevelGridMetrics&          a_gridMetrics,
                             const LayoutData<FluxBox>& a_unitNormals,
                             const Real                 a_time,
                             const int                  a_level) const
{
  // References to pre and post-shock states in an array
  std::reference_wrapper<CRDState> states[] =
    {
      CRDState::get(m_idxStatePreShock),
      CRDState::get(m_idxStatePostShock),
    };

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box4Dom = grow(box, 4);
      box4Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box4Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box4Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box4Dom, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box1Dom));
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      // Pointwise values of W (required on valid +4, for turbInitialize, except
      // at physical boundaries)
      FABSTACKTEMP(Wc, box4Dom, CRDparam::g_CRDPhysics->numPrimitive());

      MD_BOXLOOP(box4Dom, i)
        {
          const int idxState =
            (m_beta >= 0.0) !=                           // T: left running
            (XFab[MD_IX(i, 1)] > m_m*XFab[MD_IX(i, 0)] + m_b);  // T: above
          for (int c = 0, c_end = CRDparam::g_CRDPhysics->numPrimitive();
               c != c_end; ++c)
            {
              Wc[MD_IX(i, c)] = states[idxState](c);
            }
        }

      // Set conservative state
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);

      // Get fourth-order averages
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
}

  


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the exterior (or farfield) primitive state at flow BC
/** This override is used to impose an oblique shock in front of the
 *  strut.  This discontinuity can extend to the back side where
 *  supersonic outflow obviates the need for BC.
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
CNSIBCShockShock::setImposedBCprimState(
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
  if (a_bcIdx.m_dir == 1 && a_bcIdx.m_side == Side::Hi)
    {
      // References to pre and post-shock states in an array
      std::reference_wrapper<CRDState> states[] =
        {
          CRDState::get(m_idxStatePreShock),
          CRDState::get(m_idxStatePostShock),
        };

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
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          const int idxState =
            (m_beta >= 0.0) !=                           // T: left running
            (XFab[MD_IX(i, 1)] > m_m*XFab[MD_IX(i, 0)] + m_b);  // T: above
          for (int c = 0, c_end = CRDparam::g_CRDPhysics->numPrimitive();
               c != c_end; ++c)
            {
              a_Wface[MD_IX(i, c)] = states[idxState](c);
            }
        }
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
CNSIBCShockShock::readBCInfo()
{
  m_readInput = true;
}
