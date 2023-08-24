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
 * \file CNSIBCSphere.cpp
 *
 * \brief Member functions for CNSIBCSphere
 *
 *//*+*************************************************************************/

#include <algorithm>

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
 
//----- Internal -----//

#include "CNSIBCSphere.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "CRDutil.H"
#include "PatchMappedFunc.H"


/*******************************************************************************
 *
 * Class CNSIBCGeneralized: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCSphere::CNSIBCSphere(const bool a_readInputs)
  :
  CNSIBCGeneralized(),
  m_radiusInner(1.),
  m_radiusOuter(2.),
  m_rotationRateRad(0.),
  m_omega(RealVect::Zero),
  m_specialInit(false)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCSphere::~CNSIBCSphere()
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
CNSIBCSphere::IBCName() const
{
  return "Flow over a sphere";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCSphere::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
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
CNSIBCSphere::initialize(LevelData<FArrayBox>&      a_U,
                         LevelGridMetrics&          a_gridMetrics,
                         const LayoutData<FluxBox>& a_unitNormals,
                         const Real                 a_time,
                         const int                  a_level) const
{
  // Set the initial values
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const CRDState& state = CRDState::get(m_idxStateInit);
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box4Dom = grow(box, 4);
      box4Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box4Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box4Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box4Dom, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      // Pointwise values of W (required on valid +4, for turbInitialize, except
      // at physical boundaries)
      FABSTACKTEMP(Wc, box4Dom, numWVar);
      Wc.setVal(state.density(), rhoIndx);
      // Set velocity based on potential flow solution
      const Real vel_inf = state.velocity()[0];
      const Real srcMu = -0.5*m_radiusInner*m_radiusInner*m_radiusInner*vel_inf;
      MD_BOXLOOP(box4Dom, i)
        {
          // Set location vector
          RealVect loc(D_DECL(XFab[MD_IX(i, 0)],
                              XFab[MD_IX(i, 1)],
                              XFab[MD_IX(i, 2)]));
          // Radial location vector with respect to the x-axis
          RealVect radLoc = loc;
          // Set x-component to zero (we only care about distance normal to x)
          radLoc[0] = 0.;
          // Get the length of the quasi-radial vector
          Real r = radLoc.vectorLength();
          // Create a normalized quasi-radial vector
          RealVect nrmdRadLoc = radLoc/r;
          Real L = r*r + loc[0]*loc[0];
          Real u_r = 3.*srcMu*r*loc[0]/std::pow(L, (5./2.));
          Real u_z = vel_inf
            + (srcMu/std::pow(L, (3./2.)))*((3.*loc[0]*loc[0]/L) - 1.);
          D_TERM(
            Wc[MD_IX(i, WvelIndx)] = u_z;,
            Wc[MD_IX(i, WvelIndx + 1)] = nrmdRadLoc[1]*u_r;,
            Wc[MD_IX(i, WvelIndx + 2)] = nrmdRadLoc[2]*u_r;);
        }
      if (m_specialInit)
        {
          // Only overwrite the first layer of cells
          Box firstLayerBox = box4Dom;
          firstLayerBox.setBig(1, 0);
          // We assume the center of the sphere is always at (0,0,0)
          RealVect rotationOrigin = RealVect::Zero;

#if (CH_SPACEDIM == 3)
          MD_BOXLOOP(firstLayerBox, i)
            {
              // We know the cell-center is the end point of the vector
              RealVect surfaceLoc(D_DECL(XFab[MD_IX(i, 0)],
                                         XFab[MD_IX(i, 1)],
                                         XFab[MD_IX(i, 2)]));
              RealVect rVect = surfaceLoc - rotationOrigin;

              RealVect surfaceVel = RealVect::Zero;
              D_TERM(,,
                     surfaceVel[0] = m_omega[1]*rVect[2]-m_omega[2]*rVect[1];
                     surfaceVel[1] = m_omega[2]*rVect[0]-m_omega[0]*rVect[2];
                     surfaceVel[2] = m_omega[0]*rVect[1]-m_omega[1]*rVect[0];);

              D_TERM(Wc[MD_IX(i, WvelIndx)]   = surfaceVel[0];,
                     Wc[MD_IX(i, WvelIndx+1)] = surfaceVel[1];,
                     Wc[MD_IX(i, WvelIndx+2)] = surfaceVel[2];);
            }
#elif (CH_SPACEDIM == 2)
          MD_BOXLOOP(firstLayerBox, i)
            {
              // We know the cell-center is the end point of the vector
              RealVect surfaceLoc(D_DECL(XFab[MD_IX(i, 0)],
                                         XFab[MD_IX(i, 1)],
                                         XFab[MD_IX(i, 2)]));
              RealVect rVect = surfaceLoc - rotationOrigin;

              RealVect surfaceVel = RealVect::Zero;
              // For 2D, m_omega = [0, 0, m_rotationRateRad]
              D_TERM(,
                     surfaceVel[0] = -m_rotationRateRad*rVect[1];
                     surfaceVel[1] =  m_rotationRateRad*rVect[0];,);

              D_TERM(Wc[MD_IX(i, WvelIndx)]   = surfaceVel[0];,
                     Wc[MD_IX(i, WvelIndx+1)] = surfaceVel[1];,);
            }
#endif
        }
      Wc.setVal(state.pressure(), presIndx);
      Wc.setVal(state.temperature(), tempIndx);
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      if (!m_specialInit)
        {
          fourthOrderAverageCell(UFab, blockDomain, box1Dom);
        }
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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCSphere::addSourceTerm(FArrayBox&           a_sourceFab,
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
  return;
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
CNSIBCSphere::setBndryWallVelocity(FArrayBox&           a_wallVelocity,
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
  // Set boundary velocity based on where we're at on the sphere

  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));

  // Grow the evaluation box by 1 in tangential directions for convolution
  Box grownBox = a_boundaryFaceBox;
  for (const int dir : EachDir)
    {
      if (dir != a_dir)
        {
          grownBox.grow(dir, 1);
        }
    }

  // Get physical coordinates
  FABSTACKTEMP(XiFab, grownBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, grownBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getFaceCoordinates(grownBox, XiFab, XFab, a_dir, blockCoordSys);

  // We assume the center of the sphere is always at (0,0,0)
  RealVect rotationOrigin = RealVect::Zero;

  // Pointwise velocity field on the faces
  FABSTACKTEMP(velFab, grownBox, SpaceDim);

#if (CH_SPACEDIM == 3)
  MD_BOXLOOP(grownBox, i)
    {
      // We know that the point on the sphere is the end point of the vector
      RealVect surfaceLoc(D_DECL(XFab[MD_IX(i, 0)],
                                 XFab[MD_IX(i, 1)],
                                 XFab[MD_IX(i, 2)]));
      RealVect rVect = surfaceLoc - rotationOrigin;

      RealVect surfaceVel = RealVect::Zero;
      D_TERM(,,
             surfaceVel[0] = m_omega[1]*rVect[2] - m_omega[2]*rVect[1];
             surfaceVel[1] = m_omega[2]*rVect[0] - m_omega[0]*rVect[2];
             surfaceVel[2] = m_omega[0]*rVect[1] - m_omega[1]*rVect[0];);

      D_TERM(velFab[MD_IX(i, 0)] = surfaceVel[0];,
             velFab[MD_IX(i, 1)] = surfaceVel[1];,
             velFab[MD_IX(i, 2)] = surfaceVel[2];);
    }
#elif (CH_SPACEDIM == 2)
  MD_BOXLOOP(grownBox, i)
    {
      // We know that the point on the sphere is the end point of the vector
      RealVect surfaceLoc(D_DECL(XFab[MD_IX(i, 0)],
                                 XFab[MD_IX(i, 1)],
                                 XFab[MD_IX(i, 2)]));
      RealVect rVect = surfaceLoc - rotationOrigin;

      RealVect surfaceVel = RealVect::Zero;
      // For 2D, it's always true that m_omega = [0, 0, m_rotationRateRad]
      D_TERM(,
             surfaceVel[0] = -m_rotationRateRad*rVect[1];
             surfaceVel[1] =  m_rotationRateRad*rVect[0];,);

      D_TERM(velFab[MD_IX(i, 0)] = surfaceVel[0];,
             velFab[MD_IX(i, 1)] = surfaceVel[1];,
             velFab[MD_IX(i, 2)] = surfaceVel[2];);
    }
#endif
  // Convolve the surface velocity
  int order = 4;
  CRDutil::convolveFace(a_wallVelocity, velFab,
                        a_boundaryFaceBox, a_domain,
                        Interval(0,SpaceDim-1),
                        a_dir, order, true, false, false);
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
CNSIBCSphere::readBCInfo()
{
  ParmParse ppCOORD("coordsys");

  ppCOORD.query("inner_radius", m_radiusInner);
  ppCOORD.query("outer_radius", m_radiusOuter);

  ParmParse ppIBC("ibc");

  Real rotationRateRad = 0.;
  ppIBC.query("rotational_rate_rad_per_sec", rotationRateRad);
  if (rotationRateRad == 0.)
    {
      Real rotationRateRev = 0.;
      ppIBC.query("rotational_rate_rev_per_sec", rotationRateRev);
      rotationRateRad = rotationRateRev*2.*PI;
    }
  m_rotationRateRad = rotationRateRad;

  std::vector<Real> rotationAxis(SpaceDim);
  ppIBC.queryarr("rotational_axis", rotationAxis, 0, SpaceDim);
  RealVect omega = RealVect::Zero;
  SpaceDimArray<Real, Real>::loadFromArray(omega.dataPtr(),
                                           &rotationAxis.front());
  m_omega = omega*rotationRateRad;

  ppIBC.query("special_initialization", m_specialInit);

  m_readInput = true;
}
