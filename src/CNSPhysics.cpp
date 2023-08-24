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
 * \file AMRLevelCNSPhysics.cpp
 *
 * \brief Member functions for CNSPhysics
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "LoHiCenter.H"
#include "CellToEdge.H"
#include "EdgeToCellF_F.H"
#include "GodunovUtilitiesF_F.H"
// Artificial viscosity for Euler equations
#include "MOLPhysicsMappedArtViscF_F.H"
#include "PolytropicPhysicsF_F.H"
#include "MaxScaledAcousticSpeedF_F.H"
#include "FourthOrderUtil.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "VEx.H"
#include "CNSPhysics.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "PatchMappedFunc.H"
#include "CNSIBC.H"
#include "PolytropicPhysicsF_F.H"
#include "PolytropicPhysics_vex.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class CNSPhysics: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

CNSPhysics::CNSPhysics()
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSPhysics::~CNSPhysics()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Compute the maximum wave speed on mapped grids
/** \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_box   Box to find 1/dt_convective
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[out] a_invDtFab
 *                      Sum of inverse of time step sizes
 *  \param[in]  a_cellAvgW
 *                      Average prim state in cells to find max wavespeed
 *  \param[in]  a_N     Grid metrics on the box
 *  \param[in]  a_J     Metrics Jacobian on the box
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_stabilityCons
 *                      Inertial time step stability constraint
 *  \param[in]  a_dxVect
 *                      RealVect of grid spacing
 *  \param[out] a_minConvDt
 *                      Minimum convective time step
 *  \param[out] a_minConvDtCell
 *                      Cell with minimum convective time step
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::getMaxWaveSpeed(const ProblemDomain&    a_problemDomain,
                            const Box&              a_box,
                            const Box&              a_disjointBox,
                            FArrayBox&              a_invDtFab,
                            const FArrayBox&        a_cellAvgW,
                            const FluxBox&          a_N,
                            const FArrayBox&        a_J,
                            const LevelGridMetrics& a_levelGridMetrics,
                            const Real&             a_stabilityCons,
                            const RealVect&         a_dxVect,
                            Real&                   a_minConvDt,
                            IntVect&                a_minConvDtCell) const
{
  const Box box1 = grow(a_box, 1);
  const Box box1Dom = box1 & a_problemDomain;

  // Convert to primitive state
  //**FIXME Better to use a 1-sided op at domain boundaries for
  //**CellToEdge instead of extrapolating
  FABSTACKTEMP(cellAvgW, box1, numNativePrimitive());
  cellAvgW.copy(a_cellAvgW, box1Dom, 0, box1Dom, 0, numNativePrimitive());

  // Second-order extrapolation at domain boundaries
  secondOrderCellExtrapAtDomainBdry(cellAvgW,
                                    a_box,
                                    a_problemDomain);

  // Get face values (second-order)
  FLUXBOXSTACKTEMP(faceAvgW, a_box, numConservative());
  CellToEdge(cellAvgW, faceAvgW);

  // Velocity in computational space, given as component 'd' on face 'd'
  FLUXBOXSTACKTEMP(faceAvgVelCSpc, a_box, 1);
  Interval velIntv(velocityInterval().begin(), velocityInterval().begin());
  a_levelGridMetrics.getCoordSys(a_disjointBox)
    ->computeMetricTermProductAverage(
      faceAvgVelCSpc,       // Result: Velocity in CSpc
      faceAvgW,             // Velocity in PSpc
      a_N,
      SpaceDim,
      faceAvgW,             // Unused for 2nd order
      a_box,
      false,                // 2nd order
      Interval(0, 0),       // Interval for faceAvgVelCSpc
      velIntv,              // Interval of Vel. in PSpc
      -1);                  // "F" (or velocity) is space contiguous but this
                            // is just a vector so do not adjust starting
                            // location

  // Face ratios (magnitude of N)
  FLUXBOXSTACKTEMP(faceNMag, a_box, 1);
  a_levelGridMetrics.getCoordSys(a_disjointBox)->magnitudeN(faceNMag, a_N,
                                                            a_box);
  
  // Reusing cellAvgW as temp space for sum of directional wavespeeds
  cellAvgW.setVal(0., 0);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box faceBox(a_box);
      faceBox.surroundingNodes(dir);
      FArrayBox faceSound(faceBox, 1);
      soundSpeed(faceSound, faceAvgW[dir], faceBox);

      // From MaxScaledAcousticSpeedF.ChF, scale the speed of sound
      // by the face area ratios and add it to |fluid speed| to get the
      // maximum wavespeed.
      FORT_MAXSCALEDACOUSTICSPEED(CHF_FRA1(faceAvgVelCSpc[dir], 0),
                                  CHF_CONST_FRA1(faceSound, 0),
                                  CHF_CONST_FRA1(faceNMag[dir], 0),
                                  CHF_BOX(faceBox));

      // Add an average of the max wave speed from the faces to the
      // cell.
      FORT_EDGETOINCREMENTCELL(CHF_CONST_FRA1(faceAvgVelCSpc[dir], 0),
                               CHF_FRA1(cellAvgW, 0),
                               CHF_BOX(a_box),
                               CHF_CONST_INT(dir));
    }

  // Divide by <J> to get the true computational wave speed
  cellAvgW.divide(a_J, a_box, 0, 0, 1);

  Real cfl = CRDparam::g_cfl;
  MD_ARRAY_RESTRICT(arrInvDt, a_invDtFab);
  MD_ARRAY_RESTRICT(arrVelSum, cellAvgW);
  MD_BOXLOOP(a_box, i)
    {
      Real maxVelSum = arrVelSum[MD_IX(i,0)];
#ifndef NDEBUG
      if (maxVelSum < 0.)
        {
          CRD::msg.setTerminateOnError(false);
          CRD::msg << "Failure in box " << a_box << " at cell " << MD_GETIV(i)
                   << CRD::error;
          CRD::msg.setTerminateOnError(true);
          CH_assert(maxVelSum >= 0.);
        }
#endif
      if (maxVelSum > 1.E-15)
        {
          Real convDt = a_stabilityCons*a_dxVect[0]/maxVelSum;
          Real invConvDt = 1./(cfl*convDt);
          arrInvDt[MD_IX(i,0)] += invConvDt;
          if (convDt < a_minConvDt)
            {
              a_minConvDt = convDt;
              a_minConvDtCell = MD_GETIV(i);
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Compute the maximum wave speed on mapped grids
/** \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_box   Box to find 1/dt_convective
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[out] a_invDtFab
 *                      Sum of inverse of time step sizes
 *  \param[in]  a_WfacePntFxb
 *                      Face-centered primitive variables
 *  \param[in]  a_N     Grid metrics on the box
 *  \param[in]  a_J     Metrics Jacobian on the box
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_stabilityCons
 *                      Inertial time step stability constraint
 *  \param[in]  a_dxVect
 *                      RealVect of grid spacing
 *  \param[out] a_minConvDt
 *                      Minimum convective time step
 *  \param[out] a_minConvDtCell
 *                      Cell with minimum convective time step
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::getMaxWaveSpeedEvalRHS(
  const ProblemDomain&    a_problemDomain,
  const Box&              a_box,
  const Box&              a_disjointBox,
  FArrayBox&              a_invDtFab,
  const FluxBox&          a_WfacePntFxb,
  const FluxBox&          a_N,
  const FArrayBox&        a_J,
  const LevelGridMetrics& a_levelGridMetrics,
  const Real&             a_stabilityCons,
  const RealVect&         a_dxVect,
  Real&                   a_minConvDt,
  IntVect&                a_minConvDtCell) const
{
  // Velocity in computational space, given as component 'd' on face 'd'
  FLUXBOXSTACKTEMP(faceAvgVelCSpc, a_box, 1);
  Interval velIntv(velocityInterval().begin(), velocityInterval().begin());
  a_levelGridMetrics.getCoordSys(a_disjointBox)
    ->computeMetricTermProductAverage(
      faceAvgVelCSpc,       // Result: Velocity in CSpc
      a_WfacePntFxb,        // Velocity in PSpc
      a_N,
      SpaceDim,
      a_WfacePntFxb,        // Unused for 2nd order
      a_box,
      false,                // 2nd order
      Interval(0, 0),       // Interval for faceAvgVelCSpc
      velIntv,              // Interval of Vel. in PSpc
      -1);                  // "F" (or velocity) is space contiguous but this
                            // is just a vector so do not adjust starting
                            // location

  // Face ratios (magnitude of N)
  FLUXBOXSTACKTEMP(faceNMag, a_box, 1);
  a_levelGridMetrics.getCoordSys(a_disjointBox)->magnitudeN(faceNMag, a_N,
                                                            a_box);
  
  // FAB for storing the mapped wave speeds
  FABSTACKTEMP(waveSpeedFab, a_box, 1);
  waveSpeedFab.setVal(0.);
  for (auto dir : EachDir)
    {
      Box faceBox(a_box);
      faceBox.surroundingNodes(dir);
      FArrayBox faceSound(faceBox, 1);
      soundSpeed(faceSound, a_WfacePntFxb[dir], faceBox);

      // From MaxScaledAcousticSpeedF.ChF, scale the speed of sound
      // by the face area ratios and add it to |fluid speed| to get the
      // maximum wavespeed.
      FORT_MAXSCALEDACOUSTICSPEED(CHF_FRA1(faceAvgVelCSpc[dir], 0),
                                  CHF_CONST_FRA1(faceSound, 0),
                                  CHF_CONST_FRA1(faceNMag[dir], 0),
                                  CHF_BOX(faceBox));

      // Add an average of the max wave speed from the faces to the
      // cell.
      FORT_EDGETOINCREMENTCELL(CHF_CONST_FRA1(faceAvgVelCSpc[dir], 0),
                               CHF_FRA1(waveSpeedFab, 0),
                               CHF_BOX(a_box),
                               CHF_CONST_INT(dir));
    }

  // Divide by <J> to get the true computational wave speed
  waveSpeedFab.divide(a_J, a_box, 0, 0, 1);

  Real cfl = CRDparam::g_cfl;
  MD_ARRAY_RESTRICT(arrInvDt, a_invDtFab);
  MD_ARRAY_RESTRICT(arrVelSum, waveSpeedFab);
  MD_BOXLOOP(a_box, i)
    {
      Real maxVelSum = arrVelSum[MD_IX(i,0)];
      CH_assert(maxVelSum >= 0.);
      if (maxVelSum > 1.E-15)
        {
          Real convDt = a_stabilityCons*a_dxVect[0]/maxVelSum;
          Real invConvDt = 1./(cfl*convDt);
          arrInvDt[MD_IX(i,0)] += invConvDt;
          if (convDt < a_minConvDt)
            {
              a_minConvDt = convDt;
              a_minConvDtCell = MD_GETIV(i);
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Compute the time step for the elliptic component, this operation
//  is only done once, at the beginning. The rest of the operations are
//  done in evalRHS
/**
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_box   Box to find 1/dt_diffusive
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[out] a_invDtFab
 *                      Sum of inverse of time step sizes
 *  \param[in]  a_cellAvgW
 *                      Cell primitive variables (2nd order)
 *  \param[in]  a_lambda
 *                      Constant value of lambda = 5.3333
 *  \param[in]  a_NTJ   Grid metrics on faces
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dxVect
 *                      RealVect of grid spacing
 *  \param[out] a_minDiffDt
 *                      Minimum diffusive time step
 *  \param[out] a_minDiffDtCell
 *                      Cell with minimum diffusive time step
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::ellipticDt(const ProblemDomain&    a_problemDomain,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       FArrayBox&              a_invDtFab,
                       const FArrayBox&        a_cellAvgW,
                       const Real&             a_lambda,
                       const FluxBox&          a_NTJ,
                       const LevelGridMetrics& a_levelGridMetrics,
                       const RealVect&         a_dxVect,
                       Real&                   a_minDiffDt,
                       IntVect&                a_minDiffDtCell) const
{
  CH_TIME("CNSPhysics::ellipticDt");
  // Calculate the dynamic viscosity
  FABSTACKTEMP(muFab, a_box, 1);
  FABSTACKTEMP(kappaFab, a_box, 1);
  calcCoeffKappaMu(a_box,
                   muFab,
                   kappaFab,
                   a_cellAvgW);
  FABSTACKTEMP(muTurb, a_box, 1);
  muTurb.setVal(0.);
  FABSTACKTEMP(kappaTurb, a_box, 1);
  if (CRDparam::g_turbModelType)
    {
      // FIXME: Must solve for strain rate should be solved for
      FABSTACKTEMP(dummyFab, a_box, SpaceDim*SpaceDim);
      dummyFab.setVal(0.);
      CRDparam::g_CRDPhysics->calcCoeffKappatMut(a_box,
                                                 muTurb,
                                                 kappaTurb,
                                                 muFab,
                                                 kappaFab,
                                                 a_cellAvgW,
                                                 dummyFab,
                                                 a_levelGridMetrics,
                                                 0); // dummy direction
    }
  muFab.plus(muTurb);
  solveViscDt(a_box, a_disjointBox, a_invDtFab, a_cellAvgW,
              muFab, a_lambda, a_NTJ, a_dxVect, a_minDiffDt, a_minDiffDtCell);
}

/*--------------------------------------------------------------------*/
//  Compute the time step for the elliptic component during evalRHS
/**
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_box   Box to find 1/dt_diffusive
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[out] a_invDtFab
 *                      Sum of inverse of time step sizes
 *  \param[in]  a_WcellPntFab
 *                      Cell-centered primitive values
 *  \param[in]  a_muFab Dynamic viscosity values averaged from adjacent faces
 *  \param[in]  a_lambda
 *                      Constant value of lambda = 5.3333
 *  \param[in]  a_NTJ   Grid metrics on faces
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dxVect
 *                      RealVect of grid spacing
 *  \param[out] a_minDiffDt
 *                      Minimum diffusive time step
 *  \param[out] a_minDiffDtCell
 *                      Cell with minimum diffusive time step
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::ellipticDtEvalRHS(const ProblemDomain&    a_problemDomain,
                              const Box&              a_box,
                              const Box&              a_disjointBox,
                              FArrayBox&              a_invDtFab,
                              const FArrayBox&        a_WcellPntFab,
                              const FArrayBox&        a_muFab,
                              const Real&             a_lambda,
                              const FluxBox&          a_NTJ,
                              const LevelGridMetrics& a_levelGridMetrics,
                              const RealVect&         a_dxVect,
                              Real&                   a_minDiffDt,
                              IntVect&                a_minDiffDtCell) const
{
  CH_TIME("CNSPhysics::ellipticDtEvalRHS");
  CH_assert(a_WcellPntFab.contains(a_box));
  CH_assert(a_muFab.contains(a_box));
  solveViscDt(a_box, a_disjointBox, a_invDtFab, a_WcellPntFab,
              a_muFab, a_lambda, a_NTJ, a_dxVect, a_minDiffDt, a_minDiffDtCell);
}

/*--------------------------------------------------------------------*/
//  Compute the speed of sound
/** \param[out] a_speed Sound speed in cells
 *  \param[in]  a_W     Native primitive state
 *  \param[in]  a_box   Box to compute sound speed in
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::soundSpeed(FArrayBox&       a_speed,
                       const FArrayBox& a_W,
                       const Box&       a_box) const
{
  constexpr int cRho  = CRDPhysics::densityIndex();
  constexpr int cPres = CRDPhysics::pressureIndex();
  CH_assert(a_W.contains(a_box));
  CH_assert(a_speed.contains(a_box));

  const Real gamma = CRDparam::g_gamma;
  MD_BOXLOOP(a_box, i)
    {
      Real rho  = a_W[MD_IX(i, cRho)];
      if (std::isnan(rho) || std::isinf(rho))
        {
          rho = 0.;
        }
      rho  = std::max(rho,  CRDparam::g_smallr);
      Real pres = a_W[MD_IX(i, cPres)];
      if (std::isnan(pres) || std::isinf(pres))
        {
          pres = 0.;
        }
      pres = std::max(pres, CRDparam::g_smallp);
      a_speed[MD_IX(i, 0)] = std::sqrt(gamma*pres/rho);
    }
}

/*--------------------------------------------------------------------*/
//  Compute a flux from primitive variable values on a face
/** \param[out] a_flux  Flux on the faces
 *  \param[in]  a_WFace Primitive state on the face
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_box   Box on which to compute flux
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::getFlux(FArrayBox&       a_flux,
                    const FArrayBox& a_WFace,
                    const int&       a_dir,
                    const Box&       a_box) const
{
  CH_assert(a_flux.contains(a_box));
  CH_assert(a_WFace.contains(a_box));

  FORT_GETFLUXF(CHF_FRA(a_flux),
                CHF_CONST_FRA(a_WFace),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->addTurbConvectiveFlux(a_flux, a_WFace, a_dir, a_box);
    }
  if (CRDparam::g_numTransportScalars > 0)
    {
      MD_ARRAY_RESTRICT(arrFlux, a_flux);
      MD_ARRAY_RESTRICT(arrWFace, a_WFace);
      const int cVel = vectorFluxInterval().begin() + a_dir;
      MD_BOXLOOP(a_box, i)
        {
          for (int cTr = transportConsInterval().begin(),
                 cTrEnd = transportConsInterval().end() + 1;
               cTr != cTrEnd; ++cTr)
            {
              arrFlux[MD_IX(i, cTr)] =
                arrWFace[MD_IX(i, cVel)]*arrWFace[MD_IX(i, cTr)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the solution to the Riemann problem
/** The CRD version of the Riemann problem does NOT apply boundary
 *  conditions.  This is different from most Chombo versions.
 *  \param[out] a_WStar Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_box   Face-centered box on which to compute a_WStar
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::riemann(FArrayBox&       a_WStar,
                    const FArrayBox& a_WLeft,
                    const FArrayBox& a_WRight,
                    const int&       a_dir,
                    const Box&       a_box) const
{
  CH_assert(a_WStar.contains(a_box));

  // Get the numbers of relevant variables
  const int numPrim = numFluxes();

  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  // Riemann solver computes Wgdnv all edges that are not on the physical
  // boundary.
#if (defined CRD_USE_VEX) &&                                            \
  (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  riemann_vex(a_WStar, shiftWLeft, shiftWRight, a_dir, a_box);
#else
  riemannSolution(a_WStar,
                  shiftWLeft,
                  shiftWRight,
                  a_dir,
                  a_box);
#endif

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

/*--------------------------------------------------------------------*/
//  Compute the solution to the Riemann problem
/** The CRD version of the Riemann problem does NOT apply boundary
 *  conditions.  This is different from most Chombo versions.
 * \param[out]  a_WavgFace
 *                      Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_unitNormalBasisFab
 *                      Unit normal basis Fab
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_boundaryFaceBox
 *                      Face-centered box on which to compute a_WavgFace
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::riemannBC(FArrayBox&            a_WavgFace,
                      FArrayBox&            a_WLeft,
                      FArrayBox&            a_WRight,
                      const FArrayBox&      a_unitNormalBasisFab,
                      const int             a_dir,
                      const Side::LoHiSide& a_side,
                      const Box&            a_boundaryFaceBox) const
{
  // Forward transform, the following is an alias
  FArrayBox velLeftFab(velocityInterval(), a_WLeft);
  PatchMappedFunc::forwardTransform(velLeftFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  // Forward transform, the following is an alias
  FArrayBox velRightFab(velocityInterval(), a_WRight);
  PatchMappedFunc::forwardTransform(velRightFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  riemannSolution(a_WavgFace,
                  a_WLeft,
                  a_WRight,
                  a_dir,
                  a_boundaryFaceBox);
  // Reverse transform, the following is an alias
  FArrayBox velStarFab(velocityInterval(),
                       a_WavgFace);
  PatchMappedFunc::reverseTransform(velStarFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
}

/*--------------------------------------------------------------------*/
//  Compute primitive variables from conserved variables
/** \param[out] a_W     Computed primitive state
 *  \param[in]  a_U     Conservative state
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::consToPrim(FArrayBox&       a_W,
                       const FArrayBox& a_U,
                       const Box&       a_box,
                       const FArrayBox& a_Wold) const
{
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_CONSTOPRIMF(CHF_FRA(a_W),
                   CHF_CONST_FRA(a_U),
                   CHF_BOX(a_box));
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->turbConsToPrim(a_W, a_U, a_box);
    }
  if (CRDparam::g_numTransportScalars > 0)
    {
      a_W.copy(a_box,
               transportPrimInterval(),
               a_box,
               a_U,
               transportConsInterval());
    }
}

/*--------------------------------------------------------------------*/
//  Compute conservative variables from primitive variables
/** \param[out] a_U     Computed conservative state
 *  \param[in]  a_W     Primitive state
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::primToCons(FArrayBox&       a_U,
                       const FArrayBox& a_W,
                       const Box&       a_box) const
{
  // CH_assert(false);  // I don't think this is ever used.  If so, I would like
  //                    // to know when
  // We're using this now in CNSIBCBluffBodyCombustion in addSourceTerm
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  FORT_PRIMTOCONSF(CHF_FRA(a_U),
                   CHF_CONST_FRA(a_W),
                   CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Set full average primitive state from native primitive state
/** The native state is the minimal state (e.g., pressure and density)
 *  from which other values are derived (.e.g, T = p/(rho*R)).  For
 *  CNS, the full state adds temperature.
 *  \param[out] a_Wx    Location to store extra primitive state
 *  \param[in]  a_compWxBeg
 *                      Location in 'a_Wx' where extra components (T)
 *                      begin
 *  \param[in]  a_W     Native primitive state (rho, vel, pres)
 *  \param[in]  a_box   Where to compute
 *  
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::extraPrimitiveState(FArrayBox&       a_Wx,
                                const int        a_compWxBeg,
                                const FArrayBox& a_Wp,
                                const Box&       a_box) const
{
  CH_assert(a_Wx.box().contains(a_box));
  CH_assert(a_Wp.box().contains(a_box));

  FORT_PRIMTEMPERATURE(CHF_FRA1(a_Wx, a_compWxBeg),
                       CHF_CONST_FRA(a_Wp),
                       CHF_BOX(a_box),
                       CHF_CONST_REAL(CRDparam::g_R));
}

/*--------------------------------------------------------------------*/
//  Set full average primitive state from native primitive state
//  in the same FAB
/** The native state is the minimal state (e.g., pressure and density)
 *  from which other values are derived (.e.g, T = p/(rho*R)).  For
 *  CNS, the full state adds temperature.
 *  \param[out] a_W     Primitive state variables to solve for extra
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::extraPrimitiveState(FArrayBox& a_W,
                                const Box& a_box) const
{
  CH_assert(a_W.box().contains(a_box));

  FORT_PRIMTEMPERATURE(CHF_FRA1(a_W, temperatureIndex()),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(a_box),
                       CHF_CONST_REAL(CRDparam::g_R));
}
                                  
/*--------------------------------------------------------------------*/
//  Compute the temperature from primary primitive variables
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with temperature
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::temperature(FArrayBox& a_W,
                        const Box& a_box) const
{
  CH_assert(a_W.box().contains(a_box));

  FORT_PRIMTEMPERATURE(CHF_FRA1(a_W, temperatureIndex()),
                       CHF_CONST_FRA(a_W),
                       CHF_BOX(a_box),
                       CHF_CONST_REAL(CRDparam::g_R));
}

/*--------------------------------------------------------------------*/
//  Compute the pressure from the temperature and density
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with pressure
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::pressure(FArrayBox& a_W,
                     const Box& a_box) const
{
  CH_assert(a_W.box().contains(a_box));
  const int tComp = temperatureIndex();
  const int pComp = pressureIndex();
  FORT_PRIMPRESSURE(CHF_FRA1(a_W, pComp),
                    CHF_CONST_FRA(a_W),
                    CHF_BOX(a_box),
                    CHF_CONST_INT(tComp),
                    CHF_CONST_REAL(CRDparam::g_R));
}

/*--------------------------------------------------------------------*/
//  Compute the density from the temperature and pressure
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with pressure
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::density(FArrayBox& a_W,
                    const Box& a_box) const
{
  CH_assert(a_W.box().contains(a_box));
  const int tComp = temperatureIndex();
  const int pComp = pressureIndex();
  const int rComp = densityIndex();
  FORT_PRIMDENSITY(CHF_FRA1(a_W, rComp),
                   CHF_CONST_FRA(a_W),
                   CHF_BOX(a_box),
                   CHF_CONST_INT(tComp),
                   CHF_CONST_INT(pComp),
                   CHF_CONST_REAL(CRDparam::g_R));
}

/*--------------------------------------------------------------------*/
//  Compute the artificial viscosity contribution to the flux
/** \param[in]  a_box   Cell box.  Flux needs to be computed on
 *                      surrounding faces.  Must be <= disjoint.
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[out] a_NtF   Flux due to artificial viscosity
 *  \param[in]  a_U     Solution
 *  \param[in]  a_N     Grid metrics on faces
 *  \param[in]  a_J     Metrics Jacobian in cells
 *  \param[in]  a_unitNormals
 *                      Unit normals on faces
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::artVisc(const Box&         a_box,
                    const BlockDomain& a_domain,
                    FluxBox&           a_NtF,
                    const FArrayBox&   a_U,
                    const FArrayBox&   a_WOld,
                    const FluxBox&     a_N,
                    const FArrayBox&   a_J,
                    const FluxBox&     a_unitNormals,
                    const RealVect&    a_dx,
                    LevelGridMetrics&  a_gridMetrics,
                    const Real         a_time,
                    const int          a_level) const
{
  CH_assert(a_NtF.box().contains(a_box));

  Box bx1 = grow(a_box, 1);
  Box bx1inDomain(bx1);
  bx1inDomain &= a_domain;
  CH_assert(a_U.box().contains(bx1inDomain));

  // Cell-centered boxes providing marking low and high side of the domain
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  // Get the primitive variables from the conserved variables (as needed).
  int numPrim = numPrimitive();

  FABSTACKTEMP(W, bx1inDomain, numPrim);
  consToPrim(W, a_U, bx1inDomain, a_WOld);
  extraPrimitiveState(W, bx1inDomain);
  FArrayBox vel(velocityInterval(), W);  // Alias

//--Precompute gradients of vel & U in all cells.  Note that the dyads have
//--numComp*SpaceDim for all gradient directions.  gradVel is needed for
//--divergence of the velocity and gradU for the artificial viscosity.

  const int numVelComp = SpaceDim;
  // Components stored as (velocity direction, gradient direction) in Fortran
  // ordering
  const int velCompStride = 1;
  // CHArray<Real, SpaceDim+1, ArRangeCol> gradVel(numVelComp*SpaceDim,
  //                                               bx1inDomain);
  CHARRAYNBSTACKTEMP(gradVel, Real, SpaceDim+1, ArRangeCol,
                     numVelComp*SpaceDim, bx1inDomain);
  gradVel = -1.;  // Not strictly required but avoids computation with
                  // uninitialized values

  const int numUComp = numConservative();
  // Components stored as (gradient direction, conserved variable) in Fortran
  // ordering
  const int UCompStride = SpaceDim;
  // CHArray<Real, SpaceDim+1, ArRangeCol> gradU(SpaceDim*numUComp,
  //                                             bx1inDomain);
  CHARRAYNBSTACKTEMP(gradU, Real, SpaceDim+1, ArRangeCol,
                     numUComp*SpaceDim, bx1inDomain);
  gradU = -1.;    // Not strictly required but avoids computation with
                  // uninitialized values

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // LHCbox needs to be grown by 1 in dir and not intersected with the
      // domain to satisfy the requirements of loHiCenter.
      // In the other directions, it needs to be grown by 1 and intersected with
      // the domain.  This is so we have a centered understanding of gradients
      // tangential to a face at the edge of a box (but not the edge of the
      // domain).  Because we only need tangential gradients outside a_box, the
      // velocity and U are still only required within 1 ghost cell of a_box.
      Box LHCbox(bx1);
      LHCbox.grow(dir, -1);
      LHCbox &= a_domain;
      LHCbox.grow(dir, 1);
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, LHCbox,
                 a_domain, dir);

      const int velCompBegin = dir*numVelComp;
      FORT_CELLGRADDIR(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradVel),
        CHF_CONST_FRA(vel),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi),
        CHF_BOX(centerBox),
        CHF_CONST_INT(dir),
        CHF_CONST_INT(numVelComp),
        CHF_CONST_INT(velCompBegin),
        CHF_CONST_INT(velCompStride),
        CHF_CONST_REAL(a_dx[0]));  //**FIXME Make vector

      const int UCompBegin = dir;
      FORT_CELLGRADDIR(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradU),
        CHF_CONST_FRA(a_U),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi),
        CHF_BOX(centerBox),
        CHF_CONST_INT(dir),
        CHF_CONST_INT(numUComp),
        CHF_CONST_INT(UCompBegin),
        CHF_CONST_INT(UCompStride),
        CHF_CONST_REAL(a_dx[0]));  //**FIXME Make vector
    }

//--Compute cell-centered c^2

  FABSTACKTEMP(c2, bx1inDomain, 1);
  calcGamma(bx1inDomain, c2, W);
  c2.mult(W, bx1inDomain, bulkModulusIndex(), 0, 1);
  c2.divide(W, bx1inDomain, densityIndex(), 0, 1);

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box faceBox0 = a_box;
      faceBox0.surroundingNodes(dir);

//--Store a contiguous N

      const int numNComp = SpaceDim*SpaceDim;
      const int zeroVal = 0;
      CHArray<Real, SpaceDim+1, ArRangeCol> Nctg(numNComp, faceBox0);
      FORT_REVERSEFABCOMPONENTSTRIDE(CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                                     CHF_CONST_FRA(a_N[dir]),
                                     CHF_BOX(faceBox0),
                                     CHF_CONST_INT(zeroVal),
                                     CHF_CONST_INT(zeroVal),
                                     CHF_CONST_INT(numNComp));

      FABSTACKTEMP(csqA, bx1inDomain, 1);
      FABSTACKTEMP(csqB, bx1inDomain, 1);
      csqA.copy(c2);

      FArrayBox* csqIptr = &csqA;
      FArrayBox* csqOptr = &csqB;

      // Compute min of csq in transverse direction(s).
      Box csqBox(bx1);
      for (int trDir = 0; trDir != SpaceDim; ++trDir)
        {
          if (trDir != dir)
            {
              FArrayBox& csqI = *csqIptr;
              FArrayBox& csqO = *csqOptr;
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         csqBox, a_domain, trDir);

              FORT_MIN3PTSF(CHF_FRA1(csqO, 0),
                            CHF_CONST_FRA1(csqI, 0),
                            CHF_CONST_INT(trDir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
              csqBox.grow(trDir, -1);
              FArrayBox* csqTmp = csqIptr;
              csqIptr = csqOptr;
              csqOptr = csqTmp;
            }
        }
      FArrayBox& csq = *csqIptr;

//--Face boxes for this direction

      Box bx1dir(a_box);
      bx1dir.grow(dir, 1);
      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, bx1dir,
                     a_domain, dir);

//--Compute the divergence of velocity

      FArrayBox divVel(faceBox0, 1);
      FORT_MAPPEDDIVVEL(CHF_FRA1(divVel, 0),
                        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradVel),
                        CHF_CONST_FRA(vel),
                        CHF_CONST_FRA1(a_J, 0),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox),
                        CHF_CONST_INT(dir),
                        CHF_CONST_REAL(a_dx[0]));  //**FIXME Make vector

//--Compute the physical cell spacing across the faces

      FArrayBox dxFace(faceBox0, 1);
      FORT_PHYSICALCELLSPACINGONFACE(CHF_FRA1(dxFace, 0),
                                     CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                                     CHF_CONST_FRA1(a_J, 0),
                                     CHF_BOX(loBox),
                                     CHF_CONST_INT(hasLo),
                                     CHF_BOX(hiBox),
                                     CHF_CONST_INT(hasHi),
                                     CHF_BOX(centerBox),
                                     CHF_CONST_INT(dir),
                                     CHF_CONST_REAL(a_dx[0]));  //**FIXME Make vector

//--Compute the flux due to artificial viscosity on the faces of the cells in
//--*computational* space (i.e., they have already been multiplied by a row of
//--N^T).  Only interior faces are affected.

      // Need to set boundary values to zero in case they are not modified
      if (hasLo)
        {
          a_NtF[dir].setVal(0., loBox, 0, a_NtF.nComp());
        }
      if (hasHi)
        {
          a_NtF[dir].setVal(0., hiBox, 0, a_NtF.nComp());
        }

      const int hasLoHiFalse = 0;
      const Real alpha = CRDparam::g_artificialViscosityCoef;
      const Real beta = CRDparam::g_artificialViscosity4thOCoef;
      FORT_MAPPEDARTVISC(CHF_FRA(a_NtF[dir]),
                         CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, Nctg),
                         CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, gradU),
                         CHF_CONST_FRA(a_U),
                         CHF_CONST_FRA1(divVel, 0),
                         CHF_CONST_FRA1(csq, 0),
                         CHF_CONST_FRA1(a_J, 0),
                         CHF_CONST_FRA1(dxFace, 0),
                         CHF_CONST_REAL(alpha),
                         CHF_CONST_REAL(beta),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLoHiFalse),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasLoHiFalse),
                         CHF_BOX(centerBox),
                         CHF_CONST_INT(dir),
                         CHF_CONST_REAL(a_dx[0]));  //**FIXME Make vector

//--Change fluxes due to artificial viscosity on the boundary faces

      CRDparam::g_CNSIBC->artViscBC(
        a_NtF[dir],
        Nctg,
        a_U,
        a_unitNormals[dir],
        divVel,
        csq,
        dxFace,
        velocityInterval(),  // Expected to be the same as momentum
        alpha,
        beta,
        loBox,
        hasLo,
        hiBox,
        hasHi,
        dir,
        a_box,
        a_gridMetrics,
        a_time,
        a_level);
    }
}

/*--------------------------------------------------------------------*/
//  Solve for the thermal conductivity and dynamic viscosity and returns
//  the maximum dynamic viscosity
/** \param[in]  a_box  Cell box.  Flux needs to be computed on
 *                     surrounding faces.
 *  \param[out] a_muFab
 *                     FAB of dynamic viscosity
 *  \param[out] a_kappaFab
 *                     FAB of thermal conductivity
 *  \param[in]  a_WfacePntFab
 *                     FAB containing the face-centered primitive variables
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::calcCoeffKappaMu(const Box&              a_box,
                             FArrayBox&              a_muFab,
                             FArrayBox&              a_kappaFab,
                             const FArrayBox&        a_WfacePntFab) const
{
  const int comp = 0;
  a_muFab.setVal(CRDparam::g_mu, a_box, comp);
  a_kappaFab.setVal(CRDparam::g_K, a_box, comp);
}

/*--------------------------------------------------------------------*/
//  Initialize the flow field
/** \param[out] a_U    Cell-centered conservative variables
 *  \param[in]  a_W    Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                     Current disjointBox
 *  \param[in]  a_box  Box to initialize over 
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::initialize(FArrayBox&              a_U,
                       const FArrayBox&        a_W,
                       const LevelGridMetrics& a_gridMetrics,
                       const FluxBox&          a_unitNormals,
                       const DataIndex&        a_didx,
                       const Box&              a_disjointBox,
                       const Box&              a_box) const
{
  const int rhoIndx = densityIndex();
  const int velIndx = velocityInterval().begin();
  const int presIndx = pressureIndex();
  const int tempIndx = temperatureIndex();
  const Real gamma = CRDparam::g_gamma;
  const Real Rval = CRDparam::g_R;
  FORT_CNSPHYSICSINIT(CHF_FRA(a_U),
                      CHF_CONST_FRA(a_W),
                      CHF_BOX(a_box),
                      CHF_CONST_INT(rhoIndx),
                      CHF_CONST_INT(velIndx),
                      CHF_CONST_INT(presIndx),
                      CHF_CONST_INT(tempIndx),
                      CHF_CONST_REAL(gamma),
                      CHF_CONST_REAL(Rval));
  //**FIXME this isn't always called.  See CNSIBCShockTube.cpp
  // Initialize the turbulent variables
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->turbInitialize(a_U,
                                  a_W,
                                  a_gridMetrics,
                                  a_unitNormals,
                                  a_didx,
                                  a_disjointBox,
                                  a_box);
    }
}

/*--------------------------------------------------------------------*/
//  Initialize wall-model state after restart
/** \param[out] a_U    Cell-centered conservative variables
 *  \param[in]  a_W    Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                     Current disjointBox
 *  \param[in]  a_box  Box to initialize over
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::initWallModelAfterRestart(
  FArrayBox&              a_JU,
  const FluxBox&          a_unitNormals,
  const LevelGridMetrics& a_gridMetrics,
  const DataIndex&        a_didx,
  const Box&              a_disjointBox) const
{
  // Get 2nd-order <U>
  FABSTACKTEMP(SecondOrderU, a_disjointBox, a_JU.nComp());
  const FArrayBox& JFab = a_gridMetrics.m_J[a_didx];
  for (int comp = 0; comp != a_JU.nComp(); ++comp)
    {
      MD_BOXLOOP(a_disjointBox, i)
        {
          SecondOrderU[MD_IX(i, comp)] = a_JU[MD_IX(i, comp)]/JFab[MD_IX(i, 0)];
        }
    }
  // Get 2nd-order <W>
  FABSTACKTEMP(SecondOrderW,
               a_disjointBox,
               CRDparam::g_CRDPhysics->numPrimitive());
  CNSPhysics::consToPrim(SecondOrderW,
                         SecondOrderU,
                         a_disjointBox,
                         SecondOrderW);
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->turbInitialize(SecondOrderU,
                                  SecondOrderW,
                                  a_gridMetrics,
                                  a_unitNormals,
                                  a_didx,
                                  a_disjointBox,
                                  a_disjointBox);
    }
  // Get 2nd-order <JU> from <U>
  const int turbCompBegin = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int numTurbComp = CRDparam::g_CRDPhysics->turbConsInterval().size();
  for (int comp = 0; comp != numTurbComp; ++comp)
    {
      const int turbComp = turbCompBegin + comp;
      MD_BOXLOOP(a_disjointBox, i)
        {
          a_JU[MD_IX(i, turbComp)] =
            SecondOrderU[MD_IX(i, turbComp)]*JFab[MD_IX(i, 0)];
        }
    }
}

/*--------------------------------------------------------------------*/
//  Puts the values for gamma in a_gamma
/** \param[in]  a_box  Box to put values for gamma in
 *  \param[out] a_gamma
 *                     FAB containing gamma values
 *  \param[in]  a_W    Primitive variables
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::calcGamma(const Box&       a_box,
                      FArrayBox&       a_gamma,
                      const FArrayBox& a_W) const
{
  a_gamma.setVal(CRDparam::g_gamma, a_box, 0);
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/**
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::speciesDiffHeatFlux(const Box&       a_box,
                                FArrayBox&       a_JnfacePntFab,
                                FArrayBox&       a_energyFab,
                                const FArrayBox& a_muFab,
                                const FArrayBox& a_kappaFab,
                                const FArrayBox& a_WfacePntFab) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::speciesCorrection(const Box&           a_box,
                              FArrayBox&           a_JU,
                              const ProblemDomain& a_domain,
                              const RealVect&      a_dx,
                              const Real           a_tolNeg) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

Real
CNSPhysics::speciesGasConstant(const int a_speciesNum) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
  return 1.;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

Real
CNSPhysics::speciesMolarMass(const int a_speciesNum) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
  return 1.;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

int
CNSPhysics::addReactionSource(const Box&       a_box,
                              FArrayBox&       a_RCTcellPntFab,
                              FArrayBox&       a_invDtFab,
                              const FArrayBox& a_WcellPntFab,
                              const Real       a_time,
                              const int        a_level,
                              Real&            a_minChemDt,
                              IntVect&         a_minChemDtCell) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
  return -1;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

int
CNSPhysics::ARSwithDiagonalFluxCorrection(
  const Box&       a_box,
  FArrayBox&       a_RCTcellPntFab,
  FArrayBox&       a_invDtFab,
  const FArrayBox& a_WcellPntFab,
  const FArrayBox& a_UcellPntFab,
  const Real       a_time,
  const int        a_level,
  Real&            a_minChemDt,
  IntVect&         a_minChemDtCell) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
  return -1;
}

/*--------------------------------------------------------------------*/
//  This shouldn't be used for CNS code
/** 
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::computeReactionJacobian(const Box& a_box,
                                    FArrayBox& a_rxnJacobianFab,
                                    const FArrayBox& a_WcellPntFab,
                                    const Real       a_dt) const
{
  CRD::msg << "CNSPhysics::CNSPhysics: Error using combustion only function "
           << "in CNS code " << CRD::error;
}

/*--------------------------------------------------------------------*/
//  Return the name of the primitive variable in the argument
//  component
/**
 *//*-----------------------------------------------------------------*/

const char *const
CNSPhysics::primStateName(const int a_iComp) const
{
  if (a_iComp < SpaceDim + 2)
    {
      static constexpr const char* NSname[] =
        {
          "density",
          D_DECL(
            "x-velocity",
            "y-velocity",
            "z-velocity"),
          "pressure"
        };
      return NSname[a_iComp];
    }
  if (turbConsInterval().contains(a_iComp))
    {
      const int relativeComp = a_iComp - turbConsInterval().begin();
      return m_turbModel->turbStateName(relativeComp);
    }
  if (transportConsInterval().contains(a_iComp))
    {
      return "transport";  //**FIXME at some point when transport is used
    }
  if (a_iComp == temperatureIndex())
    {
      return "temperature";
    }
  return nullptr;
}

/*--------------------------------------------------------------------*/
//  Return the name of the conservative variable in the argument
//  component
/**
 *//*-----------------------------------------------------------------*/

const char *const
CNSPhysics::consvStateName(const int a_iComp) const
{
  if (a_iComp < SpaceDim + 2)
    {
      static constexpr const char* NSname[] =
        {
          "density",
          D_DECL(
            "x-momentum",
            "y-momentum",
            "z-momentum"),
          "energy-density"
        };
      return NSname[a_iComp];
    }
  if (turbConsInterval().contains(a_iComp))
    {
      const int relativeComp = a_iComp - turbConsInterval().begin();
      return m_turbModel->turbStateName(relativeComp);
    }
  if (transportConsInterval().contains(a_iComp))
    {
      return "transport";  //**FIXME at some point when transport is used
    }
  return nullptr;
}

/*--------------------------------------------------------------------*/
//  Return name of the physics described by the class
/** \return             Name of physics
 *//*-----------------------------------------------------------------*/

const char *const
CNSPhysics::physicsName() const
{
  return "compressible Navier-Stokes";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the physics to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::writePhysicsInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CRD::msg << "Gas law\nideal" << CRD::var;
  CRD::msg << "Specific heat ratio\n" << CRDparam::g_gamma << CRD::var;
  CRD::msg << "Specific gas constant\n" << CRDparam::g_R  << " J/kg-K"
           << CRD::var;
  CRD::msg << "Gravitational constant\n" << CRDparam::g_grav << " m/s^2"
           << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Expressions for VisIt
/**
 *//*-----------------------------------------------------------------*/

#ifdef CH_USE_HDF5

void
CNSPhysics::expressions(HDF5HeaderData& a_holder) const
{
  std::string sgamma = "zonal_constant(Mesh,"
    + std::to_string(CRDparam::g_gamma) + ")";
  a_holder.m_string["scalar gamma"] = sgamma;
  std::string sR = "zonal_constant(Mesh,"
    + std::to_string(CRDparam::g_R) + ")";
  a_holder.m_string["scalar R"] = sR;
  a_holder.m_string["vector velocity"] = "momentum/density";
  a_holder.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_holder.m_string["scalar pressure"] =
    "(gamma-1)*(<energy-density>-kinetic_energy*density)";
  a_holder.m_string["scalar soundspeed"] = "sqrt(gamma*(pressure/density))";
  a_holder.m_string["scalar log10entropy"] =
    "log10(pressure) - gamma*log10(<density>)";
  a_holder.m_string["scalar machnumber"] =
    "magnitude(<velocity>)/<soundspeed>";
  a_holder.m_string["scalar temperature"] = "pressure/(R*density)";
  D_TERM(
    a_holder.m_string["scalar x-velocity"] = "<x-momentum>/density";,
    a_holder.m_string["scalar y-velocity"] = "<y-momentum>/density";,
    a_holder.m_string["scalar z-velocity"] = "<z-momentum>/density";)
  // visit can only calculate cartesian vorticity
  if (SpaceDim == 2)
    {
      a_holder.m_string["scalar cart_vorticity"] = "curl(velocity)";
    }
  else if (SpaceDim == 3)
    {
      a_holder.m_string["vector cart_vorticity"] = "curl(velocity)";
      a_holder.m_string["scalar cart_vorticity_mag"] = "magnitude(vorticityCart)";
    }
}

/*--------------------------------------------------------------------*/
//  Set the output LevelData
/** \param[out] a_outputLD The level data to be outputted
 *  \param[in]  a_U        Level data containing the conservative variables
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::outputLevelData(
  LevelData<FArrayBox>&       a_outputLD,
  const LevelData<FArrayBox>& a_U,
  const LevelData<FArrayBox>& a_WOld,
  const LevelGridMetrics&     a_levelGridMetrics) const
{
  const Interval outInterval(0, numOutputVar() - 1);
  a_U.copyTo(outInterval, a_outputLD, outInterval);
}

#endif

/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Solve the Riemann problem
/** \param[out] a_WStar Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_box   Face-centered box on which to compute a_WStar    
 *//*-----------------------------------------------------------------*/

template <typename T>
static inline T sqr(T x)
{
  return x*x;
}

template <typename T>
static inline T cube(T x)
{
  return x*x*x;
}

void
CNSPhysics::riemannSolution(FArrayBox&       a_WStar,
                            const FArrayBox& a_WLeft,
                            const FArrayBox& a_WRight,
                            const int&       a_dir,
                            const Box&       a_box) const
{
  CH_TIMELEAF("riemannSolution");
  const Real small = 1.E-6;
  const Real smallr = 1.4E-6;
  const Real smallp = 7.E-7;
  const Real gamma = CRDparam::g_gamma;
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  const int velIndx = velocityInterval().begin();
  D_TERM(
    const int uNormIndx = velIndx + a_dir;,
    const int uTan1Indx = velIndx + (a_dir+1)%SpaceDim;,
    const int uTan2Indx = velIndx + (a_dir+2)%SpaceDim;);

  bool haveScalars = false;
  int scalarBeg = 0;
  int scalarEnd = 0;
  if (numTurbVar() + numTransportVar() > 0)
    {
      haveScalars = true;
      scalarBeg = turbPrimInterval().begin();
      scalarEnd = transportPrimInterval().end() + 1;
    }

  constexpr bool doExact = false;
  if (doExact)
    {
      const Real gL   = gamma;
      const Real gR   = gamma;
      const Real gLm1 = gL - 1.;
      const Real gRm1 = gR - 1.;

      // double dml, dmr, vm, pm, aml, amr;
      // double msl, pml, dpmldum, msr, pmr, dpmrdum;
      // double vsl, vhl, vtl, vsr, vhr, vtr;

      MD_BOXLOOP(a_box, i)
        {
          // Extract left and right density
          const Real rhoL = std::max(smallr, a_WLeft [MD_IX(i, rhoIndx)]);
          const Real rhoR = std::max(smallr, a_WRight[MD_IX(i, rhoIndx)]);
          // Extract left and right pressure
          const Real pL = std::max(smallp, a_WLeft [MD_IX(i, presIndx)]);
          const Real pR = std::max(smallp, a_WRight[MD_IX(i, presIndx)]);
          // Extract left and right normal and tangential velocity
          D_TERM(const Real uL     = a_WLeft [MD_IX(i, uNormIndx)];,
                 const Real uLTan1 = a_WLeft [MD_IX(i, uTan1Indx)];,
                 const Real uLTan2 = a_WLeft [MD_IX(i, uTan2Indx)];)
          D_TERM(const Real uR     = a_WRight[MD_IX(i, uNormIndx)];,
                 const Real uRTan1 = a_WRight[MD_IX(i, uTan1Indx)];,
                 const Real uRTan2 = a_WRight[MD_IX(i, uTan2Indx)];)

          /* Determine the left and right state sound speeds. */

          const Real al = std::sqrt(gL*pL/rhoL);
          const Real ar = std::sqrt(gR*pR/rhoR);

          /* Compute the left and right state Riemann invariants. */

          const Real CL = uL + 2*al/gLm1;
          const Real CR = uR - 2*ar/gRm1;

          if (CL - CR <= 0.)
            {
              a_WStar[MD_IX(i, rhoIndx)]          = 0.;
              D_TERM(a_WStar[MD_IX(i, uNormIndx)] = 0.;,
                     a_WStar[MD_IX(i, uTan1Indx)] = 0.;,
                     a_WStar[MD_IX(i, uTan2Indx)] = 0.;)
              a_WStar[MD_IX(i, presIndx)]         = 0.;
              return;
            }

          /* Make an initial estimate of the intermediate state flow
             velocity to begin the Newton-Raphson iterative solution
             procedure.  The initial guess state velocity is made
             based on isentropic flow theory. */

          const Real Z  = (ar/al)*std::pow(pL/pR, 0.5*gLm1/gL);
          Real um       = (CL*Z + CR)/(1. + Z);

          /* In the case that two rarefaction waves are present,
             then an exact solution has been found and the iterative
             procedure is not required.  Check for this. */

          if (um >= uL && um <= uR)
            {
              Real dm, pm;
              if (um >= 0.)
                {
                  D_TERM((void)0;,
                         a_WStar[MD_IX(i, uTan1Indx)] = uLTan1;,
                         a_WStar[MD_IX(i, uTan2Indx)] = uLTan2;)
                  const Real aml = al - 0.5*gLm1*(um - uL);
                  pm             = pL*std::pow(aml/al, 2*gL/gLm1);
                  const Real vhl = uL - al;
                  const Real vtl = um - aml;
                  if (vhl >= 0.)
                    {
                      dm = rhoL;
                      um = uL;
                      pm = pL;
                    }
                  else if (vtl <= 0.)
                    {
                      dm = gL*pm/sqr(aml);
                    }
                  else
                    {
                      um = (gLm1*uL + 2*al)/(gL + 1.);
                      pm = pL*std::pow(um/al, 2*gL/gLm1);
                      dm = gL*pm/sqr(um);
                    }
                }
              else
                {
                  D_TERM((void)0;,
                         a_WStar[MD_IX(i, uTan1Indx)] = uRTan1;,
                         a_WStar[MD_IX(i, uTan2Indx)] = uRTan2;)
                  const Real amr = ar + 0.5*gRm1*(um - uR);
                  pm             = pR*std::pow(amr/ar, 2*gR/gRm1);
                  const Real vhr = uR + ar;
                  const Real vtr = um + amr;
                  if (vhr <= 0.)
                    {
                      dm = rhoR;
                      um = uR;
                      pm = pR;
                    }
                  else if (vtr >= 0.)
                    {
                      dm = gR*pm/sqr(amr);
                    }
                  else
                    {
                      um = (gRm1*uR - 2*ar)/(gR + 1.);
                      pm = pR*std::pow(-um/ar, 2*gR/gRm1);
                      dm = gR*pm/sqr(um);
                    }
                }
              a_WStar[MD_IX(i, rhoIndx)]   = dm;
              a_WStar[MD_IX(i, uNormIndx)] = um;
              a_WStar[MD_IX(i, presIndx)]  = pm;
              continue;
            }
    
          /* Perform the Newton-Raphson iterative procedure and solve for
             the velocity in the intermediate state.  During this iterative
             process the pressure in the intermediate state is also found. */

          int number_of_iterations = 0;
          constexpr int max_iterations = 10;
          constexpr Real TOLER = 1.E-7;

          Real msl, msr, pml, pmr, aml, amr;
          while (true)
            {
              /* Update the iteration counter. */
              number_of_iterations = number_of_iterations + 1;
 
              /* Determine solution changes for left wave. */
              Real dpmldum;
              if (um < uL)
                {
                  msl = (gL + 1.)*(um - uL)/(4*al);
                  msl = msl - std::sqrt(1. + sqr(msl));
                  pml = pL*(1. + gL*(um - uL)*msl/al);
                  dpmldum = 2*gL*pL*cube(msl)/(al*(1. + sqr(msl)));
                }
              else
                {
                  aml = al - 0.5*gLm1*(um - uL);
                  pml = pL*std::pow((aml/al), 2*gL/gLm1);
                  dpmldum = -gL*pml/aml;
                }
	  
              /* Determine solution changes for right wave. */
              Real dpmrdum;
              if ( um > uR )
                {
                  msr = (gR + 1.)*(um - uR)/(4*ar);
                  msr = msr + std::sqrt(1. + sqr(msr));
                  pmr = pR*(1. + gR*(um - uR)*msr/ar);
                  dpmrdum = 2*gR*pR*cube(msr)/(ar*(1. + sqr(msr)));
                }
              else
                {
                  amr = ar + 0.5*gRm1*(um - uR);
                  pmr = pR*std::pow((amr/ar), 2*gR/gRm1);
                  dpmrdum = gR*pmr/amr;
                }

              /* Check for convergence (i.e., pml=pmr). */

              if (std::abs(1. - pml/pmr) <= TOLER) break;
              if (number_of_iterations > max_iterations)
                {
                  pout() << "\n ERROR: convergence problem in Riemann solver: " 
                         << "n = " << number_of_iterations
                         << " tol = " << std::abs(1. - pml/pmr);
                  break;
                }

              /* Compute next estimate for the intermediate
                 state velocity, um. */

              um = um - (pml - pmr)/(dpmldum - dpmrdum);
            }

          Real dm;
          Real pm = 0.5*(pml + pmr);

          /* Return the intermediate state solution. */

          if (um >= 0.)
            {
              D_TERM((void)0;,
                     a_WStar[MD_IX(i, uTan1Indx)] = uLTan1;,
                     a_WStar[MD_IX(i, uTan2Indx)] = uLTan2;)
                if (um < uL)
                  {
                    aml = al*std::sqrt(((gL + 1.) + gLm1*pm/pL)/
                                       ((gL + 1.) + gLm1*pL/pm));
                    const Real vsl = uL + msl*al;
                    if (vsl >= 0.)
                      {
                        dm = rhoL;
                        um = uL;
                        pm = pL;
                      }
                    else
                      {
                        dm = gL*pm/sqr(aml);
                      }
                  }
                else  // um >= uL
                  {
                    const Real vhl = uL - al;
                    const Real vtl = um - aml;
                    if (vhl >= 0.)
                      {
                        dm = rhoL;
                        um = uL;
                        pm = pL;
                      }
                    else if (vtl <= 0.)
                      {
                        dm = gL*pm/sqr(aml);
                      }
                    else
                      {
                        um = (gLm1*uL + 2*al)/(gL + 1.);
                        pm = pL*std::pow(um/al, 2*gL/gLm1);
                        dm = gL*pm/sqr(um);
                      }
                  }
            }
          else  // um < 0.
            {
              D_TERM((void)0;,
                     a_WStar[MD_IX(i, uTan1Indx)] = uRTan1;,
                     a_WStar[MD_IX(i, uTan2Indx)] = uRTan2;)
                if (um > uR)
                  {
                    amr = ar*std::sqrt(((gR + 1.) + gRm1*pm/pR)/
                                       ((gR + 1.) + gRm1*pR/pm));
                    const Real vsr = uR + msr*ar;
                    if (vsr <= 0.)
                      {
                        dm = rhoR;
                        um = uR;
                        pm = pR;
                      }
                    else
                      {
                        dm = gR*pm/sqr(amr);
                      }
                  }
                else  // um <= uR
                  {
                    const Real vhr = uR + ar;
                    const Real vtr = um + amr;
                    if (vhr <= 0.)
                      {
                        dm = rhoR;
                        um = uR;
                        pm = pR;
                      }
                    else if (vtr >= 0.)
                      {
                        dm = gR*pm/sqr(amr);
                      }
                    else 
                      {
                        um = (gRm1*uR - 2*ar)/(gR + 1.);
                        pm = pR*std::pow(-um/ar, 2*gR/gRm1);
                        dm = gR*pm/sqr(um);
                      }
                  }
            }
          a_WStar[MD_IX(i, rhoIndx)]   = dm;
          a_WStar[MD_IX(i, uNormIndx)] = um;
          a_WStar[MD_IX(i, presIndx)]  = pm;
        }  // Loop over cells
    }

  else
    {

  MD_ARRAY_RESTRICT(arrWL, a_WLeft);
  MD_ARRAY_RESTRICT(arrWR, a_WRight);
  MD_ARRAY_RESTRICT(arrWStar, a_WStar);
  MD_BOXLOOP(a_box, i)
    {
      // Extract left and right density
      Real rhoL = std::max(smallr, arrWL[MD_IX(i, rhoIndx)]);
      Real rhoR = std::max(smallr, arrWR[MD_IX(i, rhoIndx)]);
      // Extract left and right pressure
      Real pL = std::max(smallp, arrWL[MD_IX(i, presIndx)]);
      Real pR = std::max(smallp, arrWR[MD_IX(i, presIndx)]);
      // Extract left and right normal velocity
      Real uL = arrWL[MD_IX(i, uNormIndx)];
      Real uR = arrWR[MD_IX(i, uNormIndx)];
      Real cL = std::sqrt(gamma*pL/rhoL);
      Real cR = std::sqrt(gamma*pR/rhoR);
      Real wL = rhoL*cL;
      Real wR = rhoR*cR;
      Real pStar = (wR*pL + wL*pR + wL*wR*(uL - uR))/(wL + wR);
      Real uStar = (wL*uL + wR*uR + pL - pR)/(wL + wR);
      Real rho0, p0, u0, c0, sgnm;
      if (uStar > 0.)
        {
          rho0 = rhoL;
          p0 = pL;
          D_TERM(u0 = uL;,
                 arrWStar[MD_IX(i, uTan1Indx)] = arrWL[MD_IX(i, uTan1Indx)];,
                 arrWStar[MD_IX(i, uTan2Indx)] = arrWL[MD_IX(i, uTan2Indx)];);
          if (haveScalars)
            {
              for (int sComp = scalarBeg; sComp != scalarEnd; ++sComp)
                {
                  arrWStar[MD_IX(i, sComp)] = arrWL[MD_IX(i, sComp)];
                }
            }
          c0 = cL;
          sgnm = 1.;
        }
      else
        {
          rho0 = rhoR;
          p0 = pR;
          D_TERM(u0 = uR;,
                 arrWStar[MD_IX(i, uTan1Indx)] = arrWR[MD_IX(i, uTan1Indx)];,
                 arrWStar[MD_IX(i, uTan2Indx)] = arrWR[MD_IX(i, uTan2Indx)];);
          if (haveScalars)
            {
              for (int sComp = scalarBeg; sComp != scalarEnd; ++sComp)
                {
                  arrWStar[MD_IX(i, sComp)] = arrWR[MD_IX(i, sComp)];
                }
            }
          c0 = cR;
          sgnm = -1.;
        }
      Real rStar = rho0 + (pStar - p0)/(c0*c0);
      rStar = std::max(rStar, smallr);
      Real cStar = std::sqrt(std::abs(gamma*pStar/rStar));
      Real wStar = 0.5*(cStar*rStar + c0*rho0);
      Real spout = c0 - sgnm*u0;
      Real spin = cStar - sgnm*uStar;
      Real uShock = wStar/rStar - sgnm*uStar;
      if (pStar > p0)
        {
          spout = uShock;
          spin = uShock;
        }

      Real frac = ((1. + (spout + spin)/std::max(spout-spin,small))/2.);
      frac = std::max(0., std::min(1., frac));
      arrWStar[MD_IX(i, rhoIndx)] = rho0 + frac*(rStar - rho0);
      arrWStar[MD_IX(i, presIndx)] = p0 + frac*(pStar - p0);
      arrWStar[MD_IX(i, uNormIndx)] = u0 + frac*(uStar - u0);

      if (spout <= 0.)
        {
          arrWStar[MD_IX(i, rhoIndx)] = std::max(rho0, smallr);
          arrWStar[MD_IX(i, presIndx)] = std::max(p0, smallp);
          arrWStar[MD_IX(i, uNormIndx)] = u0;
        }

      if (spin > 0.)
        {
          arrWStar[MD_IX(i, rhoIndx)] = std::max(rStar, smallr);
          arrWStar[MD_IX(i, presIndx)] = std::max(pStar, smallp);
          arrWStar[MD_IX(i, uNormIndx)] = uStar;
        }
    }
    }
}

/*--------------------------------------------------------------------*/
//  Solve for the elliptic time step
/** \param[in]  a_box   Box to find 1/dt_diffusive
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[out] a_invDtFab
 *                      Sum of inverse of time step sizes
 *  \param[in]  a_WFab  Primitive variables
 *  \param[in]  a_muFab Dynamic viscosity
 *  \param[in]  a_lambda
 *                      Constant value of lambda = 5.3333
 *  \param[in]  a_NTJ   Grid metrics on faces
 *  \param[in]  a_dxVect
 *                      RealVect of grid spacing
 *  \param[out] a_minDiffDt
 *                      Minimum diffusive time step
 *  \param[out] a_minDiffDtCell
 *                      Cell with minimum diffusive time step
 *//*-----------------------------------------------------------------*/

void
CNSPhysics::solveViscDt(const Box&       a_box,
                        const Box&       a_disjointBox,
                        FArrayBox&       a_invDtFab,
                        const FArrayBox& a_WFab,
                        const FArrayBox& a_muFab,
                        const Real&      a_lambda,
                        const FluxBox&   a_NTJ,
                        const RealVect&  a_dxVect,
                        Real&            a_minDiffDt,
                        IntVect&         a_minDiffDtCell) const
{
  // Solve for |N^T/J|^2, which is the magnitude of the matrix
  FABSTACKTEMP(NTColMagCell, a_box, SpaceDim);
  NTColMagCell.setVal(0.);
  {
    FLUXBOXSTACKTEMP(NTColMagFace, a_box, 1);
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        Box solveBox = surroundingNodes(a_box, dir);
        FArrayBox& NTColMag = NTColMagFace[dir];
        const FArrayBox& NTJFab = a_NTJ[dir];
        CH_assert(NTJFab.contains(solveBox));
        FORT_NTMAGDIFFDT(CHF_FRA1(NTColMag, 0),
                         CHF_CONST_FRA(NTJFab),
                         CHF_BOX(solveBox));
        FORT_EDGETOINCREMENTCELL(CHF_CONST_FRA1(NTColMag, 0),
                                 CHF_FRA1(NTColMagCell, dir),
                                 CHF_BOX(a_box),
                                 CHF_CONST_INT(dir));
      }
  }
  const Real factor = a_lambda*SpaceDim;
  MD_ARRAY_RESTRICT(arrW, a_WFab);
  MD_ARRAY_RESTRICT(arrMu, a_muFab);
  MD_ARRAY_RESTRICT(arrNTColMag, NTColMagCell);
  MD_ARRAY_RESTRICT(arrInvDt, a_invDtFab);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, densityIndex())];
      Real muVal = arrMu[MD_IX(i, 0)];
      Real NTMag = D_SELECT(arrNTColMag[MD_IX(i,0)],
                            std::max(arrNTColMag[MD_IX(i,0)],
                                     arrNTColMag[MD_IX(i,1)]),
                            std::max(arrNTColMag[MD_IX(i,0)],
                                     std::max(arrNTColMag[MD_IX(i,1)],
                                              arrNTColMag[MD_IX(i,2)])));
      Real diffDt = 2.5*a_dxVect[0]*a_dxVect[0]*rho/(muVal*factor*NTMag*NTMag);
      if (diffDt > 1.E-15)
        {
          Real invDiffDt = 1./diffDt;
          arrInvDt[MD_IX(i,0)] += invDiffDt;
          if (diffDt < a_minDiffDt)
            {
              a_minDiffDt = diffDt;
              a_minDiffDtCell = MD_GETIV(i);
            }
        }
    }
}

