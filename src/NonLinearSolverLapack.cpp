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
 * \file NonLinearSolverLapack.cpp
 *
 * \brief Member functions for NonLinearSolverLapack
 *
 *//*+*************************************************************************/

#include "NonLinearSolverLapack.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"
#include "CRDPhysics.H"
#include "ARKUtil.H"
#include "NonLinearSolver.H"

NonLinearSolverLapack::NonLinearSolverLapack()
{}

/*--------------------------------------------------------------------*/
//  Solves the nonlinear problem that arises from doing an implicit ARK term
/**
 *  \param[out]  a_newSoln
 *                      Solution of stage
 *  \param[in]   a_prevStageSoln
 *                      Solution calculated at the previous stage
 *  \param[in]   a_prevTimeSoln
 *                      Solution calculated at the previous time step
 *  \param[in]   a_rhs  Contains the rhs of the Ax=b problem (so... its b)
 *  \param[in]   a_stage
 *                      Stage index [0-5] for ARK4
 *  \param[in]   a_time
 *                      Time at start of stage
 *  \param[in]   a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]   a_fnFluxRegister
 *                      Flux register between this and a finer level
 *  \param[in]   a_crFluxRegister
 *                      Flux register between this and a coarser level
 *  \param[in]   a_dt   Time step on this level.  This is for the
 *                      complete time step, not just the stage.
 *  \param[in]   a_timeOld
 *                      Time at start of time step on this level
 *  \param[in]   a_timeCoarseOld
 *                      Time at start of time step on coarser level
 *  \param[in]   a_timeCoarseNew
 *                      Time at end of time step on coarser level
 *  \param[in]   a_subcycleParams
 *                      Parameters from subcyling algorithm (see
 *                      AMRLevel.H for more info)
 *  \param[in]   a_levelCNSOp
 *                      Operators for computing stiff and nonstiff terms
 *  \param[in/out]
 *               a_WOld Primitive solution from previous stage
 *//*-----------------------------------------------------------------*/
void
NonLinearSolverLapack::solve(SOLN&              a_newSoln,
                             const SOLN&        a_prevStageSoln,
                             const SOLN&        a_prevTimeSoln,
                             const RHS&         a_rhs,
                             int                a_stage,
                             Real               a_time,
                             Real               a_stageWeight,
                             LevelFluxRegister& a_finerFluxRegister,
                             LevelFluxRegister& a_coarserFluxRegister,
                             const Real         a_dt,
                             const Real         a_timeOld,
                             const Real         a_timeCoarseOld,
                             const Real         a_timeCoarseNew,
                             SubcycleParams     a_subcycleParams,
                             LevelCNSOp*        a_levelCNSOp,
                             SOLN&              a_WOld)
{
  CH_TIME("NonLinearSolverLapack::solve");
  CRD::msg << CRD::fv3 << "NonLinearSolverLapack::solve" << CRD::end;
  CH_assert(a_newSoln.getBoxes() == a_prevStageSoln.getBoxes());
  CH_assert(a_prevTimeSoln.getBoxes() == a_prevStageSoln.getBoxes());
  CH_assert(a_rhs.getBoxes() == a_newSoln.getBoxes());
  CH_assert(a_stage >= 0 && a_stage < 6);

  const int densityIdx = CRDparam::g_CRDPhysics->densityIndex();
  const int energyIdx = CRDparam::g_CRDPhysics->energyFluxIndex();
  const Interval& consSpeciesIval = CRDparam::g_CRDPhysics->speciesConsInterval();

  const DisjointBoxLayout& dbl = a_newSoln.disjointBoxLayout();

  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = dbl[dit];
      FArrayBox& JUNew = a_newSoln[dit];
      const FArrayBox& JinitGuess = a_prevStageSoln[dit];
      const FArrayBox& JUOld = a_prevTimeSoln[dit];
      const FArrayBox& JRHS = a_rhs[dit];
      FArrayBox& WOldFab = a_WOld[dit];

      // Set our intiial guess
      Interval NSEComponents(densityIdx,energyIdx);
      // For the continuity, momentum, and energy equations the final solution can
      // be calculated exactly using the explicit update portion stored in JRHS
      JUNew.copy(JUOld, NSEComponents.begin(), NSEComponents.begin(),  NSEComponents.size());
      JUNew.plus(JRHS, NSEComponents.begin(), NSEComponents.begin(), NSEComponents.size());

      // For the species transport equations, the initial gues is set from
      // JinitGuess which is set by the configuration of ARK4
      JUNew.copy(JinitGuess, consSpeciesIval.begin(), consSpeciesIval.begin(), consSpeciesIval.size());
      // Avoid negative species in the initial guess
      CRDparam::g_CRDPhysics->speciesCorrection(disjointBox,
                                                JUNew,
                                                ProblemDomain{},  // Dummy
                                                RealVect_unit,    // Dummy
                                                0.);

      MD_BOXLOOP(disjointBox, i)
        {
          IntVect cell(D_DECL(i0,i1,i2));
          nonlinearSolve(JUNew,
                         JUOld,
                         JinitGuess,
                         JRHS,
                         cell,
                         dit,
                         a_timeOld,
                         a_dt,
                         a_levelCNSOp,
                         WOldFab);
        }
    }

}

/*--------------------------------------------------------------------*/
//  Solves the nonlinear problem that arises from doing an implicit ARK term
//  for a single cell
/**
 *  \param[out]  a_JUNew
 *                      New solution computed from the nonlinear solve
 *  \param[in]   a_JUOld
 *                      Previous time step solution
 *  \param[in]   a_JinitGUess
 *                      Initial guess to nonlinear solver
 *  \param[in]   a_Jrhs Explicit term computed from previous stages
 *  \param[in]   a_cell Cell which to do the nonlinear solve of
 *  \param[in]   a_dit  DataIterator corresponding to the working cell
 *  \param[in]   a_time Current solution time
 *  \param[in]   a_dt   Time step on this level.  This is for the
 *                      complete time step, not just the stage.
 *  \param[in]   a_levelCNSOp
 *                      Operators for computing stiff and nonstiff terms
 *  \param[in/out]
 *               a_WOld Primitive solution from previous stage
 *//*-----------------------------------------------------------------*/
void
NonLinearSolverLapack::nonlinearSolve(FArrayBox&          a_JUNew,
                                      const FArrayBox&    a_JUOld,
                                      const FArrayBox&    a_JinitGuess,
                                      const FArrayBox&    a_Jrhs,
                                      const IntVect&      a_cell,
                                      const DataIterator& a_dit,
                                      Real                a_time,
                                      Real                a_dt,
                                      LevelCNSOp*         a_levelCNSOp,
                                      FArrayBox&          a_WOld) const
{
  CH_TIME("NonLinearSolverLapack::nonlinearSolve");
  Box cellBox = Box(a_cell, a_cell);
  //  Set up the various components of our nonlinear solver
  ARKData arkData;
  arkData.m_levelCNSOp = a_levelCNSOp;
  arkData.m_cellBox = cellBox;
  arkData.m_time = a_time;
  arkData.m_dt = a_dt;
  arkData.m_gamma = 1./4.;
  arkData.m_tolerance = CRDparam::g_ARKNonlinearTol;
  arkData.m_dit = &a_dit;
  arkData.m_WOld = &a_WOld;
  arkData.m_explicitTerm = &a_Jrhs;
  arkData.m_prevTimeSolution = &a_JUOld;
  ARKJacobian jacFunc(arkData);
  ARKRHS RHSFunc(arkData);
  LinearSolverLapack linearSolver;
  StepSizeComputerUnit stepSizeComputer;
  ARKConvergenceChecker arkConvergenceChecker(arkData);

  NonLinearSolverNewton solver;
  solver.setMaxIters(40);
  solver.setJacobianFunction(jacFunc);
  solver.setRHSFunction(RHSFunc);
  solver.setLinearSolver(linearSolver);
  solver.setStepSizeComputer(stepSizeComputer);
  solver.setConvergenceChecker(arkConvergenceChecker);

  solver.solve(a_JUNew, cellBox);
}

/*--------------------------------------------------------------------*/
//  Constructor for the Jacobian computation
/**
 *  \param[in]   a_arkData
 *                      ARK solver data needed to compute the Jacobian
 *//*-----------------------------------------------------------------*/
ARKJacobian::ARKJacobian(const ARKData& a_arkData)
{
  m_arkData = &a_arkData;
}

/*--------------------------------------------------------------------*/
//  Computes the ARK Jacobian (I-dt*gamma*dS/dU) given the curent JU solution
/**
 *  \param[in]   a_currentSolution
 *                      The current JU solution
 *  \param[out]  a_matrix
 *                      The computed Jacobian
 *  \param[in]   a_box  The box over which to compute the Jacobian
 *//*-----------------------------------------------------------------*/
void
ARKJacobian::operator()(const FArrayBox& a_currentSolution,
                        FArrayBox&       a_matrix,
                        const Box&       a_box) const
{
  CH_TIME("ARKJacobian::operator()");

  const int numComps = CRDparam::g_numSpecies;
  LevelCNSOp* levelCNSOp = m_arkData->m_levelCNSOp;
  const Box& cellBox = m_arkData->m_cellBox;
  Real time = m_arkData->m_time;
  Real dt = m_arkData->m_dt;
  Real gamma = m_arkData->m_gamma;
  const DataIterator* dit = m_arkData->m_dit;
  FArrayBox* WOld = m_arkData->m_WOld;
  CH_assert(cellBox.smallEnd() == cellBox.bigEnd());
  levelCNSOp->calcRxnJacobian(a_matrix,
                              a_currentSolution,
                              cellBox,
                              time,
                              dt,
                              *dit,
                              *WOld);

  IntVect cell = cellBox.smallEnd();
  for(int row = 0; row != numComps;  ++row)
    {
      for (int col = 0; col != numComps; ++col)
        {
          int linIndx = row*numComps+col;
          a_matrix(cell, linIndx) *= -dt*gamma;
          if(row == col)
            {
              a_matrix(cell, linIndx) += 1.0;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Constructor for the computatin of the RHS
/**
 *  \param[in]   a_arkData
 *                      ARK solver data used to copmute the RHS
 *//*-----------------------------------------------------------------*/
ARKRHS::ARKRHS(const ARKData& a_arkData)
{
  m_arkData = &a_arkData;
}

/*--------------------------------------------------------------------*/
//  Computes the RHS for the linear solve
/**
 *  \param[in]   a_currentSolution
 *                      The current JU solution
 *  \param[out]  a_result
 *                      The RHS
 *  \param[in]   a_box  The box over which to compute the RHS
 *//*-----------------------------------------------------------------*/
void
ARKRHS::operator()(const FArrayBox& a_currentSolution,
                   FArrayBox&       a_result,
                   const Box&       a_box) const
{
  CH_TIME("ARKRHS::operator()");
  const int numConservative = a_currentSolution.nComp();
  LevelCNSOp* levelCNSOp = m_arkData->m_levelCNSOp;
  const Box& cellBox = m_arkData->m_cellBox;
  Real time = m_arkData->m_time;
  Real dt = m_arkData->m_dt;
  Real gamma = m_arkData->m_gamma;
  const DataIterator* dit = m_arkData->m_dit;
  FArrayBox* WOld = m_arkData->m_WOld;
  const FArrayBox* explicitTerm = m_arkData->m_explicitTerm;
  const FArrayBox* prevTimeSolution = m_arkData->m_prevTimeSolution;
  CH_assert(cellBox.smallEnd() == cellBox.bigEnd());
  levelCNSOp->evalStiff(a_result,
                        a_currentSolution,
                        cellBox,
                        time,
                        *dit,
                        *WOld);

  const IntVect cell = cellBox.smallEnd();

  // Set up RHS
  for(int row = 0; row != numConservative;  ++row)
    {
      a_result(cell, row) *= dt*gamma;
      a_result(cell, row) += -(a_currentSolution(cell, row)
                               - prevTimeSolution->operator()(cell, row))
        + explicitTerm->operator()(cell, row);
    }
}

/*--------------------------------------------------------------------*/
//  Constructor the ARK convergence checker class
/**
 *  \param[in]   a_arkData
 *                      ARK solver data needed for checking convergence
 *//*-----------------------------------------------------------------*/
ARKConvergenceChecker::ARKConvergenceChecker(const ARKData& a_arkData)
{
  m_arkData = &a_arkData;
}

/*--------------------------------------------------------------------*/
// Checks the convergence of the nonlinear solver
/**
 *  \param[in]   a_convData
 *                      Current nonlinear solver iteration data
 *//*-----------------------------------------------------------------*/

bool
ARKConvergenceChecker::converged(const ConvergenceData& a_convData) const
{
  const FArrayBox* residual = a_convData.m_residual;
  const FArrayBox* displacement = a_convData.m_displacement;
  const Real ref1Tol = m_arkData->m_tolerance;
  const Box& cellBox = m_arkData->m_cellBox;
  const Interval dispIval = displacement->interval();
  const Interval resIval = residual->interval();

  // Use L-infinity norms to test for convergence.
  Real displacementL2Norm = displacement->norm(cellBox,
                                               0,
                                               dispIval.begin(),
                                               dispIval.size());

  Real residualL2Norm = residual->norm(cellBox,
                                       0,
                                       resIval.begin(),
                                       resIval.size());

  const Real scaledTol = a_convData.m_refScale*ref1Tol;
  return (residualL2Norm < scaledTol || displacementL2Norm < scaledTol);
}

/*--------------------------------------------------------------------*/
//  Divides the dividend by the divisor and stores the result in the quotient.
//  If the divisor has a zero value, then division is not performed in that cell
//  and the quotient for that cell is just set to the dividend.
/** \param[out] a_quotient
 *                      Result of the division
 *  \param[in]  a_dividend
 *                      The dividend of the division
 *  \param[in]  a_divisor
 *                      The divisor of the division that may contain zero values
 *  \param[in]  a_box   Box over which the division should be performed
 *  \param[in]  a_srcComp
 *                      Starting component of the dividend and divisor
 *  \param[in]  a_destComp
 *                      Starting component of the quotient
 *  \param[in]  a_numComp
 *                      Number of components to do the division on
 *  \param[in]  a_epsilon
 *                      How close to zero the divisor must be to not performe division
 *//*-----------------------------------------------------------------*/
void
ARKConvergenceChecker::divideFABsContainingZeros(FArrayBox&       a_quotient,
                                                 const FArrayBox& a_dividend,
                                                 const FArrayBox& a_divisor,
                                                 const Box&       a_box,
                                                 int              a_srcComp,
                                                 int              a_destComp,
                                                 int              a_numComp,
                                                 Real             a_epsilon) const
{
  MD_ARRAY_RESTRICT(arrQuotient,  a_quotient);
  MD_ARRAY_RESTRICT(arrDividend,  a_dividend);
  MD_ARRAY_RESTRICT(arrDivisor,  a_divisor);
  for(int comp = 0; comp != a_numComp; ++comp)
    {
      const int srcComp = comp + a_srcComp;
      const int destComp = comp + a_destComp;
      MD_BOXLOOP(a_box, cell)
        {
          Real absDivisor = Abs(a_divisor[MD_IX(cell, srcComp)]);
          if(absDivisor >= a_epsilon)
            {
              a_quotient[MD_IX(cell,destComp)] = Abs(a_dividend[MD_IX(cell, srcComp)])/absDivisor;
            }
          else
            {
              a_quotient[MD_IX(cell, destComp)] = 0.;
            }
        }
    }
}
