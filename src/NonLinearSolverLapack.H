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
 * \file NonLinearSolver.H
 *
 * \brief Nonlinear solver interface for use with semi-implicit time marching
 *        methods
 *
 *//*+*************************************************************************/

#ifndef _NONLINEARSOVLERLAPACK_H_
#define _NONLINEARSOVLERLAPACK_H_

//----- Internal -----//
#include "NonLinearSolver.H"
#include "LevelCNSOp.H"

//----- Chombo -----//
#include "LevelFluxRegister.H"
#include "LevelGridMetrics.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"

class NonLinearSolverLapack : public NonLinearSolverI
{
public:
  typedef LevelData<FArrayBox> RHS;
  typedef LevelData<FArrayBox> SOLN;

  //! Default constructor
  NonLinearSolverLapack();

  //! Destructor
  virtual ~NonLinearSolverLapack() noexcept {};

  /// Solve the nonlinear problem that arises from doing the implicit ARK term
  virtual void solve(SOLN&              a_newSoln,
                     const SOLN&        a_prevStageSoln,
                     const SOLN&        a_prevTimeSoln,
                     const RHS&         a_rhs,
                     int                a_stage,
                     Real               a_time,
                     Real               a_stageweight,
                     LevelFluxRegister& a_finerFluxRegister,
                     LevelFluxRegister& a_coarserFluxRegister,
                     const Real         a_dt,
                     const Real         a_timeOld,
                     const Real         a_timeCoarseOld,
                     const Real         a_timeCoarseNew,
                     SubcycleParams     a_subcycleParams,
                     LevelCNSOp*        a_levelCNSOp,
                     SOLN&              a_WOld) override;


protected:

  /// Performs the nonlinear solve on a single cell
  virtual void nonlinearSolve(FArrayBox&          a_JUNew,
                              const FArrayBox&    a_JUOld,
                              const FArrayBox&    a_JinitGuess,
                              const FArrayBox&    a_Jrhs,
                              const IntVect&      a_cell,
                              const DataIterator& a_dit,
                              Real                a_time,
                              Real                a_dt,
                              LevelCNSOp*         a_levelCNSOp,
                              FArrayBox&          a_WOld) const;

private:
  //! Copy constructor forbidden
  NonLinearSolverLapack(const NonLinearSolverLapack &other) = delete;

  //! Move constructor forbidden
  NonLinearSolverLapack(NonLinearSolverLapack &&other) noexcept = delete;

  //! Copy assignment operator forbidden
  NonLinearSolverLapack& operator=(const NonLinearSolverLapack &other) = delete;

  //! Move assignment operator forbidden
  NonLinearSolverLapack& operator=(NonLinearSolverLapack &&other) noexcept = delete;
};

struct ARKData
{
  LevelCNSOp*   m_levelCNSOp;         ///< Operator for computing the stiff and nonstiff terms
  Box           m_cellBox;            ///< Box containing only the current cell
  Real          m_time;               ///< Current solution time
  Real          m_dt;                 ///< Time step size
  Real          m_gamma;              ///< ARK gamma from the Butcher Tableau
  Real          m_tolerance;          ///< ARK nonlinear convergence tolerance
  const DataIterator* m_dit;          ///< DataIterator corresponding to the current cell
  FArrayBox*    m_WOld;               ///< Previous stage primitive values
  const FArrayBox*    m_explicitTerm; ///< Explicitly computed term from previous stages
  const FArrayBox*    m_prevTimeSolution;
                                      ///< The previous time step's solution

};

class ARKJacobian : public JacobianFunctionI
{
public:
  /// Constructor, requires passing some ARK solver data
  ARKJacobian(const ARKData& a_arkData);

  /// Computes the Jacobian of ARK's linear solver
  void operator()(const FArrayBox&  a_currentSolution,
                  FArrayBox&        a_matrix,
                  const Box&        a_box) const override;
private:
  const ARKData* m_arkData;         ///< ARK solver data used to compute the Jacobian
};

class ARKRHS : public RHSFunctionI
{
public:
  /// Constructor, requires passing some ARK solver data
  ARKRHS(const ARKData& a_arkData);

  /// Computes the RHS of ARK's linear solver
  void operator()(const FABVector& a_currentSolution,
                  FABMatrix&       a_result,
                  const Box&       a_box) const override;
private:
  const ARKData* m_arkData;          ///< ARK solver data used to compute the RHS
};

class ARKConvergenceChecker : public ConvergenceCheckerI
{
public:
  /// Constructor, requires passing some ARK solver data
  ARKConvergenceChecker(const ARKData& a_arkData);

  /// Checks for convergence of the nonlinear ARK solver
  bool converged(const ConvergenceData& a_convData) const override;
private:

  /// Does division where some values of the divisor might be zero
  void divideFABsContainingZeros(FArrayBox&       a_quotient,
                                 const FArrayBox& a_dividend,
                                 const FArrayBox& a_divisor,
                                 const Box&       a_box,
                                 int              a_srcComp = 0,
                                 int              a_destComp = 0,
                                 int              a_numComp = 1,
                                 Real             a_epsilon = 1.e-14) const;

  const ARKData* m_arkData;           ///< ARK solver data used to check convergence
};

#endif // _NONLINEARSOLVERLAPACK_H_
