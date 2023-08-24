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

#ifndef _NONLINEARSOVLER_H_
#define _NONLINEARSOVLER_H_

//----- Chombo -----//
#include "FArrayBox.H"
#include "LevelFluxRegister.H"
#include "LevelGridMetrics.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"

//----- Internal -----//
class LevelCNSOp;

class NonLinearSolverI
{
public:
  typedef LevelData<FArrayBox> RHS;
  typedef LevelData<FArrayBox> SOLN;
  //! Default constructor
  NonLinearSolverI() {};

  //! Destructor
  virtual ~NonLinearSolverI() noexcept {};

  //! Solves a nonlinear problem for ARK
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
                     SOLN&              a_WOld) = 0;

protected:
private:
  //! Copy constructor forbidden
  NonLinearSolverI(const NonLinearSolverI &other) = delete;

  //! Move constructor forbidden
  NonLinearSolverI(NonLinearSolverI &&other) noexcept = delete;

  //! Copy assignment operator forbidden
  NonLinearSolverI& operator=(const NonLinearSolverI &other) = delete;

  //! Move assignment operator forbidden
  NonLinearSolverI& operator=(NonLinearSolverI &&other) noexcept = delete;
};

using FABVector = FArrayBox;
using FABMatrix = FArrayBox;

//! An interface to a Jacobian computation class
class JacobianFunctionI
{
public:
  //! Default constructor forbidden
  JacobianFunctionI() = default;
  //! Default destructor forbidden
  virtual ~JacobianFunctionI() = default;
  //! Default copy constructor forbidden
  JacobianFunctionI(const JacobianFunctionI&) = delete;
  //! Default move constructor forbidden
  JacobianFunctionI(JacobianFunctionI&&) = delete;

  //! Computes the Jacobian given a solution
  virtual void operator()(const FABVector& a_currentSolution,
                          FABMatrix&       a_matrix,
                          const Box&       a_box) const = 0;
};

class RHSFunctionI
{
public:
        RHSFunctionI() = default;
        virtual ~RHSFunctionI() = default;
        RHSFunctionI(const RHSFunctionI&) = delete;
        RHSFunctionI(RHSFunctionI&&) = delete;

        virtual void operator()(const FABVector& a_currentSolution,
                                FABVector&       a_result,
                                const Box&         a_box) const = 0;
};

//! An interface to a solving a linear problem
class LinearSolverI
{
public:
  //! Default constructor forbidden
  LinearSolverI() = default;
  //! Default destructor forbidden
  virtual ~LinearSolverI() = default;
  //! Default copy constructor forbidden
  LinearSolverI(const LinearSolverI&) = delete;
  //! Default move constructor forbidden
  LinearSolverI(LinearSolverI&&) = delete;

  //! Solves a linear problem of the form Ax=b
  virtual void solve(const FArrayBox& a_matrix,
                     FArrayBox&       a_result,
                     FArrayBox&       a_residual,
                     const FArrayBox& a_rhs,
                     const Box&       a_box) const = 0;
};

//! Implementation of a linear solver using Lapack
class LinearSolverLapack : public LinearSolverI
{
public:
  //! Solves a linear problem of the form Ax=b
  void solve(const FArrayBox& a_matrix,
             FArrayBox&       a_result,
             FArrayBox&       a_residual,
             const FArrayBox& a_rhs,
             const Box&       a_box) const override;
};

//! A struct that data to compute the step size
struct StepSizeData
{
  Real m_currentNorm;               ///< The current iteration's norm
  Real m_previousNorm;              ///< The previous iteration's norm
};

//! A class to compute the step size
class StepSizeComputerI
{
public:
  //! Default constructor forbidden
  StepSizeComputerI() = default;
  //! Default destructor forbidden
  virtual ~StepSizeComputerI() = default;
  //! Default copy constructir forbidden
  StepSizeComputerI(const StepSizeComputerI&) = delete;
  //! Default move constructor forbidden
  StepSizeComputerI(StepSizeComputerI&&) = delete;

  //! Computes the nonlinear step size
  virtual Real computeStepSize(const StepSizeData& a_data) const = 0;
};

//! An implementation of the step size commpute that just returns 1
class StepSizeComputerUnit : public StepSizeComputerI
{
public:
  //! Computes the step size (but its just one)
  Real computeStepSize(const StepSizeData& a_data) const override
  { return 1.; }
};

//! A struct to contain data needed to compute convergence conditions
struct ConvergenceData
{
  FArrayBox *m_residual;              ///< The residual of the nonlinear problem
  FArrayBox *m_displacement;          ///< The displacement from nonlinear
                                      ///< problem
  FArrayBox *m_prevSolution;          ///< The previous solution values
  Real m_refScale;                    ///< Reference scale of the species mass
                                      ///< J*rho*c_n (used for scaling
                                      ///< tolerances)
  Real m_currentNorm;                 ///< The current norm
};

//! An interface class for comuting convergence conditions
class ConvergenceCheckerI
{
public:
  //! Default constructor forbidden
  ConvergenceCheckerI() = default;
  //! Default destructor forbidden
  virtual ~ConvergenceCheckerI() = default;
  //! Default copy constructor forbidden
  ConvergenceCheckerI(const ConvergenceCheckerI&) = delete;
  //! Default move constructor forbidden
  ConvergenceCheckerI(ConvergenceCheckerI&&) = delete;

  //! Computes convergence of the nonlinear rpoblem
  virtual bool converged(const ConvergenceData& a_data) const = 0;
};

//! An implementation of the nonlinear solver using the Newton method
class NonLinearSolverNewton
{
public:
  //! Sets the maximum number of Newton iterations
  void setMaxIters(int a_maxIters)
  { m_maxIters = a_maxIters; };

  //! Sets the specific Jacobian implementation to use
  void setJacobianFunction(const JacobianFunctionI& a_jacFunc)
  { m_jacFunc = &a_jacFunc; };

  //! Sets the specific RHS implementation to use
  void setRHSFunction(const RHSFunctionI& a_rhsFunc)
  { m_rhsFunc = &a_rhsFunc; };

  //! Sets the specific linear solver impelementation to use
  void setLinearSolver(const LinearSolverI& a_linearSolver)
  { m_linearSolver = &a_linearSolver; };

  //! Sets specific step size method to use
  void setStepSizeComputer(const StepSizeComputerI& a_stepSizeComputer)
  { m_stepSizeComputer = &a_stepSizeComputer; };

  //! Sets specific convergence checking method to use
  void setConvergenceChecker(const ConvergenceCheckerI& a_convergenceChecker)
  { m_convergenceChecker = &a_convergenceChecker; };

  //! Solves the nonlinear problem
  void solve(FABVector& a_solution,
             const Box& a_box);

private:
  Real m_maxIters;                   ///< Maximum number of Newton iterations to perform
  const JacobianFunctionI* m_jacFunc;
                                     ///< The Jacobian computation class to use
  const RHSFunctionI* m_rhsFunc;     ///< The RHS computation class to use
  const LinearSolverI* m_linearSolver;
                                     ///< The linear solver class to use
  const StepSizeComputerI* m_stepSizeComputer;
                                     ///< The step size computation class to use
  const ConvergenceCheckerI* m_convergenceChecker;
                                     ///< The convergence checking computation class to use
};


// Useful free functions

//! Converts an FArrayBox to a CHVector
void FabToMatrix(CHMatrix&        a_matrix,
                 const FArrayBox& a_matrixFAB,
                 const IntVect&   a_cell,
                 int              a_numRows,
                 int              a_numCols);

//! Converts a CHMatrix to an FArrayBox
void MatrixToFab(FArrayBox&      a_matrixFAB,
                 const CHMatrix& a_matrix,
                 const IntVect&  a_cell,
                 int             a_numRows,
                 int             a_numCols);

//! Converts an FArrayBox to a CHVector
void FabToVector(CHVector&        a_vector,
                 const FArrayBox& a_vectorFAB,
                 const IntVect&   a_cell,
                 int              a_numRows);

//! Converts a CHVector to an FArrayBox
void VectorToFab(FArrayBox&      a_vectorFAB,
                 const CHVector& a_vector,
                 const IntVect&  a_cell,
                 int             a_numRows);

#endif // _NONLINEARSOLVER_H_
