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
 * \file NonLinearSolver.cpp
 *
 * \brief A generic nonlinear solver for Chord.
 *
 *//*+*************************************************************************/

#include "NonLinearSolver.H"
#include <iostream>

//----- Chombo -----//

#include "CHMatrixOps.H"

//----- Internal -----//

#include "CRDmsg.H"
#include "CRDparam.H"
#include "ThermPhysics.H"
#include "DataTemp.H"

/*--------------------------------------------------------------------*/
//  Solves the nonlinear problem
/** \param[in/out]
 *              a_solution
 *                      Input is the initial guess, transforms into the solution
 *  \param[in]  a_box   Box over which the nonlinear problem should be solved
 *//*-----------------------------------------------------------------*/

void NonLinearSolverNewton::solve(FABVector& a_solution,
                                  const Box& a_box)
{
  CH_TIME("NonLinearSolverNewton::solve");
  // The fraction of species which are considered not relevant.  Consider a
  // sum of N species (sorted small to large) sumcj = \Sum_{j=1}^N cj.  The
  // largest value k where \Sum_{j=1}^k cj < c_cjRelv*sumcj defines species
  // not considered relevant.  This number is defined with reference to 1.0.
  constexpr Real c_fracIrrelv = 1.E-6;
  // The maximum amount a relevant species that is already negative is allowed
  // to change in the negative direction.  This number is defined with
  // reference to 1.0.
  // If this is set greater than CRDparam::g_ARKNonlinearTol, it may stall
  // convergence via displacement (since both are relative to Jrho)
  const Real c_minNegDelta = CRDparam::g_ARKNonlinearTol;
  CH_assert(a_box.numPts() == 1);
  const IntVect iv = a_box.smallEnd();
  int numIters = 0;

  const int numComps = a_solution.nComp();
  const int numSpecies = CRDparam::g_numSpecies;
  const int consSpecStart =
    CRDparam::g_CRDPhysics->speciesConsInterval().begin();
  const int rhoIdx = CRDparam::g_CRDPhysics->densityIndex();
  // const int Tidx = CRDparam::g_CRDPhysics->energyFluxIndex(); // Unused
  FABMatrix jacobianMatrix(a_box, (numSpecies)*(numSpecies));
  FABVector rhsVector(a_box, numComps);
  FABVector imRHS(a_box,numSpecies); // im = Only implicit terms
  FABVector imDeltaVector(a_box, numSpecies);
  FABVector deltaVector(a_box, numComps);
  FABVector imResidual(a_box, numSpecies);
  FABVector residual(a_box, numComps);
  StepSizeData stepSizeData;
  const Real JrhoFracIrrelv = a_solution[MD_IV(iv, rhoIdx)]*c_minNegDelta;

  // Vector of mass fractions
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);

  m_rhsFunc->operator()(a_solution, rhsVector, a_box);
  while (true)
    {
      ++numIters;
      // Find current solution delta
      m_jacFunc->operator()(a_solution, jacobianMatrix, a_box);

      for(int sp = 0; sp != numSpecies; ++sp)
        {
          int consSp = consSpecStart + sp;
          imRHS(a_box.smallEnd(), sp) = rhsVector(a_box.smallEnd(), consSp);
        }
      m_linearSolver->solve(jacobianMatrix, imDeltaVector, imResidual, imRHS,
                            a_box);

      // Update step size data
      stepSizeData.m_previousNorm = stepSizeData.m_currentNorm;
      Real stepSize = m_stepSizeComputer->computeStepSize(stepSizeData);
      // Expecting step size = 1.0 for following code

//--Step-size reduction based on reaction rates

      /* The general idea here is that the step size should be reduced to avoid
         negative species.  The challenge is that species already near zero
         could drastically limit the step, even if the changes are a result of
         errors.  So we limit the effect to "relevant" species which sum to a
         mass fraction > c_fracIrrelv (more or less, see code below).
         "Irrelevant" species are not allowed to affect the step size.
       */

      for (int sp = 0; sp != numSpecies; ++sp)
        {
          // Still working with JU here and throughout
          spec[sp] = a_solution[MD_IV(iv, consSpecStart + sp)];
        }
      // This only performs a sort as we are worried about FP error in sums
      // across a range of scales.  May not be necessary.
      static_cast<const ThermPhysics *const>(CRDparam::g_CRDPhysics)
        ->normalizePrimSpecies(
          CRDPhysics::NormalizeTypeNone,
          false,
          spec.data(),
          1,
          sortcnj.data());
      Real sumJrhocj = 0.;  // Note: sum |J rho c_j| != J rho unless all c_j > 0
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = sortcnj[j];
          sumJrhocj += std::abs(spec[sp]);
        }
      // Now go in reverse order, large to small.  Subtract from sumJrhocj and
      // when the sum falls below the threshold, the remaining species are
      // "irrelevant"
      const Real threshold = c_fracIrrelv*sumJrhocj;
      for (int j = numSpecies; j--;)
        {
          const int sp = sortcnj[j];
          const Real JcnNew = spec[sp] + imDeltaVector[MD_IV(iv, sp)];
          if (sumJrhocj > threshold)  // Relevant
            {
              sumJrhocj -= std::abs(spec[sp]);
              if (JcnNew < 0. && spec[sp] >= 0.)  // Went from + to -
                {
                  // U^(i)_{k+1} = U^(i)_k + stepSize*imDeltaVector  is bad
                  //           0 = U^(i)_k + stepSize*imDeltaVector  so use this
                  stepSize =
                    std::min(stepSize,
                             -spec[sp]/imDeltaVector[MD_IV(iv, sp)]);
                }
              else if (spec[sp] < 0. && imDeltaVector[MD_IV(iv, sp)] < 0.)
                // This particular species is considered significant.  But it
                // starts as negative and is becoming even more negative.
                // Note:  All species started out as >= 0.  To get here, it had
                // to start as zero or approach zero even with the step size
                // limiting given above.  In other words it had to first
                // approach a value near zero considered irrelevant.
                {
                  if (imDeltaVector[MD_IV(iv, sp)] < -JrhoFracIrrelv)
                    {
                      // This is not a fix to step size, rather a hard fix to
                      // the update itself.  We allow any species to change by
                      // an amount considered irrelevant for any single step.
                      imDeltaVector[MD_IV(iv, sp)] = -JrhoFracIrrelv;
                    }
                }
            }
        }
      CH_assert(stepSize > 0.);

//--Step-size reduction based on iteration count

      // Usually by 10 iterations, something is wrong.  Most often, we see the
      // solution oscillating between two exact states.  Halving the step size
      // breaks this cycle.  This is followed by further reductions.

      //**Note: If it is not converging even with 0.1, the system is probably
      //**not well condititioned.  Consider preconditioning of some form.
      if (numIters > 30)
        {
          stepSize = std::min(stepSize, 0.1);
        }
      else if (numIters > 20)
        {
          stepSize = std::min(stepSize, 0.25);
        }
      else if (numIters > 10)
        {
          stepSize = std::min(stepSize, 0.5);
        }

//--Actual update

      a_solution.plus(imDeltaVector, stepSize, 0, consSpecStart, numSpecies);
      // Any negative species should be either a result of floating-point error
      // or for any species below the threshold used for adjustments
      // for (int c = consSpecStart, c_end = consSpecStart + numSpecies;
      //      c != c_end; ++c)
      //   {
      //     a_solution[MD_IV(iv, c)] = std::max(0., a_solution[MD_IV(iv, c)]);
      //   }

      // Adjust for the convergence check
      imDeltaVector *= stepSize;

      // Test halting criteria
      ConvergenceData convergenceData;
      convergenceData.m_refScale = a_solution(iv, rhoIdx);
      convergenceData.m_displacement = &imDeltaVector;
      m_rhsFunc->operator()(a_solution, rhsVector, a_box);
      // Alias to only access species
      FArrayBox rhsVectorSpec(Interval(consSpecStart,
                                       consSpecStart + numSpecies - 1),
                              rhsVector);
      convergenceData.m_residual = &rhsVectorSpec;

      bool converged = m_convergenceChecker->converged(convergenceData);
      if (converged)
        break;
      if (numIters >= m_maxIters)
        {
          CRD::msg << "Reactions: NonLinearSolver max iterations reached: "
                   << m_maxIters << " on level " << CRDparam::g_level
                   << " in cell " << iv << CRD::warn;
          Box cellBox(iv, iv);
          // const Interval dispIval = imDeltaVector.interval(); // Unused
          // const Interval resIval = rhsVector.interval(); // Unused
          Real displacementInfNorm = imDeltaVector.norm(cellBox,
                                                        0,
                                                        0,
                                                        numSpecies);
          Real residualInfNorm = rhsVector.norm(cellBox,
                                                0,
                                                consSpecStart,
                                                numSpecies);
          // The way this is structured, there is no access to the tolerance
          // through the class data
          CRD::msg << "  residual    : " << residualInfNorm << " >= "
                   << a_solution(iv, rhoIdx)*CRDparam::g_ARKNonlinearTol
                   << CRD::warn;
          CRD::msg << "  displacement: " << displacementInfNorm << " >= "
                   << a_solution(iv, rhoIdx)*CRDparam::g_ARKNonlinearTol
                   << CRD::warn;
          break;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Solves the linear Ax=b problem that arises from the linearization of the nonlinear problem
/** \param[in]  a_matrix
 *                      The A matrix in the Ax=b problem
 *  \param[out] a_result
 *                      Solution to the linear problem, x in Ax=b
 *  \param[out] a_residual
 *                      Residual after solution, that is A*x
 *  \param[in]  a_rhs   RHS of the linear problem, the b vector in Ax=b
 *  \param[in]  a_box   Box over which the nonlinear problem should be solved
 *//*-----------------------------------------------------------------*/
void
LinearSolverLapack::solve(const FArrayBox& a_matrix,
                          FArrayBox&       a_result,
                          FArrayBox&       a_residual,
                          const FArrayBox& a_rhs,
                          const Box&       a_box) const
{
  CH_TIME("LinearSolverLapack::solve");
    const int numComponents = a_result.nComp();
    CHMatrix CHJac(numComponents,numComponents);
    CHVector CHRHS(numComponents);
    CHVector CHResidual(numComponents);
    CHArray<int,1> ipiv;

    CH_assert(a_box.smallEnd() == a_box.bigEnd());
    IntVect cell = a_box.smallEnd();

    FabToMatrix(CHJac, a_matrix, cell, numComponents, numComponents);
    FabToVector(CHRHS, a_rhs,  cell,  numComponents);

    // Solve Ax=b
    // Modifies CHRHS into Delta U
    CHgesv(CHJac, CHRHS, ipiv);

    // Compute the residual
    CHgemv(CHJac, CHRHS, CHResidual);

    VectorToFab(a_result, CHRHS, cell, numComponents);
    VectorToFab(a_residual, CHResidual, cell, numComponents);
}

/*--------------------------------------------------------------------*/
//  Converts an FArrayBox into a CHMatrix
/** \param[out] a_matrix
 *                      The output matrix of type CHMatrix
 *  \param[in] a_matrixFAB
 *                      FArrayBox to be converted into a CHMatrix
 *  \param[in]  a_cell  Which cell of the FArrayBox to convert into a matrix
 *  \param[in]  a_numRows
 *                      The number of rows in the matrix
 *  \param[in]  a_numCols
 *                      The number of columns in the marix
 *//*-----------------------------------------------------------------*/
void
FabToMatrix(CHMatrix&        a_matrix,
            const FArrayBox& a_matrixFAB,
            const IntVect&   a_cell,
            int              a_numRows,
            int              a_numCols)
{
  for(int row = 0; row != a_numRows; ++row)
    {
      for(int col = 0; col != a_numCols; ++col)
        {
          int linIndx = row*a_numCols+col;
          a_matrix(row,col) = a_matrixFAB(a_cell, linIndx);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Converts an CHMatrix into a FArrayBox
/** \param[out] a_matrixFAB
 *                      The output FArrayBox
 *  \param[in] a_matrix
 *                      CHMatrix to be converted into a FArrayBox
 *  \param[in]  a_cell  Which cell of the FArrayBox to store the matrix values
 *  \param[in]  a_numRows
 *                      The number of rows in the matrix
 *  \param[in]  a_numCols
 *                      The number of columns in the marix
 *//*-----------------------------------------------------------------*/
void
MatrixToFab(FArrayBox&      a_matrixFAB,
            const CHMatrix& a_matrix,
            const IntVect&  a_cell,
            int             a_numRows,
            int             a_numCols)
{
  for(int row = 0; row != a_numRows; ++row)
    {
      for(int col = 0; col != a_numCols; ++col)
        {
          int linIndx = row*a_numCols+col;
          a_matrixFAB(a_cell, linIndx) = a_matrix(row,col);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Converts an FArrayBox into a CHVector
/** \param[out] a_vector
 *                      The output CHVector
 *  \param[in] a_vectorFAB
 *                      The FArrayBox to be converted into a CHVector
 *  \param[in]  a_cell  Which cell of the FArrayBox to convert to a vector
 *  \param[in]  a_numElements
 *                      The number of elements in the vector
 *//*-----------------------------------------------------------------*/
void
FabToVector(CHVector&        a_vector,
            const FArrayBox& a_vectorFAB,
            const IntVect&   a_cell,
            int              a_numElements)
{
  for(int row = 0; row != a_numElements; ++row)
    {
      a_vector(row) = a_vectorFAB(a_cell, row);
    }
}

/*--------------------------------------------------------------------*/
//  Converts a CHVector into an FArrayBox
/** \param[out] a_vector
 *                      The output FArrayBox
 *  \param[in] a_vectorFAB
 *                      The CHVector to be converted into a FArrayBox
 *  \param[in]  a_cell  Which cell of the FArrayBox store the CHVector into
 *  \param[in]  a_numElements
 *                      The number of elements in the vector
 *//*-----------------------------------------------------------------*/
void
VectorToFab(FArrayBox&      a_vectorFAB,
            const CHVector& a_vector,
            const IntVect&  a_cell,
            int             a_numElements)
{
  for(int row = 0; row != a_numElements; ++row)
    {
      a_vectorFAB(a_cell, row) = a_vector(row);
    }
}

