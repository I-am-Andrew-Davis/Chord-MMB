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
 * \file DCFlattening.H
 *
 * \brief Flattens cell and face deconvolution and convolution
 *
 *//*+*************************************************************************/

#ifndef _DCFLATTENING_H_
#define _DCFLATTENING_H_

//----- Chombo Library -----//

#include "LGintegrator.H"

/*******************************************************************************
 */
///  A simple flattener for deconvolution and convolution operations
/**
 *   The parameters for this flattener are set in the cpp file
 *//*+*************************************************************************/

class DCFlattening
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  DCFlattening(const Real& a_nlvTol,
               const Real& a_fcorTol);

  /// Destructor
  ~DCFlattening();

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  DCFlattening(const DCFlattening&);

//--Assignment not permitted

  DCFlattening &operator=(const DCFlattening&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Directly flatten cell and face values
  void directDCFlatten(const Box&       a_box,
                       FArrayBox&       a_HOFab,
                       const FArrayBox& a_LOFab,
                       const int        a_numComp,
                       const bool       a_useExtraChecks = true) const;

  // Apply FCOR based on differences between HO and LO left and right states
  void ppmFCOR(FArrayBox&           a_Wminus,
               FArrayBox&           a_Wplus,
               const FArrayBox&     a_W,
               const FArrayBox&     a_WLO,
               const FArrayBox&     a_WminusLO,
               const FArrayBox&     a_WplusLO,
               const int&           a_numSlopes,
               const int&           a_dir,
               const Box&           a_box,
               const ProblemDomain& a_domain) const;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Real m_nlvTol;                      ///< Tolerance for DC flatten
  Real m_fcorTol;                     ///< Tolerance for mass fractions
};

#endif  /* ! defined _DCFLATTENING_H_ */
