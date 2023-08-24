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
 * \file CNSIBCCombustionTest.H
 *
 * \brief Initial and boundary conditions for the Rayleigh-Taylor instability
 *        problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCCOMBUSTIONTEST_H_
#define _CNSIBCCOMBUSTIONTEST_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  
/**
 *//*+*************************************************************************/

class CNSIBCCombustionTest : public CNSIBCCombustionReference
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCCombustionTest();

  /// Destructor
  virtual ~CNSIBCCombustionTest();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCCombustionTest(const CNSIBCCombustionTest&);

//--Assignment not permitted

  CNSIBCCombustionTest &operator=(const CNSIBCCombustionTest&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;

//--Patch operators

  /// Add source term
  virtual void addSourceTerm(FArrayBox&           a_sourceFab,
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
                             const Real           a_globalHelicity) const;

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;
  
/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<Real> m_upMassFraction; ///< Species concentration for upper part
  std::vector<Real> m_loMassFraction; ///< Species concentration for lower part
  Real m_upRho;                       ///< Upper density
  Real m_loRho;                       ///< Lower density
  Real m_grav;                        ///< Gravity force
  Real m_pressure;                    ///< Initial pressure throughout domain
};

#endif /* ! defined _CNSIBCCOMBUSTIONTEST_H_ */
