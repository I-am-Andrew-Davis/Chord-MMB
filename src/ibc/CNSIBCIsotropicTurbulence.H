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
 * \file CNSIBCIsotropicTurbulence.H
 *
 * \brief Initial and boundary conditions for homogeneous isotropic turbulence
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCISOTROPICTURBULENCE_H_
#define _CNSIBCISOTROPICTURBULENCE_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for Taylor-Green test case.
/**
 *//*+*************************************************************************/

class CNSIBCIsotropicTurbulence : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCIsotropicTurbulence();

  /// Destructor
  virtual ~CNSIBCIsotropicTurbulence();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCIsotropicTurbulence(const CNSIBCIsotropicTurbulence&);

//--Assignment not permitted

  CNSIBCIsotropicTurbulence &operator=(const CNSIBCIsotropicTurbulence&);


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

  /// Add source term (turned on with "PhysicsSource")
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
                             const Real           a_GlobalKE,
                             const Real           a_globalHelicity) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

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
  Real     m_P0;                      ///< Mean thermodynamic pressure
  Real     m_rho0;                    ///< mean density
  Real     m_T0;
  Real     m_Einit;
};

#endif /* ! defined _CNSIBCISOTROPICTURBULENCE_H_ */
