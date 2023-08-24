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
 * \file CNSIBCTaylorGreen.H
 *
 * \brief Initial and boundary conditions for Taylor-Green test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTAYLORGREEN_H_
#define _CNSIBCTAYLORGREEN_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for Taylor-Green test case.
/**
 *//*+*************************************************************************/

class CNSIBCTaylorGreen : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCTaylorGreen();

  /// Destructor
  virtual ~CNSIBCTaylorGreen();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTaylorGreen(const CNSIBCTaylorGreen&);

//--Assignment not permitted

  CNSIBCTaylorGreen &operator=(const CNSIBCTaylorGreen&);


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

  ///  Output pnt-values and domain-sums specific to simulation
  virtual void inSituSumPntProcessing(
    const LevelData<FluxBox>&   a_faceAvgPlotData,
    const LevelData<FluxBox>&   a_WfaceAvgFxb,
    const LevelData<FArrayBox>& a_WcellAvgFab,
    const LevelData<FArrayBox>& a_WcellPntFab,
    const LayoutData<FluxBox>&  a_unitNormalsFxb,
    const LevelGridMetrics&     a_gridMetrics,
    const Real                  a_time,
    const int                   a_stage,
    const int                   a_level) const;

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

  int      m_turbState;               ///< Turbulent fluctuation state
  int      m_turbForcing;             ///< 1 = using forcing
  RealVect m_numVort;                 ///< Number of vortices in each direction
  bool     m_zVelIC;                  ///< Z velocity IC: 0 = sine, 1 = cosine
  bool     m_uniformRhoIC;            ///< Density IC: 0 = varying, 1 = constant
  Real     m_coeffG;                  ///< Proportional controller magnitude for
                                      ///< mean turbulent energy forcing. Higher
                                      ///< magnitude forces mean energy to more
                                      ///< closely follow target energy
};

#endif /* ! defined _CNSIBCTAYLORGREEN_H_ */
