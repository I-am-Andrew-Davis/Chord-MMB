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
 * \file CNSIBCTemperatureDiffusion.H
 *
 * \brief Initial and boundary conditions for temperature diffusion test
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTEMPERATUREDIFFUSION_H_
#define _CNSIBCTEMPERATUREDIFFUSION_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions for a problem with a heated wall
/**
 *//*+*************************************************************************/

class CNSIBCTemperatureDiffusion : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCTemperatureDiffusion();

  /// Destructor
  virtual ~CNSIBCTemperatureDiffusion();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTemperatureDiffusion(const CNSIBCTemperatureDiffusion&);

//--Assignment not permitted

  CNSIBCTemperatureDiffusion &operator=(const CNSIBCTemperatureDiffusion&);


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

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells
  virtual int exactSol(FArrayBox&              a_Ux,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       const LevelGridMetrics& a_gridMetrics,
                       const FluxBox&          a_unitNormals,
                       const DataIndex&        a_didx,
                       const Real              a_time,
                       const int               a_level) const;

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
  RealVect m_trueSourceLoc;                    ///< Location of true source term
  RealVect m_sourceLoc;                        ///< Location of source term
  Real m_sourceSigma;                          ///< Source "spread"
  Real m_sourceAlpha;                          ///< Source amplitude
  Real m_timeSourceEnd;                        ///< Time to end source
  Real m_tHi;                                  ///< Initial hi temp zone temperature
  Real m_tLo;                                  ///< Initial lo temp zone temperature
  Real m_initRadius;                           ///< Initial hi temp region radius size
};

#endif /* ! defined _CNSIBCTEMPERATUREDIFFUSION_H_ */
