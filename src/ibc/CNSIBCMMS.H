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
 * \file CNSIBCMMS.H
 *
 * \brief Initial and boundary conditions for MMS test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCMMS_H_
#define _CNSIBCMMS_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"

/*******************************************************************************
 */
///  BC and initial conditions for method-of-manufactured-solutions case
/**
 *//*+*************************************************************************/

class CNSIBCMMS : public CNSIBCGeneralized
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCMMS();

  /// Destructor
  virtual ~CNSIBCMMS();


/*==============================================================================
 * Deleted constructors and destructors
 *============================================================================*/

  /// Copy not permitted
  CNSIBCMMS(const CNSIBCMMS&) = delete;

  /// Assignment not permitted
  CNSIBCMMS &operator=(const CNSIBCMMS&) = delete;


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
 * Protected member functions
 *============================================================================*/

protected:
  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const override;
  
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

  Real     m_rA;                      ///< Density-wave amplitude
  RealVect m_uA;                      ///< Velocity-wave amplitude
  Real     m_pA;                      ///< Pressure-wave amplitude
  Real     m_TA;                      ///< Temperature-wave amplitude
  RealVect m_sOmegaRho;               ///< Density spatial-frequency
  Real     m_tOmegaRho;               ///< Density temporal-frequency
  RealVect m_sPhiRho;                 ///< Density spatial phase-shift
  Real     m_tPhiRho;                 ///< Density temporal phase-shift
  RealVect m_sOmegaU;                 ///< X-vel spatial-frequency
  Real     m_tOmegaU;                 ///< X-vel temporal-frequency
  RealVect m_sPhiU;                   ///< X-vel spatial phase-shift
  Real     m_tPhiU;                   ///< X-vel temporal phase-shift
  RealVect m_sOmegaV;                 ///< Y-vel spatial-frequency
  Real     m_tOmegaV;                 ///< Y-vel temporal-frequency
  RealVect m_sPhiV;                   ///< Y-vel spatial phase-shift
  Real     m_tPhiV;                   ///< Y-vel temporal phase-shift
  RealVect m_sOmegaW;                 ///< Z-vel spatial-frequency
  Real     m_tOmegaW;                 ///< Z-vel temporal-frequency
  RealVect m_sPhiW;                   ///< Z-vel spatial phase-shift
  Real     m_tPhiW;                   ///< Z-vel temporal phase-shift
  RealVect m_sOmegaP;                 ///< Pressure spatial-frequency
  Real     m_tOmegaP;                 ///< Pressure temporal-frequency
  RealVect m_sPhiP;                   ///< Pressure spatial phase-shift
  Real     m_tPhiP;                   ///< Pressure temporal phase-shift
  RealVect m_sOmegaT;                 ///< Temperature spatial-frequency
  Real     m_tOmegaT;                 ///< Temperature temporal-frequency
  RealVect m_sPhiT;                   ///< Temperature spatial phase-shift
  Real     m_tPhiT;                   ///< Temperature temporal phase-shift
  RealVect m_charLength;              ///< Characteristic length of frequency 1
  Real     m_charTime;                ///< Characteristic time of frequency 1
  int      m_bcTestCase;              ///< Boundary condition test case number
  bool     m_exactInit;               ///< True = exact cell-averaged IC
  bool     m_sixthOrderInit;          ///< True = sixth-order IC
};

#endif /* ! defined _CNSIBCSINGLEBLOCKMMS_H_ */
