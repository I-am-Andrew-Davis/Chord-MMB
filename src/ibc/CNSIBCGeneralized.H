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
 * \file CNSIBCGeneralized.H
 *
 * \brief Initial and boundary conditions for a generalized case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCGENERALIZED_H_
#define _CNSIBCGENERALIZED_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  BC and initial conditions for a generalized problem
/**
 *//*+*************************************************************************/

class CNSIBCGeneralized : public CNSIBCCombustionReference
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCGeneralized(const bool a_readInputs = true);

  /// Destructor
  virtual ~CNSIBCGeneralized();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCGeneralized(const CNSIBCGeneralized&) = delete;

//--Assignment not permitted

  CNSIBCGeneralized &operator=(const CNSIBCGeneralized&) = delete;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const override;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const override;

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const override;

//--Patch operators

  /// Add a source term to the flow
  virtual void addSourceTerm(
    FArrayBox&           a_sourceFab,
    FArrayBox&           a_invDtFab,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_UcellAvg,
    const FluxBox&       a_WfaceAvgFxb,
    const BlockDomain&   a_domain,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const Real           a_stageWeight,
    const int            a_level,
    const Box&           a_disjointBox,
    const Box&           a_solveBox,
    const DataIndex&     a_dataIndx,
    const Real           a_globalKE,
    const Real           a_globalHelicity) const override;
  
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

  /// Set the primitive state at wall BC
  virtual void setWallBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_boundarySlipVelocity,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const override;

  /// Set characteristic inlet profiles for velocity, temperature, species
  virtual void setCNSCBCProfiles(
    FArrayBox&           a_BCProfile,
    const Box&           a_boundaryBox,
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

  int m_idxStateInit;                 ///< State for initialization
  RealVect m_bodyForce;               ///< Body force in source term


  // FIXME These are obsolete, stop using them!
  std::vector<Real> m_initMassFraction;
                                      ///< Initial species concentration
  Real m_initRho = -1;                     ///< Initial density
  Real m_initT = -1;                       ///< Initial temperature
  Real m_initP = -1;                       ///< Initial pressure
  RealVect m_initVel;                 ///< Initial velocity of flow field
};

#endif /* ! defined _CNSIBCGENERALIZED_H_ */
