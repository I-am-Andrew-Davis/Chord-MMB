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
 * \file CNSIBCSphere.H
 *
 * \brief Initial and boundary conditions for a sphere case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSPHERE_H_
#define _CNSIBCSPHERE_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions for a sphere problem
/**
 *//*+*************************************************************************/

class CNSIBCSphere : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCSphere(const bool a_readInputs = true);

  /// Destructor
  virtual ~CNSIBCSphere();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCSphere(const CNSIBCSphere&) = delete;

//--Assignment not permitted

  CNSIBCSphere &operator=(const CNSIBCSphere&) = delete;


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

  /// Set wall velocity on specific boundary -- user can specialize this
  virtual void setBndryWallVelocity(
    FArrayBox&           a_wallVelocity,
    const FArrayBox&     a_unitNormalBasisFab,
    const Box&           a_boundaryFaceBox,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const ProblemDomain& a_domain,
    const Real           a_time,
    const int            a_level,
    const int            a_dir,
    const Side::LoHiSide a_side,
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

  Real     m_radiusInner;             ///< Inner radius of sphere (surface)
  Real     m_radiusOuter;             ///< Outer radius (edge of domain)
  Real     m_rotationRateRad;         ///< Rotation rate in radians per second
  RealVect m_omega;                   ///< Vector of rotational axis
  bool     m_specialInit;             ///< True: debugging initialization
};

#endif /* ! defined _CNSIBCSPHERE_H_ */
