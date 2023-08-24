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
 * \file CNSIBCShockShock.H
 *
 * \brief Initial and boundary conditions for shock-shock intersection cases
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSHOCKSHOCK_H_
#define _CNSIBCSHOCKSHOCK_H_

//----- Internal -----//

#include "CNSIBC.H"


/*******************************************************************************
 */
///  BC and initial conditions for a sphere problem
/**
 *//*+*************************************************************************/

class CNSIBCShockShock : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCShockShock();

  /// Destructor
  virtual ~CNSIBCShockShock();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCShockShock(const CNSIBCShockShock&) = delete;

//--Assignment not permitted

  CNSIBCShockShock &operator=(const CNSIBCShockShock&) = delete;


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


/*==============================================================================
 * Protected member functions
 *============================================================================*/

  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const BoundaryIndex&          a_bcIdx,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const BCInfo&                 a_domT) const override;


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

  Real m_beta;                        ///< Oblique shock angle
  Real m_m;                           ///< Oblique shock slope
  Real m_b;                           ///< Oblique shock intercept
  int m_idxStatePreShock;             ///< Index for state in front of shock
                                      ///< (aka "initial" state)
  int m_idxStatePostShock;            ///< Index for state behind the shock
  

};

#endif /* ! defined _CNSIBCSPHERE_H_ */
