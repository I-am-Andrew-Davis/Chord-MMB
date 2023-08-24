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
 * \file CNSIBCWedgeDetonation.H
 *
 * \brief Initial and boundary conditions for shock reflecting off ramp
 *        with multispecies
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCWEDGEDETONATION_H_
#define _CNSIBCWEDGEDETONATION_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for wedge detonation case
/**
 *
 *//*+*************************************************************************/

class CNSIBCWedgeDetonation : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCWedgeDetonation();

  /// Destructor
  virtual ~CNSIBCWedgeDetonation();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCWedgeDetonation(const CNSIBCWedgeDetonation&);

//--Assignment not permitted

  CNSIBCWedgeDetonation &operator=(const CNSIBCWedgeDetonation&);


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
  Real m_alpha;                       ///< Angle of the ramp (rad)
  Real m_xLead;                       ///< Lead length before ramp
  Real m_xRamp;                       ///< Length of ramp (projected to x-axis)
  int m_wallTagCells;                 ///< Number of wall normal cells to tag
  bool m_refAlongWall;                ///< T - Refine along entire wall
  BCInfo m_rampBC;                    ///< BC for ramp
};

#endif /* ! defined _CNSIBCWEDGEDETONATION_H_ */
