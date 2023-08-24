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
 * \file CNSIBCMachReflection.H
 *
 * \brief Initial and boundary conditions for shock reflecting off ramp
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCMACHREFLECTION_H_
#define _CNSIBCMACHREFLECTION_H_

//----- Internal -----//

#include "CNSIBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for Mach reflection problem
/**
 *   +--------------------------------------------+
 *   |        H                                   |
 *   |        H                                   |
 *   |        H                                   j
 *   |        H    c0                          __/
 *   |  p1    H    p0 = 1                   __/
 *   |  r1    H    r0 = gamma            __/
 *   |  u1    H    u0 = 0             __/
 *   |        H                    __/
 *   |        H--> shock Mach # __/
 *   |        H              __/    alpha = 30.0
 *   +----------------------'
 *            ^             ^- x = 0
 *            +- xBase = -0.5
 *   <- lead length = 0.8 -><-- ramp length = 4 -->
 *
 *   The flow can be either viscous or inviscid.
 *
 *//*+*************************************************************************/

class CNSIBCMachReflection : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCMachReflection();

  /// Destructor
  virtual ~CNSIBCMachReflection();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCMachReflection(const CNSIBCMachReflection&);

//--Assignment not permitted

  CNSIBCMachReflection &operator=(const CNSIBCMachReflection&);


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

  Real m_alpha;                       ///< Angle of the ramp (rad)
  Real m_xBase;                       ///< Starting x-location of shock wave
  Real m_c0;                          ///< Sound speed in front of shock
  Real m_Mach;                        ///< Mach number of the shock wave
  int m_wallOrder;                    ///< Order of extrapolation for the wall
                                      ///< Default is 4, can be changed to 1
  int m_idxStatePreShock;             ///< Index for state in front of shock
  int m_idxStatePostShock;            ///< Index for state behind the shock
};

#endif /* ! defined _CNSIBCMACHREFLECTION_H_ */
