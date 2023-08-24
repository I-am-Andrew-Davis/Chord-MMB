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
 * \file CNSIBCSpecMachReflection.H
 *
 * \brief Initial and boundary conditions for shock reflecting off ramp
 *        with multispecies
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSPECMACHREFLECTION_H_
#define _CNSIBCSPECMACHREFLECTION_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


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

class CNSIBCSpecMachReflection : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCSpecMachReflection();

  /// Destructor
  virtual ~CNSIBCSpecMachReflection();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCSpecMachReflection(const CNSIBCSpecMachReflection&);

//--Assignment not permitted

  CNSIBCSpecMachReflection &operator=(const CNSIBCSpecMachReflection&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// Set the primitive state at wall BC
  virtual void setWallBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_boundarySlipVelocity,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

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
  int m_lowXBound;                    ///< Sets the lower x boundary
                                      ///< 0 - Dirichlet
                                      ///< 1 - Extrapolated
  int m_hiYBound;                     ///< Sets the upper y boundary
                                      ///< 0 - Dirichlet
                                      ///< 1 - Extrapolated
                                      ///< 2 - Farfield
  Real m_alpha;                       ///< Angle of the ramp (rad)
  Real m_xLead;                       ///< Lead length before ramp
  Real m_xRamp;                       ///< Length of ramp (projected to x-axis)
  Real m_xBase;                       ///< Starting x-location of shock wave
  Real m_xloclow;                     ///< Location of start of no slip wall
  Real m_p0;                          ///< Pressure in front of shock
  Real m_r0;                          ///< Density in front of shock
  Real m_c0;                          ///< Sound speed in front of shock
  Real m_Mach;                        ///< Mach number of the shock wave
  Real m_p1;                          ///< Pressure behind the shock
  Real m_r1;                          ///< Density behind the shock
  Real m_u1;                          ///< Velocity behind the shock
  Real m_threshold;                   ///< Threshold of relative density
                                      ///< gradient for tagging cells
  int m_wallOrder;                    ///< Order of extrapolation for the wall
                                      ///< Default is 4, can be changed to 1
  int m_tagComp;                      ///< Set the component to tag
  bool m_tagWall;                     ///< T - Tag to refine near wall
                                      ///< F - Only tag based on gradients
  std::vector<Real> m_upstreamMF;     ///< Mass fractions upstream of shock
  std::vector<Real> m_downstreamMF;   ///< Mass fractions downstream of shock
};

#endif /* ! defined _CNSIBCSPECMACHREFLECTION_H_ */
