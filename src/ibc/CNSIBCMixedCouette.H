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
 * \file CNSIBCMixedCouette.H
 *
 * \brief Initial and boundary conditions for transient Couette flow with mixing
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCMIXEDCOUETTE_H_
#define _CNSIBCMIXEDCOUETTE_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for transient Couette flow
///  with mixing
/**
 *//*+*************************************************************************/

class CNSIBCMixedCouette : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  /*  Note that Re*nu/L is the velocity (V) according to the Reynolds number.
   *  For SpaceDim > 2, the wall moves in all tangential directions, set by
   *  Vi = sqrt(V^2/(SpaceDim-1))
   */
  CNSIBCMixedCouette(const RealVect a_wallVelocity = RealVect::Zero,
                         const int  a_wallNormalDir = 1);

  /// Destructor
  virtual ~CNSIBCMixedCouette();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCMixedCouette(const CNSIBCMixedCouette&);

//--Assignment not permitted

  CNSIBCMixedCouette &operator=(const CNSIBCMixedCouette&);


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

  RealVect m_wallVelocity;            ///< Velocity of the wall
  int m_wallNormalDir;                ///< Direction normal to the moving wall
  int m_movWall;                      ///< 0 - low wall moves, 1 - high wall
  std::vector<Real> m_initMassFraction;
                                      ///< Initial species concentration
  std::string m_outputFileName;       ///< Name of file with exact solution
  Real m_initVel;                     ///< For setting an initial velocity
  bool m_sonicTunnel;                 ///< True if sonic tunnel case,
                                      ///< False if mixed couette flow
};

#endif /* ! defined _CNSIBCMIXEDCOUETTE_H_ */
