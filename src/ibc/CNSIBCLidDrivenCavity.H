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
 * \file CNSIBCLidDrivenCavity.H
 *
 * \brief Initial and boundary conditions for species mixing test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCLIDDRIVENCAVITY_H_
#define _CNSIBCLIDDRIVENCAVITY_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for species mixing test case.
/**
 *          Moving wall
 *          --------->
 *    _________________________
 *   |            !            |
 *   |            !            |
 *   |            !            |
 *   |            !            |
 *   |  Left      !  Right     |
 *   |  side      !  side      |
 *   |  species   !  species   |
 *   |            !            |
 *   |            !            |
 *   |            !            |
 *   |            !            |
 *   +-------------------------+
 *
 *//*+*************************************************************************/

class CNSIBCLidDrivenCavity : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCLidDrivenCavity();

  /// Destructor
  virtual ~CNSIBCLidDrivenCavity();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCLidDrivenCavity(const CNSIBCLidDrivenCavity&);

//--Assignment not permitted

  CNSIBCLidDrivenCavity &operator=(const CNSIBCLidDrivenCavity&);


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

protected:
  
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
    const BCInfo&        a_domT) const;

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

  int m_wallDir;                      ///< Direction of wall
  int m_wallNormDir;                  ///< Direction normal to the moving wall
  std::vector<Real> m_leftMassFraction;
                                      ///< Species concentration on left side
  std::vector<Real> m_rightMassFraction;
                                      ///< Species concentration for right side
  Real m_wallVel;                     ///< Wall inlet velocity
  int m_wallOrder;                    ///< Order of extrapolation to use at wall
                                      ///< Default is 1
};

#endif /* ! defined _CNSIBCLIDDRIVENCAVITY_H_ */
