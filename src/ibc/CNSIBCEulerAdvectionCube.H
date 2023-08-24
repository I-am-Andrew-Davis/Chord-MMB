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
 * \file CNSIBCEulerAdvectionCube.H
 *
 * \brief Initial and boundary conditions for Euler advection in a cube
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCEULERADVECTIONCUBE_H_
#define _CNSIBCEULERADVECTIONCUBE_H_

//----- Internal -----//

#include "CNSIBC.H"


/*******************************************************************************
 */
/// BC and initial conditions on the mesh for Euler advection in a cube
/** This class implements advection of a Gaussian profile in a periodic domain.
 *
 *//*+*************************************************************************/

class CNSIBCEulerAdvectionCube : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCEulerAdvectionCube();

  /// Destructor
  virtual ~CNSIBCEulerAdvectionCube();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCEulerAdvectionCube(const CNSIBCEulerAdvectionCube&);

//--Assignment not permitted

  CNSIBCEulerAdvectionCube &operator=(const CNSIBCEulerAdvectionCube&);


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

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells (defaults to a no-op)
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

  RealVect m_velocity;                ///< Velocity field
  RealVect m_center;                  ///< Center of the Gaussian at time t=0
                                      ///< in dimensions of the grid
  Real m_deltaRho;                    ///< Change in density from center of
                                      ///< Gaussian to ambient
  Real m_sizeGaussian;                ///< Size of the Gaussian given as one
                                      ///< standard deviation
  Real m_radMax;                      ///< Smoothing is applied outside this
                                      ///< radius
  bool m_useSharpProfiles;            ///< Use "ball and Jacks" instead of guassian
  Real m_threshold;                   ///< Threshold of relative density
                                      ///< gradient for tagging cells
};

#endif /* ! defined _CNSIBCEULERADVECTIONCUBE_H_ */
