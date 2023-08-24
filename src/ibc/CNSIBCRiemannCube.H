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
 * \file CNSIBCRiemannCube.H
 *
 * \brief Initial and boundary conditions for general 3D Riemann problems
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCRIEMANNCUBE_H_
#define _CNSIBCRIEMANNCUBE_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions for general 3D Riemann problems
/**
 *//*+*************************************************************************/

class CNSIBCRiemannCube : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCRiemannCube();

  /// Destructor
  virtual ~CNSIBCRiemannCube();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCRiemannCube(const CNSIBCRiemannCube&);

//--Assignment not permitted

  CNSIBCRiemannCube &operator=(const CNSIBCRiemannCube&);


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

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;
  
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

  std::vector<Real> m_density;        ///< Initial density
  std::vector<Real> m_temperature;    ///< Initial temperature
  std::vector<Real> m_pressure;       ///< Initial pressure
  std::vector<Real> m_velocity;       ///< Initial velocity of flow field
  int               m_case;           ///< Case number for 2D Riemann problems
};

#endif /* ! defined _CNSIBCRIEMANNCUBE_H_ */
