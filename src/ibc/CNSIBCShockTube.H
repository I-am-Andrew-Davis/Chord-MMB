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
 * \file CNSIBCShockBox.H
 *
 * \brief Initial and boundary conditions for Sod's shock tube problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSHOCKTUBE_H_
#define _CNSIBCSHOCKTUBE_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for Sod's shock tube problem
/**
 *                                 Diaphragm
 *   +---------------------------------+---------------------------------+
 *   | Left: rho      = 1              H Right: rho      = 0.125         |
 *   |       pressure = 1              H        pressure = 0.1           |
 *   +---------------------------------+---------------------------------|
 *                                Velocity = 0
 *
 *   The flow can be either viscous or inviscid.  Time is 0.2
 *
 *//*+*************************************************************************/

class CNSIBCShockTube : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCShockTube();

  /// Destructor
  virtual ~CNSIBCShockTube();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCShockTube(const CNSIBCShockTube&);

//--Assignment not permitted

  CNSIBCShockTube &operator=(const CNSIBCShockTube&);


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
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

 
/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Real m_rhoL;                        ///< Initial density on left
  Real m_pressureL;                   ///< Initial pressure on left
  Real m_rhoR;                        ///< Initial density on right
  Real m_pressureR;                   ///< Initial pressure on right
  Real m_threshold;                   ///< Threshold of relative density
                                      ///< gradient for tagging cells

  int m_tubeDir;                      ///< Direction of the tube itself
};

#endif /* ! defined _CNSIBCShockTube_H_ */
