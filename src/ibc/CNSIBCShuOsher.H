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
 * \file CNSIBCShuOsher.H
 *
 * \brief Initial and boundary conditions for the Shu-Osher problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSHUOSHER_H_
#define _CNSIBCSHUOSHER_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions for the Shu-Osher problem
/**
 *
 *//***************************************************************************/

class CNSIBCShuOsher : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCShuOsher();

  /// Destructor
  virtual ~CNSIBCShuOsher();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCShuOsher(const CNSIBCShuOsher&);

//--Assignment not permitted

  CNSIBCShuOsher &operator=(const CNSIBCShuOsher&);


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

  Real m_shockLoc;                    ///< Location of the shock
  Real m_rhoR;                        ///< Initial density on right
  Real m_presR;                       ///< Initial pressure on right
  Real m_tempR;                       ///< Initial temperature on right
  Real m_velR;                        ///< Initial velocity on the right
  Real m_rhoAmp;                      ///< Amplitude of density perturbations
  Real m_rhoFreq;                     ///< Frequency of density perturbations
};

#endif /* ! defined _CNSIBCSHUOSHER_H_ */
