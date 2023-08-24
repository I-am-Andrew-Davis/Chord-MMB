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
 * \file CNSIBCObliqueWave.H
 *
 * \brief Initial and boundary conditions for steady oblique wave
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCOBLIQUEWAVE_H_
#define _CNSIBCOBLIQUEWAVE_H_

//----- Internal -----//

#include "CNSIBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for Mach reflection problem
/**
 *   +--------------------------------------------+
 *   |                                            |
 *   |                                            |
 *   |                                            j
 *   |      Initialized to state "initial"     __/
 *   |                                      __/
 *   |                                   __/
 *   |                                __/
 *   |                             __/
 *   |                          __/
 *   |                       __/    alpha = 30.0
 *   +----------------------'
 *                          ^- x = 0
 *   <- lead length = 0.8 -><-- ramp length = 4 -->
 *
 *   The flow can be either viscous or inviscid.
 *
 *//*+*************************************************************************/

class CNSIBCObliqueWave : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCObliqueWave();

  /// Destructor
  virtual ~CNSIBCObliqueWave();

  // Copy not permitted
  CNSIBCObliqueWave(const CNSIBCObliqueWave&) = delete;

  // Assignment not permitted
  CNSIBCObliqueWave &operator=(const CNSIBCObliqueWave&) = delete;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const override;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const override;


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
  int m_wallOrder;                    ///< Order of extrapolation for the wall
                                      ///< Default is 4, can be changed to 1
};

#endif /* ! defined _CNSIBCOBLIQUEWAVE_H_ */
