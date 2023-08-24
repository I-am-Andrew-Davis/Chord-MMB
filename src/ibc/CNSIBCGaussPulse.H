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
 * \file CNSIBCGaussPulse.H
 *
 * \brief Initial and boundary conditions for a Gaussian pulse problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCGAUSSPULSE_H_
#define _CNSIBCGAUSSPULSE_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for Gaussian pulse problem
/**
 *
 *//*+*************************************************************************/

class CNSIBCGaussPulse : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCGaussPulse();

  /// Destructor
  virtual ~CNSIBCGaussPulse();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCGaussPulse(const CNSIBCGaussPulse&);

//--Assignment not permitted

  CNSIBCGaussPulse &operator=(const CNSIBCGaussPulse&);


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

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;
  
/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

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

  Real m_deltaRho;                    ///< Change in density for pulse
  RealVect m_center;                  ///< Center of pulse
};

#endif /* ! defined _CNSIBCGAUSSPULSE_H_ */
