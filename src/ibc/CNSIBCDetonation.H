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
 * \file CNSIBCDetonation.H
 *
 * \brief Initial and boundary conditions for detonation problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCDETONATION_H_
#define _CNSIBCDETONATION_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for species mixing test case.
/**
 *                             Periodic
 *         +===============================================
 *         |Hi T__                                        )
 *         |       \                                      )
 *         |         \                                    )
 * Wall    |           \                                  )  Open end
 *         |             \                                )
 *         |               \ __ Low T                     )
 *         +===============================================
 *         |---grad length---|
 *
 *//*+*************************************************************************/

class CNSIBCDetonation : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCDetonation();

  /// Destructor
  virtual ~CNSIBCDetonation();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCDetonation(const CNSIBCDetonation&);

//--Assignment not permitted

  CNSIBCDetonation &operator=(const CNSIBCDetonation&);


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

  Real m_highT;                       ///< High temperature along gradient
  Real m_lowT;                        ///< Low temperature along gradient
  Real m_gradLength;                  ///< Length of the gradient of temperature
};

#endif /* ! defined _CNSIBCDETONATION_H_ */
