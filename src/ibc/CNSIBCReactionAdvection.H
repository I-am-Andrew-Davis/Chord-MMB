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
 * \file CNSIBCReactionAdvection.H
 *
 * \brief Initial and boundary conditions for H2-O2 reaction front
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCREACTIONADVECTION_H_
#define _CNSIBCREACTIONADVECTION_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  Initialization and BC for H2-O2 reaction front
/**
 *//*+*************************************************************************/

class CNSIBCReactionAdvection : public CNSIBCCombustionReference
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCReactionAdvection();

  /// Destructor
  virtual ~CNSIBCReactionAdvection();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCReactionAdvection(const CNSIBCReactionAdvection&);

//--Assignment not permitted

  CNSIBCReactionAdvection &operator=(const CNSIBCReactionAdvection&);


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
  
  Real m_C;                           ///< Constant in initialization equation
  Real m_x0;                          ///< Location of the center
  Real m_L;                           ///< Width of hydrogen plateau
  Real m_o2T;                         ///< Initial temperature of O2
  Real m_h2T;                         ///< Initial temperature of H2
  Real m_pressure;                    ///< Initial pressure
  RealVect m_initVel;                 ///< Initial velocity
  int m_O2comp;                       ///< O2 comp number
  int m_H2comp;                       ///< H2 comp number
};

#endif /* ! defined _CNSIBCREACTIONADVECTION_H_ */
