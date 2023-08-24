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
 * \file CNSIBCShear.H
 *
 * \brief Initial and boundary conditions for shear problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSHEAR_H_
#define _CNSIBCSHEAR_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for the shear problem
/**
 *
 *//*+*************************************************************************/

class CNSIBCShear : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCShear();

  /// Destructor
  virtual ~CNSIBCShear();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCShear(const CNSIBCShear&);

//--Assignment not permitted

  CNSIBCShear &operator=(const CNSIBCShear&);


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

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);
  
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
  RealVect m_velConst;                ///< Velocity multiplier
};

#endif /* ! defined _CNSIBCSHEARPROBLEM_H_ */
