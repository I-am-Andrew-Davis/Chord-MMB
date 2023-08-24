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
 * \file CNSIBCSpecShockBox.H
 *
 * \brief Initial and boundary conditions for multispecies shock box problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSPECSHOCKBOX_H_
#define _CNSIBCSPECSHOCKBOX_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  BC and initial conditions for the shock box problem with species transport
/**
 *
 *//*+*************************************************************************/

class CNSIBCSpecShockBox : public CNSIBCCombustionReference
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCSpecShockBox();

  /// Destructor
  virtual ~CNSIBCSpecShockBox();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCSpecShockBox(const CNSIBCSpecShockBox&);

//--Assignment not permitted

  CNSIBCSpecShockBox &operator=(const CNSIBCSpecShockBox&);


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

  Real m_rhoL;                        ///< Initial density, lower
  Real m_presL;                       ///< Initial pressure, lower
  Real m_tL;                          ///< Iniital temperature, lower
  Real m_rhoU;                        ///< Initial density, upper
  Real m_presU;                       ///< Initial pressure, upper
  Real m_tU;                          ///< Initial temperature, upper
  std::vector<Real> m_lowerMassFraction;
                                      ///< Mass fractions for lower square
  std::vector<Real> m_upperMassFraction;
                                      ///< Mass fractions for rest of domain
};

#endif /* ! defined _CNSIBCSPECSHOCKBOX_H_ */
