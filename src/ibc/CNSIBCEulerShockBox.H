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
 * \file CNSIBCEulerShockBox.H
 *
 * \brief Initial and boundary conditions for Euler shock box problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCEULERSHOCKBOX_H_
#define _CNSIBCEULERSHOCKBOX_H_

//----- Internal -----//

#include "CNSIBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the lattice.
/**  This class serves as reference for setting initial and boundary conditions
 *   on the lattice.  It is expected that more complicated setups derive from
 *   this class.
 *
 *//*+*************************************************************************/

class CNSIBCEulerShockBox : public CNSIBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCEulerShockBox();

  /// Destructor
  virtual ~CNSIBCEulerShockBox();

  // Prohibit all others
  CNSIBCEulerShockBox(const CNSIBCEulerShockBox&) = delete;
  CNSIBCEulerShockBox(CNSIBCEulerShockBox&&) = delete;
  CNSIBCEulerShockBox& operator=(const CNSIBCEulerShockBox&) = delete;
  CNSIBCEulerShockBox& operator=(CNSIBCEulerShockBox&&) = delete;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Allocations of new physics states.  Customize for each derivative class.
  static std::vector<CRDPhysics*> allocatePhysics();

  /// Return a name describing the IBC
  virtual const char *const IBCName() const override;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const override;

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const override;


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

  RealVect m_bxLo;                    ///< Lower physical coordinate of high
                                      ///< pressure box
  RealVect m_bxHi;                    ///< Upper physical coordinate of high
                                      ///< pressure box
};

#endif /* ! defined _CNSIBCEulerShockBox_H_ */
