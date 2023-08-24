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
 * \file CNSIBCSpecShock.H
 *
 * \brief Initial and boundary conditions for thermally perfect, multispecies
 *        Sod's shock tube problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSPECSHOCK_H_
#define _CNSIBCSPECSHOCK_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for Sod's shock tube problem
/**
 *                                 Diaphragm
 *   +---------------------------------+---------------------------------+
 *   | Left: rho      = 3.6            H Right: rho      = 1.2           |
 *   |       pressure = 3.6E5          H        pressure = 1.2E5         |
 *   +---------------------------------+---------------------------------|
 *                                Velocity = 0
 *
 *   The flow can be either viscous or inviscid.  Time is 0.2
 *
 *//*+*************************************************************************/

class CNSIBCSpecShock : public CNSIBCCombustionReference
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCSpecShock();

  /// Destructor
  virtual ~CNSIBCSpecShock();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCSpecShock(const CNSIBCSpecShock&);

//--Assignment not permitted

  CNSIBCSpecShock &operator=(const CNSIBCSpecShock&);


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

  Real m_rhoL;                        ///< Initial density on left
  Real m_presL;                       ///< Initial pressure on left
  Real m_tempL;                       ///< Initial temperature on left
  Real m_rhoR;                        ///< Initial density on right
  Real m_presR;                       ///< Initial pressure on right
  Real m_tempR;                       ///< Initial temperature on right
  Real m_velL;                        ///< Initial velocity on the left
  Real m_velR;                        ///< Initial velocity on the right
  Real m_Cr;                          ///< Smearing parameter
  std::vector<Real> m_leftMassFraction;
                                      ///< Left species concentrations
  std::vector<Real> m_rightMassFraction;
                                      ///< Right species concentrations
  Real m_wallVel;                     ///< Wall inlet velocity

  int m_tubeDir;                      ///< Direction of the tube itself
};

#endif /* ! defined _CNSIBCSSPECSHOCK_H_ */
