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
 * \file CNSIBCCylinderFlow.H
 *
 * \brief Initial and boundary conditions for flow over a cylinder
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCCYLINDERFLOW_H_
#define _CNSIBCCYLINDERFLOW_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for cylinder flow.
/**
 *//*+*************************************************************************/

class CNSIBCCylinderFlow : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  /* 
   */
  CNSIBCCylinderFlow();

  /// Destructor
  virtual ~CNSIBCCylinderFlow();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCCylinderFlow(const CNSIBCCylinderFlow&);

//--Assignment not permitted

  CNSIBCCylinderFlow &operator=(const CNSIBCCylinderFlow&);


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

  std::vector<Real> m_threshold;      ///< Tag for vorticity per level    
  Real m_tagPerc;                     ///< Percent to tag in radial direction
  int m_tagDir;                       ///< The radial direction
  bool m_tagSide;                     ///< Hi or lo side for tagging
  
};

#endif /* ! defined _CNSIBCCYLINDERFLOW_H_ */
