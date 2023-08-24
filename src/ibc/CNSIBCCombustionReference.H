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
 * \file CNSIBCReferenceCubeBC.H
 *
 * \brief IBC implementing BC methods for ThermPhysics cases
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCCOMBUSTIONREFERENCE_H_
#define _CNSIBCCOMBUSTIONREFERENCE_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "CHArray.H"

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"
#include "TagLevel.H"

//---- Debug ----//
#include "CRDhardconfig.H"

/*******************************************************************************
 */
///  BC and initial conditions on the lattice.
/**  This class serves as reference for setting initial and boundary conditions
 *   for cases involving the ThermPhysics class. Note that pure function 
 *   'initialize' is not implemented so this class is not sufficient as an IBC.
 *
 *//*+*************************************************************************/

class CNSIBCCombustionReference : public CNSIBC
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor sets BC at domain extents
  CNSIBCCombustionReference();

  /// Destructor
  virtual ~CNSIBCCombustionReference();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCCombustionReference(const CNSIBCCombustionReference&);

//--Assignment not permitted

  CNSIBCCombustionReference &operator=(const CNSIBCCombustionReference&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Allocations of new physics states.  Customize for each derivative class.
  static std::vector<CRDPhysics*> allocatePhysics();

  /// Return a name describing the IBC
  virtual const char *const IBCName() const  override;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize) override;

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  void setTagMethodLevel(
    const int        a_tagBufferSize,
    TagLevel*        a_tagLevel,
    std::vector<Box> a_restrictTagBoxes = std::vector<Box>());

  /// Print info about refinement tagging
  void writeTagInfo() const;
  
/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<string> m_tagVar;       ///< Vector of strings of the components
                                      ///< to tag
  std::vector<Real> m_threshold;      ///< Vector of threshold for refinment
                                      ///< tagging, also counts as hi threshold
  std::vector<Real> m_loThreshold;    ///< Vector of lower thresholds for
                                      ///< value based tagging
  std::vector<int> m_maxLevel;        ///< Vector of max levels each refinement
                                      ///< method can tag up to
  std::vector<int> m_tagType;         ///< Set tag method:
                                      ///< 0 - Gradient based
                                      ///< 1 - Value based
                                      ///< 2 - Vorticity based
                                      ///< 3 - VMS-value based
                                      ///< 4 - VMS-vorticity based
};

#endif /* ! defined _CNSIBCCOMBUSTIONREFERENCE_H_ */
