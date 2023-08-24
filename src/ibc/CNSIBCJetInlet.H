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
 * \file CNSIBCJetInlet.H
 *
 * \brief Initial and boundary conditions for jet inlet
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCJETINLET_H_
#define _CNSIBCJETINLET_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for jet inlet
/**
 *//*+*************************************************************************/

class CNSIBCJetInlet : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCJetInlet();

  /// Destructor
  virtual ~CNSIBCJetInlet();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCJetInlet(const CNSIBCJetInlet&);

//--Assignment not permitted

  CNSIBCJetInlet &operator=(const CNSIBCJetInlet&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;
  
/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set boxes for inflow or outflow conditions in mixed boundaries
  virtual int setMixedBC(const Box&           a_boundaryFaceBox,
                         const BoundaryIndex& a_bcIdx,
                         const Box&           a_disjointBox,
                         LevelGridMetrics&    a_gridMetrics,
                         const Real           a_time,
                         const int            a_level,
                         Vector<Box>&         a_boxVect,
                         Vector<BCInfo>&      a_domainBC) const override;

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

  bool m_flatplateSetup;              ///< False- channel flow, true- flatplate
  Box m_slipWall;                     ///< Box of slip walls
  Box m_noSlipWall;                   ///< Box of no-slip walls
  Real m_slipWallFrac;                ///< Percent of the domain length
                                      ///< that is slip wall
};

#endif /* ! defined _CNSIBCJETINLET_H_ */
