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
 * \file CNSIBCGeneralizedSingleBlock.H
 *
 * \brief Initial and boundary conditions for a generalized case on
 *        a single block domain
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCGENERALIZEDSINGLEBLOCK_H_
#define _CNSIBCGENERALIZEDSINGLEBLOCK_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions a generalized case on a single block domain
/**
 *   This is just a simple interface to the true generalized IBC for when
 *   we don't want to bother with multi-block cases
 *//*+*************************************************************************/

class CNSIBCGeneralizedSingleBlock : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCGeneralizedSingleBlock();

  /// Destructor
  virtual ~CNSIBCGeneralizedSingleBlock();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCGeneralizedSingleBlock(const CNSIBCGeneralizedSingleBlock&);

//--Assignment not permitted

  CNSIBCGeneralizedSingleBlock &operator=(const CNSIBCGeneralizedSingleBlock&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

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
  
};

#endif /* ! defined _CNSIBCGENERALIZEDSINGLEBLOCK_H_ */
