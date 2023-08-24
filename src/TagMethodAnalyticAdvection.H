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
 * \file TagMethodAnalyticAdvection.H
 *
 * \brief Tag cells in a refinement of a box on the coarsest level
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODANALYTICADVECTION_H_
#define _TAGMETHODANALYTICADVECTION_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodAnalyticAdvection
/**
 *   TagMethod class that tags all cells in a refinment of a box given on the
 *   coarsest (base) level
 *
 *//*+*************************************************************************/

class TagMethodAnalyticAdvection : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodAnalyticAdvection(Real a_toggleCycle = -1);

  /// Destructor
  virtual ~TagMethodAnalyticAdvection();

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells in baseBox
  virtual void operator()(IntVectSet&               a_tags,
                          const Box&                a_box,
                          const DataIndex&          a_didx,
                          const MappedLevelData&    a_data,
                          LevelGridMetrics&         a_gridMetrics,
                          const MultiBlockCoordSys& a_MBCoordSys,
                          const Real                a_time,
                          const int                 a_level);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  Real m_toggleCycle;
};

#endif  /* ! defined _TAGMETHODANALYTICADVECTION_H_ */
