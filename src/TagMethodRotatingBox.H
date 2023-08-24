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
 * \file TagMethodRotatingBox.H
 *
 * \brief Tag cells in a refinement of a box on the coarsest level
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODROTATINGBOX_H_
#define _TAGMETHODROTATINGBOX_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodRotatingBox
/**
 *   TagMethod class that tags all cells in a refinment of a box given on the
 *   coarsest (base) level
 *
 *//*+*************************************************************************/

class TagMethodRotatingBox : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodRotatingBox(const Real& a_radius,
                       const Real& a_cycleTime,
                       const int   a_boxSize);

  /// Destructor
  virtual ~TagMethodRotatingBox();

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

  Real m_radius;             ///< Radius of rotation
  Real m_cycleTime;          ///< Time to compleate one rotation
  int m_boxSize;             ///< Size of box to refine for this level
};

#endif  /* ! defined _TAGMETHODROTATINGBOX_H_ */
