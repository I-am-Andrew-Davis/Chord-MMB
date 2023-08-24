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
 * \file TagMethod.H
 *
 * \brief Abstract base object for tagging methods
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHOD_H_
#define _TAGMETHOD_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class IntVectSet;
class DataIndex;
class MappedLevelData;
class LevelGridMetrics;
class MultiBlockCoordSys;
#include "NamespaceFooter.H"

//--Forward declarations

class TagMethod;


/*******************************************************************************
 */
///  TagMethod
/**
 *   Abstract base class for tag method functors.  All methods must derive
 *   from this class.
 *
 *//*+*************************************************************************/

class TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  TagMethod()
    { }

  /// Destructor
  virtual ~TagMethod()
    { }

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Apply some tagging method to the cells in a box.
  /** \param[out] a_tags  IntVectSet of cells to refine in 'a_box'
   *  \param[in]  a_box   Disjoint box in which cells should be marked
   *                      for refinement
   *  \param[in]  a_didx  Data index of the current box
   *  \param[in]  a_data  Data on the level, both mapped and physical.
   *  \param[in]  a_gridMetrics
   *                      Grid metrics for the level
   *  \param[out] a_gridMetrics
   *                      Only caches should be modified
   *  \param[in]  a_MBCoordSys
   *                      Coordinate system describing the mapping
   *  \param[in]  a_time  Current solution time
   *  \param[in]  a_level Index of the AMR level
   */
  virtual void operator()(IntVectSet&               a_tags,
                          const Box&                a_box,
                          const DataIndex&          a_didx,
                          const MappedLevelData&    a_data,
                          LevelGridMetrics&         a_gridMetrics,
                          const MultiBlockCoordSys& a_MBCoordSys,
                          const Real                a_time,
                          const int                 a_level) = 0;
};

#endif  /* ! defined _TAGMETHOD_H_ */
