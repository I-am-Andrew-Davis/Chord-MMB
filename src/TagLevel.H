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
 * \file TagLevel.H
 *
 * \brief Implements tagging on a level
 *
 *//*+*************************************************************************/

#ifndef _TAGLEVEL_H_
#define _TAGLEVEL_H_

//----- Standard Library -----//

#include <list>

//----- Chombo Library -----//

#include "REAL.H"
#include "RefCountedPtr.H"
#include "UsingNamespace.H"
#include "NamespaceVar.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class IntVectSet;
class DisjointBoxLayout;
class MappedLevelData;
class LevelGridMetrics;
class MultiBlockCoordSys;
#include "NamespaceFooter.H"

//--Forward declarations

class TagMethod;


/*******************************************************************************
 */
///  TagLevel
/**
 *   This class creates TagLevel objects.
 *   There are a large number of ways to modify the behavior of tagging
 *   <ol>
 *     <li> A single default TagLevel object is used for all levels with
 *          the same TagMethods.
 *     <li> Default TagLevel objects are used for all levels but with
 *          different TagMethods.  You will need to derive from
 *          TagLevelFactory and redefine new_tagLevel to create TagLevels
 *          with different methods on different levels.
 *     <li> You can derive from this class to create a customized TagLevel
 *          object to use on all levels by the default TagLevelFactory.
 *     <li> You can derive from this class to create multiple customized
 *          TagLevel objects and derive from TagLevelFactory assign those
 *          different objects to different levels.
 *   </ol>
 *   Start by trying to simply use the same TagMethod objects on all levels
 *   and add complexity as needed.
 *   
 *   \note
 *   <ul>
 *     <li> WARNING - This class is not fully thought out yet and may change
 *                    significantly to incorporate mapped and multi-block
 *                    grids.
 *   </ul>
 *
 *//*+*************************************************************************/

class TagLevel
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  TagLevel();

  /// Destructor
  virtual ~TagLevel();

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Set the tag buffer size (temporary until HashIVS is ready)
  void FIXMEsetTagBufferSize(const int a_tagBufferSize)
    {
      m_tagBufferSize = a_tagBufferSize;
    }

  /// Prepend a tag method to the list
  void prependTagMethod(TagMethod *const a_tagMethod);

  /// Append a tag method to the list
  void appendTagMethod(TagMethod *const a_tagMethod);

  /// Tag all cells on a level
  virtual void tagCells(IntVectSet&               a_tags,
                        const MappedLevelData&    a_data,
                        LevelGridMetrics&         a_gridMetrics,
                        const MultiBlockCoordSys& a_MBCoordSys,
                        const Real                a_time,
                        const int                 a_level);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::list<RefCountedPtr<TagMethod> > m_tagMethodList;
                                      ///< List of functors that tag cells in
                                      ///< a box
  int m_tagBufferSize;                ///< Needs to be known to size the DIVS
                                      ///< correctly.  FIXME Once HashIVS is
                                      ///< available, this can be removed
};

#endif  /* ! defined _TAGLEVEL_H_ */
