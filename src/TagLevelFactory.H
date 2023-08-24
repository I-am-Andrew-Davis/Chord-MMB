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
 * \file TagLevelFactory.H
 *
 * \brief Factory object for cell tagging strategies
 *
 *//*+*************************************************************************/

#ifndef _TAGLEVELFACTORY_H_
#define _TAGLEVELFACTORY_H_

//----- Chombo Library -----//

#include "RefCountedPtr.H"
#include "UsingNamespace.H"

//--Forward declarations

class TagLevel;


/*******************************************************************************
 */
///  TagLevelFactory
/**
 *   This class creates TagLevel objects.  The default implementation allows
 *   for a vector of tag level objects and a sequential map of which AMR levels
 *   they belong to.
 *
 *//*+*************************************************************************/

class TagLevelFactory
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Construct by reading the input file
  TagLevelFactory(const int a_tagBufferSize);

  /// Constructor for single TagLevel object
  TagLevelFactory(TagLevel *const a_tagLevel);

  /// Constructor for multiple TagLevel objects
  TagLevelFactory(const std::vector<TagLevel*> a_tagLevelVec,
                  const std::vector<int>       a_levelMapVec);

  /// Destructor
  virtual ~TagLevelFactory();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  TagLevelFactory(const TagLevelFactory&);

//--Assignment not permitted

  TagLevelFactory &operator=(const TagLevelFactory&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Create a new TagLevel object
  virtual const RefCountedPtr<TagLevel>& new_tagLevel(
    const int a_level,
    const int a_fnRefRatio) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Read details of a tag method and add to appropriate level(s)
  virtual void readTagMethod(ParmParse&                a_ppTAG,
                             const char*               a_tagLevelName,
                             const int                 a_level,
                             const int                 a_idxMethod,
                             const char*               a_methodName,
                             std::map<int, TagLevel*>& a_tagLevels);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<RefCountedPtr<TagLevel> > m_tagLevelVec;
                                      ///< A vector of TagLevel objects.  Each
                                      ///< is used for tagging on some levels
  std::vector<int> m_levelMapVec;     ///< A vector of AMR level indices.
                                      ///< If an AMR level index (idx) is in
                                      ///< the range, (m_levelMapVec[i] <= idx
                                      ///< < m_levelMapVec[i+1]), then
                                      ///< m_tagLevelVec[i] is used to tag cells
                                      ///< on that level
};

#endif
