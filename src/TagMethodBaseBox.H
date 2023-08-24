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
 * \file TagMethodBaseBox.H
 *
 * \brief Tag cells in a refinement of a box on the coarsest level
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODBASEBOX_H_
#define _TAGMETHODBASEBOX_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodBaseBox
/**
 *   TagMethod class that tags all cells in a refinment of a box given on the
 *   coarsest (base) level
 *
 *//*+*************************************************************************/

class TagMethodBaseBox : public TagMethod
{


/*==============================================================================
 * Types
 *============================================================================*/

public:

  enum Mode
  {
    ModeAddTagsInsideBox,
    ModeExcludeTagsOutsideBox
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodBaseBox(const Box& a_tagBaseBox,
                   const Mode a_mode = ModeAddTagsInsideBox);

  /// Destructor
  virtual ~TagMethodBaseBox();

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

  Box m_tagBaseBox;                   ///< Tag any cells in this box.  This box
                                      ///< is defined on the coarsest (base)
                                      ///< level.
  Mode m_mode;                        ///< Mode of operation
};

#endif  /* ! defined _TAGMETHODBASEBOX_H_ */
