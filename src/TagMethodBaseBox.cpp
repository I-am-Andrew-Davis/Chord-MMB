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
 * \file TagMethodBaseBox.cpp
 *
 * \brief Member functions for TagMethodBaseBox
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "Box.H"
#include "DataIndex.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"

//----- Internal -----//

#include "TagMethodBaseBox.H"
#include "CRDparam.H"


/*******************************************************************************
 *
 * Class TagMethodBaseBox: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_tagBaseBox
 *                      Tag any cells in this box.  This box is
 *                      defined on the coarsest (base) level.  The
 *                      tagBaseBox must be contained in the problem
 *                      domain.
 *//*-----------------------------------------------------------------*/

TagMethodBaseBox::TagMethodBaseBox(const Box& a_tagBaseBox,
                                   const Mode a_mode)
  :
  m_tagBaseBox(a_tagBaseBox),
  m_mode(a_mode)
{
  CH_assert(!a_tagBaseBox.isEmpty());
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodBaseBox::~TagMethodBaseBox()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Tag all cells in baseBox
/** If the baseBox is face-centered, it is first refined to the
 *  current level, grown by 1 in face-centered directions, and then
 *  enclosed cells are found.  This way, one can specify a thin
 *  layer of cells at the resolution of the current level.
 *  \param[out] a_tags  IntVectSet of cells to refine local to
 *                      'a_box'.
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
 *//*-----------------------------------------------------------------*/

void
TagMethodBaseBox::operator()(IntVectSet&               a_tags,
                             const Box&                a_box,
                             const DataIndex&          a_didx,
                             const MappedLevelData&    a_data,
                             LevelGridMetrics&         a_gridMetrics,
                             const MultiBlockCoordSys& a_MBCoordSys,
                             const Real                a_time,
                             const int                 a_level)
{
  Box tagBox(m_tagBaseBox);
  tagBox.refine(CRDparam::g_refFromBase[a_level]);
  if (!tagBox.cellCentered())
    {
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (tagBox.type(dir) == IndexType::NODE)
            {
              tagBox.grow(dir, 1);
            }
        }
      tagBox.enclosedCells();
    }

  switch (m_mode)
    {
    case ModeAddTagsInsideBox:
      tagBox &= a_box;
      a_tags |= tagBox;
      break;
    case ModeExcludeTagsOutsideBox:
      if (tagBox.intersectsNotEmpty(a_box))
        {
          IntVectSet removeTags(a_box);
          removeTags -= tagBox;
          a_tags -= removeTags;
        }
      else
        {
          a_tags -= a_box;
        }
      break;
    }
}
