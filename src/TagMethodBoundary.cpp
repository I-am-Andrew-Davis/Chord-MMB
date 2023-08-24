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
 * \file TagMethodBoundary.cpp
 *
 * \brief Member functions for TagMethodBoundary
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

#include "TagMethodBoundary.H"


/*******************************************************************************
 *
 * Class TagMethodBoundary: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_Xlo   Low physical coordinates of the box
 *  \param[in]  a_Xhi   High physical coordinates of the box
 *  \param[in]  a_fastExlusion
 *                      T - performs a quick test to see if a box
 *                          should be further analysed for having
 *                          cells in the tag box.  For very complex
 *                          geometry, where it is difficult to find
 *                          the minbox, set to false.
 *//*-----------------------------------------------------------------*/

TagMethodBoundary::TagMethodBoundary(const int            a_idxBlk,
                                     const int            a_dir,
                                     const Side::LoHiSide a_side,
                                     const int            a_normalGrowth)
  :
  m_idxBlk(a_idxBlk),
  m_dir(a_dir),
  m_side(a_side),
  m_normalGrowth(a_normalGrowth)
{
  CH_assert(a_normalGrowth > 0);
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Tag all cells in a box in physical space
/** If m_fastExlcusion is T, a minbox representing the physical space
 *  of a_box is estimated and intersected with the tag box to exclude
 *  boxes from detailed analysis by sampling only corners and face
 *  centers.
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
TagMethodBoundary::operator()(IntVectSet&               a_tags,
                              const Box&                a_box,
                              const DataIndex&          a_didx,
                              const MappedLevelData&    a_data,
                              LevelGridMetrics&         a_gridMetrics,
                              const MultiBlockCoordSys& a_MBCoordSys,
                              const Real                a_time,
                              const int                 a_level)
{
  const int idxBlk = a_gridMetrics.getBoxes().blockIndex(a_didx);
  CH_assert(idxBlk != -1);
  if (idxBlk == m_idxBlk)
    {
      Box blockBox = a_MBCoordSys.mappingBlocks()[idxBlk];
      blockBox.growDir(m_dir, m_side, -1);
      if (!blockBox.contains(a_box))
        {
          Box tagBox = adjCellBox(a_box, m_dir, m_side, -m_normalGrowth);
          tagBox &= a_box;
          a_tags |= tagBox;
        }
    }
}
