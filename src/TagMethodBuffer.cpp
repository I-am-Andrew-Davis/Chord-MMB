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
 * \file TagMethodBuffer.cpp
 *
 * \brief Member functions for TagMethodBuffer
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "Box.H"
#include "DataIndex.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"
#include "AdjacentBlock.H"

//----- Internal -----//

#include "TagMethodBuffer.H"


/*******************************************************************************
 *
 * Class TagMethodBuffer: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_tagBufferSize
 *                      Number of extra cells to tag around a cell
 *                      already marked for refinement (i.e., buffer)
 *//*-----------------------------------------------------------------*/

TagMethodBuffer::TagMethodBuffer(const int a_tagBufferSize)
  :
  m_tagBufferSize(a_tagBufferSize)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodBuffer::~TagMethodBuffer()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Tag all cells that are a given distance from already tagged cells
/** \param[out] a_tags  IntVectSet of cells to refine local to
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
TagMethodBuffer::operator()(IntVectSet&               a_tags,
                            const Box&                a_box,
                            const DataIndex&          a_didx,
                            const MappedLevelData&    a_data,
                            LevelGridMetrics&         a_gridMetrics,
                            const MultiBlockCoordSys& a_MBCoordSys,
                            const Real                a_time,
                            const int                 a_level)
{
  // Get the block index from the layout and get the latter from LGM
  const DisjointBoxLayout& boxes = a_gridMetrics.getBoxes();
  const int idxBlk = boxes.blockIndex(a_didx);
  CH_assert(idxBlk != -1);
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(idxBlk);
  a_tags.grow(m_tagBufferSize);
  // Intersect with problem domain.
  Box tagBox = grow(a_box, m_tagBufferSize);
  tagBox &= blockDomain;
  a_tags &= tagBox;
  if (a_gridMetrics.isMultiBlock())
    {
      // If the tags extend across a block boundary, we need to transform them
      // to the other block
      a_tags.recalcMinBox();
      const Box minBox = a_tags.minBox();
      const Box blockBox = blockDomain.domainBox();
      if (blockBox.contains(minBox)) return;
      IntVectSet MBtags = a_tags;
      MBtags &= blockBox;
      // Walk through all blocks that may intersect with minBox and transform
      // tags to that location.
      AdjacentBlock::forEachApply(
        a_MBCoordSys,
        minBox,
        idxBlk,
        [&a_MBCoordSys, &a_tags, &MBtags]
         // The parts of minBox that intersect the neighbor block and in that
         // that space (for complex edges and corners, it is not quite that
         // simple---see AdjacentBlock for more info).
        (const Box& a_box,
         // The index of the neighbor box
         const int a_idxBlk,
         // The total transformation to the neighbor block
         const IndicesTransformation& a_trfm)
          {
            // Transform back to the original space
            Box box = a_trfm.transformBack(a_box);
            IntVectSet tags = a_tags;
            tags &= box;
            for (IVSIterator ivsIt(tags); ivsIt.ok(); ++ivsIt)
              {
                IntVect tag = a_trfm.transformFwd(ivsIt());
                MBtags |= tag;
              }
          });
      // All tags in a_tags should now be in MBtags
      // CH_assert(a_tags.numPts() == MBtags.numPts());
      /* In simple cases, the above is true.  However, tags in a_tags are not
         cropped across edges or corner, only across faces.  One could use
         AdjacentBlock to properly crop all tags outside of physical boundaries.
         However, the transformation of tags does the right thing.  I.e., MBtags
         is correct and neglects tags the go across a physical boundary.  So
         making the cropping correct for a_tags only facilitates the above
         assertion at some expense.  Probably not worth it.
      */
      a_tags = MBtags;
    }
}
