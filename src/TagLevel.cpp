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
 * \file TagLevel.cpp
 *
 * \brief Member functions for TagLevel
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"

//----- Internal -----//

#include "TagLevel.H"
#include "TagMethod.H"


/*******************************************************************************
 *
 * Class TagLevel: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor
/**
 *//*-----------------------------------------------------------------*/

TagLevel::TagLevel()
  :
  m_tagBufferSize(2)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagLevel::~TagLevel()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Prepend a tag method to the list
/** \param[in]  a_tagMethod
 *                      A method for tagging cells.  Allocate with
 *                      new and do not delete.
 *//*-----------------------------------------------------------------*/

void
TagLevel::prependTagMethod(TagMethod *const a_tagMethod)
{
  m_tagMethodList.push_front(RefCountedPtr<TagMethod>(a_tagMethod));
}

/*--------------------------------------------------------------------*/
//  Append a tag method to the list
/** \param[in]  a_tagMethod
 *                      A method for tagging cells.  Allocate with
 *                      new and do not delete.
 *//*-----------------------------------------------------------------*/

void
TagLevel::appendTagMethod(TagMethod *const a_tagMethod)
{
  m_tagMethodList.push_back(RefCountedPtr<TagMethod>(a_tagMethod));
}
  
/*--------------------------------------------------------------------*/
//  Tag all cells on a level
/** In the default implementation by apply tagging methods to a domain
 *  \param[out] a_tags  Global (to process) IntVectSet of cells to
 *                      refine on this level.
 *  \param[in]  a_boxes Boxes on the level
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
TagLevel::tagCells(IntVectSet&               a_tags,
                   const MappedLevelData&    a_data,
                   LevelGridMetrics&         a_gridMetrics,
                   const MultiBlockCoordSys& a_MBCoordSys,
                   const Real                a_time,
                   const int                 a_level)
{
  // Retrive boxes from the LGM and make sure data is compatible
  const DisjointBoxLayout& boxes = a_gridMetrics.getBoxes();
  CH_assert(boxes.compatible(a_data.rawU().getBoxes()));
  typedef std::list<RefCountedPtr<TagMethod> >::iterator IterType;
  // For each box
  for (DataIterator dit(boxes); dit.ok(); ++dit)
    {
      const Box box = boxes[dit];
      //**FIXME Switch to HashIVS
      IntVectSet localTagsInBox(grow(box, m_tagBufferSize));
      localTagsInBox.makeEmptyBits();
      const IterType tagItEnd = m_tagMethodList.end();
      // Loop through and apply the tag methods.
      for (IterType tagIt = m_tagMethodList.begin(); tagIt != tagItEnd; ++tagIt)
        {
          (*tagIt)->operator()(localTagsInBox,
                               box,
                               dit(),
                               a_data,
                               a_gridMetrics,
                               a_MBCoordSys,
                               a_time,
                               a_level);
        }
      a_tags |= localTagsInBox;
    }
}
