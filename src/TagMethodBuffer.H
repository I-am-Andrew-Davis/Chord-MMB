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
 * \file TagMethodBuffer.H
 *
 * \brief Applies buffers to the tagged cells.
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODBUFFER_H_
#define _TAGMETHODBUFFER_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodBuffer
/**
 *   TagMethod class that applies buffers to tagged cells and crops to domains
 *   as necessary
 *
 *//*+*************************************************************************/

class TagMethodBuffer : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodBuffer(const int a_tagBufferSize);

  /// Destructor
  virtual ~TagMethodBuffer();

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells that are a given distance from already tagged cells
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

  int m_tagBufferSize;                ///< Number of extra cells to tag around
                                      ///< a cell already marked for refinement
};

#endif  /* ! defined _TAGMETHODBUFFER_H_ */
