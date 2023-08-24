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
 * \file TagMethodBoundary.H
 *
 * \brief Tag cells on the boundary of a block
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODBOUNDARY_H_
#define _TAGMETHODBOUNDARY_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodPhysBox
/**
 *   TagMethod class that tags all cells in a physical box.
 *
 *//*+*************************************************************************/

class TagMethodBoundary : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodBoundary(const int            a_idxBlk,
                    const int            a_dir,
                    const Side::LoHiSide a_side,
                    const int            a_normalGrowth = 1);

  /// Destructor
  virtual ~TagMethodBoundary() = default;

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells in a box in physical space
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

  const int m_idxBlk;                 ///< Block
  const int m_dir;                    ///< Direction and
  const Side::LoHiSide m_side;        ///< side of the block
  const int m_normalGrowth;           ///< Growth normal to the surface
};

#endif  /* ! defined _TAGMETHODBOUNDARY_H_ */
