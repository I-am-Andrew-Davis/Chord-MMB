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
 * \file TagMethodTBL.H
 *
 * \brief Tag cells in a turbulent boundary-layer
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODTBL_H_
#define _TAGMETHODTBL_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodTBL
/**
 *   TagMethod class that tags all cells in a turbulent boundary-layer.
 *
 *//*+*************************************************************************/

class TagMethodTBL : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodTBL(const int            a_idxBlk,
               const int            a_dir,
               const Side::LoHiSide a_side,
               const int            a_streamwiseDir,
               const Real           a_delta,
               const Real           a_ReynoldsNumber);

  /// Destructor
  virtual ~TagMethodTBL() = default;

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells in turbulent boundary-layer
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
  const int m_streamwiseDir;          ///< Streamwise direction of flow
  const Real m_delta;                 ///< Initial boundary-layer thicknes
  const Real m_Re;                    ///< Initial Reynolds number
};

#endif  /* ! defined _TAGMETHODTBL_H_ */
