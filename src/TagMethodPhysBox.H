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
 * \file TagMethodPhysBox.H
 *
 * \brief Tag cells in a physical box
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODPHYSBOX_H_
#define _TAGMETHODPHYSBOX_H_

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

class TagMethodPhysBox : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodPhysBox(const RealVect& a_Xlo,
                   const RealVect& a_Xhi,
                   const bool      a_fastExclusion = true);

  /// Destructor
  virtual ~TagMethodPhysBox() = default;

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
 * Protected member functions
 *============================================================================*/

protected:

  /// Is a physical-space coordinate contained in the tag box?
  bool contains(const RealVect& a_X) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  const RealVect m_Xlo;               ///< Lower corner of box
  const RealVect m_Xhi;               ///< Upper corner of box
  const bool m_fastExclusion;         ///< Perform a quick test to see if any
                                      ///< tags might be found in a_box
};

/*--------------------------------------------------------------------*/
//  Is a physical-space coordinate contained in the tag box?
/** \param[in] a_X      Physical space coordinate
 *//*-----------------------------------------------------------------*/

inline bool
TagMethodPhysBox::contains(const RealVect& a_X) const
{
  return a_X >= m_Xlo && a_X <= m_Xhi;
}

#endif  /* ! defined _TAGMETHODPHYSBOX_H_ */
