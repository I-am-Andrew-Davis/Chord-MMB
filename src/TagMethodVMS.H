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
 * \file TagMethodVMS.H
 *
 * \brief Implements variational-multiscale-based tagging
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODVMS_H_
#define _TAGMETHODVMS_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodVMS
/**
 *   TagMethod class that implements variational-multiscale-based tagging
 *
 *   \note
 *   <ul>
 *     <li> Consider improvements for metrics
 *   </ul>
 *   \note
 *   <ul>
 *     <li> This should, ideally, only be called once for multiple components.
 *     <li> This tagging method should be reworked considerably.
 *   <\ul>
 *
 *//*+*************************************************************************/

class TagMethodVMS : public TagMethod
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodVMS(const int        a_comp,
               const Real       a_loThreshold,
               const Real       a_upThreshold,
               const bool       a_usePrimitive    = true,
               const bool       a_insideThreshold = false,
               const bool       a_useVorticity    = false,
               const int        a_maxLevel = -1);

  /// Destructor
  virtual ~TagMethodVMS();

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells with values that exceed a threshold
  virtual void operator()(IntVectSet&               a_tags,
                          const Box&                a_box,
                          const DataIndex&          a_didx,
                          const MappedLevelData&    a_data,
                          LevelGridMetrics&         a_gridMetrics,
                          const MultiBlockCoordSys& a_MBCoordSys,
                          const Real                a_time,
                          const int                 a_level);

  /// Separate scales of variables of interest
  void scaleSep(FArrayBox&                a_coarserData,
                const FArrayBox&          a_origData,
                const Box&                a_box,
                LevelGridMetrics&         a_gridMetrics,
                const MultiBlockCoordSys& a_MBCoordSys);

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  int  m_comp;                        ///< Component to consider
  Real m_loThreshold;                 ///< Lower threshold for refinement
  Real m_upThreshold;                 ///< Upper threshold for refinement
  bool m_usePrimitive;                ///< Compute gradients of primitive
                                      ///< variables
  bool m_insideThreshold;             ///< True - Tag values that fall between
                                      ///< the thresholds
  bool m_useVorticity;                ///< Tag based on vorticity
  int  m_maxLevel;                    ///< Maximum level allowable for
                                      ///< refinement, set to -1 for all levels
};

#endif  /* ! defined _TAGMETHODVMS_H_ */
