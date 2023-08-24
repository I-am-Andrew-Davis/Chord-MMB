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
 * \file TagMethodGradient.H
 *
 * \brief Implements gradient-based tagging
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODGRADIENT_H_
#define _TAGMETHODGRADIENT_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodGradient
/**
 *   TagMethod class that implements gradient-based tagging
 *
 *   \note
 *   <ul>
 *     <li> Consider improvements for metrics
 *   </ul>
 *
 *//*+*************************************************************************/

class TagMethodGradient : public TagMethod
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodGradient(const int        a_comp,
                    const Real       a_threshold,
                    const bool       a_usePrimitive,
                    const Real       a_loVal       = 0.,
                    std::vector<Box> a_restrictBox = std::vector<Box>());

  /// Destructor
  virtual ~TagMethodGradient();

  // Use synthesized copy and assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Tag all cells with gradients that exceed a threshold
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

  int m_comp;                         ///< Component to consider
  Real m_threshold;                   ///< Threshold for refinement
  Real m_maxGradMag;                  ///< Maximum gradient found on the box.
                                      ///< Returned to the user for help with
                                      ///< setting the threshold
  bool m_usePrimitive;                ///< Compute gradients of primitive
                                      ///< variables
  Real m_loVal;                       ///< If value is below this, do not tag
                                      ///< even if gradient is high enough
  std::vector<Box> m_restrictBox;     ///< Boxes where refinement is restricted
};

#endif  /* ! defined _TAGMETHODGRADIENT_H_ */
