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
 * \file TagMethodValue.H
 *
 * \brief Implements value-based tagging
 *
 *//*+*************************************************************************/

#ifndef _TAGMETHODVALUE_H_
#define _TAGMETHODVALUE_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagMethod.H"


/*******************************************************************************
 */
///  TagMethodValue
/**
 *   TagMethod class that implements value-based tagging
 *
 *   \note
 *   <ul>
 *     <li> Consider improvements for metrics
 *   </ul>
 *
 *//*+*************************************************************************/

class TagMethodValue : public TagMethod
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

enum ValueType : int
{
  ValueTypeConservative = 0,
  ValueTypePrimitive,
  ValueTypeVorticity,
  ValueTypeHydrocarbonFlame
};


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TagMethodValue(const int        a_comp,
                 const Real       a_loThreshold,
                 const Real       a_upThreshold,
                 const bool       a_insideThreshold = true,
                 const ValueType  a_valueType       = ValueTypeConservative,
                 std::vector<Box> a_restrictBox     = std::vector<Box>());

  /// Destructor
  virtual ~TagMethodValue();

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


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  int m_comp;                         ///< Component to consider
  Real m_loThreshold;                 ///< Lower threshold for refinement
  Real m_upThreshold;                 ///< Upper threshold for refinement
  bool m_insideThreshold;             ///< True - Tag values that fall between
                                      ///< the thresholds
  ValueType m_valueType;              ///< Enum on type of value
  std::vector<Box> m_restrictBox;     ///< Boxes where refinement is restricted
                                      ///< (DEPRECIATED)
  int m_cOH;                          ///< Consv. variable component of OH for
                                      ///< tagging hydrocarbon flames
  int m_cCH2O;                        ///< Consv. variable component of CH2O for
                                      ///< tagging hydrocarbon flames
};

#endif  /* ! defined _TAGMETHODVALUE_H_ */
