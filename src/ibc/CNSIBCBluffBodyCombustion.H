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
 * \file CNSIBCBluffBodyCombustion.H
 *
 * \brief Initial and boundary conditions for Bluff-Body Combustion test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCBLUFFBODYCOMBUSTION_H_
#define _CNSIBCBLUFFBODYCOMBUSTION_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for Bluff-Body combustion test case
/**
 *//*+*************************************************************************/

class CNSIBCBluffBodyCombustion : public CNSIBCGeneralized
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCBluffBodyCombustion(const bool a_readInputs = true);

  /// Destructor
  virtual ~CNSIBCBluffBodyCombustion();

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCBluffBodyCombustion(const CNSIBCBluffBodyCombustion&) = delete;

//--Assignment not permitted

  CNSIBCBluffBodyCombustion &operator=(const CNSIBCBluffBodyCombustion&) = delete;
  
/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const override;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const override;

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const override;

  /// Add source term
  virtual void addSourceTerm(FArrayBox&           a_sourceFab,
                             FArrayBox&           a_invDtFab,
                             const FArrayBox&     a_Wcell,
                             const FArrayBox&     a_UcellAvg,
                             const FluxBox&       a_WfaceAvgFxb,
                             const ProblemDomain& a_domain,
                             LevelGridMetrics&    a_gridMetrics,
                             const Real           a_time,
                             const Real           a_stageWeight,
                             const int            a_level,
                             const Box&           a_disjointBox,
                             const Box&           a_solveBox,
                             const DataIndex&     a_dataIndx,
                             const Real           a_GlobalKE,
                             const Real           a_globalHelicity) const;
  
/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:
  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const override;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  int m_idxStateRegion2;              ///< State for initialization
  RealVect m_region2_lo;              ///< Lo-bounded physical location of region-2
  RealVect m_region2_hi;              ///< Hi-bounded physical location of region-2
};

#endif /* ! defined _CNSIBCBLUFFBODYCOMBUSTION_H_ */


