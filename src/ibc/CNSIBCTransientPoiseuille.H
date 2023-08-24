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
 * \file CNSIBCTransientPoiseuille.H
 *
 * \brief Initial and boundary conditions for transient Poiseuille flow
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTRANSIENTPOISEUILLE_H_
#define _CNSIBCTRANSIENTPOISEUILLE_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for transient Poiseuille flow.
/**
 *//*+*************************************************************************/

class CNSIBCTransientPoiseuille : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  /*  Note that Re*nu/L is the velocity (V) according to the Reynolds number.
   *  For SpaceDim > 2, the wall moves in all tangential directions, set by
   *  Vi = sqrt(V^2/(SpaceDim-1))
   */
  CNSIBCTransientPoiseuille(const int a_wallNormalDir = 1);

  /// Destructor
  virtual ~CNSIBCTransientPoiseuille();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTransientPoiseuille(const CNSIBCTransientPoiseuille&);

//--Assignment not permitted

  CNSIBCTransientPoiseuille &operator=(const CNSIBCTransientPoiseuille&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;

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
                             const Real           a_globalKE,
                             const Real           a_globalHelicity) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells
  virtual int exactSol(FArrayBox&              a_Ux,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       const LevelGridMetrics& a_gridMetrics,
                       const FluxBox&          a_unitNormals,
                       const DataIndex&        a_didx,
                       const Real              a_time,
                       const int               a_level) const;

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
  int m_wallNormalDir;                ///< Direction normal to the moving wall
  int m_forceDir;                     ///< Direction the force is applied
  Real m_forceVal;                    ///< Force added to the flow as a source
  Real m_U;                           ///< Bulk velocity
  Real m_refPerc;                     ///< Percent of the domain to refine near
                                      ///< the walls
};

#endif /* ! defined _CNSIBCTRANSIENTPOISEUILLE_H_ */
