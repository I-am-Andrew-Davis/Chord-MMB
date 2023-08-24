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
 * \file CNSIBCSpatiallyEvolvingShear.H
 *
 * \brief Initial and boundary conditions for mixing layer test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSPATIALLYEVOLVINGSHEAR_H_
#define _CNSIBCSPATIALLYEVOLVINGSHEAR_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for mixing layer test case
/**
 *//*+*************************************************************************/

class CNSIBCSpatiallyEvolvingShear : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCSpatiallyEvolvingShear();

  /// Destructor
  virtual ~CNSIBCSpatiallyEvolvingShear();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCSpatiallyEvolvingShear(const CNSIBCSpatiallyEvolvingShear&);

//--Assignment not permitted

  CNSIBCSpatiallyEvolvingShear &operator=(const CNSIBCSpatiallyEvolvingShear&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

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
                             const Real           a_GlobalKE,
                             const Real           a_globalHelicity) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:
  
  /// Set the imposed BC for primitive variables
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
    const BCInfo&        a_domT) const;

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

  std::vector<IntVect>  m_perturbHi;  ///< Perturbation modes (wavenumbers)
  std::vector<IntVect>  m_perturbLo;  ///< Perturbation modes (wavenumbers)
  std::vector<RealVect> m_phaseHi;    ///< Phases of upper perturbation modes
  std::vector<RealVect> m_phaseLo;    ///< Phases of lower perturbation modes
  std::vector<RealVect> m_modeAmpHi;  ///< Mode amplitude as percent of
                                      ///< single-layer (upper layer)
  std::vector<RealVect> m_modeAmpLo;  ///< Mode amplitude as percent of
                                      ///< single-layer (lower layer)
  Real                  m_layerShift; ///< Shear layer offset from center-line
                                      ///< as percent of sing.-layer thickness
  Real                  m_momThick;   ///< Momentum thickness of shear layer
  Real                  m_rhoThick;   ///< Density transition thickness
  Real                  m_specThick;  ///< Species transition thickness
  Real                  m_loRho;      ///< Density of high layer
  RealVect              m_loU;        ///< Velocity of high layer
  Real                  m_loT;        ///< Temperature of high layer
  bool                  m_singleShear;///< T = single shear, F = double
  bool                  m_reaction;   ///< T = reacting, F = non-reacting
  int                   m_idxStateInletLow;
                                      ///< Lower inlet state
  int                   m_idxStateInletHigh;
                                      ///< Upper inlet state
  std::vector<Real> m_hiMassFractions;///< Upper species concentration
  std::vector<Real> m_loMassFractions;///< Lower Species concentration

  // Methane-Oxygen(air) Reaction case
  int                   m_O2comp;       ///< O2 comp number
  int                   m_CH4comp;      ///< CH4 comp number
  int                   m_N2comp;       ///< N2 comp number
  // Note: gcc ~6.3 won't tolerate emplacement of an FArrayBox so use a pointer
  // (using FArrayBox directly as follows does work on gcc ~8.3)
  // auto ins = m_hackRelax.insert(
  //   { a_boundaryFaceBox, FArrayBox(a_boundaryFaceBox, 1) });
  // This would make destruction very ugly for AMR...
  mutable std::unordered_map<Box, FArrayBox*, CH_Hash::google_CityHash<Box>>
  m_hackRelax;
};

#endif /* ! defined _CNSIBCSPATIALLYEVOLVINGSHEAR_H_ */
