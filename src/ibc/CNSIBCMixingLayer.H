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
 * \file CNSIBCMixingLayer.H
 *
 * \brief Initial and boundary conditions for mixing layer test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCMIXINGLAYER_H_
#define _CNSIBCMIXINGLAYER_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for mixing layer test case
/**
 *//*+*************************************************************************/

class CNSIBCMixingLayer : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCMixingLayer();

  /// Destructor
  virtual ~CNSIBCMixingLayer();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCMixingLayer(const CNSIBCMixingLayer&);

//--Assignment not permitted

  CNSIBCMixingLayer &operator=(const CNSIBCMixingLayer&);


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
                             const Real           a_globalKE,
                             const Real           a_globalHelicity) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells
  virtual int exactSol(
    FArrayBox&        a_Ux,
    const Box&        a_box,
    const Box&        a_disjointBox,
    LevelGridMetrics& a_gridMetrics,
    const Real        a_time,
    const int         a_level) const;

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
                                      ///< as percent of single-layer thickness
  Real                  m_momThick;   ///< Momentum thickness of shear layer
  Real                  m_rhoThick;   ///< Density transition thickness
  Real                  m_specThick;  ///< Species transition thickness
  Real                  m_loRho;      ///< Density of high layer
  RealVect              m_loU;        ///< Velocity of high layer
  Real                  m_loT;        ///< Temperature of high layer
  bool                  m_singleShear;
                                      ///< True == single shear, false == double
};

#endif /* ! defined _CNSIBCMIXINGLAYER_H_ */
