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
 * \file PatchCNSOp.H
 *
 * \brief Patch operator for compressible Navier-Stokes
 *
 *//*+*************************************************************************/

#ifndef _PATCHCNSOP_H_
#define _PATCHCNSOP_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "Inertial4thOrderOp.H"
#include "ViscousTensor4thOrderOp.H"

//--Forward declarations

#include "CRDhardconfig.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class ProblemDomain;
using BlockDomain = ProblemDomain;
class FArrayBox;
class FluxBox;
class DataIndex;
class LevelFluxRegister;
class LevelGridMetrics;
#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"
namespace stc
{

template <typename T,
          array_size_type N,
          template<typename, array_size_type> class Impl>
class VectorImpl;

}
#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"
using RealVect = stc::VectorImpl<Real, SpaceDim, stc::DefaultImpl>;
#include "NamespaceFooter.H"


/*******************************************************************************
 */
///  Patch operator for compressible Navier-Stokes
/**
 *//*+*************************************************************************/

class PatchCNSOp
{


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  PatchCNSOp(LevelGridMetrics& a_levelGridMetrics);

  // Use synthesized destructor

  // Copy, move, and assignment not permitted
  PatchCNSOp(const PatchCNSOp&) = delete;
  PatchCNSOp(PatchCNSOp&&) = delete;
  PatchCNSOp& operator=(const PatchCNSOp&) = delete;
  PatchCNSOp& operator=(PatchCNSOp&&) = delete;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:                 

  /// Compute primitive cell point and average states
  void computeWcell(FArrayBox&         a_UcellPntFab,
                    FArrayBox&         a_WcellPntFab,
                    FArrayBox&         a_WcellAvgFab,
                    FArrayBox&         a_WfromUavgFab,
                    const DataIndex&   a_dataIndx,
                    const Box&         a_box,
                    const BlockDomain& a_domain,
                    const FArrayBox&   a_UcellAvgFab,
                    FArrayBox&         a_WOld) const;

  /// Compute the primitive cell-centered state from conservative state
  void computeWstate(FArrayBox&         a_WcellPntFab,
                     const FArrayBox&   a_UcellAvgFab,
                     FArrayBox&         a_WOld,
                     const Box&         a_box,
                     const BlockDomain& a_domain,
                     bool               a_cellAvgInput = true,
                     bool               a_fourthOrder = true) const;

  /// Compute the face-average primitive state
  //  This is the heart of the canonical inertial solve with limited
  //  reconstruction
  void computeWfaceAvg(const Box&         a_disjointBox,
                       const Box&         a_box,
                       const BlockDomain& a_domain,
                       FluxBox&           a_WfaceAvgFxb,
                       const FArrayBox&   a_WcellAvgFab,
                       const FArrayBox&   a_WfromUavgFab,
                       FArrayBox&         a_flattening,
                       const FluxBox&     a_unitNormalFxb,
                       const DataIndex&   a_dataIndx,
                       const Real         a_time,
                       const Real         a_dt,
                       const Real         a_stageWeight,
                       const bool         a_setFlattening,
                       const int          a_level) const;

  /// Compute the inertial (hyperbolic) flux
  void addInertialFlux(const Box&                          a_disjointBox,
                       const Box&                          a_box,
                       const BlockDomain&                  a_domain,
                       FluxBox&                            a_flux,
                       FluxBox&                            a_fluxFromWavg,
                       FluxBox&                            a_WfaceAvgFxb,
                       FluxBox&                            a_WfacePntFxb,
                       FArrayBox&                          a_WcellAvgFab,
                       const FArrayBox&                    a_WcellPntFab,
                       FluxBox&                            a_faceAvgPlotFxb,
                       FArrayBox&                          a_flattening,
                       stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellData,
                       stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
                       const FArrayBox&                    a_U,
                       const FluxBox&                      a_unitNormalFxb,
                       const DataIndex&                    a_dataIndx,
                       const Real                          a_time,
                       const Real                          a_dt,
                       const Real                          a_prevDt,
                       const Real                          a_stageWeight,
                       const bool                          a_setFlattening,
                       const int                           a_level) const;

  /// Compute the viscous (elliptic) flux
  void addViscousFlux(const Box&         a_box,
                      const BlockDomain& a_domain,
                      FluxBox&           a_FluxfaceAvgFxb,
                      FluxBox&           a_FluxfacePntFxb,
                      FArrayBox&         a_invDtFab,
                      FArrayBox&         a_turbSourceAvgFab,
                      FArrayBox&         a_WcellAvgFab,
                      FArrayBox&         a_WcellPntFab,
                      FArrayBox&         a_timeAvgFab,
                      FluxBox&           a_faceAvgPlotFxb,
                      FluxBox&           a_WfaceAvgFxb,
                      FluxBox&           a_WfacePntFxb,
                      FluxBox&           a_facePntVelGradFxb,
                      FluxBox&           a_stressFluxFxb,
                      const FluxBox&     a_facePntDeltaC,
                      const FluxBox&     a_faceCoord,
                      const FluxBox&     a_unitNormalsFxb,
                      const DataIndex&   a_dataIndx,
                      const Real         a_dt,
                      const Real         a_time,
                      const Real         a_timeAfterFilterInit,
                      const int          a_level,
                      Real&              a_minDiffDt,
                      IntVect&           a_minDiffDtCell) const;
                         
  /// Compute flux due to artificial viscosity and add to a_JUnew
  void addArtificialViscosity(const Box&         a_box,
                              const BlockDomain& a_domain,
                              FArrayBox&         a_JUnew,
                              FluxBox&           a_NtF,
                              const FArrayBox&   a_Uold,
                              const FArrayBox&   a_WOld,
                              const FluxBox&     a_N,
                              const FArrayBox&   a_J,
                              const FluxBox&     a_unitNormals,
                              const RealVect&    a_dx,
                              const Real         a_dt,
                              const Real         a_time,
                              const int          a_level) const;

  /// Compute correction to the species equations and add to a_JUnew
  void speciesCorrection(const Box&         a_box,
                         const BlockDomain& a_domain,
                         FArrayBox&         a_JU,
                         const RealVect&    a_dx) const;

  /// Compute the divergence of a flux
  void fluxDivergence(const Box&      a_box,
                      FArrayBox&      a_rhs,
                      const FluxBox&  a_flux,
                      const RealVect& a_dx) const;

  /// Update flux registers
  void updateFluxRegisters(LevelFluxRegister& a_crFluxRegister,
                           LevelFluxRegister& a_fnFluxRegister,
                           const FluxBox&     a_flux,
                           const Box&         a_workingBox,
                           const DataIndex&   a_dataIdx,
                           const Real         a_weight,
                           const bool         a_hasCoarserGrid,
                           const bool         a_hasFinerGrid) const;

  /// Add source term
  void addSourceTerm(FArrayBox&         a_sourceFab,
                     FArrayBox&         a_invDtFab,
                     const FArrayBox&   a_Wcell,
                     const FArrayBox&   a_UcellAvg,
                     const FluxBox&     a_WfaceAvgFxb,
                     const Real         a_time,
                     const Real         a_stageWeight,
                     const int          a_level,
                     const Box&         a_disjointBox,
                     const BlockDomain& a_domain,
                     const Box&         a_solveBox,
                     const DataIndex&   a_dataIndx,
                     const Real         a_globalKE,
                     const Real         a_globalHelicity) const;

  /// Solve for the mapped, cell-averaged source term
  void solveJSfromS(FArrayBox&         a_cellAvgJS,
                    FArrayBox&         a_cellPntS,
                    const Box&         a_boxJS,
                    const Box&         a_boxS,
                    const DataIndex&   a_dataIndx,
                    const BlockDomain& a_domain,
                    const bool         a_fourthOrder) const;

  /// Compute the max norm of JU-\hat{JU}
  Real computeUminusUhatMaxNorm(const FArrayBox& a_JU,
                                const FArrayBox& a_JUhat,
                                const Box&       a_box,
                                const Interval&  a_interval) const;

//--Static members

  /// Things to do before calling the Riemann solver
  static void preRiemann(FArrayBox&     a_WLeft,
                         FArrayBox&     a_WRight,
                         const FluxBox& a_unitNormalFxb,
                         const int      a_dir,
                         const Box&     a_box);

  /// Things to do after calling the Riemann solver
  static void postRiemann(FArrayBox&     a_Wface,
                          const FluxBox& a_unitNormalFxb,
                          const int      a_dir,
                          const Box&     a_box);

  /// Forward basis transform on mapped grids
  static void forwardBasisTransform(const Box&     a_box,
                                    FArrayBox&     a_W,
                                    const FluxBox& a_unitNormalFxb,
                                    const int      a_dir);

  /// Reverse basis transform on mapped grids
  static void reverseBasisTransform(const Box&     a_box,
                                    FArrayBox&     a_W,
                                    const FluxBox& a_unitNormalFxb,
                                    const int      a_dir);

  /// Compute the primitive cell-centered state from conservative state
  static void computeWpntCell(FArrayBox&         a_WcellPntFab,
                              const FArrayBox&   a_UcellAvgFab,
                              const FArrayBox&   a_WOld,
                              const Box&         a_box,
                              const BlockDomain& a_domain,
                              bool               a_cellAvgInput = true,
                              bool               a_fourthOrder = true);

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

//--References data defined in AMRLevelCNS.  These are all set during
//--construction and are fully defined by AMRLevelCNS

  LevelGridMetrics& m_levelGridMetrics;
                                      ///< Grid metrics for the level.

//--New data for this class

  Inertial4thOrderOp m_inertial4thOrderOp;
                                      ///< Operator for inertial (hyperbolic)
                                      ///< flux in compressible Navier-Stokes

  ViscousTensor4thOrderOp m_viscous4thOrderOp;
                                      ///< Operator for viscous (elliptic)
                                      ///< flux in compressible Navier-Stokes
};

#endif  /* ! defined _PATCHCNSOP_H_ */
