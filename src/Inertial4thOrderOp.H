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
 * \file Inertial4thOrderOp.H
 *
 * \brief Patch operator for inertial flux in compressible Navier-Stokes
 *
 *//*+*************************************************************************/

#ifndef _INERTIAL4THORDEROP_H_
#define _INERTIAL4THORDEROP_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CNSIBC.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class ProblemDomain;
using BlockDomain = ProblemDomain;
class FArrayBox;
class FluxBox;
class LevelGridMetrics;
#include "NamespaceFooter.H"


/*******************************************************************************
 */
///  Patch operator for inertial flux in compressible Navier-Stokes
/**
 *//*+*************************************************************************/

class Inertial4thOrderOp
{


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  Inertial4thOrderOp(LevelGridMetrics& a_levelGridMetrics);

  // Use synthesized destructor

  // Copy, move, and assignment not permitted
  Inertial4thOrderOp(const Inertial4thOrderOp&) = default;
  Inertial4thOrderOp(Inertial4thOrderOp&&) = default;
  Inertial4thOrderOp& operator=(const Inertial4thOrderOp&) = default;
  Inertial4thOrderOp& operator=(Inertial4thOrderOp&&) = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Compute averages of W on faces
  void computeFaceAverage(const Box&         a_box,
                          const BlockDomain& a_domain,
                          FArrayBox&         a_flattening,
                          FluxBox&           a_WfaceAvgFxb,
                          const FArrayBox&   a_WcellAvgFab,
                          const FArrayBox&   a_WfromUavgFab,
                          const FluxBox&     a_unitNormalFxb,
                          const Real         a_time,
                          const bool         a_setFlattening,
                          const int          a_level) const;

  /// Compute the inertial (hyperbolic) flux
  void flux(const Box&                          a_disjointBox,
            const Box&                          a_box,
            const BlockDomain&                  a_domain,
            FluxBox&                            a_fluxFaceAvgFxb,
            FluxBox&                            a_fluxFromWavgFxb,
            FluxBox&                            a_WfaceAvgFxb,
            FluxBox&                            a_WfacePntFxb,
            FArrayBox&                          a_WcellAvgFab,
            const FArrayBox&                    a_WcellPntFab,
            FluxBox&                            a_faceAvgPlotFxb,
            FArrayBox&                          a_flattening,
            stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellFab,
            stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
            const FArrayBox&                    a_UcellAvgFab,
            const FluxBox&                      a_unitNormalFxb,
            const DataIndex&                    a_dataIndx,
            const Real                          a_time,
            const Real                          a_dt,
            const Real                          a_prevDt,
            const Real                          a_stageWeight,
            const bool                          a_setFlattening,
            const int                           a_level) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Compute fluxes on the faces from primitive state W on faces
  void computeAllFluxes(const Box&         a_disjointBox,
                        const BlockDomain& a_domain,
                        FluxBox&           a_fluxFaceFxb,
                        const FluxBox&     a_WfaceFxb,
                        const int          a_numTanGhost) const;

  /// Compute the full average and point primitive state on all the faces
  void computeFullFaceState(
    const Box&                          a_cellBoxWfacePnt,
    const Box&                          a_cellBoxWpFaceAvg,
    const Box&                          a_boxWcellBdry,
    const Box&                          a_disjointBox,
    const BlockDomain&                  a_domain,
    FArrayBox&                          a_WcellAvgFab,
    stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellFab,
    stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
    FluxBox&                            a_WfaceAvgFxb,
    FluxBox&                            a_WfacePntFxb,
    const FluxBox&                      a_unitNormalFxb,
    const DataIndex&                    a_dataIndx,
    const Real                          a_time,
    const Real                          a_prevDt,
    const int                           a_level) const;

  /// Compute boundary conditions and ghost cells
  void setBndryConditions(
    const Box&            a_boundaryFaceBox,
    const Box&            a_boundaryFaceGhostBox,
    const Box&            a_disjointBox,
    const BlockDomain&    a_domain,
    FArrayBox&            a_WfaceAvgBndryDirFab,
    FArrayBox&            a_WfaceAvgDirFab,
    FArrayBox&            a_WfacePntDirFab,
    FArrayBox&            a_WcellAvgFab,
    const FArrayBox&      a_bndryFaceAvgSlipVelDirFab,
    FArrayBox&            a_bndryCellFab,
    BCInfo&               a_bcInfo,
    Box&                  a_bcBox,
    const int             a_dir,
    const Side::LoHiSide& a_side,
    const FArrayBox&      a_unitNormalDirFab,
    const DataIndex&      a_dataIndx,
    const Real            a_time,
    const Real            a_prevDt,
    const int             a_level) const;

  /// Use boundary conditions to fill ghost cells
  void extrapolateDomainGhostCells(
    const Box&            a_boundaryFaceBox,
    const Box&            a_boundaryFaceGhostBox,
    const Box&            a_disjointBox,
    const BlockDomain&    a_domain,
    FArrayBox&            a_WcellAvgFab,
    const FArrayBox&      a_WfaceAvgBndryDirFab,
    const FArrayBox&      a_WfaceAvgDirFab,
    const BCInfo&         a_bcInfo,
    const Box&            a_bcBox,
    const int             a_dir,
    const Side::LoHiSide& a_side,
    const FArrayBox&      a_unitNormalDirFab,
    const DataIndex&      a_dataIndx,
    const Real            a_time,
    const int             a_level) const;

  /// Solves for the second-order W values
  void solveSecondOrderVals(FArrayBox&         a_WfaceAvgDirFab,
                            const FArrayBox&   a_WcellAvgFab,
                            const int          a_numWcomp,
                            const int          a_dir,
                            const Box&         a_box,
                            const BlockDomain& a_domain) const;

  /// Solves for the fifth-order W- and W+ values
  void solveFifthOrderVals(FArrayBox&         a_WMinus,
                           FArrayBox&         a_WPlus,
                           const FArrayBox&   a_WcellAvgFab,
                           const int          a_numSlopes,
                           const int          a_dir,
                           const Box&         a_box,
                           const BlockDomain& a_domain) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

//--References data defined in AMRLevelCNS.  These are all set during
//--construction and are fully defined by AMRLevelCNS

  LevelGridMetrics& m_levelGridMetrics;
                                      ///< Grid metrics for the level
};

#endif  /* ! defined _INERTIAL4THORDEROP_H_ */
