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
 * \file ViscousTensor4thOrderOp.H
 *
 * \brief Class encompassing viscous fourth-order tensor operations
 *
 *//*+*************************************************************************/

#ifndef _VISCOUSTENSOR4THORDEROP_H_
#define _VISCOUSTENSOR4THORDEROP_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

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
///  Patch operator for viscous flux in compressible Navier-Stokes
/**
 *//*+*************************************************************************/

class ViscousTensor4thOrderOp
{


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  ViscousTensor4thOrderOp(LevelGridMetrics& a_levelGridMetrics);

  // Use synthesized destructor

  // Copy, move, and assignment not permitted
  ViscousTensor4thOrderOp(const ViscousTensor4thOrderOp&) = default;
  ViscousTensor4thOrderOp(ViscousTensor4thOrderOp&&) = default;
  ViscousTensor4thOrderOp& operator=(const ViscousTensor4thOrderOp&) = default;
  ViscousTensor4thOrderOp& operator=(ViscousTensor4thOrderOp&&) = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Add in the viscous flux
  void flux(const Box&         a_box,
            const BlockDomain& a_domain,
            FluxBox&           a_fluxfaceAvgFxb,
            FluxBox&           a_fluxfacePntFxb,
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

  /// Returns component of the viscous stress tensor from row and column
  int static tensorIdxRowOrder(const int a_row, const int a_col)
    {
      return a_row*SpaceDim + a_col;
    }

  /// Returns component of the viscous stress tensor from row and column
  int static gradIdx(const int a_comp, const int a_dir)
    {
      return a_comp*SpaceDim + a_dir;
    }

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Compute the viscous stress tensor
  void computeVST(const Box&         a_box,
                  const BlockDomain& a_domain,
                  const Box&         a_domain5Pp2BCbox,
                  const Box&         a_domain5Pm1BCbox,
                  FluxBox&           a_VSTfacePntFxb,
                  FluxBox&           a_stressRateTensorFxb,
                  FluxBox&           a_NGradUfacePntFxb,
                  FArrayBox&         a_GradWcellAvgFab,
                  FArrayBox&         a_WcellAvgFab,
                  FluxBox&           a_WfaceAvgFxb,
                  const FluxBox&     a_unitNormalsFxb,
                  const DataIndex&   a_dataIndx,
                  const Real         a_time,
                  const int          a_level) const;

  /// Compute the cell-avgeraged velocity gradient in the interior
  void computeInteriorCellAvgVelGrad(FArrayBox&         a_GradUcellAvgFab,
                                     const FArrayBox&   a_UcellAvgFab,
                                     const FluxBox&     a_WfaceAvgFxb,
                                     const Box&         a_box,
                                     const BlockDomain& a_domain,
                                     const RealVect&    a_dxVect) const;

  /// Compute the face-avgeraged velocity gradient
  void computeInteriorFaceAvgVelGrad(FluxBox&           a_GradUfaceAvgFxb,
                                     const FArrayBox&   a_GradUcellAvgFab,
                                     const FArrayBox&   a_UcellAvgFab,
                                     const Box&         a_box,
                                     const BlockDomain& a_domain,
                                     const RealVect&    a_dxVect) const;

  /// Compute the face-averaged <N^T/J Grad Phi>
  void computeNGradPhi(const Box&         a_box,
                       const BlockDomain& a_domain,
                       const Box&         a_domain5Pp2BCbox,
                       const Box&         a_domain5Pm1BCbox,
                       FluxBox&           a_NGradPhifacePntFxb,
                       FArrayBox&         a_GradWcellAvgFab,
                       FArrayBox&         a_WcellAvgFab,
                       FluxBox&           a_WfaceAvgFxb,
                       const DataIndex&   a_dataIndx,
                       const int          a_indexPrim,
                       const int          a_vecStart,
                       const bool         a_rowMult,
                       const Real         a_time,
                       const int          a_level,
                       const bool         a_onlyCellGrad = false) const;

  /// Compute the cell-averaged temperature gradient in the interior
  void computeInteriorCellAvgPhiGrad(FArrayBox&         a_GradPhicellAvgFab,
                                     const FArrayBox&   a_PhicellAvgFab,
                                     const FluxBox&     a_PhifaceAvgFxb,
                                     const Box&         a_box,
                                     const BlockDomain& a_domain,
                                     const RealVect&    a_dxVect) const;

  /// Compute the face-averaged temperature gradient
  void computeInteriorFaceAvgPhiGrad(FluxBox&           a_GradPhifaceAvgFxb,
                                     const FArrayBox&   a_GradPhicellAvgFab,
                                     const FArrayBox&   a_PhicellAvgFab,
                                     const Box&         a_box,
                                     const BlockDomain& a_domain,
                                     const RealVect&    a_dxVect) const;

  /// Compute the 4th order product of a vector and matrix
  void computeFOMatVecProd(const Box&     a_box,
                           FluxBox&       a_outputFxb,
                           const FluxBox& a_MatfaceAvgFxb,
                           const FluxBox& a_VecfaceAvgFxb,
                           const int      a_vecStart,
                           const bool     a_rowMult) const;

  /// Compute a face-averaged matrix multiplication
  void computeMatrixProduct(FluxBox&       a_C,
                            const FluxBox& a_A,
                            const FluxBox& a_B,
                            const int      a_varStart) const;
  
  /// Retrieve the component from a flux dyad
  int retrieveFluxDyadComp(const int a_var,
                           const int a_space) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  enum : int
  {
    c_tanGhost = 2                    ///< Additional number of ghost cells 
                                      ///< to m_pBndryGhost to expand for
                                      ///< the cell-averaged velocity
                                      ///< gradients that are tangential 
                                      ///< to the periodic boundary
  };

//--References data defined in AMRLevelCNS.  These are all set during
//--construction and are fully defined by AMRLevelCNS

  LevelGridMetrics& m_levelGridMetrics;
                                      ///< Grid metrics for the level
};

#endif  /* ! defined _VISCOUSTENSOR4THORDEROP_H_ */
