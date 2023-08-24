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
 * \file LES.H
 *
 * \brief Class for Large Eddy Simulation
 *
 *//*+*************************************************************************/

#ifndef _LES_H_
#define _LES_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "TurbModeling.H"

//--Forward declarations

/*******************************************************************************
 */
///  LES derived from turbulence modeling base class
/**
 *   This class is an implementation of Large Eddy Simulation
 *
 *//*+*************************************************************************/

class LES : public TurbModeling
{
  
/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  LES()
    { }

  /// Destructor
  virtual ~LES()
    { }

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  LES(const LES&);

//--Assignment not permitted

  LES &operator=(const LES&);

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Calculate momentum flux subgrid-scale correction
  virtual void sgsModelFlux(FluxBox&                a_sgsMomentumFxb,
                            FluxBox&                a_sgsEnergyFxb,
                            FluxBox&                a_faceAvgPlotFxb,
                            const FluxBox&          a_WfacePntFxb,
                            const FArrayBox&        a_WcellAvgFab,
                            const FluxBox&          a_NGradUfacePntFxb,
                            const FluxBox&          a_strainRateTensorFxb,
                            const FluxBox&          a_NGradTfacePntFxb,
                            const FluxBox&          a_facePntDeltaC,
                            const FluxBox&          a_faceCoord,
                            const FluxBox&          a_unitNormalsDirFxb,
                            const LevelGridMetrics& a_gridMetrics,
                            const DataIndex&        a_dataIndx,
                            const ProblemDomain&    a_domain,
                            const Box&              a_box) const;

  ///  Compute cell-avg sgs kinetic energy estimate with minimal input
  virtual void cellAvgSGSKineticEnergy(FArrayBox&              a_U,
                                       FArrayBox&              a_JU,
                                       const LevelGridMetrics& a_gridMetrics,
                                       const DataIndex&        a_dataIndx,
                                       const bool              a_fourthOrder,
                                       const Box&              a_box) const;

  /// Compute cell-avg sgs kinetic energy estimate on coarsened coarsest level
  virtual void crsLevCellAvgSGSKineticEnergy(
    LevelData<FArrayBox>&       a_SGSKE,
    const LevelData<FArrayBox>& a_JW,
    const LevelData<FArrayBox>& a_JGradW,
    const LevelData<FArrayBox>& a_cellPntJ,
    const LevelData<FArrayBox>& a_crsXFab,
    const LevelData<FArrayBox>& a_crsDeltaC,
    const ProblemDomain&        a_probDom) const;

  /// Compute sgs kinetic energy estimate for consToPrim correction
  virtual void consToPrimCorrectionLES(FArrayBox& a_WcellFab,
                                       const Box& a_box) const;

  /// Store SGS kinetic energy from current time to act as source term
  virtual void storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
                                     LevelData<FArrayBox>& a_JU,
                                     const bool            a_store) const;

  /// Set wall-model data from finer level
  virtual void setWallModelFromFiner(FArrayBox&              a_JU,
                                     const ProblemDomain&    a_blockDomain,
                                     const Box&              a_box,
                                     const LevelGridMetrics& a_gridMetrics,
                                     const DataIndex&        a_dataIndex) const;

  /// Compute the slip-velocity at no-slip walls for turb model
  virtual void turbModelSlipVel(FArrayBox&              a_bndryVelFab,
                                FArrayBox&              a_WcellAvgFab,
                                const FluxBox&          a_WfaceAvgFxb,
                                const LevelGridMetrics& a_gridMetrics,
                                const ProblemDomain&    a_domain,
                                const DataIndex&        a_dataIndx,
                                const FArrayBox&        a_unitNormalsDirFab,
                                const FArrayBox&        a_bndryNtJFab,
                                const Box&              a_boundaryFaceBox,
                                const Box&              a_boundaryFaceGhostBox,
                                const Box&              a_disjointBox,
                                const int               a_dir,
                                const Side::LoHiSide&   a_side) const;

  /// Compute source term for eta_0 in order to advance eta_0 in time
  virtual void updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                             FArrayBox&              a_WcellAvgFab,
                             const FArrayBox&        a_JUFab,
                             const FluxBox&          a_facePntVelGradFxb,
                             const FluxBox&          a_WfaceAvgFxb,
                             const FluxBox&          a_fluxFxb,
                             FluxBox&                a_stressFluxFxb,
                             const FArrayBox&        a_RHSfab,
                             const FluxBox&          a_unitNormals,
                             const ProblemDomain&    a_domain,
                             const LevelGridMetrics& a_gridMetrics,
                             const DataIndex&        a_dataIndx,
                             const Box&              a_disjointBox,
                             const Box&              a_box,
                             const Real              a_dt) const;

  // We should have functions that provide
  // (1) filtering capability
  // (2) general eddy-viscosity capability
  // (3) wall models
  // Where should corrections for BCs go?

};

/*******************************************************************************
 *
 * Class LES: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Calculate the momentum flux subgrid-scale correction
/** \param[out] a_sgsMomentumFxb
 *                      Subgrid-scale momentum fluxes on faces
 *  \param[out] a_sgsEnergyFxb
 *                      Subgrid-scale energy fluxes on faces
 *  \param[in]  a_WfacePntFxb
 *                      Face-centered primitive variables
 *  \param[in]  a_NGradUfacePntFxb
 *                      Face-centered physical-space velocity gradients
 *  \param[in]  a_strainRateTensorFxb
 *                      Face-centered, physical-space 
 *                      strain-rate tensor
 *  \param[in]  a_NGradTfacePntFxb
 *                      Face-centered, physical-space
 *                      temperature gradient
 *  \param[in]  a_facePntDeltaC
 *                      Face-centered cell-cutoff length for LES
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dataIndx
 *                      Data index for the grid information
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Cell-centered box with faces over which the
 *                      subgrid-scale flux must be computed
 *//*-----------------------------------------------------------------*/

inline void
LES::sgsModelFlux(FluxBox&                a_sgsMomentumFxb,
                  FluxBox&                a_sgsEnergyFxb,
                  FluxBox&                a_faceAvgPlotFxb,
                  const FluxBox&          a_WfacePntFxb,
                  const FArrayBox&        a_WcellAvgFab,
                  const FluxBox&          a_NGradUfacePntFxb,
                  const FluxBox&          a_strainRateTensorFxb,
                  const FluxBox&          a_NGradTfacePntFxb,
                  const FluxBox&          a_facePntDeltaC,
                  const FluxBox&          a_faceCoord,
                  const FluxBox&          a_unitNormalsDirFxb,
                  const LevelGridMetrics& a_gridMetrics,
                  const DataIndex&        a_dataIndx,
                  const ProblemDomain&    a_domain,
                  const Box&              a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Compute cell-avg sgs kinetic energy estimate with minimal input
/** \param[out] a_U     Cell-averaged conservative state
 *  \param[out] a_JU    Mapped, cell-averaged conservative state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for this level
 *  \param[out] a_dataIndx
 *                      Data index for box within this level
 *  \param[in]  a_fourthOrder
 *                      If true, SGS KE is computed to 4th-order
 *  \param[in]  a_box   Centered box over which subgrid-scale
 *                      kinetic energy is calculated
 *//*-----------------------------------------------------------------*/

inline void
LES::cellAvgSGSKineticEnergy(FArrayBox&              a_U,
                             FArrayBox&              a_JU,
                             const LevelGridMetrics& a_gridMetrics,
                             const DataIndex&        a_dataIndx,
                             const bool              a_fourthOrder,
                             const Box&              a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Compute cell-avg sgs kinetic energy estimate with minimal input
/** \param[out] a_SGSKE Cell-averaged SGS KE estimate
 *  \param[in]  a_JW    Mapped, cell-averaged conservative state
 *  \param[in]  a_JGradW
 *                      Mapped, cell-averaged velocity gradient
 *  \param[in]  a_cellPntJ
 *                      Cell-centered physical cell volume
 *  \param[in]  a_crsXFab
 *                      Cell-centered physical space coordinates
 *  \param[in]  a_crsDeltaC
 *                      Cell-centered grid-scale cutoff-length
 *  \param[out] a_probDom
 *                      Problem domain
 *//*-----------------------------------------------------------------*/

inline void
LES::crsLevCellAvgSGSKineticEnergy(
  LevelData<FArrayBox>&       a_SGSKE,
  const LevelData<FArrayBox>& a_JW,
  const LevelData<FArrayBox>& a_JGradW,
  const LevelData<FArrayBox>& a_cellPntJ,
  const LevelData<FArrayBox>& a_crsXFab,
  const LevelData<FArrayBox>& a_crsDeltaC,
  const ProblemDomain&        a_probDom) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Compute sgs kinetic energy correction estimate for consToPrim
/** \param[out] a_WcellFab
 *                      Cell primitive variables with pressure
 *                      corrected by sgs KE estimate
 *  \param[in]  a_box   Cell-centered box over which pressure is
 *                      corrected after consToPrim calculation
 *//*-----------------------------------------------------------------*/

inline void
LES::consToPrimCorrectionLES(FArrayBox& a_WcellFab,
                             const Box& a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Store SGS kinetic energy from current time to act as source term
/** \param[out] a_sgsKE
 *                     Cell-averaged, mapped SGS kinetic energy
 *  \param[in]  a_JU   Mapped conservative variables
 *  \param[in]  a_store
 *                     Store the SGS kinetic energy; otherwise, update
 *//*-----------------------------------------------------------------*/

inline void
LES::storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
                           LevelData<FArrayBox>& a_JU,
                           const bool            a_store) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Set wall-model data from finer level
/** \param[in]  a_JU    \<\JU\> state on this level
 *  \param[in]  a_blockDomain
 *                      Domain for this block on this level
 *  \param[in]  a_box   Current disjoint box
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the current level
 *  \param[in]  a_dataIndex
 *                      DataIndex for the current box
 *//*-----------------------------------------------------------------*/

inline void
LES::setWallModelFromFiner(FArrayBox&              a_JU,
                           const ProblemDomain&    a_blockDomain,
                           const Box&              a_box,
                           const LevelGridMetrics& a_gridMetrics,
                           const DataIndex&        a_dataIndex) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Applies the wall model to the faces and ghost cells
/** \param[out] a_bndryVelFab
 *                      Slip velocities on the boundary
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Domain of the current block
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_bndryNtJFab
 *                      Face-centered NtJ
 *  \param[in]  a_cellFaceAvgBdryBox
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_cellBoxWcellBdry
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_disjointBox
 *                      Disjoint box for grid metric functions
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *//*-----------------------------------------------------------------*/

inline void
LES::turbModelSlipVel(FArrayBox&              a_bndryVelFab,
                      FArrayBox&              a_WcellAvgFab,
                      const FluxBox&          a_WfaceAvgFxb,
                      const LevelGridMetrics& a_gridMetrics,
                      const ProblemDomain&    a_domain,
                      const DataIndex&        a_dataIndx,
                      const FArrayBox&        a_unitNormalsDirFab,
                      const FArrayBox&        a_bndryNtJFab,
                      const Box&              a_boundaryFaceBox,
                      const Box&              a_boundaryFaceGhostBox,
                      const Box&              a_disjointBox,
                      const int               a_dir,
                      const Side::LoHiSide&   a_side) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Compute source term for eta_0 in order to advance eta_0 in time
/** \param[out] a_turbSourceAvgFab
 *                      Cell-averaged turbulent source terms
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_facePntVelGradFxb
 *                      The face-centered velocity gradient
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_RHSfab
 *                      Cell-averaged conservative state update
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_domain
 *                      Domain of this block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics on this level
 *  \param[in]  a_dataIndx
 *                      Index of current box in DisjointBoxLayout
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_box   Cell-centered box over which eta_0 source-term
 *                      is computed
 *//*-----------------------------------------------------------------*/

inline void
LES::updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                   FArrayBox&              a_WcellAvgFab,
                   const FArrayBox&        a_JUFab,
                   const FluxBox&          a_facePntVelGradFxb,
                   const FluxBox&          a_WfaceAvgFxb,
                   const FluxBox&          a_fluxFxb,
                   FluxBox&                a_stressFluxFxb,
                   const FArrayBox&        a_RHSfab,
                   const FluxBox&          a_unitNormals,
                   const ProblemDomain&    a_domain,
                   const LevelGridMetrics& a_gridMetrics,
                   const DataIndex&        a_dataIndx,
                   const Box&              a_disjointBox,
                   const Box&              a_box,
                   const Real              a_dt) const
{
  return;
}

#endif  /* ! defined _LES_H_ */
