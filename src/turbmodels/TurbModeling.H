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
 * \file TurbModeling.H
 *
 * \brief Turbulence modeling base class
 *
 *//*+*************************************************************************/

#ifndef _TURBMODELING_H_
#define _TURBMODELING_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"
#include "ViscousTensor4thOrderOp.H"

//--Forward declarations

class Box;
class FArrayBox;
struct BCInfo;


/*******************************************************************************
 */
///  Base class for turbulence modeling
/**
 *   This class is a Base class for turbulence models
 *
 *//*+*************************************************************************/

class TurbModeling
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  TurbModeling()
    { }

  /// Destructor
  virtual ~TurbModeling()
    { }

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  TurbModeling(const TurbModeling&);

//--Assignment not permitted

  TurbModeling &operator=(const TurbModeling&);

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--Variables

  /// Return the number of turbulent variables
  virtual int numTurbVar() const;

//--Misc

  /// Return the name the variable names, index a_iComp only relates to
  /// turbulence variables, not any other variables
  virtual const char *const turbStateName(const int a_iComp) const;

  /// Output turbulence model info
  virtual void writeTurbModelInfo() const;

  /// Output turbulence model name
  virtual const char *const modelName() const;

//--Algorithms (on a patch)

  /// Calculate turbulent thermal conductivity and dynamic viscosity
  virtual void calcCoeffKappatMut(const Box&              a_box,
                                  FArrayBox&              a_mutFab,
                                  FArrayBox&              a_kappatFab,
                                  const FArrayBox&        a_muFab,
                                  const FArrayBox&        a_kappaFab,
                                  const FArrayBox&        a_WfacePntFab,
                                  const FArrayBox&        a_strainRateFab,
                                  const LevelGridMetrics& a_gridMetrics,
                                  const int               a_dir) const;

  /// Calculate the turbulent convective fluxes
  virtual void addTurbConvectiveFlux(FArrayBox&       a_flux,
                                     const FArrayBox& a_WFace,
                                     const int&       a_dir,
                                     const Box&       a_box) const;

  /// Calculate the turbulent diffusive fluxes
  virtual void calcTurbDiffFlux(const Box&       a_box,
                                FArrayBox&       a_NGradNufacePntFab,
                                const FArrayBox& a_muFab,
                                const FArrayBox& a_kappaFab,
                                const FArrayBox& a_WfacePntFab,
                                const int        a_dir,
                                const int        a_wTurbComp,
                                const int        a_tComp) const;

  /// Calculate the turbulent source terms
  virtual void calcTurbSourceTerms(
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_turbSourcePntFab,
    const FArrayBox&        a_GradWcellPntFab,
    const FArrayBox&        a_WcellPntFab,
    const DataIndex&        a_dataIndx,
    const LevelGridMetrics& a_gridMetrics) const;

  /// Initialize turbulent variables
  virtual void turbInitialize(FArrayBox&              a_U,
                              const FArrayBox&        a_W,
                              const LevelGridMetrics& a_gridMetrics,
                              const FluxBox&          a_unitNormals,
                              const DataIndex&        a_didx,
                              const Box&              a_disjointBox,
                              const Box&              a_box) const;

  /// Solve primitive turbulence variables from conservative variables
  virtual void turbConsToPrim(FArrayBox&       a_W,
                              const FArrayBox& a_U,
                              const Box&       a_box) const;

  /// Set turbulence values at boundaries
  virtual void setTurbulentBC(
    const Box&                    a_boundaryFaceBox,
    FArrayBox&                    a_Wface,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const LevelGridMetrics&       a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// Apply wall model
  virtual void applyWallModel(
    const Box&              a_boundaryFaceBox,
    const Box&              a_boundaryFaceGhostBox,
    const Box&              a_disjointBox,
    FArrayBox&              a_WfaceBdryFab,
    FArrayBox&              a_WfaceAvgDirFab,
    FArrayBox&              a_WcellAvgFab,
    const FArrayBox&        a_unitNormalBasisFab,
    const BCInfo&           a_bcType,
    const Box&              a_bcBox,
    const LevelGridMetrics& a_gridMetrics,
    const DataIndex&        a_dataIndx,
    const int               a_dir,
    const Side::LoHiSide&   a_side,
    const Real              a_time,
    const int               a_level) const;

  /// Calculate the subgrid fluxes
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

  /// Compute cell-avg sgs kinetic energy estimate with minimal input
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

  /// Modify velocity gradient to incorporate wall-shear-stress from model
  virtual void modeledNoSlipWallVelocityGradient(
    FluxBox&                a_NGradUfacePntFxb,
    const FluxBox&          a_unitNormalsFxb,
    const FArrayBox&        a_WcellAvgFab,
    const LevelGridMetrics& a_gridMetrics,
    const ProblemDomain&    a_domain,
    const Box&              a_box,
    const Box&              a_disjointBox) const;
};

/*******************************************************************************
 *
 * Class TurbModeling: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of turbulence variables
/** \return             Number of turbulence variables
 *//*-----------------------------------------------------------------*/

inline int
TurbModeling::numTurbVar() const
{
  return 0;
}

/*--------------------------------------------------------------------*/
//  Return the name of the turbulence variables, default is nothing
/** \return             Name of the a_iComp variable
 *//*-----------------------------------------------------------------*/

inline const char *const
TurbModeling::turbStateName(const int a_iComp) const
{
  return " ";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the turbulence model
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::writeTurbModelInfo() const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Write the name of the turbulence model, should not be called
/**
 *//*-----------------------------------------------------------------*/

inline const char *const
TurbModeling::modelName() const
{
  return " ";
}

/*--------------------------------------------------------------------*/
//  Returns nothing since no turbulence model is being used
/** \param[in]  a_box  Cell box.  Flux needs to be computed on
 *                     surrounding faces.
 *  \param[out] a_mutFab
 *                     FAB of turbulent dynamic viscosity
 *  \param[out] a_kappatFab
 *                     FAB of turbulent thermal conductivity
 *  \param[in]  a_muFab
 *                     FAB of dynamic viscosity
 *  \param[in]  a_kappaFab
 *                     FAB of thermal conductivity
 *  \param[in]  a_WfacePntFab
 *                     FAB containing the face-centered primitive variables
 *  \param[in]  a_strainRateFab
 *                     FAB of face-centered strain rate tensor
 *                     This has SpaceDim^2 components when being used and
 *                     only 1 component when serving as a placeholder
 *                     Retrieve the components using 
 *                     ViscousTensor4thOrderOp::tensorIdxRowOrder(row, col)
 *  \param[in]  a_gridMetrics
 *                     Level grid metrics
 *  \param[in]  a_dir  Current direction of solution
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::calcCoeffKappatMut(const Box&              a_box,
                                 FArrayBox&              a_mutFab,
                                 FArrayBox&              a_kappatFab,
                                 const FArrayBox&        a_muFab,
                                 const FArrayBox&        a_kappaFab,
                                 const FArrayBox&        a_WfacePntFab,
                                 const FArrayBox&        a_strainRateFab,
                                 const LevelGridMetrics& a_gridMetrics,
                                 const int               a_dir) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Compute the intertial flux from primitive variable values on a face
/** \param[out] a_flux  Flux on the faces
 *  \param[in]  a_WFace Primitive state on the face
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_box   Box on which to compute flux
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::addTurbConvectiveFlux(FArrayBox&       a_flux,
                                    const FArrayBox& a_WFace,
                                    const int&       a_dir,
                                    const Box&       a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Calculates the final diffusive flux for a single turbulent state
/** \param[in]  a_box  Box of faces to solve over
 *  \param[in]  a_NGradNufacePntFab
 *                     FAB of the face-centered mapped gradients 
 *                     of turbulent variables
 *  \param[out] a_NGradNufacePntFab
 *                     FAB of gradients modified accordingly
 *  \param[in]  a_muFab
 *                     FAB of the dynamic viscosity on each face
 *  \param[in]  a_kappaFab
 *                     FAB of the thermal conductivity values
 *  \param[in]  a_WfacePntFab
 *                     FAB containing face-centered primitive variables
 *  \param[in]  a_dir  Face direction
 *  \param[in]  a_wTurbComp
 *                     Primitive component of the current turbulent state
 *  \param[in]  a_tComp
 *                     Gradient component for use in a_NGradNufacePntFab
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::calcTurbDiffFlux(const Box&       a_box,
                               FArrayBox&       a_NGradNufacePntFab,
                               const FArrayBox& a_muFab,
                               const FArrayBox& a_kappaFab,
                               const FArrayBox& a_WfacePntFab,
                               const int        a_dir,
                               const int        a_wTurbComp,
                               const int        a_tComp) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Calculates the turbulent source terms
/** \param[in]  a_box  Box of cells to solve cell-centered source
 *  \param[in]  a_disjointBox
 *                     Disjoint box
 *  \param[out] a_turbSourcePntFab
 *                     FAB of cell-centered turbulent source terms
 *  \param[in]  a_GradWcellPntFab
 *                     FAB of cell-centered gradients
 *  \param[in]  a_WcellPntFab
 *                     FAB containing cell-centered primitive variables
 *  \param[in]  a_dataIndx
 *                     Current data index
 *  \param[in]  a_gridMetrics
 *                     Level grid metrics
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::calcTurbSourceTerms(const Box&              a_box,
                                  const Box&              a_disjointBox,
                                  FArrayBox&              a_turbSourcePntFab,
                                  const FArrayBox&        a_GradWcellPntFab,
                                  const FArrayBox&        a_WcellPntFab,
                                  const DataIndex&        a_dataIndx,
                                  const LevelGridMetrics& a_gridMetrics) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Initialize the turbulent variables in the flow field
/** \param[out] a_U    Cell-centered conservative variables
 *  \param[in]  a_W    Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                     Current disjointBox
 *  \param[in]  a_box  Box to initialize over
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::turbInitialize(FArrayBox&              a_U,
                             const FArrayBox&        a_W,
                             const LevelGridMetrics& a_gridMetrics,
                             const FluxBox&          a_unitNormals,
                             const DataIndex&        a_didx,
                             const Box&              a_disjointBox,
                             const Box&              a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Find primitive turbulent variables from conservative variables
/** \param[out] a_W    Primitive variables
 *  \param[in]  a_U    Conservative variables
 *  \param[in]  a_box  Box to solve over
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::turbConsToPrim(FArrayBox&       a_W,
                             const FArrayBox& a_U,
                             const Box&       a_box) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Set the primitive state for turbulent variables
/**
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state corrected for presence of wall
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::setTurbulentBC(
  const Box&                    a_boundaryFaceBox,
  FArrayBox&                    a_Wface,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const LevelGridMetrics&       a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Applies the wall model to the faces and ghost cells
/** \param[in]  a_boundaryFaceBox
 *                      Box on which to apply the wall model
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box on which to apply the wall model for ghost cells
 *  \param[in]  a_disjointBox
 *                      Disjoint box for grid metric functions
 *  \param[in]  a_WfaceBdryFab
 *                      Face-averaged primitive state including on boundaries
 *  \param[out] a_WfaceBdryFab
 *                      Face-averaged primitive state modified with wall model
 *  \param[in]  a_WfaceAvgDirFab
 *                      Face-averaged primitive state including on boundaries
 *  \param[out] a_WfaceAvgDirFab
 *                      Boundary face state out to a_boundaryFaceBox with
 *                      turbulent variables modified by wall model
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive with extrapolated ghost cells
 *  \param[out] a_WcellAvgFab
 *                      Cell-averaged primitive modified with wall model
 *  \param[in]  a_unitNormalFab
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_bcTypes
 *                      Vector of domain types
 *  \param[in]  a_bcBoxes
 *                      Vector of corresponding boxes for each domain type
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of the domain
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::applyWallModel(
  const Box&              a_boundaryFaceBox,
  const Box&              a_boundaryFaceGhostBox,
  const Box&              a_disjointBox,
  FArrayBox&              a_WfaceBdryFab,
  FArrayBox&              a_WfaceAvgDirFab,
  FArrayBox&              a_WcellAvgFab,
  const FArrayBox&        a_unitNormalBasisFab,
  const BCInfo&           a_bcType,
  const Box&              a_bcBox,
  const LevelGridMetrics& a_gridMetrics,
  const DataIndex&        a_dataIndx,
  const int               a_dir,
  const Side::LoHiSide&   a_side,
  const Real              a_time,
  const int               a_level) const
{
  return;
}

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
TurbModeling::sgsModelFlux(FluxBox&                a_sgsMomentumFxb,
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
TurbModeling::cellAvgSGSKineticEnergy(FArrayBox&              a_U,
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
TurbModeling::crsLevCellAvgSGSKineticEnergy(
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
TurbModeling::consToPrimCorrectionLES(FArrayBox& a_WcellFab,
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
TurbModeling::storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
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
TurbModeling::setWallModelFromFiner(FArrayBox&              a_JU,
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
TurbModeling::turbModelSlipVel(FArrayBox&              a_bndryVelFab,
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
TurbModeling::updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
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

/*--------------------------------------------------------------------*/
//  Modify velocity gradient to incorporate wall-shear-stress from model
/** \param[out] a_NGradUfacePntFxb
 *                      Physical-space velocity-gradient
 *  \param[in]  a_unitNormalsFxb     
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Box over which wall velocity-gradients
 *                      are modified based on wall-model
 *  \param[in]  a_disjointBox
 *                      Disjoint box of this dataIndex
 *//*-----------------------------------------------------------------*/

inline void
TurbModeling::modeledNoSlipWallVelocityGradient(
  FluxBox&                a_NGradUfacePntFxb,
  const FluxBox&          a_unitNormalsFxb,
  const FArrayBox&        a_WcellAvgFab,
  const LevelGridMetrics& a_gridMetrics,
  const ProblemDomain&    a_domain,
  const Box&              a_box,
  const Box&              a_disjointBox) const
{
  return;
}

#endif  /* ! defined _TURBMODELING_H_ */
