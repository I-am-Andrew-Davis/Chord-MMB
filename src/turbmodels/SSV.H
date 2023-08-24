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
 * \file SSV.H
 *
 * \brief Class for Strained-spiral vortex (stretched-vortex) SGS model
 *
 *//*+*************************************************************************/

#ifndef _SSV_H_
#define _SSV_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "LevelGridMetrics.H"
#include "CHMatrixOps.H"

//----- Internal -----//

#include "LES.H"

//--Forward declarations

/*******************************************************************************
 */
///  SSV derived from LES base class
/**
 *   This class is an implementation of the stretched-vortex SGS model
 *   The implementation is based on:
 *   (1)  A physical-space version of the stretched-vortex subgrid-stress model
 *        for large-eddy simulation (Voelkl, Pullin, Chan 2000)
 *   (2)  A vortex-based model for the subgrid flux of a passive scalar
 *        (Pullin 2000)
 *   (3)  Subgrid-scale modeling for large-eddy simulations of compressible
 *        turbulence (Kosovic, Pullin, Samtaney 2002)
 *   (4)  Large-eddy simulation and wall modelling of turbulent channel flow
 *        (Chung and Pullin 2009)
 *   (5)  Large-eddy simulation of separation and reattachment of a flat plate
 *        turbulent boundary layer (Cheng, Pullin, Samtaney 2015)
 *   (6)  Large-eddy simulation of flow over a cylinder with ReD from
 *        3.9x10^3 to 8.5x10^5: a skin-friction perspective (Cheng et al. 2017)
 *   (7)  Wall-modelled large-eddy simulation of turbulent flow past airfoils
 *        (Gao, Zhang, Cheng, Samtaney 2019)
 *
 *//*+*************************************************************************/

class SSV : public LES
{
  
/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  SSV();

  /// Destructor
  ~SSV();

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  SSV(const SSV&);

//--Assignment not permitted

  SSV &operator=(const SSV&);

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

  /// Calculate strain-rate tensor eigenvalues and eigenvectors
  void eigenDecomp(FArrayBox&       a_eigenVect,
                   const FArrayBox& a_strainRateTensorFab,
                   const Box&       a_box) const;

  /// Compute the orientation of sgs-vortices
  void sgsVortexVector(FArrayBox&              a_orientationFxb,
                       const FArrayBox&        a_NGradUfacePntFxb,
                       const FArrayBox&        a_eigenVect,
                       const LevelGridMetrics& a_gridMetrics,
                       const ProblemDomain&    a_domain,
                       const Box&              a_box) const;

  /// Compute kappa_C
  void computeKappaC(FArrayBox&       a_kappaC,
                     const FArrayBox& a_WfacePntFab,
                     const FArrayBox& a_eigenVect,
                     const FArrayBox& a_strainRateTensorFab,
                     const FArrayBox& a_deltaC,
                     const Box&       a_box) const;

  /// Compute the Kolmogorov prefactor for the SGS KE calculation
  void prefactorKolmogorov(FArrayBox&           a_KP,
                           const FArrayBox&     a_WfacePntFab,
                           const FArrayBox&     a_eigenVect,
                           const FArrayBox&     a_kappaC,
                           const FArrayBox&     a_XFab,
                           const FArrayBox&     a_strainRateTensor,
                           const FArrayBox&     a_deltaC,
                           const ProblemDomain& a_domain,
                           const Box&           a_box,
                           const int            a_localCoarsening=false) const;

  /// Compute cell-avg sgs kinetic energy estimate with minimal input
  void cellAvgSGSKineticEnergy(FArrayBox&              a_U,
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

  /// Compute sgs kinetic energy estimate
  void sgsKineticEnergy(FArrayBox&           a_sgsKE,
                        FArrayBox&           a_eigVect,
                        const FArrayBox&     a_NGradUfacePntFab,
                        const FArrayBox&     a_srtFab,
                        const FArrayBox&     a_WfacePntFab,
                        const FArrayBox&     a_deltaC,
                        const FArrayBox&     a_XFab,
                        const ProblemDomain& a_domain,
                        const Box&           a_box,
                        const int            a_localCoarsening = false) const;

  /// Compute sgs kinetic energy estimate for consToPrim correction
  virtual void consToPrimCorrectionLES(FArrayBox& a_WcellFab,
                                       const Box& a_box) const;

  /// Applies the wall model to the faces and ghost cells
  virtual void applyWallModel(const Box&              a_boundaryFaceBox,
                              const Box&              a_boundaryFaceGhostBox,
                              const Box&              a_disjointBox,
                              FArrayBox&              a_WfaceBdryFab,
                              FArrayBox&              a_WfaceAvgDirFab,
                              FArrayBox&              a_WcellAvgFab,
                              const FArrayBox&        a_unitNormalBasisFab,
                              const BCInfo&           a_bcTypes,
                              const Box&              a_bcBoxes,
                              const LevelGridMetrics& a_gridMetrics,
                              const DataIndex&        a_dataIndx,
                              const int               a_dir,
                              const Side::LoHiSide&   a_side,
                              const Real              a_time,
                              const int               a_level) const;

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

  /// Compute the slip-velocity for the turbulence model
  void computeSlipVel(FArrayBox&              a_bndryVelFab,
                      FArrayBox&              a_WcellAvgFab,
                      const FluxBox&          a_WfaceAvgFxb,
                      const LevelGridMetrics& a_gridMetrics,
                      const ProblemDomain&    a_domain,
                      const DataIndex&        a_dataIndx,
                      const FArrayBox&        a_unitNormalsDirFab,
                      const FArrayBox&        a_bndryNtJFab,
                      const Box&              a_cellFaceAvgBdryBox,
                      const Box&              a_cellBoxWcellBdry,
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

  /// Compute source term for eta_0 in order to advance eta_0 in time
  virtual void updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                             FArrayBox&              a_WcellAvgFab,
                             const FArrayBox&        a_facePntVelGradFab,
                             const FArrayBox&        a_WfaceAvgFab,
                             const FArrayBox&        a_RHSfab,
                             const FArrayBox&        a_unitNormals,
                             const ProblemDomain&    a_domain,
                             const LevelGridMetrics& a_gridMetrics,
                             const DataIndex&        a_dataIndx,
                             const Side::LoHiSide&   a_side,
                             const int               a_dir,
                             const Box&              a_disjointBox,
                             const Box&              a_box,
                             const Box&              a_boundaryFaceBox,
                             const Real              a_dt) const;

  /// Compute source term for eta_0 in order to advance eta_0 in time
  virtual void newUpdateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                                FArrayBox&              a_WcellAvgFab,
                                const FArrayBox&        a_JUFab,
                                const FArrayBox&        a_facePntVelGradFab,
                                const FArrayBox&        a_WfaceAvgFab,
                                const FluxBox&          a_WfaceAvgFxb,
                                FluxBox&                a_stressFluxFxb,
                                const FArrayBox&        a_RHSfab,
                                const FArrayBox&        a_unitNormals,
                                const FluxBox&          a_unitNormalsFxb,
                                const ProblemDomain&    a_domain,
                                const LevelGridMetrics& a_gridMetrics,
                                const DataIndex&        a_dataIndx,
                                const Side::LoHiSide&   a_side,
                                const int               a_dir,
                                const Box&              a_disjointBox,
                                const Box&              a_box,
                                const Box&              a_boundaryFaceBox,
                                const Real              a_dt) const;

  /// Compute source term for eta_0 in order to advance eta_0 in time
  virtual void updateEtaZeroLocal(FArrayBox&              a_turbSourceAvgFab,
                                  FArrayBox&              a_WcellAvgFab,
                                  const FArrayBox&        a_JUFab,
                                  const FArrayBox&        a_facePntVelGradFab,
                                  const FArrayBox&        a_WfaceAvgFab,
                                  const FluxBox&          a_WfaceAvgFxb,
                                  const FluxBox&          a_fluxFxb,
                                  FluxBox&                a_stressFluxFxb,
                                  const FArrayBox&        a_RHSfab,
                                  const FArrayBox&        a_unitNormals,
                                  const FluxBox&          a_unitNormalsFxb,
                                  const ProblemDomain&    a_domain,
                                  const LevelGridMetrics& a_gridMetrics,
                                  const DataIndex&        a_dataIndx,
                                  const Side::LoHiSide&   a_side,
                                  const int               a_dir,
                                  const Box&              a_disjointBox,
                                  const Box&              a_box,
                                  const Box&              a_boundaryFaceBox,
                                  const Real              a_dt) const;

  /// Calculate the turbulent convective fluxes
  virtual void addTurbConvectiveFlux(FArrayBox&       a_flux,
                                     const FArrayBox& a_WFace,
                                     const int&       a_dir,
                                     const Box&       a_box) const;

  /// Initialize turbulent variables
  virtual void turbInitialize(FArrayBox&              a_U,
                              const FArrayBox&        a_W,
                              const LevelGridMetrics& a_gridMetrics,
                              const FluxBox&          a_unitNormals,
                              const DataIndex&        a_didx,
                              const Box&              a_disjointBox,
                              const Box&              a_box) const;

  /// Find primitive turbulent variables from conservative variables
  virtual void turbConsToPrim(FArrayBox&       a_W,
                              const FArrayBox& a_U,
                              const Box&       a_box) const;

  /// Calculate the turbulent source terms
  virtual void calcTurbSourceTerms(
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_turbSourcePntFab,
    const FArrayBox&        a_GradWcellPntFab,
    const FArrayBox&        a_WcellPntFab,
    const DataIndex&        a_dataIndx,
    const LevelGridMetrics& a_gridMetrics) const;

  /// Store SGS kinetic energy from current time to act as source term
  void storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
                             LevelData<FArrayBox>& a_JU,
                             const bool            a_store) const;

  /// Set wall-model data from finer level
  virtual void setWallModelFromFiner(FArrayBox&              a_JU,
                                     const ProblemDomain&    a_blockDomain,
                                     const Box&              a_box,
                                     const LevelGridMetrics& a_gridMetrics,
                                     const DataIndex&        a_dataIndex) const;

  /// Estimate eta_0 based on near-wall state
  virtual void etaEst(FArrayBox&              a_etaFab,
                      const FArrayBox&        a_WcellAvgFab,
                      const LevelGridMetrics& a_gridMetrics,
                      const FArrayBox&        a_unitNormals,
                      const Box&              a_box,
                      const Box&              a_disjointBox,
                      const int               a_dir,
                      const Side::LoHiSide&   a_side) const;

  ///  Check if a box has a no-slip wall boundary
  virtual void hasNoSlipWall(int&                    a_hasNoSlipWall,
                             const LevelGridMetrics& a_gridMetrics,
                             const Box&              a_disjointBox,
                             const int               a_dir,
                             const Side::LoHiSide&   a_side) const;

  /// Modify velocity gradient to incorporate wall-shear-stress from model
  virtual void modeledNoSlipWallVelocityGradient(
    FluxBox&                a_NGradUfacePntFxb,
    const FluxBox&          a_unitNormalsFxb,
    const FArrayBox&        a_WcellAvgFab,
    const LevelGridMetrics& a_gridMetrics,
    const ProblemDomain&    a_domain,
    const Box&              a_box,
    const Box&              a_disjointBox) const;

  /// Returns component of the viscous stress tensor from row and column
  int static srtIdx(const int a_row, const int a_col)
    {
      return a_row*SpaceDim + a_col;
    }

/*--------------------------------------------------------------------*/
//  Non-linear function for estimating wall friction velocity
/**
 *  \param[in]  a_y     Wall-normal evaluation distance
 *  \param[in]  a_nu    Kinematic viscosity at wall
 *  \param[in]  a_u     van Driest velocity at evaluation point
 *//*-----------------------------------------------------------------*/

  struct UTauMuskerFunc
  {
    UTauMuskerFunc(const Real a_y,
                   const Real a_nu,
                   const Real a_u)
      :
      m_y(a_y),
      m_nu(a_nu),
      m_u(a_u)
      {}
    Real operator()(const Real& a_uTau) const
      {
        Real y_plus = m_y*a_uTau/m_nu;
        Real c_0 = 5.424*std::atan2((2.*y_plus - 8.15), 16.7);
        Real c_1 = std::pow((y_plus + 10.6), 9.6);
        Real c_2 = (y_plus*y_plus - 8.15*y_plus + 86);
        Real c_3 = std::log10(c_1/(c_2*c_2));
        return a_uTau*(c_0 + c_3 - 3.52) - m_u;
      }
    const Real m_y;
    const Real m_nu;
    const Real m_u;
  };

/*--------------------------------------------------------------------*/
//  Non-linear function for estimating wall friction velocity
/**
 *  \param[in]  a_y     Wall-normal evaluation distance
 *  \param[in]  a_nu    Kinematic viscosity at wall
 *  \param[in]  a_u     van Driest velocity at evaluation point
 *//*-----------------------------------------------------------------*/

  struct UTauSpaldingFunc
  {
    UTauSpaldingFunc(const Real a_y,
                     const Real a_nu,
                     const Real a_u)
      :
      m_y(a_y),
      m_nu(a_nu),
      m_u(a_u)
      {}
    Real operator()(const Real& a_uTau) const
      {
        Real y_plus = m_y*a_uTau/m_nu;
        Real u_plus = std::abs(m_u)/a_uTau;
        Real b = 0.4*u_plus;
        Real b2 = b*b;
        Real b3 = b*b2;
        Real b4 = b*b3;
        Real c_0 = std::exp(b);
        Real c_1 = u_plus + 0.1108*(c_0 - 1. - b - 0.5*b2 - (b3/6.) - (b4/24.));
        return y_plus - c_1;
      }
    const Real m_y;
    const Real m_nu;
    const Real m_u;
  };

/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  Real m_vortexProportion;            ///< Weighting for two direction models
  int  m_vortexModel;                 ///< SGS vortex model
                                      ///< 0 = normal model
                                      ///< 1 = \vec{e} \times \vec{\omega}
  int  m_sgsKineticEnergyApprox;      ///< Method for approximating SGS KE
                                      ///< 0 = Voelkl, Chung, Pullin's method
                                      ///< 1 = Mattner's method - faster
  int  m_useWallModel;                ///< The SV wall-model developed by
                                      ///< Chung and Pullin (2009) and others
                                      ///< 1 = use wall-model, 0 = don't use it
  int  m_useOnlyWallModel;            ///< Whether or not to use the interior
                                      ///< model with the wall model
                                      ///< 1 = no interior model, 0 = interior
  int  m_bndryVSTCorrection;          ///< 1 = correction to wall-model update
                                      ///< based on estimated wall-shear-stress
  Real m_virtualWallHeight;           ///< Height of virtual-wall as a fraction
                                      ///< of the first cell wall-normal height
  int  m_enforceWallShearStress;      ///< 1 = modify wall-face velocity-derivs
                                      ///< in ViscousOp to enforce
                                      ///< wall-model-predicted shear-stress
  int  m_suppressWallModel;           ///< 1 = don't use wall model even if
                                      ///< wall-model components are around
  int  m_sgsEnergyFlux;               ///< 1 = use SGS energy flux model
                                      ///< 0 = don't use SGS energy flux model
  int  m_constantK1;                  ///< 1 = set K1 to 0.41 in wall-model
  Real m_etaMaxBound;                 ///< Max bound on wall-modeled eta_0
};

/*******************************************************************************
 *
 * Class SSV: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of turbulence variables
/** 
 *//*-----------------------------------------------------------------*/

inline int
SSV::numTurbVar() const
{
  int turbVarNum = 1;
  if (m_useWallModel)
    {
      turbVarNum += SpaceDim;
      // eta_0, SpaceDim = num faces
    }
  return turbVarNum;
}

/*--------------------------------------------------------------------*/
//  Return the name of the turbulence variables, default is nothing
/** \return             Name of the a_iComp variable
 *//*-----------------------------------------------------------------*/

inline const char *const
SSV::turbStateName(const int a_iComp) const
{
  static constexpr const char* turbName[] =
    {
      "sgs kinetic energy",
      D_DECL("eta_0-x-face",  "eta_0-y-face",  "eta_0-z-face")
    };
  if (!m_useWallModel)
    {
      return "sgs kinetic energy";
    }
  else
    {
      return turbName[a_iComp];
    }
}

/*--------------------------------------------------------------------*/
//  Write any information related to the turbulence model
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

inline void
SSV::writeTurbModelInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "SGS vortex orientation: principal strain direction\n"
           << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Write the name of the turbulence model
/**
 *//*-----------------------------------------------------------------*/

inline const char *const
SSV::modelName() const
{
  return "LES: Stretched-vortex (strained-spiral vortex) SGS model";
}

#endif  /* ! defined _SSV_H_ */
