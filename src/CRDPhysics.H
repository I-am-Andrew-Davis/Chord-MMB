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
 * \file CRDPhysics.H
 *
 * \brief MOL physics abstract class for physics in Chord applications
 *
 *//*+*************************************************************************/

#ifndef _CRDPHYSICS_H_
#define _CRDPHYSICS_H_

//----- Chombo Library -----//

#include "CH_HDF5.H"
#include "RefCountedPtr.H"
#include "UsingNamespace.H"
#include "NamespaceVar.H"

//----- Internal -----//

#include "TurbModeling.H"
#include "LES.H"
#include "CRDparam.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class FArrayBox;
class Fluxbox;
class ProblemDomain;
using BlockDomain = ProblemDomain;
class LevelGridMetrics;
#include "NamespaceFooter.H"

struct BCinfo;


/*******************************************************************************
 */
///  MOL Physics for compressible Navier-Stokes on mapped grids
/**
 *   This class is a reimplementation of MOLMappedPhysics and isolates
 *   physics-related algorithms.  A main difference is that all functions are
 *   const.
 *
 *//*+*************************************************************************/

class CRDPhysics
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

  /// How to normalize species mass fractions
  enum SpeciesNormalizeType
  {
    NormalizeTypeNone,                ///< No normalization is applied
    NormalizeTypeRedistribute,        ///< Redistribute all fractions
    NormalizeTypeLumpInert            ///< Lump errors into an inert species
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  CRDPhysics()
    :
    m_turbModel(NULL)
    { }

  /// Destructor
  virtual ~CRDPhysics()
    { }


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CRDPhysics(const CRDPhysics&);

//--Assignment not permitted

  CRDPhysics &operator=(const CRDPhysics&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--Numbers

  /// Number of primitive variables
  virtual int numPrimitive() const = 0;

  /// Number of primitive variables required to calculate the flux
  virtual int numNativePrimitive() const = 0;

  /// Number of conserved variables (number of state variables)
  virtual int numConservative() const = 0;

  /// Number of flux variables
  virtual int numFluxes() const = 0;

  /// Number of output variables for plot files (not checkpoint files)
  virtual int numOutputVar() const = 0;

  /// Number of turbulent variables
  virtual int numTurbVar() const
    { if (m_turbModel) return m_turbModel->numTurbVar();
      else return 0;}

  /// Number of transport variables
  virtual int numTransportVar() const
    { return CRDparam::g_numTransportScalars; }

//--Indices

  /* Try to preserve the constexpr ones if possible.  If not, make virtual
   * and uncomment in all derived classes.
   */

  /// Component index within the primitive variables for density
  // virtual int densityIndex() const = 0;
  static constexpr int densityIndex()
    { return 0; }

  /// Interval within the primitive variables corresponding to the velocities
  // virtual Interval velocityInterval() const = 0;
  static constexpr Interval velocityInterval()
    { return Interval(1, SpaceDim); }

  /// Component index within the primitive variables for pressure
  // virtual int pressureIndex() const = 0;
  static constexpr int pressureIndex()
    { return 1 + SpaceDim; }

  /// Component index within the primitive variables for bulk modulus
  // virtual int bulkModulusIndex() const = 0;
  static constexpr int bulkModulusIndex()
    { return pressureIndex(); }

  /// Component index within the primitive variables for temperature
  virtual int temperatureIndex() const = 0;

  /// Interval within the flux variables corresponding to vector flux
  // virtual Interval vectorFluxInterval() const = 0;
  static constexpr Interval vectorFluxInterval()
    { return Interval(1, SpaceDim); }

  /// Component index within the flux variables for energy
  // virtual int energyFluxIndex() const = 0;
  static constexpr int energyFluxIndex()
    { return 1 + SpaceDim; }

  /// Interval of native primitive state
  virtual Interval nativePrimInterval() const = 0;

  /// Interval of extra primitive state
  virtual Interval extraPrimInterval() const = 0;

  /// Interval in primitive variables of the species concentration
  virtual Interval speciesPrimInterval() const = 0;

  /// Interval in conservative variables of the species concentration
  virtual Interval speciesConsInterval() const = 0;

  /// Interval of primitive turbulence variables
  virtual Interval turbPrimInterval() const = 0;

  /// Interval of conservative turbulence variables
  virtual Interval turbConsInterval() const = 0;

  /// Interval of primitive scalar transport variables
  virtual Interval transportPrimInterval() const = 0;

  /// Interval of conservative scalar transport variables
  virtual Interval transportConsInterval() const = 0;

  /// Component index of species named 'a_specName' in primitive variables
  virtual int speciesPrimIndex(const std::string& a_specName) const
    { return -1; }

  /// Component index of species named 'a_specName' in conservative variables
  virtual int speciesConsIndex(const std::string& a_specName) const
    { return -1; }

//--Algorithms (on a patch)

  /// Compute the maximum wave speed on mapped grids 
  virtual void getMaxWaveSpeed(
    const ProblemDomain&    a_problemDomain,
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_invDtFab,
    const FArrayBox&        a_cellAvgW,
    const FluxBox&          a_N,
    const FArrayBox&        a_J,
    const LevelGridMetrics& a_levelGridMetrics,
    const Real&             a_stabilityCons,
    const RealVect&         a_dxVect,
    Real&                   a_minConvDt,
    IntVect&                a_minConvDtCell) const = 0;

  /// Compute the maximum wave speed on mapped grids during evalRHS
  virtual void getMaxWaveSpeedEvalRHS(
    const ProblemDomain&    a_problemDomain,
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_invDtFab,
    const FluxBox&          a_WfacePntFxb,
    const FluxBox&          a_N,
    const FArrayBox&        a_J,
    const LevelGridMetrics& a_levelGridMetrics,
    const Real&             a_stabilityCons,
    const RealVect&         a_dxVect,
    Real&                   a_minConvDt,
    IntVect&                a_minConvDtCell) const = 0;

  /// Compute the time step on mapped grids for the elliptic component
  virtual void ellipticDt(
    const ProblemDomain&    a_problemDomain,
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_invDtFab,
    const FArrayBox&        a_cellAvgW,
    const Real&             a_lambda,
    const FluxBox&          a_NTJ,
    const LevelGridMetrics& a_levelGridMetrics,
    const RealVect&         a_dxVect,
    Real&                   a_minDiffDt,
    IntVect&                a_minDiffDtCell) const = 0;

  /// Compute the time step on mapped grids for the elliptic component
  virtual void ellipticDtEvalRHS(
    const ProblemDomain&    a_problemDomain,
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_invDtFab,
    const FArrayBox&        a_WcellPntFab,
    const FArrayBox&        a_muFab,
    const Real&             a_lambda,
    const FluxBox&          a_NTJ,
    const LevelGridMetrics& a_levelGridMetrics,
    const RealVect&         a_dxVect,
    Real&                   a_minDiffDt,
    IntVect&                a_minDiffDtCell) const = 0;

  /// Compute the speed of sound
  virtual void soundSpeed(FArrayBox&       a_speed,
                          const FArrayBox& a_W,
                          const Box&       a_box) const = 0;

  /// Compute a flux from primitive variable values on a face
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_whalf,
                       const int&       a_dir,
                       const Box&       a_box) const = 0;

  /// Compute the solution to the Riemann problem.
  virtual void riemann(FArrayBox&       a_WStar,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const int&       a_dir,
                       const Box&       a_box) const = 0;

  /// Solve Riemann problem for inflow or outflow
  virtual void riemannBC(FArrayBox&            a_WavgFace,
                         FArrayBox&            a_WLeft,
                         FArrayBox&            a_WRight,
                         const FArrayBox&      a_unitNormalBasisFab,
                         const int             a_dir,
                         const Side::LoHiSide& a_side,
                         const Box&            a_boundaryFaceBox) const = 0;

  /// Compute primitive variables from conserved variables.
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box,
                          const FArrayBox& a_WOld) const = 0;

  /// Compute conservative variables from primitive variables
  virtual void primToCons(FArrayBox&       a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box) const = 0;

  /// Set extra primitive state from native primitive state
  virtual void extraPrimitiveState(FArrayBox&       a_Wx,
                                   const int        a_compWxBeg,
                                   const FArrayBox& a_Wp,
                                   const Box&       a_box) const = 0;

  /// Solve for the extra primitive state in a_W
  virtual void extraPrimitiveState(FArrayBox& a_W,
                                   const Box& a_box) const = 0;

  /// Compute the temperature from primary primitive variables
  virtual void temperature(FArrayBox& a_W,
                           const Box& a_box) const = 0;

  /// Compute the pressure from the temperature and density
  virtual void pressure(FArrayBox& a_W,
                        const Box& a_box) const = 0;

  /// Compute the density from the temperature and pressure
  virtual void density(FArrayBox& a_W,
                       const Box& a_box) const = 0;

  /// Compute the artificial viscosity contribution to the flux
  virtual void artVisc(const Box&         a_box,
                       const BlockDomain& a_domain,
                       FluxBox&           a_NtF,
                       const FArrayBox&   a_U,
                       const FArrayBox&   a_WOld,
                       const FluxBox&     a_N,
                       const FArrayBox&   a_J,
                       const FluxBox&     a_unitNormals,
                       const RealVect&    a_dx,
                       LevelGridMetrics&  a_gridMetrics,
                       const Real         a_time,
                       const int          a_level) const = 0;

  /// Solves for the linearly related primitive variables
  virtual void intermediateConsToPrim(FArrayBox&       a_W,
                                      const FArrayBox& a_U,
                                      const Box&       a_box,
                                      const FArrayBox& a_WOld) const
    {
      consToPrim(a_W, a_U, a_box, a_WOld);
    }

  /// Calculate the intermediate temperature and the derivative
  virtual void intermediateTandDiffF(FArrayBox&       a_Tbar,
                                     FArrayBox&       a_dFdTbar,
                                     const FArrayBox& a_WcellPntFab,
                                     const FArrayBox& a_engCellAvgFab,
                                     const Box&       a_box) const
    { }

  /// Calculate thermal conductivity and dynamic viscosity
  virtual void calcCoeffKappaMu(const Box&       a_box,
                                FArrayBox&       a_muFab,
                                FArrayBox&       a_kappaFab,
                                const FArrayBox& a_WfacePntFab) const = 0;

  /// Solves for J_n and adds h_n*J_n to the energy
  virtual void speciesDiffHeatFlux(const Box&       a_box,
                                   FArrayBox&       a_JnfacePntFab,
                                   FArrayBox&       a_energyFab,
                                   const FArrayBox& a_muFab,
                                   const FArrayBox& a_kappaFab,
                                   const FArrayBox& a_WfacePntFab) const = 0;

  /// Normalize the species mass fractions (primitive state)
  //  The default implementation is a no-op
  virtual void normalizePrimSpecies(
    const SpeciesNormalizeType a_normalizeType,
    const bool                 a_bound,
    const bool                 a_sortcnj,
    const Box&                 a_box,
    FArrayBox&                 a_W) const
    { }

  /// Species correction
  virtual void speciesCorrection(
    const Box&           a_box,
    FArrayBox&           a_JU,
    const ProblemDomain& a_domain,
    const RealVect&      a_dx,
    const Real           a_tolNeg = 0.01) const = 0;

  /// Return the species gas constants
  virtual Real speciesGasConstant(const int a_speciesNum) const = 0;

  /// Return the species molar mass
  virtual Real speciesMolarMass(const int a_speciesNum) const = 0;

  /// Solve for the reaction rate source term
  virtual int addReactionSource(const Box&       a_box,
                                FArrayBox&       a_RCTcellPntFab,
                                FArrayBox&       a_invDtFab,
                                const FArrayBox& a_WcellPntFab,
                                const Real       a_time,
                                const int        a_level,
                                Real&            a_minChemDt,
                                IntVect&         a_minChemDtCell) const = 0;

  /// Solve for the reaction rate source term with diagonal flux corrections
  virtual int ARSwithDiagonalFluxCorrection(
    const Box&       a_box,
    FArrayBox&       a_RCTcellPntFab,
    FArrayBox&       a_invDtFab,
    const FArrayBox& a_WcellPntFab,
    const FArrayBox& a_UcellPntFab,
    const Real       a_time,
    const int        a_level,
    Real&            a_minChemDt,
    IntVect&         a_minChemDtCell) const = 0;

  /// Compute the reaction source term jacobian, dS/dW*dW/dU
  virtual void computeReactionJacobian(const Box&       a_box,
                                       FArrayBox&       a_rxnJacobianFab,
                                       const FArrayBox& a_WcellPntFab,
                                       const Real       a_dt) const = 0;

  /// Initialize the flow field for a given density, velocity, pressure
  /// (or temperature) and species concentration if applicable
  virtual void initialize(FArrayBox&              a_U,
                          const FArrayBox&        a_W,
                          const LevelGridMetrics& a_gridMetrics,
                          const FluxBox&          a_unitNormals,
                          const DataIndex&        a_didx,
                          const Box&              a_disjointBox,
                          const Box&              a_box) const = 0;

  /// Initialize wall-model state after restart
  virtual void initWallModelAfterRestart(
    FArrayBox&              a_JU,
    const FluxBox&          a_unitNormals,
    const LevelGridMetrics& a_gridMetrics,
    const DataIndex&        a_didx,
    const Box&              a_disjointBox) const = 0;

  /// Return the specific gas constant
  virtual Real Rgas(const Real*          a_cn = nullptr,
                    const std::ptrdiff_t a_cnStride = 1) const = 0;

  /// Return the ratio of specific heats
  virtual Real gamma(const Real           a_T = 1.0,
                     const Real*          a_cn = nullptr,
                     const std::ptrdiff_t a_cnStride = 1) const = 0;

  /// Calculate values of gamma from primitive variables
  virtual void calcGamma(const Box&       a_box,
                         FArrayBox&       a_gamma,
                         const FArrayBox& a_W) const = 0;

  /// Return the specific heat at constant pressure
  virtual Real cp(const Real           a_T = 1.0,
                  const Real*          a_cn = nullptr,
                  const std::ptrdiff_t a_cnStride = 1) const = 0;

  /// Return the sensible enthalpy
  virtual Real enthalpy(const Real           a_T,
                        const Real*          a_cn = nullptr,
                        const std::ptrdiff_t a_cnStride = 1) const = 0;

//--Misc

  /// Return the name of the primitive variable in the argument component
  virtual const char *const primStateName(const int a_iComp) const = 0;

  /// Return the name of the conservative variable in the argument component
  virtual const char *const consvStateName(const int a_iComp) const = 0;

  /// Return name of the physics described by the class
  virtual const char *const physicsName() const = 0;

  /// Write any information related to the physics to output
  /*  Terminate with empty line
   */
  virtual void writePhysicsInfo() const = 0;

#ifdef CH_USE_HDF5
  /// Set the output LevelData
  virtual void outputLevelData(
    LevelData<FArrayBox>&       a_outputLD,
    const LevelData<FArrayBox>& a_U,
    const LevelData<FArrayBox>& a_WOld,
    const LevelGridMetrics&     a_levelGridMetrics) const = 0;

  /// Expressions for VisIt
  virtual void expressions(HDF5HeaderData& a_holder) const = 0;
#endif

//--Turbulence modeling

  /// Define the turbulence model pointer
  void defineTurbModeling(const TurbModeling *const a_turbModel);

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

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  RefCountedPtr<const TurbModeling> m_turbModel;
  bool m_turbModelDefined;            ///< Check if turbulence model is defined

};


/*******************************************************************************
 *
 * Class CRDPhysics: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Assign the pointer for the turbulence model
/**
 *//*-----------------------------------------------------------------*/

inline void
CRDPhysics::defineTurbModeling(const TurbModeling *const a_turbModel)
{
  m_turbModel = RefCountedPtr<const TurbModeling>(a_turbModel);
  m_turbModelDefined = true;
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
CRDPhysics::calcTurbDiffFlux(const Box&              a_box,
                             FArrayBox&              a_NGradNufacePntFab,
                             const FArrayBox&        a_muFab,
                             const FArrayBox&        a_kappaFab,
                             const FArrayBox&        a_WfacePntFab,
                             const int               a_dir,
                             const int               a_wTurbComp,
                             const int               a_tComp) const
{
  CH_assert(m_turbModelDefined);
  m_turbModel->calcTurbDiffFlux(a_box,
                                a_NGradNufacePntFab,
                                a_muFab,
                                a_kappaFab,
                                a_WfacePntFab,
                                a_dir,
                                a_wTurbComp,
                                a_tComp);
}

/*--------------------------------------------------------------------*/
//  Calculates the final diffusive flux for a single turbulent state
/** \param[in]  a_box  Box of cells to solve over
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
CRDPhysics::calcTurbSourceTerms(
  const Box&              a_box,
  const Box&              a_disjointBox,
  FArrayBox&              a_turbSourcePntFab,
  const FArrayBox&        a_GradWcellPntFab,
  const FArrayBox&        a_WcellPntFab,
  const DataIndex&        a_dataIndx,
  const LevelGridMetrics& a_gridMetrics) const
{
  CH_assert(m_turbModelDefined);
  m_turbModel->calcTurbSourceTerms(a_box,
                                   a_disjointBox,
                                   a_turbSourcePntFab,
                                   a_GradWcellPntFab,
                                   a_WcellPntFab,
                                   a_dataIndx,
                                   a_gridMetrics);
}

/*--------------------------------------------------------------------*/
//  Solve for the turbulent thermal conductivity and dynamic viscosity
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
CRDPhysics::calcCoeffKappatMut(const Box&              a_box,
                               FArrayBox&              a_mutFab,
                               FArrayBox&              a_kappatFab,
                               const FArrayBox&        a_muFab,
                               const FArrayBox&        a_kappaFab,
                               const FArrayBox&        a_WfacePntFab,
                               const FArrayBox&        a_strainRateFab,
                               const LevelGridMetrics& a_gridMetrics,
                               const int               a_dir) const
{
  m_turbModel->calcCoeffKappatMut(a_box,
                                  a_mutFab,
                                  a_kappatFab,
                                  a_muFab,
                                  a_kappaFab,
                                  a_WfacePntFab,
                                  a_strainRateFab,
                                  a_gridMetrics,
                                  a_dir);
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
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

inline void
CRDPhysics::setTurbulentBC(
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
  m_turbModel->setTurbulentBC(a_boundaryFaceBox,
                              a_Wface,
                              a_Wcell,
                              a_unitNormalBasisFab,
                              a_dir,
                              a_side,
                              a_gridMetrics,
                              a_time,
                              a_level,
                              a_domT);
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
 *                      Boundary face values out to a_boundaryFaceGhostBox
 *  \param[out] a_WfaceAvgDirFab
 *                      Boundary face state out to a_boundaryFaceBox with
 *                      turbulent variables modified by wall model
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive with extrapolated ghost cells
 *  \param[out] a_WcellAvgFab
 *                      Cell-averaged ghost cells modified with wall model
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
CRDPhysics::applyWallModel(
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
  m_turbModel->applyWallModel(a_boundaryFaceBox,
                              a_boundaryFaceGhostBox,
                              a_disjointBox,
                              a_WfaceBdryFab,
                              a_WfaceAvgDirFab,
                              a_WcellAvgFab,
                              a_unitNormalBasisFab,
                              a_bcType,
                              a_bcBox,
                              a_gridMetrics,
                              a_dataIndx,
                              a_dir,
                              a_side,
                              a_time,
                              a_level);
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
CRDPhysics::sgsModelFlux(FluxBox&                a_sgsMomentumFxb,
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
  m_turbModel->sgsModelFlux(a_sgsMomentumFxb,
                            a_sgsEnergyFxb,
                            a_faceAvgPlotFxb,
                            a_WfacePntFxb,
                            a_WcellAvgFab,
                            a_NGradUfacePntFxb,
                            a_strainRateTensorFxb,
                            a_NGradTfacePntFxb,
                            a_facePntDeltaC,
                            a_faceCoord,
                            a_unitNormalsDirFxb,
                            a_gridMetrics,
                            a_dataIndx,
                            a_domain,
                            a_box);
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
CRDPhysics::cellAvgSGSKineticEnergy(FArrayBox&              a_U,
                                    FArrayBox&              a_JU,
                                    const LevelGridMetrics& a_gridMetrics,
                                    const DataIndex&        a_dataIndx,
                                    const bool              a_fourthOrder,
                                    const Box&              a_box) const
{
  m_turbModel->cellAvgSGSKineticEnergy(a_U,
                                       a_JU,
                                       a_gridMetrics,
                                       a_dataIndx,
                                       a_fourthOrder,
                                       a_box);
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
CRDPhysics::crsLevCellAvgSGSKineticEnergy(
  LevelData<FArrayBox>&       a_SGSKE,
  const LevelData<FArrayBox>& a_JW,
  const LevelData<FArrayBox>& a_JGradW,
  const LevelData<FArrayBox>& a_cellPntJ,
  const LevelData<FArrayBox>& a_crsXFab,
  const LevelData<FArrayBox>& a_crsDeltaC,
  const ProblemDomain&        a_probDom) const
{
  m_turbModel->crsLevCellAvgSGSKineticEnergy(a_SGSKE,
                                             a_JW,
                                             a_JGradW,
                                             a_cellPntJ,
                                             a_crsXFab,
                                             a_crsDeltaC,
                                             a_probDom);
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
CRDPhysics::consToPrimCorrectionLES(FArrayBox& a_WcellFab,
                                    const Box& a_box) const
{
  m_turbModel->consToPrimCorrectionLES(a_WcellFab,
                                       a_box);
}

/*--------------------------------------------------------------------*/
//  Store SGS kinetic energy from current time to act as source term
/** \param[out] a_sgsKE Cell-averaged, mapped SGS kinetic energy
 *  \param[in]  a_JU    Mapped conservative variables
 *  \param[in]  a_store Store the SGS kinetic energy; otherwise, update
 *//*-----------------------------------------------------------------*/

inline void
CRDPhysics::storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
                                  LevelData<FArrayBox>& a_JU,
                                  const bool            a_store) const
{
  m_turbModel->storeSGSKineticEnergy(a_sgsKE, a_JU, a_store);
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
CRDPhysics::setWallModelFromFiner(FArrayBox&              a_JU,
                                  const ProblemDomain&    a_blockDomain,
                                  const Box&              a_box,
                                  const LevelGridMetrics& a_gridMetrics,
                                  const DataIndex&        a_dataIndex) const
{
  m_turbModel->setWallModelFromFiner(a_JU,
                                     a_blockDomain,
                                     a_box,
                                     a_gridMetrics,
                                     a_dataIndex);
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
CRDPhysics::turbModelSlipVel(FArrayBox&              a_bndryVelFab,
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
  m_turbModel->turbModelSlipVel(a_bndryVelFab,
                                a_WcellAvgFab,
                                a_WfaceAvgFxb,
                                a_gridMetrics,
                                a_domain,
                                a_dataIndx,
                                a_unitNormalsDirFab,
                                a_bndryNtJFab,
                                a_boundaryFaceBox,
                                a_boundaryFaceGhostBox,
                                a_disjointBox,
                                a_dir,
                                a_side);
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
CRDPhysics::updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
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
  m_turbModel->updateEtaZero(a_turbSourceAvgFab,
                             a_WcellAvgFab,
                             a_JUFab,
                             a_facePntVelGradFxb,
                             a_WfaceAvgFxb,
                             a_fluxFxb,
                             a_stressFluxFxb,
                             a_RHSfab,
                             a_unitNormals,
                             a_domain,
                             a_gridMetrics,
                             a_dataIndx,
                             a_disjointBox,
                             a_box,
                             a_dt);
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
CRDPhysics::modeledNoSlipWallVelocityGradient(
  FluxBox&                a_NGradUfacePntFxb,
  const FluxBox&          a_unitNormalsFxb,
  const FArrayBox&        a_WcellAvgFab,
  const LevelGridMetrics& a_gridMetrics,
  const ProblemDomain&    a_domain,
  const Box&              a_box,
  const Box&              a_disjointBox) const
{
  m_turbModel->modeledNoSlipWallVelocityGradient(a_NGradUfacePntFxb,
                                                 a_unitNormalsFxb,
                                                 a_WcellAvgFab,
                                                 a_gridMetrics,
                                                 a_domain,
                                                 a_box,
                                                 a_disjointBox);
}

#endif  /* ! defined _CRDPHYSICS_H_ */
