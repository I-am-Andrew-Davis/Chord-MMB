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
 * \file CNSPhysics.H
 *
 * \brief MOL physics for compressible Navier-Stokes equations
 *
 *//*+*************************************************************************/

#ifndef _CNSPHYSICS_H_
#define _CNSPHYSICS_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDPhysics.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class FArrayBox;
#include "NamespaceFooter.H"


/*******************************************************************************
 */
///  MOL Physics for compressible Navier-Stokes on mapped grids
/**
 *   This class is a reimplementation of MOLMappedPhysics and isolates
 *   physics-related algorithms.  A main difference is that all functions are
 *   const.
 *
 *//*+*************************************************************************/

class CNSPhysics : public CRDPhysics
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  CNSPhysics();

  /// Destructor
  ~CNSPhysics();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSPhysics(const CNSPhysics&);

//--Assignment not permitted

  CNSPhysics &operator=(const CNSPhysics&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:


//--Numbers

  /// Number of primitive variables
  virtual int numPrimitive() const;

  /// Number of primitive variables required to calculate the flux
  virtual int numNativePrimitive() const;

  /// Number of conserved variables (number of state variables)
  virtual int numConservative() const;

  /// Number of flux variables
  virtual int numFluxes() const;

  /// Number of output variables for plot files (not checkpoint files)
  virtual int numOutputVar() const;

//--Indices

  /// Component index within the primitive variables for density
  // virtual int densityIndex() const
  //   { return 0; }

  /// Interval within the primitive variables corresponding to the velocities
  // virtual Interval velocityInterval() const
  //   { return Interval(1, SpaceDim); }

  /// Component index within the primitive variables for pressure
  // virtual int pressureIndex() const
  //   { return 1 + SpaceDim; }

  /// Component index within the primitive variables for bulk modulus
  // virtual int bulkModulusIndex() const
  //   { return pressureIndex(); }

  /// Component index within the primitive variables for temperature
  virtual int temperatureIndex() const;

  /// Interval within the flux variables corresponding to vector flux
  // virtual Interval vectorFluxInterval() const
  //   { return Interval(1, SpaceDim); }

  /// Component index within the flux variables for energy
  // virtual int energyFluxIndex() const
  //   { return 1 + SpaceDim; }

  /// Interval of native primitive state
  virtual Interval nativePrimInterval() const;

  /// Interval of extra primitive state
  virtual Interval extraPrimInterval() const;

  /// Interval in primitive variables of the species concentration
  virtual Interval speciesPrimInterval() const;

  /// Interval in conservative variables of the species concentration
  virtual Interval speciesConsInterval() const;

  /// Interval of primitive turbulence variables
  virtual Interval turbPrimInterval() const;

  /// Interval of conservative turbulence variables
  virtual Interval turbConsInterval() const;

  /// Interval of primitive scalar transport variables
  virtual Interval transportPrimInterval() const;

  /// Interval of conservative scalar transport variables
  virtual Interval transportConsInterval() const;

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
    IntVect&                a_minConvDtCell) const override;

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
    IntVect&                a_minConvDtCell) const override;

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
    IntVect&                a_minDiffDtCell) const override;

  /// Compute the elliptic time step during evalRHS
  virtual void ellipticDtEvalRHS(
    const ProblemDomain&    a_problemDomain,
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_invDtFab,
    const FArrayBox&        a_WfacePntFab,
    const FArrayBox&        a_muFab,
    const Real&             a_lambda,
    const FluxBox&          a_NTJ,
    const LevelGridMetrics& a_levelGridMetrics,
    const RealVect&         a_dxVect,
    Real&                   a_minDiffDt,
    IntVect&                a_minDiffDtCell) const override;

  /// Compute the speed of sound
  virtual void soundSpeed(FArrayBox&       a_speed,
                          const FArrayBox& a_W,
                          const Box&       a_box) const;

  /// Compute a flux from primitive variable values on a face
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_whalf,
                       const int&       a_dir,
                       const Box&       a_box) const;

  /// Compute the solution to the Riemann problem.
  virtual void riemann(FArrayBox&       a_WStar,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const int&       a_dir,
                       const Box&       a_box) const;

  /// Solve Riemann problem for inflow or outflow
  virtual void riemannBC(FArrayBox&            a_WavgFace,
                         FArrayBox&            a_WLeft,
                         FArrayBox&            a_WRight,
                         const FArrayBox&      a_unitNormalBasisFab,
                         const int             a_dir,
                         const Side::LoHiSide& a_side,
                         const Box&            a_boundaryFaceBox) const;

  /// Compute primitive variables from conserved variables.
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box,
                          const FArrayBox& a_WOld) const;

  /// Compute conservative variables from primitive variables
  virtual void primToCons(FArrayBox&       a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box) const;

  /// Set extra primitive state from native primitive state
  virtual void extraPrimitiveState(FArrayBox&       a_Wx,
                                   const int        a_compWxBeg,
                                   const FArrayBox& a_Wp,
                                   const Box&       a_box) const;

  /// Solve for the extra primitive state in a_W
  virtual void extraPrimitiveState(FArrayBox& a_W,
                                   const Box& a_box) const;
  
  /// Compute the temperature from native primitive variables
  virtual void temperature(FArrayBox& a_W,
                           const Box& a_box) const;

  /// Compute the pressure from the temperature
  virtual void pressure(FArrayBox& a_W,
                        const Box& a_box) const;

  /// Compute the density from the temperature and pressure
  virtual void density(FArrayBox& a_W,
                       const Box& a_box) const;

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
                       const int          a_level) const;

  /// Calculate thermal conductivity and dynamic viscosity
  virtual void calcCoeffKappaMu(const Box&       a_box,
                                FArrayBox&       a_muFab,
                                FArrayBox&       a_kappaFab,
                                const FArrayBox& a_WfacePntFab) const;

  /// This should never be used for CNS code
  virtual void speciesDiffHeatFlux(const Box&       a_box,
                                   FArrayBox&       a_JnfacePntFab,
                                   FArrayBox&       a_energyFab,
                                   const FArrayBox& a_muFab,
                                   const FArrayBox& a_kappaFab,
                                   const FArrayBox& a_WfacePntFab) const;

  /// This should never be used for CNS code
  virtual void speciesCorrection(const Box&           a_box,
                                 FArrayBox&           a_JU,
                                 const ProblemDomain& a_domain,
                                 const RealVect&      a_dx,
                                 const Real           a_tolNeg = 0.01) const;

  /// Should not be used for CNS code
  virtual Real speciesGasConstant(const int a_speciesNum) const;

  /// Should not be used for CNS code
  virtual Real speciesMolarMass(const int a_speciesNum) const;

  /// Should not be used for CNS code, computes the reaction Jacobian
  virtual void computeReactionJacobian(const Box&       a_box,
                                       FArrayBox&       a_rxnJacobianFab,
                                       const FArrayBox& a_WcellPntFab,
                                       const Real       a_dt) const;

  /// Unused for CNS code
  virtual int addReactionSource(
    const Box&       a_box,
    FArrayBox&       a_RCTcellPntFab,
    FArrayBox&       a_invDtFab,
    const FArrayBox& a_WcellPntFab,
    const Real       a_time,
    const int        a_level,
    Real&            a_minChemDt,
    IntVect&         a_minChemDtCell) const override;

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
    IntVect&         a_minChemDtCell) const override;

  /// Initialize the flow field for a given density, velocity, pressure
  /// (or temperature) and species concentration if applicable
  virtual void initialize(FArrayBox&              a_U,
                          const FArrayBox&        a_W,
                          const LevelGridMetrics& a_gridMetrics,
                          const FluxBox&          a_unitNormals,
                          const DataIndex&        a_didx,
                          const Box&              a_disjointBox,
                          const Box&              a_box) const;

  /// Initialize wall-model state after restart
  virtual void initWallModelAfterRestart(
    FArrayBox&              a_JU,
    const FluxBox&          a_unitNormals,
    const LevelGridMetrics& a_gridMetrics,
    const DataIndex&        a_didx,
    const Box&              a_disjointBox) const;

  /// Return the specific gas constant
  virtual Real Rgas(const Real*          a_cn = nullptr,
                    const std::ptrdiff_t a_cnStride = 1) const override;

  /// Return the ratio of specific heats
  virtual Real gamma(const Real           a_T = 1.0,
                     const Real*          a_cn = nullptr,
                     const std::ptrdiff_t a_cnStride = 1) const override;

  /// Calculate values of gamma from primitive variables
  virtual void calcGamma(const Box&       a_box,
                         FArrayBox&       a_gamma,
                         const FArrayBox& a_W) const override;

  /// Return the specific heat at constant  pressure
  virtual Real cp(const Real           a_T = 1.0,
                  const Real*          a_cn = nullptr,
                  const std::ptrdiff_t a_cnStride = 1) const override;

  /// Return the sensible enthalpy
  virtual Real enthalpy(const Real           a_T,
                        const Real*          a_cn = nullptr,
                        const std::ptrdiff_t a_cnStride = 1) const override;

//--Misc

  /// Return the name of the primitive variable in the argument component
  virtual const char *const primStateName(const int a_iComp) const override;

  /// Return the name of the conservative variable in the argument component
  virtual const char *const consvStateName(const int a_iComp) const override;

  /// Return name of the physics described by the class
  virtual const char *const physicsName() const override;

  /// Write any information related to the physics to output
  virtual void writePhysicsInfo() const;

#ifdef CH_USE_HDF5
  /// Set the output LevelData
  virtual void outputLevelData(
    LevelData<FArrayBox>&       a_outputLD,
    const LevelData<FArrayBox>& a_U,
    const LevelData<FArrayBox>& a_WOld,
    const LevelGridMetrics&     a_levelGridMetrics) const;

  /// Expressions for VisIt
  virtual void expressions(HDF5HeaderData& a_holder) const;
#endif

private:

  /// Solve the Riemann problem
  void riemannSolution(FArrayBox&       a_WStar,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const int&       a_dir,
                       const Box&       a_box) const;

  /// Elliptic time step calculation
  void solveViscDt(const Box&       a_box,
                   const Box&       a_disjointBox,
                   FArrayBox&       a_invDtFab,
                   const FArrayBox& a_WFab,
                   const FArrayBox& a_muFab,
                   const Real&      a_lambda,
                   const FluxBox&   a_NTJ,
                   const RealVect&  a_dxVect,
                   Real&            a_minDiffDt,
                   IntVect&         a_minDiffDtCell) const;
};


/*******************************************************************************
 *
 * Class CNSPhysics: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of primitive variables
/** This may be greater than the number of conserved variables if
 *  derived/redundant quantities are also stored for convenience.
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::numPrimitive() const
{
  return numNativePrimitive() + 1;
  //                            t
}

/*--------------------------------------------------------------------*/
//  Number of primitive variables required to calculate the flux
/** I.e., the minimal number of primitive variables
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::numNativePrimitive() const
{
  return 1 + SpaceDim + 1 + numTurbVar() + numTransportVar();
  //     rho velocity   p   turbuluence    scalar_transport
}

/*--------------------------------------------------------------------*/
//  Number of conserved variables (number of state variables)
/**
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::numConservative() const
{
  return 1 + SpaceDim + 1 + numTurbVar() + numTransportVar();
  //     rho momentum   e   turbuluence    scalar_transport
}

/*--------------------------------------------------------------------*/
//  Number of flux variables
/** This can be greater than the number of conserved variables if
 *  additional fluxes/face-centered quantities are computed.
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::numFluxes() const
{
  return numConservative();
}

/*--------------------------------------------------------------------*/
//  Number of output variables for plot files
/** 
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::numOutputVar() const
{
  int c = numConservative();
  return c;
}

/*--------------------------------------------------------------------*/
//  Component index within the primitive variables for temperature
/**
 *//*-----------------------------------------------------------------*/

inline int
CNSPhysics::temperatureIndex() const
{
  return numPrimitive() - 1;
}

/*--------------------------------------------------------------------*/
//  Interval of native primitive state
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::nativePrimInterval() const
{
  return Interval(0, numNativePrimitive() - 1);
}

/*--------------------------------------------------------------------*/
//  Interval of extra primitive state
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::extraPrimInterval() const
{
  const int cT = temperatureIndex();
  return Interval(cT, cT);
}

/*--------------------------------------------------------------------*/
//  Should not be used at all
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::speciesPrimInterval() const
{
  return Interval(-1, -1);
}

/*--------------------------------------------------------------------*/
//  Should not be used at all
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::speciesConsInterval() const
{
  return Interval(-1, -1);
}

/*--------------------------------------------------------------------*/
//  Return interval of primitive turbulence variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::turbPrimInterval() const
{
  return Interval(pressureIndex() + 1, 
                  pressureIndex() + numTurbVar());
}

/*--------------------------------------------------------------------*/
//  Return interval of conservative turbulence variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::turbConsInterval() const
{
  return Interval(energyFluxIndex() + 1, 
                  energyFluxIndex() + numTurbVar());
}

/*--------------------------------------------------------------------*/
//  Scalar transport variable interval in primitive variable array
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::transportPrimInterval() const
{
  const int offset = turbPrimInterval().end();
  return Interval(offset + 1, offset + numTransportVar());
}

/*--------------------------------------------------------------------*/
//  Scalar transport variable interval in conservative variable array
/**
 *//*-----------------------------------------------------------------*/

inline Interval
CNSPhysics::transportConsInterval() const
{
  const int offset = turbConsInterval().end();
  return Interval(offset + 1, offset + numTransportVar());
}

/*--------------------------------------------------------------------*/
//  Return the specific gas constant
/** \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific gas constant for the mixture
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed (which for CNS is 0)
 *  </ul>
 *//*-----------------------------------------------------------------*/

inline Real
CNSPhysics::Rgas(const Real *const    a_cn,
                 const std::ptrdiff_t a_cnStride) const
{
  return CRDparam::g_R;
}

/*--------------------------------------------------------------------*/
//  Return the ratio of specific heats
/** \param[in]  a_T     Temperature of the mixture
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific heat ratio for the mixture
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

inline Real
CNSPhysics::gamma(const Real           a_T,
                  const Real*          a_cn,
                  const std::ptrdiff_t a_cnStride) const
{
  return CRDparam::g_gamma;
}

/*--------------------------------------------------------------------*/
//  Solves for specific heat at constant pressure
/** \param[in]  a_T     Temperature of the mixture
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific heat at constant pressure
 *//*-----------------------------------------------------------------*/

inline Real
CNSPhysics::cp(const Real           a_T,
               const Real*          a_cn,
               const std::ptrdiff_t a_cnStride) const
{
  const Real gamma = CRDparam::g_gamma;
  return CRDparam::g_R*gamma/(gamma - 1.);
}

/*--------------------------------------------------------------------*/
//  Solves for the enthalpy value
/** \tparam Vec         Container acting like a vector
 *  \param[in]  a_T     Temperature value
 *  \param[in]  a_cn    Unused for CNS
 *  \return             Enthalpy value
 *//*-----------------------------------------------------------------*/

inline Real
CNSPhysics::enthalpy(const Real           a_T,
                     const Real*          a_cn,
                     const std::ptrdiff_t a_cnStride) const
{
  const Real gamma = CRDparam::g_gamma;
  const Real R = CRDparam::g_R;
  const Real Cp = R*gamma/(gamma - 1.);
  return a_T*Cp;
}

#endif  /* ! defined _CNSPHYSICS_H_ */
