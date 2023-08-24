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
 * \file ThermPhysics.H
 *
 * \brief MOL physics for thermally perfect compressible Navier-Stokes
 *        with species transport
 *
 *//*+*************************************************************************/

#ifndef _THERMPHYSICS_H_
#define _THERMPHYSICS_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CNSPhysics.H"
#include "CRDparam.H"
#include "CRDmsg.H"
#include "DataTemp.H"

//-- Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class FArrayBox;
#include "NamespaceFooter.H"


/*******************************************************************************
 */
///  MOL Physics for compressible Navier-Stokes on mapped grids
///  with species transport and reaction mechanisms
/**
 *   This class is a reimplementation of MOLMappedPhysics and isolates
 *   physics-related algorithms.  A main difference is that all fuctions are
 *   const.
 *
 *//*+*************************************************************************/

class ThermPhysics : public CNSPhysics
{

  friend struct C2PFunc;
  friend struct C2PDeriv;
  friend class FileThermParser;
  friend class FileThermCKParser;


/*==============================================================================
 * Types
 *============================================================================*/

  enum RiemannSolver
  {
    Exact,   // Exact solver type
    Approx,  // Approximate solver type
    Adaptive // Adaptive solver type
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  ThermPhysics();

  /// Destructor
  ~ThermPhysics() override;


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  ThermPhysics(const ThermPhysics&);

//--Assignment not permitted

  ThermPhysics &operator=(const ThermPhysics&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:


//--Numbers

  /// Number of primitive variables
  virtual int numPrimitive() const override;

  /// Number of primitive variables required to calculate the flux
  virtual int numNativePrimitive() const override;

  /// Number of conserved variables (number of state variables)
  virtual int numConservative() const override;

  /// Number of flux variables
  virtual int numFluxes() const override;

  /// Number of output variables for plot files (not checkpoint files)
  virtual int numOutputVar() const override;

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
  virtual int temperatureIndex() const override;

  /// Interval within the flux variables corresponding to vector flux
  // virtual Interval vectorFluxInterval() const
  //   { return Interval(1, SpaceDim); }

  /// Component index within the flux variables for energy
  // virtual int energyFluxIndex() const
  //   { return 1 + SpaceDim; }

  /// Interval of native primitive state
  virtual Interval nativePrimInterval() const override;

  /// Interval of extra primitive state
  virtual Interval extraPrimInterval() const override;

  /// Interval in primitive variables of the species concentration
  virtual Interval speciesPrimInterval() const override;

  /// Interval in conservative variables of the species concentration
  virtual Interval speciesConsInterval() const override;

  /// Interval of primitive turbulence variables
  virtual Interval turbPrimInterval() const override;

  /// Interval of conservative turbulence variables
  virtual Interval turbConsInterval() const override;

  /// Interval of primitive scalar transport variables
  virtual Interval transportPrimInterval() const override;

  /// Interval of conservative scalar transport variables
  virtual Interval transportConsInterval() const override;

  /// Zero-based index of a species
  int speciesZeroIndex(const std::string& a_specName) const;

  /// Component index of species named 'a_specName' in primitive variables
  virtual int speciesPrimIndex(const std::string& a_specName) const override;

  /// Component index of species named 'a_specName' in conservative variables
  virtual int speciesConsIndex(const std::string& a_specName) const override;

//--Algorithms (on a patch)

  /// Compute the speed of sound
  virtual void soundSpeed(FArrayBox&       a_speed,
                          const FArrayBox& a_W,
                          const Box&       a_box) const override;

  /// Compute a flux from primitive variable values on a face
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_whalf,
                       const int&       a_dir,
                       const Box&       a_box) const override;

  /// Compute the solution to the Riemann problem.
  virtual void riemann(FArrayBox&       a_WStar,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const int&       a_dir,
                       const Box&       a_box) const override;

  /// Solve Riemann problem for inflow or outflow
  virtual void riemannBC(FArrayBox&            a_WavgFace,
                         FArrayBox&            a_WLeft,
                         FArrayBox&            a_WRight,
                         const FArrayBox&      a_unitNormalBasisFab,
                         const int             a_dir,
                         const Side::LoHiSide& a_side,
                         const Box&            a_boundaryFaceBox) const override;

  /// Compute primitive variables from conserved variables.
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box,
                          const FArrayBox& a_WOld) const override;

  /// Compute conservative variables from primitive variables
  virtual void primToCons(FArrayBox&       a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box) const override;

  /// Set extra primitive state from native primitive state
  virtual void extraPrimitiveState(FArrayBox&       a_Wx,
                                   const int        a_compWxBeg,
                                   const FArrayBox& a_Wp,
                                   const Box&       a_box) const override;

  /// Solve for the extra primitive state in a_W
  virtual void extraPrimitiveState(FArrayBox& a_W,
                                   const Box& a_box) const override;

  /// Compute the temperature from primary primitive variables
  virtual void temperature(FArrayBox& a_W,
                           const Box& a_box) const override;

  /// Compute the pressure from the temperature and density
  virtual void pressure(FArrayBox& a_W,
                        const Box& a_box) const override;

  /// Compute the density from the temperature and pressure
  virtual void density(FArrayBox& a_W,
                       const Box& a_box) const override;

  /// Solves for the linearly related primitive variables
  virtual void intermediateConsToPrim(FArrayBox&       a_W,
                                      const FArrayBox& a_U,
                                      const Box&       a_box,
                                      const FArrayBox& a_WOld) const override;

  /// Calculate the intermediate temperature and the derivative
  virtual void intermediateTandDiffF(FArrayBox&       a_Tbar,
                                     FArrayBox&       a_dFdTbar,
                                     const FArrayBox& a_WcellPntFab,
                                     const FArrayBox& a_engCellAvgFab,
                                     const Box&       a_box) const override;

  /// Calculate thermal conductivity and dynamic viscosity and return max mux
  virtual void calcCoeffKappaMu(const Box&       a_box,
                                FArrayBox&       a_muFab,
                                FArrayBox&       a_kappaFab,
                                const FArrayBox& a_WfacePntFab) const override;

  /// Perform species diffusivity and species heat flux operations
  virtual void speciesDiffHeatFlux(const Box&       a_box,
                                   FArrayBox&       a_JnfacePntFab,
                                   FArrayBox&       a_energyFab,
                                   const FArrayBox& a_muFab,
                                   const FArrayBox& a_kappaFab,
                                   const FArrayBox& a_WfacePntFab) const override;

  /// Normalize the species mass fractions (primitive state)
  void normalizePrimSpecies(
    const SpeciesNormalizeType a_normalizeType,
    const bool                 a_bound,
    Real *const                a_cn,
    const std::ptrdiff_t       a_cnStride = 1,
    int*                       a_sortcnj = nullptr) const;

  /// Normalize the species mass fractions (primitive state)
  virtual void normalizePrimSpecies(
    const SpeciesNormalizeType a_normalizeType,
    const bool                 a_bound,
    const bool                 a_sortcnj,
    const Box&                 a_box,
    FArrayBox&                 a_W) const override;

  /// Perform species correction
  virtual void speciesCorrection(
    const Box&           a_box,
    FArrayBox&           a_JU,
    const ProblemDomain& a_domain,
    const RealVect&      a_dx,
    const Real           a_tolNeg = 0.01) const override;

  /// Returns the species gas constant
  virtual Real speciesGasConstant(const int a_speciesNum) const override;

  /// Returns the species gas constant
  virtual Real speciesMolarMass(const int a_speciesNum) const override;

  /// Solve for the reaction rate source term
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

  /// Compute the reaction source term jacobian, dS/dW*dW/dU
  virtual void computeReactionJacobian(const Box&       a_box,
                                       FArrayBox&       a_rxnJacobianFab,
                                       const FArrayBox& a_WcellPntFab,
                                       const Real       a_dt) const override;

  /// Initialize the flow field for a given density, velocity, pressure
  /// (or temperature) and species concentration if applicable
  virtual void initialize(FArrayBox&              a_U,
                          const FArrayBox&        a_W,
                          const LevelGridMetrics& a_gridMetrics,
                          const FluxBox&          a_unitNormals,
                          const DataIndex&        a_didx,
                          const Box&              a_disjointBox,
                          const Box&              a_box) const override;

  /// Return the specific gas constant
  virtual Real Rgas(const Real*          a_cn,
                    const std::ptrdiff_t a_cnStride = 1) const override;

  /// Return the ratio of specific heats
  virtual Real gamma(const Real           a_T = 1.0,
                     const Real*          a_cn = nullptr,
                     const std::ptrdiff_t a_cnStride = 1) const override;

  /// Calculate values of gamma from primitive variables for a FAB
  virtual void calcGamma(const Box&       a_box,
                         FArrayBox&       a_gamma,
                         const FArrayBox& a_W) const override;

  /// Solve for gamma using 'cp = h/T' and 'gamma = cp/(cp-Rgas)'
  virtual Real frozenGamma(const Real           a_T,
                           const Real           a_Rgas,
                           const Real*          a_cn,
                           const std::ptrdiff_t a_cnStride = 1) const;

  /// Return the specific heat at constant  pressure
  virtual Real cp(const Real           a_T = 1.0,
                  const Real*          a_cn = nullptr,
                  const std::ptrdiff_t a_cnStride = 1) const override;

  /// Return the sensible enthalpy of the mixture
  virtual Real enthalpy(const Real           a_T,
                        const Real*          a_cn = nullptr,
                        const std::ptrdiff_t a_cnStride = 1) const override;

  /// Solves for the sensible enthalpy with the heat of formation
  virtual Real enthalpy(const Real           a_T,
                        const Real*          a_cn,
                        Real*                a_hi,
                        const std::ptrdiff_t a_cnStride = 1,
                        const std::ptrdiff_t a_hiStride = 1) const;

//--Misc

  /// Return the name of the primitive variable in the argument component
  virtual const char *const primStateName(const int a_iComp) const override;

  /// Return the name of the conservative variable in the argument component
  virtual const char *const consvStateName(const int a_iComp) const override;

  /// Return name of the physics described by the class
  virtual const char *const physicsName() const;

  /// Write any information related to the physics to output
  virtual void writePhysicsInfo() const override;

#ifdef CH_USE_HDF5
  /// Set the output LevelData
  virtual void outputLevelData(
    LevelData<FArrayBox>&       a_outputLD,
    const LevelData<FArrayBox>& a_U,
    const LevelData<FArrayBox>& a_WOld,
    const LevelGridMetrics&     a_levelGridMetrics) const override;

  /// Expressions for VisIt
  virtual void expressions(HDF5HeaderData& a_holder) const;
#endif

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Non-linear solve of Riemann problem
  template <typename Func>
  void riemannSolution(const Func&      a_RStarSolve,
                       FArrayBox&       a_WStar,
                       const FArrayBox& a_shiftWLeft,
                       const FArrayBox& a_shiftWRight,
                       const int&       a_dir,
                       const Box&       a_box) const;

  /// Returns a vector of the coefficients and locations for reference tables
  void lookupRefs(const Real&        a_T,
                  std::vector<int>&  a_locs,
                  std::vector<Real>& a_P) const;

  /// Returns a vector of the coefficients and locations for reference tables
#ifdef USE_STACK
  void lookupRefs(const Real&                          a_T,
                  std::vector<int, StackAlloc<int>>&   a_locs,
                  std::vector<Real, StackAlloc<Real>>& a_P) const;
#endif

  /// Pressure dependent reaction functions
  Real TroeSolution(const Real& a_Pr,
                    const Real& a_T,
                    const int   a_curPR) const;

  Real SRISolution(const Real& a_Pr,
                   const Real& a_T,
                   const int   a_curPR) const;

  /// Find a temperature and pressure in cells where the Brent solver failed
  void findTempSpeciesFromFailedBrentSolve(
    Real&                                a_temperature,
    std::vector<Real, StackAlloc<Real>>& a_massFracs,
    std::vector<int, StackAlloc<int>>&   a_sortcnj,
    const Real&                          a_energy) const;

  /// Function to evaluate for the L-BFGS-B sovler
  Real LBFGSB_Function_Eval(const std::vector<Real>& a_currentSolution,
                            Real                     a_internalEnergy) const;

  /// Derivative of the function evaluated for the L-BFGS-B solver
  void LBFGSB_Gradient_Eval(std::vector<Real>&       a_gradVals,
                            const std::vector<Real>& a_currentSolution,
                            Real                     a_internalEnergy) const;

  /// Derivative of the function evalauted for the L-BFGS-B solver with respect to temperature
  void LBFGSB_Gradient_Eval_wrt_T(std::vector<Real>&                          a_gradVals,
                                  const std::vector<Real, StackAlloc<Real> >& a_currentMassFracs,
                                  Real                                        a_temperature,
                                  Real                                        a_internalEnergy,
                                  int                                         a_funcVal) const;

  /// Derivative of the function evalauted for the L-BFGS-B solver with respect to species mass fraction
  void LBFGSB_Gradient_Eval_wrt_cn(std::vector<Real>&                          a_gradVals,
                                   const std::vector<Real, StackAlloc<Real> >& a_currentMassFracs,
                                   Real                                        a_temperature,
                                   Real                                        a_internalEnergy,
                                   int                                         a_funcVal) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Real m_Rmol;                        ///< Universal gas constant in units
                                      ///< J/(mol K)
  std::vector<Real> m_molMass;        ///< Molar mass for each species
  std::vector<Real> m_Rn;             ///< Gas constant for each species
  RiemannSolver m_riemannSolver;      ///< Enum for Riemann solver type
  // File information
  std::string m_reactionFile;         ///< Reaction file name
  std::string m_thermoFile;           ///< Thermo coefficients file name
  std::string m_transFile;            ///< Transport coefficients file name
  int m_thermFileFormat;              ///< 0 - Thermo file format is NASA
                                      ///< 1 - Thermo file format is Chemkin
  Real m_lewisNum;                    ///< Lewis number
  Real m_schmidtNum;                  ///< Schmidt number
  // Table lookup data
  Real m_lookupLoT;                   ///< Low T value in lookup tables
  Real m_lookupHiT;                   ///< High T value in lookup tables
  Real m_lookupDelT;                  ///< Change in T in lookup tables
  Real m_midLookupT;                  ///< T value where coefficients change
  int m_lookupSize;                   ///< Number of values in lookup table
  std::vector<std::vector<Real> > m_lookupH;
                                      ///< Enthalpy lookup table
  std::vector<std::vector<Real> > m_lookupMu;
                                      ///< Dynamic viscosity lookup table
  std::vector<std::vector<Real> > m_lookupKappa;
                                      ///< Thermal conductivity lookup table
  int m_interpN;                      ///< Number of interp points to use
  std::vector<int> m_diffCL;          ///< Vector of the stencil used for
                                      ///< Lagrange polynomial interpolation
  std::vector<Real> m_interpDenoms;   ///< Denominators for the Lagrange interp
  std::vector<Real> m_H0;             ///< H(298.15) - H(0) - Delta H_F(298.15)
                                      ///< and Delta H_F(0) is added
  // Coefficients for when 200 < T < 1000
  std::vector<Real> m_hnA1L;          ///< Enthalpy A1 coeff
  std::vector<Real> m_hnA2L;          ///< Enthalpy A2 coeff
  std::vector<Real> m_hnA3L;          ///< Enthalpy A3 coeff
  std::vector<Real> m_hnA4L;          ///< Enthalpy A4 coeff
  std::vector<Real> m_hnA5L;          ///< Enthalpy A5 coeff
  std::vector<Real> m_hnA6L;          ///< Enthalpy A6 coeff
  std::vector<Real> m_hnA7L;          ///< Enthalpy A7 coeff
  std::vector<Real> m_hnB1L;          ///< Enthalpy B1 coeff
  std::vector<Real> m_hnB2L;          ///< Enthalpy B2 coeff
  // Coefficients for when 1000 < T < 6000
  std::vector<Real> m_hnA1H;          ///< Enthalpy A1 coeff
  std::vector<Real> m_hnA2H;          ///< Enthalpy A2 coeff
  std::vector<Real> m_hnA3H;          ///< Enthalpy A3 coeff
  std::vector<Real> m_hnA4H;          ///< Enthalpy A4 coeff
  std::vector<Real> m_hnA5H;          ///< Enthalpy A5 coeff
  std::vector<Real> m_hnA6H;          ///< Enthalpy A6 coeff
  std::vector<Real> m_hnA7H;          ///< Enthalpy A7 coeff
  std::vector<Real> m_hnB1H;          ///< Enthalpy B1 coeff
  std::vector<Real> m_hnB2H;          ///< Enthalpy B2 coeff
  std::vector<string> m_altSpecNames; ///< Alternate names for species
  std::unordered_map<std::string, int, CH_Hash::google_CityHash<std::string>>
  m_mapSpecNameToIndex;               ///< Map of species names to an index
                                      ///< (with first species at 0)
  int m_N2idx;                        ///< Species index (from zero) for N2, set
                                      ///< to -1 if no N2 in mechanism
  std::vector<int> m_intSeq;          ///< A sequence of integers counting from
                                      ///< 0 for all species.  This replaces an
                                      ///< indirection array when sorting of
                                      ///< species mass fractions is _not_ used.

  // Reaction data
  Real m_preAmodUnits;                ///< Factor to modify [X] to fit units
                                      ///< of A_r
  int m_numReactions;                 ///< Number of reactions
  int m_numTBReact;                   ///< Number of third body reactions and/or
                                      ///< pressure-dependent reactions
  int m_numPRReact;                   ///< Number of pressure dependent reaction
  int m_numREVReact;                  ///< Number of reactions where reverse
                                      ///< rate parameters are provided
  int m_numArbReact;                  ///< Number of arbitrary order reactions
  std::vector<int> m_numRects;        ///< List # of reactants for each reaction
  std::vector<int> m_numProds;        ///< List # of products for each reaction
  Real m_Tcutoff;                     ///< Minimum temperature for reactions
  std::vector<std::vector<Real> > m_lookupKfwd;
                                      ///< Forward reaction rate lookup table
  std::vector<std::vector<Real> > m_lookupKbkwd;
                                      ///< Backward reaction rate lookup table
  std::vector<std::vector<Real> > m_lookupKeq;
                                      ///< Equilibrium constant lookup table
  int m_maxSpecPerReact;              ///< Maximum number of species allowed
                                      ///< for a single side of a reaction
  std::vector<int> m_refNup;          ///< Ref vector for reactants
  std::vector<int> m_refNupp;         ///< Ref vector for products
  std::vector<Real> m_nup;            ///< Reactant coefficients
  std::vector<Real> m_nupp;           ///< Product coefficients
  std::vector<int> m_revReact;        ///< Denotes if the reaction is
                                      ///< reversible or irreversible
                                      ///< 1 - reversible, -1 - irreversible
                                      ///< 0 - REV react, |2| - arbitrary rate
                                      ///< Default is 1
  // Third body reactions (TB)
  std::vector<int> m_thirdBody;       ///< 1 if the reaction contains
                                      ///< third body reactions, -1 if the
                                      ///< reaction is PR but not TB, and 0
                                      ///< otherwise
  std::vector<Real> m_tbAlpha;        ///< Coefficients of third body reactions
  // Arbitrary order reactions
  std::vector<Real> m_ford;           ///< Vector of forward reaction orders
  std::vector<Real> m_rord;           ///< Vector of reverse reaction orders
  // Pressure-dependent reactions (PR)
  std::vector<int> m_presDep;         ///< Labels pressure dependent reactions
                                      ///< 0 - Not pressure dependent
                                      ///< 1 - Lindemann form pressure dependent
                                      ///< 2 - Troe form pressure dependent
                                      ///< 3 - SRI form
                                      ///< These values will be negative for
                                      ///< unimolecular and positive for
                                      ///< bimolecular
  Real m_Fc;                          ///< Lindemann fit number, default 1.
  std::vector<Real> m_betai;          ///< Temperature exponent
  std::vector<Real> m_EAR;            ///< Activation energy/R
  std::vector<Real> m_PREAR;          ///< Activation energy/R for PR
  std::vector<Real> m_PRpreAF;        ///< Forward reaction rate for PR
  std::vector<Real> m_PRbetai;        ///< Temperature exponent for PR
  // For Troe type pressure-dependent reactions
  // F_cent = (1 - alpha)*exp(-T/TS3) + alpha*exp(-T/TS) + exp(-TS2/T)
  std::vector<Real> m_TroeAlpha;      ///< Alpha value for Troe equation
  std::vector<Real> m_TroeTS3;        ///< TS3 value for Troe equation
  std::vector<Real> m_TroeTS2;        ///< TS2 value for Troe equation
  std::vector<Real> m_TroeTS;         ///< TS value for Troe equation
  // For SRI type pressure-dependent reactions
  std::vector<Real> m_SRIa;           ///< Coefficient a for SRI reaction
  std::vector<Real> m_SRIb;           ///< Coefficient b for SRI reaction
  std::vector<Real> m_SRIc;           ///< Coefficient c for SRI reaction
  std::vector<Real> m_SRId;           ///< Coefficient d for SRI reaction
  std::vector<Real> m_SRIe;           ///< Coefficient e for SRI reaction
  // FIXME: Currently, the following are NOT supported from chemkin input:
  // Landau-Teller reactions,
  // Optional rate fit expressions,
  // Other things to do:
  // Find correct kappa and mu coefficients for HO2, HCO, CH3, and CH2
private:

  void computeReactionJacobian_dW_dU(CHMatrix&                                  a_dW_dU,
                                     const Real&                                a_rhoInv,
                                     const Real&                                a_temperature,
                                     const Real&                                a_gamma,
                                     const Real&                                a_RInv,
                                     const std::vector<Real, StackAlloc<Real>>& a_hvals,
                                     const unsigned long&                       a_numSpecies,
                                     const int&                                 a_JacTindx) const;

  void computeReactionJacobian_dS_dW(CHMatrix&                                  a_dS_dW,
                                     const Real&                                a_T,
                                     const Real&                                a_rho,
                                     const std::vector<Real, StackAlloc<Real>>& a_spec,
                                     const std::vector<Real, StackAlloc<Real>>& a_mmass1,
                                     const std::vector<Real, StackAlloc<Real>>& a_mmass2,
                                     const std::vector<Real, StackAlloc<Real>>& a_Pvals,
                                     const std::vector<int, StackAlloc<int>>&   a_curL,
                                     const unsigned long&                       a_numSpecies,
                                     const unsigned long&                       a_numReactions,
                                     const int&                                 a_JacTindx) const;

};


/*******************************************************************************
 *
 * Class ThermPhysics: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of primitive variables
/** This may be greater than the number of conserved variables if
 *  derived/redundant quantities are also stored for convenience.
 *  For demonstration purposes, the variables numbers are listed for a
 *  2D case.  In the following,
 *    D is the number of dimensions
 *    Ns is the total number of species
 *    Nt is the number of turbulence variables
 *    Nr is the number of scalar transport variables
 *
 *  The conservative variables in 2D are:
 *  (0)             : density
 *  (1)             : x-momentum
 *  (D)             : y-momentum
 *  (D+1)           : energy-density
 *  (D+1+1)         : rho*c_1
 *  ...
 *  (D+1+Ns)        : rho*c_Ns
 *  (D+1+Ns+1)      : First turbulence variable
 *  ...
 *  (D+1+Ns+Nt)     : Last turbulence variable
 *  (D+1+Ns+Nt+1)   : First transport variable
 *  ...
 *  (D+1+Ns+Nt+Nr)  : Last transport variable
 *
 *  The primitive variables in 2D are:
 *  (0)             : density
 *  (1)             : x-vel
 *  (D)             : y-vel
 *  (D+1)           : pressure
 *  (D+1+1)         : c_1
 *  ...
 *  (D+1+Ns)        : c_Ns
 *  (D+1+Ns+1)      : First turbulence variable
 *  ...
 *  (D+1+Ns+Nt)     : Last turbulence variable
 *  (D+1+Ns+Nt+1)   : First transport variable
 *  ...
 *  (D+1+Ns+Nt+Nr)  : Last transport variable
 *  (D+1+Ns+Nt+Nr+1): temperature
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::numPrimitive() const
{
  return numNativePrimitive() + 1;
  //                            T
}

/*--------------------------------------------------------------------*/
//  Number of primitive variables required to calculate the flux
/** I.e., the minimal number of primitive variables
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::numNativePrimitive() const
{
  return 1 + SpaceDim + 1 + CRDparam::g_numSpecies + numTurbVar() +
  //     rho velocity   p   species_concentration    turbulence
         numTransportVar();
  //     scalar_transport
}

/*--------------------------------------------------------------------*/
//  Number of conserved variables (number of state variables)
/**
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::numConservative() const
{
  return 1 + SpaceDim + 1 + CRDparam::g_numSpecies + numTurbVar() +
  //     rho momentum   e   species_mass_fraction    turbuluence
         numTransportVar();
  //     scalar_transport
}

/*--------------------------------------------------------------------*/
//  Number of flux variables
/** This can be greater than the number of conserved variables if
 *  additional fluxes/face-centered quantities are computed.
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::numFluxes() const
{
  return numConservative();
}

/*--------------------------------------------------------------------*/
//  Number of output variables for plot files
/** We normally need two additional variables, T and c_p
 *  We can additionally include enthalpy, kappa, and mu
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::numOutputVar() const
{
  int c = numConservative() + 1;  // + T
  if (!CRDparam::g_plotDACFDCheck)
    {
      c += 2;                     // + p, c_p
    }
  
  if (CRDparam::g_plotExtraVars)
    {
      c += 1;                     // + enthalpy
      if (CRDparam::g_K < 0.)
        {
          c += 2;                 // + mu, kappa
        }
    }
  if (CRDparam::g_plotDACFDCheck)
    {
      c += 1;                     // + p
    }
  return c;
}

/*--------------------------------------------------------------------*/
//  Component index within the primitive variables for temperature
/**
 *//*-----------------------------------------------------------------*/

inline int
ThermPhysics::temperatureIndex() const
{
  return numPrimitive() - 1;
}

/*--------------------------------------------------------------------*/
//  Interval of native primitive state
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::nativePrimInterval() const
{
  return Interval(0, numNativePrimitive() - 1);
}

/*--------------------------------------------------------------------*/
//  Interval of extra primitive state
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::extraPrimInterval() const
{
  return Interval(temperatureIndex(), temperatureIndex());
}

/*--------------------------------------------------------------------*/
//  Interval of species concentration in primitive variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::speciesPrimInterval() const
{
  return Interval(pressureIndex() + 1, 
                  pressureIndex() + CRDparam::g_numSpecies);
}

/*--------------------------------------------------------------------*/
//  Interval of species concentration in conservative variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::speciesConsInterval() const
{
  return Interval(energyFluxIndex() + 1, 
                  energyFluxIndex() + CRDparam::g_numSpecies);
}

/*--------------------------------------------------------------------*/
//  Return interval of primitive turbulence variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::turbPrimInterval() const
{
  const int offset = speciesPrimInterval().end();
  return Interval(offset + 1, offset + numTurbVar());
}

/*--------------------------------------------------------------------*/
//  Return interval of conservative turbulence variables
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::turbConsInterval() const
{
  const int offset = speciesConsInterval().end();
  return Interval(offset + 1, offset + numTurbVar());
}

/*--------------------------------------------------------------------*/
//  Scalar transport variable interval in primitive variable array
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::transportPrimInterval() const
{
  const int offset = turbPrimInterval().end();
  return Interval(offset + 1, offset + numTransportVar());
}

/*--------------------------------------------------------------------*/
//  Scalar transport variable interval in conservative variable array
/**
 *//*-----------------------------------------------------------------*/

inline Interval
ThermPhysics::transportConsInterval() const
{
  const int offset = turbConsInterval().end();
  return Interval(offset + 1, offset + numTransportVar());
}

/*--------------------------------------------------------------------*/
// Return the given species gas constant
/**
 *//*-----------------------------------------------------------------*/

inline Real
ThermPhysics::speciesGasConstant(const int a_speciesNum) const
{
  CH_assert(a_speciesNum < CRDparam::g_numSpecies);
  return m_Rn[a_speciesNum];
}

/*--------------------------------------------------------------------*/
// Fill vector with molar masses
/**
 *//*-----------------------------------------------------------------*/

inline Real
ThermPhysics::speciesMolarMass(const int a_speciesNum) const
{
  CH_assert(a_speciesNum < CRDparam::g_numSpecies);
  return m_molMass[a_speciesNum];
}

#endif  /* ! defined _THERMPHYSICS_H_ */
