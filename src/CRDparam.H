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
 * \file CRDparam.H
 *
 * \brief Storage and access to parameters that are constant for a run
 *
 *//*+*************************************************************************/

#ifndef _CRDPARAM_H_
#define _CRDPARAM_H_

//----- Chombo Library -----//

#include "IntVect.H"
#include "RealVect.H"
#include "RefCountedPtr.H"
#include "Interval.H"
#include "SpectralUtil.H"
#include "MultiBlockCoordSys.H"
#include "UsingNamespace.H"

//----- Internal -----//

#ifdef CH_CTHR
#include "ThreadTeamArchitectChombo.H"
#endif

//--Forward declarations

class CRDPhysics;
class CNSIBC;
class DCFlattening;

// Using CRDparam instead of CRD::param to save colons
namespace CRDparam
{

/*--------------------------------------------------------------------*/
/// Types of simulation problems
/*--------------------------------------------------------------------*/

enum ProblemType : int
{
  ProblemUndefined = 0,                     // 0
  ProblemExternal,                          // 1
  ProblemGeneralized,                       // 2
  ProblemGeneralizedSingleBlock,            // 3
  ProblemEulerAdvectionCube,                // 4
  ProblemEulerAdvectionCylinder,            // 5
  ProblemEulerShockBox,                     // 6
  ProblemMachReflection,                    // 7
  ProblemShockTube,                         // 8
  ProblemNavierStokesTransientCouette,      // 9
  ProblemNavierStokesTransientFlatplate,    // 10
  ProblemNavierStokesTransientPoiseuille,   // 11
  ProblemCombustion,                        // 12
  ProblemLidDrivenCavity,                   // 13
  ProblemMixedCouette,                      // 14
  ProblemJetInlet,                          // 15
  ProblemSpecShock,                         // 16
  ProblemSpecShockBox,                      // 17
  ProblemFlame,                             // 18
  ProblemDetonation,                        // 19
  ProblemBurnerTest,                        // 20
  ProblemVortex,                            // 21
  ProblemReactionAdvection,                 // 22
  ProblemShockBubble,                       // 23
  ProblemRMI,                               // 24
  ProblemTemperatureDiffusion,              // 25
  ProblemSpecMachReflection,                // 26
  ProblemChannelFlow,                       // 27
  ProblemTurbulentCouette,                  // 28
  ProblemWedgeDetonation,                   // 29
  ProblemShuOsher,                          // 30
  ProblemCylinderFlow,                      // 31
  ProblemTaylorGreen,                       // 32
  ProblemIsotropicTurbulence,               // 33
  ProblemMixingLayer,                       // 34
  ProblemGaussPulse,                        // 35
  ProblemShear,                             // 36
  ProblemRiemannCube,                       // 37
  ProblemSpatiallyEvolvingShear,            // 38
  ProblemBluffBodyCombustion,               // 39
  ProblemMMS,                               // 40
  ProblemRecirculatingInletTFP,             // 41
  ProblemSphere,                            // 42
  ProblemObliqueWave,                       // 43
  ProblemChannelSeparation,                 // 44
  ProblemShockShock                         // 45
};

/*--------------------------------------------------------------------*/
/// Types of physics models
/*--------------------------------------------------------------------*/

enum PhysicsTypes
{
  PhysicsUndefined      = 0,
  PhysicsInertial       = (1<<0),  // Hyperbolic NS terms
  PhysicsViscous        = (1<<1),  // Elliptic NS terms
  PhysicsThermPerf      = (1<<2),  // Thermally perfect multi-species
  PhysicsSource         = (1<<3)
};

/*--------------------------------------------------------------------*/
/// Types of turbulence models
/*--------------------------------------------------------------------*/

enum TurbModelType
{
  TurbModelUndefined    = 0,
  TurbModelLES,  // LES model
  TurbModelSA,   // SA model
  TurbModelKO    // k-omega model
};

/*--------------------------------------------------------------------*/
/// Types of LES subgrid scale models
/*--------------------------------------------------------------------*/

enum SGSModelType
{
  SGSModelUndefined     = 0, 
  SGSModelSmagorinsky,      // Smagorinsky model
  SGSModelStretchedVortex   // Stretched-vortex model
};

/*--------------------------------------------------------------------*/
/// Types of explicit filters for LES
/*--------------------------------------------------------------------*/

enum ExplicitFilterType
{
  ExplicitFilterNone    = 0,
  ExplicitFilterSpectral   // spectral filter
};

/*--------------------------------------------------------------------*/
/// Types of forcing for Taylor-Green Vortex / HIT type problems
/*--------------------------------------------------------------------*/

enum TurbForcingType
{
  TurbForcingNone    = 0,
  TurbForcingLinear,       // Morales (2005) linear broadband forcing
  TurbForcingSpectral      // Petersen & Livescu (2010) isotropic, solenoidal,
                           // and spectrally-truncated forcing
};

/*--------------------------------------------------------------------*/
/// Types of coordinate systems
/*--------------------------------------------------------------------*/

enum CoordSysType
{
  CoordSysSingleBlockExternal,
  CoordSysSingleBlockCartesian,
  CoordSysSingleBlockWarped,
  CoordSysSingleBlockTwisted,
  CoordSysSingleBlockSchwarzChristoffelRamp,
  CoordSysSingleBlockLogSchwarzChristoffelRamp,
  CoordSysSingleBlockLogStretch,
  CoordSysSingleBlockAnnulus,
  CoordSysSingleBlockSkewed,
  CoordSysSingleBlockJoukowskiAirfoil,
  CoordSysMultiBlockExternal,
  CoordSysMultiBlockCartesian,
  CoordSysMultiBlockSphereShell,
  CoordSysExtrudedMultiBlock, //used for the ChannelSeparation case
  CoordSysMultiBlockSmoothRamp
};

/*--------------------------------------------------------------------*/
/// Time integration  method
/*--------------------------------------------------------------------*/

enum TimeIntegrationMethod : int
{
  RK2 = 0,
  RK4,
  ARK4
};

/*--------------------------------------------------------------------*/
/// Types of boundary conditions applied to a problem domain
/*--------------------------------------------------------------------*/

enum DomainBCType : int
{
  DomainBCTypeUndefined      = (0),     // For error checking
  DomainBCTypePeriodic       = (1<<0),
  DomainBCTypeDirichlet      = (1<<1),
  DomainBCTypeMixed          = (1<<2),  // Flag when more than 1 BC on boundary
  DomainBCTypeWall           = (1<<3),
  DomainBCTypeSlipWall       = (1<<4),
  DomainBCTypeMovingWall     = (1<<5),
  DomainBCTypeAdiabaticWall  = (1<<6),
  DomainBCTypeInflow         = (1<<7),  // Same as farfield but p extrapolated
  DomainBCTypeOutflow        = (1<<8),  // Same as farfield but !p extrapolated
  DomainBCTypeFarfield       = (1<<9),  // Characteristic BC
  DomainBCTypeExtrapolated   = (1<<10), // Fully extrapolated
  DomainBCTypeIsothermalWall = (1<<11), // Constant temperature wall
  DomainBCTypeJet            = (1<<12), //**FIXME: What is this?
  DomainBCTypeCNSCBCInflow   = (1<<13),
  DomainBCTypeCNSCBCOutflow  = (1<<14),
  DomainBCTypeRelaxedCBCIn   = (1<<15), // Inflow version of RelaxedCBC
  DomainBCTypeRelaxedCBCOut  = (1<<16), // Outflow version of RelaxedCBC
  DomainBCTypeRelaxedCBCFar  = (1<<17)  // Farfield version of RelaxedCBC
};
inline const char*
domainBCname(const int a_BCtype)
{
  switch (a_BCtype)
    {
    case DomainBCTypeUndefined:
      return "Undefined";
    case DomainBCTypePeriodic:
      return "Periodic";
    case DomainBCTypeDirichlet:
      return "Dirichlet";
    case DomainBCTypeMixed:
      return "Mixed";
    case DomainBCTypeWall:
      return "Wall";
    case DomainBCTypeSlipWall:
      return "SlipWall";
    case DomainBCTypeMovingWall:
      return "MovingWall";
    case DomainBCTypeAdiabaticWall:
      return "AdiabaticWall";
    case DomainBCTypeInflow:
      return "Inflow";
    case DomainBCTypeOutflow:
      return "Outflow";
    case DomainBCTypeFarfield:
      return "Farfield";
    case DomainBCTypeExtrapolated:
      return "Extrapolated";
    case DomainBCTypeIsothermalWall:
      return "IsothermalWall";
    case DomainBCTypeJet:
      return "Jet"; //???
    case DomainBCTypeCNSCBCInflow:
      return "CBCInflow";
    case DomainBCTypeCNSCBCOutflow:
      return "CBCOutflow";
    case DomainBCTypeRelaxedCBCIn:
      return "RelaxedCBCIn";
    case DomainBCTypeRelaxedCBCOut:
      return "RelaxedCBCOut";
    case DomainBCTypeRelaxedCBCFar:
      return "RelaxedCBCFar";
    default:
      return "Composite";
    };
}
const int DomainBCTypeAllFixedWall =
  (int)DomainBCTypeWall + (int)DomainBCTypeSlipWall;
// Any no-slip walls
const int DomainBCTypeSWall =
  (int)DomainBCTypeWall +
  (int)DomainBCTypeMovingWall +
  (int)DomainBCTypeAdiabaticWall +
  (int)DomainBCTypeIsothermalWall;
// All wall types
const int DomainBCTypeAllWall =
  (int)DomainBCTypeWall +
  (int)DomainBCTypeSlipWall +
  (int)DomainBCTypeMovingWall +
  (int)DomainBCTypeAdiabaticWall +
  (int)DomainBCTypeIsothermalWall;
// All characteristic boundaries
const int DomainBCTypeCNSCBC =
  (int)DomainBCTypeCNSCBCInflow + (int)DomainBCTypeCNSCBCOutflow;
// All inlet boundaries
const int DomainBCTypeInlet = (int)DomainBCTypeCNSCBCInflow
  + (int)DomainBCTypeInflow + (int)DomainBCTypeDirichlet
  + (int)DomainBCTypeJet;
// All inlet or outlet boundaries
const int DomainBCTypeInOut = (int)DomainBCTypeCNSCBCInflow
  + (int)DomainBCTypeCNSCBCOutflow + (int)DomainBCTypeDirichlet
  + (int)DomainBCTypeInflow + (int)DomainBCTypeOutflow
  + (int)DomainBCTypeJet + (int)DomainBCTypeFarfield
  + (int)DomainBCTypeExtrapolated;
// All relaxed characteristic boundaries
const int DomainBCTypeRelaxedCBC = (int)DomainBCTypeRelaxedCBCIn
  + (int)DomainBCTypeRelaxedCBCOut + (int)DomainBCTypeRelaxedCBCFar;

/*--------------------------------------------------------------------*/
/// Queries on the number of ghosts
/*--------------------------------------------------------------------*/

enum NumGhostQuery  // These are INDICES, not #ghosts.  The ones you should
                    // define in definePhysics() are noted.  That is where you
                    // specify the ghosts.
{
  NumGhostWcellExtTngBdry  = 0,  // define
  NumGhostInertialFfaceAvg = 1,  // define
  NumGhostInertialFfacePnt = 2,
  NumGhostWfaceAvgTan      = 2,
  NumGhostWfaceAvgTngBdry  = 2,
  NumGhostWfacePntTan      = 2,
  NumGhostWfacePntTngBdry  = 2,
  NumGhostWpFacePntTan     = 2,  // define
  NumGhostWfaceAvgNrm      = 3,
  NumGhostWfacePntNrm      = 3,
  NumGhostWpFacePntNrm     = 3,
  NumGhostWpFaceAvgNrm     = 3,  // define
  NumGhostFfromWpFaceAvg   = 4,
  NumGhostWpFaceAvgTan     = 4,  // define
  NumGhostWcellAvg         = 5,  // define
  NumGhostWcellPnt         = 6,
  NumGhostUcellPnt         = 6,  // define
  NumGhostWfromUcellAvg    = 7,
  NumGhostUcellAvg         = 7,  // define
};
const int NumGhostQuerySize = 8;

/*--------------------------------------------------------------------*/
/// Parameter groups to print to hdf5
/*--------------------------------------------------------------------*/

enum PlotComps //**NOTE: Most of these are not implemented yet
{
  PlotUndefined         = 0,
  PlotPrimitive         = (1<<0),  // Primitive state
  PlotConservative      = (1<<1),  // Conservative state
  PlotFluxes            = (1<<2),  // Flux (separate inertial, visc, etc.)
  PlotTurbulentComps    = (1<<3),  // Turbulence-related (<u_i*u_j>, etc.)
  PlotRHS               = (1<<4),  // Right-hand side, i.e. flux divergence
  PlotNtJ               = (1<<5),  // Mapping transformation (if needed)
  PlotFaceNormals       = (1<<6),  // Face normal transform (if needed)
  PlotVelGrad           = (1<<7),  // Velocity gradient tensor
  PlotTGrad             = (1<<8),  // Temperature gradient vector
  PlotSpeciesGrad       = (1<<9),  // Species gradients
  PlotSources           = (1<<10), // Source terms
};


/*******************************************************************************
 */
/// Global constants are stored here and accessed by const references.
/**
 *  The constants are modifiable within the class but non-modifiable when
 *  accessed by reference
 *
 *//*+*************************************************************************/

class CRDparamVar
{


/*==============================================================================
 * Private types
 *============================================================================*/

enum DefinedComponent : unsigned
{
  DefinedComponentNone     = 0,
  DefinedComponentGeometry = (1<<0)
};


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor sets defaults
  CRDparamVar();

  /// Destructor
  ~CRDparamVar();


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// Define initializing verbosity
  void defineInfo(const int  a_verbosity,
                  const bool a_verboseDt);

  /// Define initializing plot and file information
  void defineFiles(const std::string& a_plotPrefix,
                   const int          a_plotNumGhost,
                   const bool         a_plotExtraVars,
                   const bool         a_plotMappedDerivs,
                   const bool         a_plotWallDist,
                   const bool         a_plotJ,
                   const bool         a_plotJU,
                   const bool         a_plotDACFDCheck,
                   const bool         a_plotFlattening,
                   const bool         a_plotError,
                   const bool         a_plotTimeAvgTurb,
                   const bool         a_plotLoFaceCoordinates,
                   const bool         a_plotHiFaceCoordinates,
                   const bool         a_restartAddWallModel,
                   const int          a_plotLoFaceAvgTimeAvgComps,
                   const int          a_plotLoFaceAvgComps);

  /// Define initializing simulation type
  void defineSimulation(const ProblemType&      a_problemType,
                        const std::vector<int>& a_periodicity);

  /// Define coordinate system
  void defineCoordSys(const RefCountedPtr<MultiBlockCoordSys> a_coordSys);
  
  /// Define physics models
  void definePhysics(const int               a_physicsType,
                     const std::vector<int>& a_periodic);

  /// Define turbulence modeling parameters
  void defineTurbulence(
    const TurbModelType&                  a_turbModelType,
    const SGSModelType&                   a_sgsModelType,
    const RealVect&                       a_globStreamwiseDir,
    const bool                            a_turbForce,
    const bool                            a_useConsToPrimCorrection,
    const bool                            a_useSGSCoarsening,
    CRDparam::ExplicitFilterType          a_explicitFilterType,
    SpectralFilter::SpectralFilterProfile a_spectralFilterProfile,
    const std::vector<int>&               a_spectralFilterDomainResolution,
    const Real                            a_spectralFilterParam,
    const TurbForcingType&                a_turbForcingType,
    const Interval&                       a_spectralForcingInterval,
    const Real                            a_spectralForcingEps,
    const Real                            a_spectralForcingDt,
    const int                             a_sgskeCrsLevFilterRatio,
    const int                             a_wallModelCrsRatio,
    const bool                            a_enforceModeledWallShearStress,
    const Real                            a_startTimeAvgTime);

  /// Define initializing geometry
  void defineGeometry(const std::vector<int>&  a_domainBaseSize,
                      const std::vector<Real>& a_domainLength,
                      const std::vector<Real>& a_domainOrigin);

  /// Define initializing information about coordinate system type
  void defineCoordSysType(const CoordSysType a_coordSysType,
                          const RealVect&    a_physicalLength);

  /// Define initializing fluid parameters
  void defineFluid(const int  a_numStates,
                   const Real a_rho,
                   const Real a_T,
                   const Real a_speed,
                   const Real a_mu,
                   const Real a_lambda,
                   const Real a_gamma,
                   const Real a_R,
                   const Real a_grav,
                   const Real a_Re,
                   const Real a_K);

  /// Define combustion parameters
  void defineCombustion(
    const int                       a_numSpecies,
    const Real                      a_reactionStartTime,
    const Real                      a_reactionEndTime,
    const bool                      a_useSpeciesCorrection,
    const int                       a_numReactions,
    const std::vector<std::string>& a_speciesNames);

  /// Define initializing small reference thermodynamics
  void defineSmall(const Real a_smallr, const Real a_smallp);

  /// Define initializing solver
  void defineSolver(const Real a_cfl,
                    const Real a_initialDt,
                    const Real a_maxDt,
                    const Real a_initialCFL,
                    const Real a_initialTime,
                    const CRDparam::TimeIntegrationMethod
                      a_timeIntegrationMethod);

  /// Define Additive Runge Kutta parameters
  void defineARK(const bool a_additiveRK,
                 const int  a_ARKmaxLevel,
                 const bool a_additiveRKUseAnalyticJac,
                 const Real a_chemicalDtScale,
                 const int  a_initERKSteps,
                 const bool a_poutLinearSolverStats,
                 const Real a_nonlinearTol,
                 const bool a_ARKUsePIDControl,
                 const Real a_ARKPIDControlEps,
                 const bool a_ARKExtrapInitGuess);

  /// Define initializing limiting
  void defineLimiter(const int  a_faceInterpolationOrder,
                     const int  a_diffusiveDerivativeOrder,
                     const int  a_reactionOrder,
                     const bool a_usePPMlimiter,
                     const bool a_limitFaceValues,
                     const bool a_useFlattening,
                     const int  a_cellConvolveFlatten,
                     const int  a_cellDeconvolveFlatten,
                     const int  a_faceConvolveFlatten,
                     const int  a_faceDeconvolveFlatten,
                     const int  a_cellConvolveLimit,
                     const int  a_cellDeconvolveLimit,
                     const int  a_faceConvolveLimit,
                     const int  a_faceDeconvolveLimit,
                     const bool a_useFCOR,
                     const bool a_useArtificialViscosity,
                     const Real a_artificialViscosityCoef,
                     const bool a_useArtificialViscosity4thO,
                     const Real a_artificialViscosity4thOCoef,
                     const bool a_extraBoundLim,
                     const bool a_noHOchecks,
                     const bool a_wallCorrections,
                     const bool a_clipping,
                     const bool a_clippingHO,
                     const bool a_clippingPostSmooth,
                     const Real a_fifthOrderBlendingCoef);

  /// Define initializing AMR
  void defineAMR(const int               a_minBoxSize,
                 const int               a_maxBoxSize,
                 const std::vector<int>& a_refRatios,
                 const bool              a_useSubcycling);

  /// Define initializing physics class
  void defineCRDPhysics(const CRDPhysics *const a_CRDPhysics);

  /// Define initializing a vector of physics class
  void defineCRDPhysics(const int               a_numPhysics,
                        const CRDPhysics *const a_CRDPhysics[]);

  /// Define DC flattening
  void defineDCF(const DCFlattening *const a_DCF);

  /// Define initializing IBC
  void defineIBC(const CNSIBC *const a_CNSIBC);

#ifdef CH_CTHR
  /// Define ThreadTeamArchitect
  void defineThreads(ThreadTools::ThreadTeamArchitectChombo_t *const a_threads);
#endif

//--Adjustments are used to fix discrepancies, i.e., difference between
//--input file and checkpoint file or conflicts between speed and Re

  /// Because the user can enter conflicting speed and Re, allow adjustment
  void adjustFluidSpeed(const Real a_speed, const Real a_length);
  void adjustFluidRe(const Real a_Re, const Real a_length);

  /// Adjust refinement ratios to match restart file
  void adjustRefFromBase(const int a_level, const int a_val);

//--Parameters used in helper routines

  const int* numGhostList() const { return m_numGhostList; }
  int numCRDPhysics() const { return m_CRDPhysics.size(); }

//--Fetch parameters

  const int& verbosity() const { return m_verbosity; }
  const bool& verboseDt() const { return m_verboseDt; }
  const ProblemType& problemType() const { return m_problemType; }
  const CoordSysType& coordSysType() const { return m_coordSysType; }
  const RefCountedPtr<MultiBlockCoordSys>& coordSys() const
    { return m_coordSys; }
  const int& physicsModels() const { return m_physicsModels; }
  const TurbModelType& turbModelType() const { return m_turbModelType; }
  const SGSModelType& sgsModelType() const { return m_sgsModelType; }
  const IntVect& domainBaseSize() const { return m_domainBaseSize; }
  const RealVect& domainLength() const { return m_domainLength; }
  const RealVect& domainOrigin() const { return m_domainOrigin; }
  const RealVect& physicalLength() const { return m_physicalLength; }
  const IntVect& periodicity() const { return m_periodicity; }
  const bool& cartesian() const { return m_cartesian; }
  const int& numStates() const { return m_numStates; }
  const int& numTransportScalars() { return m_numTransportScalars; }
  const Real& rho() const { return m_rho; }
  const Real& T() const { return m_T; }
  const Real& K() const { return m_K; }
  const Real& speed() const { return m_speed; }
  const Real& small() const { return m_small; }
  const Real& smallr() const { return m_smallr; }
  const Real& smallp() const { return m_smallp; }
  const Real& mu() const { return m_mu; }
  const Real& lambda() const { return m_lambda; }
  const Real& gamma() const { return m_gamma; }
  const Real& R() const { return m_R; }
  const Real& grav() const { return m_grav; }
  const Real& Re() const { return m_Re; }
  const Real& cfl() const { return m_cfl; }
  const Real& initialDt() const { return m_initialDt; }
  const Real& maxDt() const { return m_maxDt; }
  const Real& initialCFL() const { return m_initialCFL; }
  const CRDparam::TimeIntegrationMethod&
              timeIntegrationMethod() const { return m_timeIntegrationMethod; }
  const bool& additiveRK() const { return m_additiveRK; }
  const int&  ARKmaxLevel() const { return m_ARKmaxLevel; }
  const bool& additiveRKUseAnalyticJac() const
    { return m_additiveRKUseAnalyticJac; }
  const Real& chemicalDtScale() const { return m_chemicalDtScale; }
  const Real& ARKNonlinearTol() const { return m_ARKNonlinearTol; }
  const int&  initERKSteps() const { return m_initERKSteps; }
  const bool& poutLinearSolverStats() const { return m_poutLinearSolverStats; }
  const bool& ARKUsePIDControl() const { return m_ARKUsePIDControl; }
  const Real& ARKPIDControlEps() const { return m_ARKPIDControlEps; }
  const bool& ARKExtrapInitGuess() const { return m_ARKExtrapInitGuess; }
  const Real& initialTime() const { return m_initialTime; }
  const int& minBoxSize() const { return m_minBoxSize; }
  const int& maxBoxSize() const { return m_maxBoxSize; }
  const std::vector<int>& refFromBase() const { return m_refFromBase; }
  const bool& useSubcycling() const { return m_useSubcycling; }
  const int& faceInterpolationOrder() const { return m_faceInterpolationOrder; }
  const int& diffusiveDerivativeOrder() const
    { return m_diffusiveDerivativeOrder; }
  const int& reactionOrder() const { return m_reactionOrder; }
  const bool& usePPMlimiter() const { return m_usePPMlimiter; }
  const bool& limitFaceValues() const { return m_limitFaceValues; }
  const bool& useFlattening() const { return m_useFlattening; }
  const int& cellConvolveFlatten() const { return m_cellConvolveFlatten; }
  const int& cellDeconvolveFlatten() const { return m_cellDeconvolveFlatten; }
  const int& faceConvolveFlatten() const { return m_faceConvolveFlatten; }
  const int& faceDeconvolveFlatten() const { return m_faceDeconvolveFlatten; }
  const int& cellConvolveLimit() const { return m_cellConvolveLimit; }
  const int& cellDeconvolveLimit() const { return m_cellDeconvolveLimit; }
  const int& faceConvolveLimit() const { return m_faceConvolveLimit; }
  const int& faceDeconvolveLimit() const { return m_faceDeconvolveLimit; }
  const bool& useFCOR() const { return m_useFCOR; }
  const bool& useArtificialViscosity() const
    { return m_useArtificialViscosity; }
  const Real& artificialViscosityCoef() const
    { return m_artificialViscosityCoef; }
  const bool& useArtificialViscosity4thO() const
    { return m_useArtificialViscosity4thO; }
  const Real& artificialViscosity4thOCoef() const
    { return m_artificialViscosity4thOCoef; }
  const bool& extraBoundLim() const { return m_extraBoundLim; }
  const bool& noHOchecks() const { return m_noHOchecks; }
  const bool& wallCorrections() const { return m_wallCorrections; }
  const bool& clipping() const { return m_clipping; }
  const bool& clippingHO() const { return m_clippingHO; }
  const bool& clippingPostSmooth() const { return m_clippingPostSmooth; }
  const Real& fifthOrderBlendingCoef() const
    { return m_fifthOrderBlendingCoef; }
  const int& numSpecies() const { return m_numSpecies; }
  const Real& reactionStartTime() const { return m_reactionStartTime; }
  const Real& reactionEndTime() const { return m_reactionEndTime; }
  const bool& useSpeciesCorrection() const { return m_useSpeciesCorrection; }
  const bool& useTurbForce() const { return m_turbForce; }
  const bool& useConsToPrimCorrection() const
    { return m_useConsToPrimCorrection; }
  const bool& useSGSCoarsening() const { return m_useSGSCoarsening; }
  const int& numReactions() const { return m_numReactions; }
  const std::vector<std::string>& speciesNames() const
    { return m_speciesNames; }
  const std::string& plotPrefix() const { return m_plotPrefix; }
  const int& plotNumGhost() const { return m_plotNumGhost; }
  const bool& plotJ() const { return m_plotJ; }
  const bool& plotJU() const { return m_plotJU; }
  const bool& plotDACFDCheck() const { return m_plotDACFDCheck; }
  const bool& plotFlattening() const { return m_plotFlattening; }
  const bool& plotExtraVars() const { return m_plotExtraVars; }
  const bool& plotMappedDerivs() const { return m_plotMappedDerivs; }
  const bool& plotWallDist() const { return m_plotWallDist; }
  const bool& plotError() const { return m_plotError; }
  const bool& plotTimeAvgTurb() const { return m_plotTimeAvgTurb; }
  const bool& plotLoFaceCoordinates() const { return m_plotLoFaceCoordinates; }
  const bool& plotHiFaceCoordinates() const { return m_plotHiFaceCoordinates; }
  const bool& restartAddWallModel() const { return m_restartAddWallModel; }
  const int& plotLoFaceAvgTimeAvgComps() const
    { return m_plotLoFaceAvgTimeAvgComps; }
  const int& plotLoFaceAvgComps() const { return m_plotLoFaceAvgComps; }
  const ExplicitFilterType& explicitFilterType() const
    { return m_explicitFilterType; }
  const SpectralFilter::SpectralFilterProfile& spectralFilterProfile() const
    { return m_spectralFilterProfile; }
  const RealVect& globStreamwiseDir() const { return m_globStreamwiseDir; }
  const IntVect& spectralFilterDomainResolution() const
    { return m_spectralFilterDomainResolution; }
  const Real& spectralFilterParam() const
    { return m_spectralFilterParam; }
  const TurbForcingType& turbForcingType() const
    { return m_turbForcingType; }
  const Interval& spectralForcingInterval() const
    { return m_spectralForcingInterval; }
  const Real& spectralForcingEps() const
    { return m_spectralForcingEps; }
  const Real& spectralForcingDt() const
    { return m_spectralForcingDt; }
  const int&  sgskeCrsLevFilterRatio() const
    { return m_sgskeCrsLevFilterRatio; }
  const int&  wallModelCrsRatio() const
    { return m_wallModelCrsRatio; }
  const bool& enforceModeledWallShearStress() const
    { return m_enforceModeledWallShearStress; }
  const Real& startTimeAvgTime() const
    { return m_startTimeAvgTime; }
  const CRDPhysics* const& CRDPhysicsOp(const int a_set = 0) const
    { return m_CRDPhysics[a_set].getRefToThePointer(); }
  const CNSIBC* const& CNSIBCOp() const
    { return m_CNSIBC.getRefToThePointer(); }
  const DCFlattening* const& DCFOp() const
    { return m_DCF.getRefToThePointer(); }
#ifdef CH_CTHR
  ThreadTools::ThreadTeamArchitectChombo_t* const& threads()
    { return m_threads.getRefToThePointer(); }
#endif

/*==============================================================================
 * Data members
 *============================================================================*/

private:

  // Storage of parameters
  int m_numGhostList[NumGhostQuerySize];
  int m_verbosity;
  bool m_verboseDt;
  ProblemType m_problemType;
  CoordSysType m_coordSysType;
  RefCountedPtr<MultiBlockCoordSys> m_coordSys;
  int m_physicsModels;
  TurbModelType m_turbModelType;
  SGSModelType m_sgsModelType;
  IntVect m_domainBaseSize;
  RealVect m_domainLength;
  RealVect m_domainOrigin;
  RealVect m_physicalLength;
  IntVect m_periodicity;
  bool m_cartesian;
  int m_numStates;
  int m_numTransportScalars;
  Real m_rho;
  Real m_T;
  Real m_K;
  Real m_speed;
  Real m_small;
  Real m_smallr;
  Real m_smallp;
  Real m_mu;
  Real m_lambda;
  Real m_gamma;
  Real m_R;
  Real m_grav;
  Real m_Re;
  Real m_cfl;
  Real m_initialDt;
  Real m_maxDt;
  Real m_initialCFL;
  CRDparam::TimeIntegrationMethod m_timeIntegrationMethod;
  bool m_additiveRK;
  int  m_ARKmaxLevel;
  bool m_additiveRKUseAnalyticJac;
  Real m_chemicalDtScale;
  Real m_ARKNonlinearTol;
  int  m_initERKSteps;
  bool m_poutLinearSolverStats;
  bool m_ARKUsePIDControl;
  Real m_ARKPIDControlEps;
  bool m_ARKExtrapInitGuess;
  Real m_initialTime;
  int  m_minBoxSize;
  int  m_maxBoxSize;
  std::vector<int> m_refFromBase;
  bool m_useSubcycling;
  int  m_faceInterpolationOrder;
  int  m_diffusiveDerivativeOrder;
  int  m_reactionOrder;
  bool m_usePPMlimiter;
  bool m_limitFaceValues;
  bool m_useFlattening;
  int  m_cellConvolveFlatten;
  int  m_cellDeconvolveFlatten;
  int  m_faceConvolveFlatten;
  int  m_faceDeconvolveFlatten;
  int  m_cellConvolveLimit;
  int  m_cellDeconvolveLimit;
  int  m_faceConvolveLimit;
  int  m_faceDeconvolveLimit;
  bool m_useFCOR;
  bool m_useArtificialViscosity;
  Real m_artificialViscosityCoef;
  bool m_useArtificialViscosity4thO;  // Unused
  Real m_artificialViscosity4thOCoef;
  bool m_extraBoundLim;
  bool m_noHOchecks;
  bool m_wallCorrections;
  bool m_clipping;
  bool m_clippingHO;
  bool m_clippingPostSmooth;
  Real m_fifthOrderBlendingCoef;
  int  m_numSpecies;
  Real m_reactionStartTime;
  Real m_reactionEndTime;
  bool m_useSpeciesCorrection;
  bool m_turbForce;
  bool m_useConsToPrimCorrection;
  bool m_useSGSCoarsening;
  int  m_numReactions;
  std::vector<std::string> m_speciesNames;
  std::string m_plotPrefix;
  int  m_plotNumGhost;
  bool m_plotJ;
  bool m_plotJU;
  bool m_plotDACFDCheck;
  bool m_plotFlattening;
  bool m_plotExtraVars;
  bool m_plotMappedDerivs;
  bool m_plotWallDist;
  bool m_plotError;
  bool m_plotTimeAvgTurb;
  bool m_plotLoFaceCoordinates;
  bool m_plotHiFaceCoordinates;
  bool m_restartAddWallModel;
  int m_plotLoFaceAvgTimeAvgComps;
  int m_plotLoFaceAvgComps;
  ExplicitFilterType m_explicitFilterType;
  SpectralFilter::SpectralFilterProfile m_spectralFilterProfile;
  RealVect m_globStreamwiseDir;
  IntVect m_spectralFilterDomainResolution;
  Real m_spectralFilterParam;
  TurbForcingType m_turbForcingType;
  Interval m_spectralForcingInterval;
  Real m_spectralForcingEps;
  Real m_spectralForcingDt;
  int  m_sgskeCrsLevFilterRatio;
  int  m_wallModelCrsRatio;
  bool m_enforceModeledWallShearStress;
  Real m_startTimeAvgTime;
  std::vector<RefCountedPtr<const CRDPhysics>> m_CRDPhysics;
  RefCountedPtr<const CNSIBC> m_CNSIBC;
  RefCountedPtr<const DCFlattening> m_DCF;
#ifdef CH_CTHR
  RefCountedPtr<ThreadTools::ThreadTeamArchitectChombo_t> m_threads;
#endif
#ifndef NDEBUG
  // Consistency flags (only checked through asserts)
  unsigned m_definedComponents;
public:
  int m_definedQueryGhosts;
#endif
};

//--References for parameters that are constant for a run

extern CRDparamVar CRDP;              ///< Location of parameters
extern const int& g_verbosity;        ///< Amount of info printed
extern const bool& g_verboseDt;       ///< T : print detailed time step info
extern const ProblemType& g_problemType;
                                      ///< Type of problem to solve
extern const CoordSysType& g_coordSysType;
                                      ///< Type of coordinate system to use
extern const RefCountedPtr<MultiBlockCoordSys>& g_coordSys;
                                      ///< The coordinate system on base level
extern const int& g_physicsModels;    ///< Types of physics models
extern const TurbModelType& g_turbModelType;
                                      ///< Type of turbulence models
extern const SGSModelType& g_sgsModelType;
                                      ///< Types of SGS models
extern const IntVect& g_domainBaseSize;
                                      ///< Number of cells in the base problem
                                      ///< domain
extern const RealVect& g_domainLength;
                                      ///< Length of the physical domain
extern const RealVect& g_domainOrigin;
			              ///< Physical location of node (0,0,0)
extern const RealVect& g_physicalLength;
                                      ///< Physical domain length for Cartesian
extern const IntVect& g_periodicity;  ///< Domain periodicity
extern const bool& g_cartesian;       ///< T - Cartesian domain.
                                      ///< Note: do not make any optimizations
                                      ///< using this without consulting
                                      ///< primary developers.
extern const int& g_numStates;        ///< Number of state variables
extern const int& g_numTransportScalars;
                                      ///< Simple transport of scalars
extern const Real& g_rho;             ///< Freestream density (kg/m^3)
extern const Real& g_T;               ///< Freestream temperature (K)
extern const Real& g_K;               ///< Thermal conductivity (W/m-K)
extern const Real& g_speed;           ///< Freestream speed (m/s)
extern const Real& g_small;           ///< Small reference
extern const Real& g_smallr;          ///< Small reference density (kg/m^3)
extern const Real& g_smallp;          ///< Small reference pressure (Pa)
extern const Real& g_mu;              ///< Freestream kinematic viscosity
                                      ///< (m^2/s)
extern const Real& g_lambda;

extern const Real& g_gamma;           ///< Specific head ratio
extern const Real& g_R;               ///< Specific gas constant
extern const Real& g_grav;            ///< Graviational constant
extern const Real& g_Re;              ///< Reynolds number (if applicable)
extern const Real& g_cfl;             ///< CFL number
extern const Real& g_initialDt;       ///< Initial time step size
extern const Real& g_maxDt;           ///< Maximum time step size
extern const Real& g_initialCFL;      ///< Initial CFL number
extern const CRDparam::TimeIntegrationMethod&
                   g_timeIntegrationMethod;
                                      ///< The time integration method to use
extern const bool& g_additiveRK;      ///< Use Additive RK instead of explicit RK
extern const int&  g_ARKmaxLevel;     ///< Max level to apply ARK time stepping
extern const bool& g_additiveRKUseAnalyticJac;
                                      ///< Uses an analytical Jacobian instead of
                                      ///< numerical Jacobian computed via
                                      ///< difference quotient.
extern const Real& g_chemicalDtScale; ///< Scaling factor for the chemical time
                                      ///< step size
extern const Real& g_ARKNonlinearTol; ///< Convergence tolerance for ARK's nonlinear solver
extern const int&  g_initERKSteps;    ///< Number of explicit RK time steps to
                                      ///< take before doing additive RK time
                                      ///< stepping
extern const bool& g_poutLinearSolverStats;
                                      ///< Print out linear solver stats
extern const bool& g_ARKUsePIDControl;///< Use the PID control for ARK time step size
extern const Real& g_ARKPIDControlEps;///< Epsilon for ARK's PID controller
extern const bool& g_ARKExtrapInitGuess;
                                      ///< Extrapolate the nonlinear solver's initial guess
                                      ///  if true, otherwise use the previous stage
                                      ///  solution as the initial guess
extern const Real& g_initialTime;     ///< Initial solution time
extern const int&  g_minBoxSize;      ///< Mininum box size
extern const int&  g_maxBoxSize;      ///< Maximum box size
extern const std::vector<int>& g_refFromBase;
                                      ///< Refinement from base grid
extern const bool& g_useSubcycling;   ///< Bool to use subcycling or not
extern const int&  g_faceInterpolationOrder;
                                      ///< 1 - First order Godunov
                                      ///< 4 - Fourth order
                                      ///< 5 - Fifth order upwind
extern const int&  g_diffusiveDerivativeOrder;
                                      ///< 2 - Second order
                                      ///< 4 - Fourth order
extern const int&  g_reactionOrder;   ///< 2 - Second order
                                      ///< 4 - Fourth order
extern const bool& g_usePPMlimiter;   ///< Selects whether to use PPM limiting
                                      ///< in cells.  Must be selected to use
                                      ///< flattening.
extern const bool& g_limitFaceValues; ///< T - Apply limiting when calculating
                                      ///<     state on faces
                                      ///< Only available when
                                      ///< g_faceInterpolationOrder == 4
extern const bool& g_useFlattening;   ///< T - Use flattening near steep
                                      ///< gradients
extern const int& g_cellConvolveFlatten;
                                      ///< 0 - Convolve cells
                                      ///< 1 - Limit cell convolution
                                      ///< 2 - Do not convolve cells
extern const int& g_cellDeconvolveFlatten;
extern const int& g_faceConvolveFlatten;
extern const int& g_faceDeconvolveFlatten;
                                      ///< Flatten the cell/face convolutions
                                      ///< and deconvolutions
                                      ///< 0 - Do no flatten
                                      ///< 1 - Apply non-linear flattening
                                      ///< 2 - No convolutions, 2nd order

extern const int& g_cellConvolveLimit;///< Limit the cell convolutions
                                      ///< 0 - Do no limit
                                      ///< 1 - Apply limiting
extern const int& g_cellDeconvolveLimit;
                                      ///< Limit the cell deconvolutions
                                      ///< 0 - Do no limit
                                      ///< 1 - Apply limiting
extern const int& g_faceConvolveLimit;///< Limit the face convolutions
                                      ///< 0 - Do no limit
                                      ///< 1 - Apply limiting
extern const int& g_faceDeconvolveLimit;
                                      ///< Limit the face deconvolutions
                                      ///< 0 - Do no limit
                                      ///< 1 - Apply limiting


extern const bool& g_useFCOR;         ///< T - Use face construction order
                                      ///< reduction
extern const bool& g_useArtificialViscosity;
                                      ///< T - Use artificial viscosity
extern const Real& g_artificialViscosityCoef;
                                      ///< T - Coefficient of artificial
                                      ///< viscosity
extern const bool& g_useArtificialViscosity4thO;
                                      ///< NOTE: This is not implemented!
                                      ///< T - Use 4thO artificial viscosity
extern const Real& g_artificialViscosity4thOCoef;
                                      ///< NOTE: This is still used for regular
                                      ///< artificial viscosity!
                                      ///< T - Strong shock threshold for
                                      ///<     fourth-order artificial viscosity
extern const bool& g_extraBoundLim;   ///< T - Skip 3rd derivative checks
                                      ///<     near boundaries
extern const bool& g_noHOchecks;      ///< T - Skip all 3rd derivative checks
extern const bool& g_wallCorrections; ///< T - Correct thermo state at walls
extern const bool& g_clipping;        ///< T - Prevent AMR interpolation overshoots
extern const bool& g_clippingHO;      ///< T - But allow smooth regions to overshoot
extern const bool& g_clippingPostSmooth;
                                      ///< T - Smooth solution over entire regrid
extern const Real& g_fifthOrderBlendingCoef;
                                      ///< Fifth-order blending coefficient
extern const int& g_numSpecies;       ///< Number of species
extern const Real& g_reactionStartTime;
                                      ///< Time to start including reactions
extern const Real& g_reactionEndTime; ///< Time to end including reactions
extern const bool& g_useSpeciesCorrection;
                                      ///< Turns species correction on/off
extern const bool& g_useTurbForce;    ///< Turns turbulent forcing on/off
extern const bool& g_useConsToPrimCorrection;
                                      ///< Turns pressure correction on/off
                                      ///< Specific to LES SGS turbulence models
extern const bool& g_useSGSCoarsening;
                                      ///< Turns sgs model coarsening on/off
                                      ///< Specific to LES SGS turbulence models
extern const int& g_numReactions;     ///< Number of reactions to be modeled
extern const std::vector<std::string>& g_speciesNames;
                                      ///< List of species names
extern const std::string& g_plotPrefix;
                                      ///< Plot prefix
extern const int& g_plotNumGhost;     ///< Number of ghost cells to write to
                                      ///< plot file
extern const bool& g_plotJ;           ///< T - Write the metrics Jacobian to
                                      ///<     the plot file
extern const bool& g_plotJU;          ///< T - Write the conservative state
                                      ///<     \<JU\> from computational space
extern const bool&  g_plotDACFDCheck; ///< T - Write the conservative state
                                      ///<     (\<U\>) and temperature to
                                      ///<     the checkpoint file
                                      ///<     (specially, only for DA application)
extern const bool& g_plotFlattening;  ///< T - Plot the flattening coefficients
extern const bool& g_plotExtraVars;   ///< T - Write out kappa and mu, expensive
extern const bool& g_plotMappedDerivs;///< T - Write mapped vorticity and
                                      ///<     similar
extern const bool& g_plotWallDist;    ///< T - Write wall distance field
extern const bool& g_plotError;       ///< T - Write error field
extern const bool& g_plotTimeAvgTurb; ///< T - Write time-avg turbulence vars
extern const bool& g_plotLoFaceCoordinates;
                                      ///< T - Write low-face coordinates
extern const bool& g_plotHiFaceCoordinates;
                                      ///< T - Write high-face coordinates
extern const bool& g_restartAddWallModel;
                                      ///< T - Restart case and add wall-model
extern const int& g_plotLoFaceAvgTimeAvgComps;
                                      ///< Component groups to plot to hdf5
extern const int& g_plotLoFaceAvgComps;
                                      ///< Component groups to plot to hdf5
extern const ExplicitFilterType& g_explicitFilterType;
                                      ///< T - Type of explicit filter to apply
extern const SpectralFilter::SpectralFilterProfile& g_spectralFilterProfile;
                                      ///< T - Profile of spectral filter
extern const RealVect& g_globStreamwiseDir;
                                      ///< Global streamwise direction for
                                      ///< SV SGS LES wall-model
extern const IntVect& g_spectralFilterDomainResolution;
                                      ///< T - Apply spectral filter to state
extern const Real& g_spectralFilterParam;
                                      ///< T - Spectral filter width
extern const TurbForcingType& g_turbForcingType;
                                      ///< type of forcing to use to sustain
                                      ///< turbulence in domain
extern const Interval& g_spectralForcingInterval;
                                      ///< Wavenumber band-pass interval of
                                      ///< spectral forcing
extern const Real& g_spectralForcingEps;  ///< Energy injection rate of spectral
                                      ///< forcing
extern const Real& g_spectralForcingDt;
                                      ///< Time interval for calculating a new
                                      ///< spectral forcing source field
extern const int&  g_sgskeCrsLevFilterRatio;
                                      ///< Ratio of coarsened-coarsest level
                                      ///< and coarsest level -- to be used with
                                      ///< the coarsened SV SGS KE estimate
extern const int& g_wallModelCrsRatio;///< Ratio of coarsening for wall-model
extern const bool& g_enforceModeledWallShearStress;
                                      ///< Enforce modeled wall shear-stress
extern const Real& g_startTimeAvgTime;
                                      ///< Start time for time-averaging data
extern const CRDPhysics* const& g_CRDPhysics;
                                      ///< Physics used in the simulation
                                      ///< (Use as a pointer to MOLPhysics)
extern const CNSIBC* const& g_CNSIBC; ///< Initial and boundary conditions used
                                      ///< in the simulation
                                      ///< (Use as a pointer to CNSIBC)
extern const DCFlattening* const& g_DCF;
                                      ///< D/C flattening method
#ifdef CH_CTHR
extern ThreadTools::ThreadTeamArchitectChombo_t* const& g_threads;
                                      ///< Creates, binds, and contains threads
#endif
inline const CRDPhysics* const& CRDPhysicsSet(const int a_set = 0)
{
  return CRDP.CRDPhysicsOp(a_set);
}

//--Parameters that are variable for a run (these should only be used for
//--debugging or reporting)

extern int g_level;                   ///< The AMR mesh level


/*******************************************************************************
 *
 * Functions returning global constants
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Return the maximum AMR level (base is 0)
/*--------------------------------------------------------------------*/

inline int
maxAMRLevel()
{
  return g_refFromBase.size() - 1;
}

/*--------------------------------------------------------------------*/
/// Return the number of  AMR levels
/*--------------------------------------------------------------------*/

inline int
numAMRLevel()
{
  return g_refFromBase.size();
}

/*--------------------------------------------------------------------*/
/// A query on the number of ghosts
/** These queries all return cell-centered boxes but often define
 *  quantities on faces.
 *  \param[in]  a_query A particular data type
 *//*-----------------------------------------------------------------*/

inline int
queryGhosts(const NumGhostQuery& a_query)
{
  CH_assert(!(CRDparam::CRDP.m_definedQueryGhosts));
  return CRDparam::CRDP.numGhostList()[(int)a_query];
}

/*--------------------------------------------------------------------*/
/// Return the number of  AMR levels
/*--------------------------------------------------------------------*/

inline int
numCRDPhysics()
{
  return CRDparam::CRDP.numCRDPhysics();
}

}  // namespace CRDparam

#endif /* ! defined _CRDPARAM_H_ */
