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
 * \file CRDparam.cpp
 *
 * \brief Member functions and definitions of global variables for CRDparam
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <numeric>

//----- Chombo Library -----//

#include "CRDPhysics.H"
#include "DCFlattening.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CNSIBC.H"

/*******************************************************************************
 *
 * Class CRDparamVar: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor sets defaults
/** \note
 *  <ul>
 *    <li> Most of these defaults are never used.  Defaults are
 *         instead assigned in 'ChordInput.cpp' where they may be
 *         overridden by users.
 *  </ul>
 *//*-----------------------------------------------------------------*/

CRDparam::CRDparamVar::CRDparamVar()
 :
  // Output
  m_verbosity(1),
  m_verboseDt(false),
  // Simulation
  m_problemType(ProblemUndefined),
  m_coordSysType(CoordSysSingleBlockCartesian),
  m_physicsModels(0),
  m_turbModelType(TurbModelUndefined),
  m_sgsModelType(SGSModelUndefined),
  // Geometry
  m_domainBaseSize(IntVect::Unit),
  m_domainLength(RealVect::Unit),
  m_domainOrigin(RealVect::Unit),
  m_physicalLength(RealVect::Unit),
  m_periodicity(IntVect::Unit),
  m_cartesian(false),
  // Fluid parameters (generally freestream ICAO STP)
  m_numStates(SpaceDim + 2),
  m_numTransportScalars(0),
  m_rho(1.2250),    // Density (kg/m^3)
  m_T(288.15),      // Temperature (K)
  m_K(2.5326E-2),   // Thermal conductivity (W/m-K)
  m_speed(0.),      // Speed (m/s)
  m_small(1.E-6),   // Small reference
  m_smallr(1.E-6),  // Small reference density (kg/m^3)
  m_smallp(1.E-6),  // Small reference pressure (Pa)
  m_mu(1.789E-5),   // Dynamic viscosity (kg/m s)
  m_lambda(-1.0/3.0),// Second kinematic viscosity coefficient
  m_gamma(1.4),     // Specific heat ratio
  m_R(287.),        // Specific gas constant (J/kg-K)
  m_grav(9.80665),  // Gravitational constant (m/s^2)
  m_Re(1.0),        // Reynolds number
  // Solver
  m_cfl(1.),
  m_initialDt(-1.),
  m_maxDt(-1.),
  m_initialCFL(1.),
  m_additiveRK(false),
  m_ARKmaxLevel(0),
  m_additiveRKUseAnalyticJac(false),
  m_chemicalDtScale(1.),
  m_ARKNonlinearTol(1.e-12),
  m_initERKSteps(0),
  m_poutLinearSolverStats(false),
  m_ARKUsePIDControl(false),
  m_ARKPIDControlEps(5.e-10),
  m_ARKExtrapInitGuess(false),
  m_initialTime(0.),
  // AMR
  m_minBoxSize(8),
  m_maxBoxSize(32),
  m_refFromBase(1, 1),
  // Limiting
  m_faceInterpolationOrder(4),
  m_diffusiveDerivativeOrder(4),
  m_reactionOrder(4),
  m_usePPMlimiter(true),
  m_limitFaceValues(false),
  m_useFlattening(true),
  m_cellConvolveFlatten(0),
  m_cellDeconvolveFlatten(0),
  m_faceConvolveFlatten(0),
  m_faceDeconvolveFlatten(0),
  m_cellConvolveLimit(0),
  m_cellDeconvolveLimit(0),
  m_faceConvolveLimit(0),
  m_faceDeconvolveLimit(0),
  m_useFCOR(false),
  m_useArtificialViscosity(true),
  m_artificialViscosityCoef(0.3),
  m_useArtificialViscosity4thO(true),
  m_artificialViscosity4thOCoef(0.42),
  m_extraBoundLim(false),
  m_noHOchecks(false),
  m_wallCorrections(true),
  m_numSpecies(0),
  m_reactionStartTime(-1.),
  m_reactionEndTime(1.E6),
  m_useSpeciesCorrection(true),
  m_numReactions(0),
  // Files & plotting
  m_plotPrefix(""),
  m_plotNumGhost(0),
  m_plotJ(true),
  m_plotJU(true),
  m_plotDACFDCheck(false),
  m_plotFlattening(false),
  m_plotExtraVars(false),
  m_plotMappedDerivs(false),
  m_plotWallDist(false),
  m_plotError(false),
  m_plotTimeAvgTurb(false),
  m_plotLoFaceCoordinates(false),
  m_plotHiFaceCoordinates(false),
  m_restartAddWallModel(false),
  m_plotLoFaceAvgTimeAvgComps(0),
  m_plotLoFaceAvgComps(0),
  // Turbulence modeling
  m_explicitFilterType(ExplicitFilterNone),
  m_spectralFilterProfile(SpectralFilter::SpectralFilterSharp),
  m_spectralFilterDomainResolution(IntVect::Unit),
  m_spectralFilterParam(1.),
  m_turbForcingType(TurbForcingNone),
  m_spectralForcingInterval(Interval(2, 4)),
  m_spectralForcingEps(0.),
  m_spectralForcingDt(1.e99),
  m_sgskeCrsLevFilterRatio(1),
  m_wallModelCrsRatio(1),
  m_enforceModeledWallShearStress(false),
  m_startTimeAvgTime(0.),
  // Singleton classes
  m_CRDPhysics(1, RefCountedPtr<const CRDPhysics>(nullptr)),
  m_CNSIBC(nullptr),
  m_DCF(nullptr)
#ifdef CH_CTHR
  ,
  m_threads(nullptr)
#endif
  // Consistency flags (only checked through asserts)
#ifndef NDEBUG
  ,
  m_definedComponents(DefinedComponentNone),
  /* Bit 0 removed when defineLimiter is called.
   * Bit 1 removed when definePhysics is called.
   * Ghosts can be queried when all bits unset.
   */
  m_definedQueryGhosts((1<<0) | (1<<1))
#endif
{
  m_speciesNames.resize(1);
  m_speciesNames.assign(1,"");
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CRDparam::CRDparamVar::~CRDparamVar()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define initializing verbosity
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineInfo(const int  a_verbosity,
                                  const bool a_verboseDt)
{
  m_verbosity = a_verbosity;
  m_verboseDt = a_verboseDt;
}

/*--------------------------------------------------------------------*/
//  Define initializing plot and file information
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineFiles(
  const std::string& a_plotPrefix,
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
  const int          a_plotLoFaceAvgComps)
{
  m_plotPrefix = a_plotPrefix;
  m_plotNumGhost = a_plotNumGhost;
  m_plotExtraVars = a_plotExtraVars;
  m_plotMappedDerivs = a_plotMappedDerivs;
  m_plotWallDist = a_plotWallDist;
  m_plotJ = a_plotJ;
  m_plotJU = a_plotJU;
  m_plotDACFDCheck = a_plotDACFDCheck;
  m_plotFlattening = a_plotFlattening;
  m_plotError = a_plotError;
  m_plotTimeAvgTurb = a_plotTimeAvgTurb;
  m_plotLoFaceCoordinates = a_plotLoFaceCoordinates;
  m_plotHiFaceCoordinates = a_plotHiFaceCoordinates;
  m_restartAddWallModel = a_restartAddWallModel;
  m_plotLoFaceAvgTimeAvgComps = a_plotLoFaceAvgTimeAvgComps;
  m_plotLoFaceAvgComps = a_plotLoFaceAvgComps;
}

/*--------------------------------------------------------------------*/
//  Define initializing simulation type
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineSimulation(const ProblemType&      a_problemType,
                                        const std::vector<int>& a_periodicity)
{
  m_problemType = a_problemType;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      m_periodicity[dir] = a_periodicity[dir];
    }
}

/*--------------------------------------------------------------------*/
//  Define the coordinate system
/*--------------------------------------------------------------------*/
void
CRDparam::CRDparamVar::defineCoordSys(
  const RefCountedPtr<MultiBlockCoordSys> a_coordSys)
{
  m_coordSys = a_coordSys;
}

/*--------------------------------------------------------------------*/
//  Define physics models
/** 'defineLimiter' must be called first
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::definePhysics(const int               a_physicsType,
                                     const std::vector<int>& a_periodic)
{
#ifndef NDEBUG
  CH_assert(!(m_definedQueryGhosts & (1<<0)));
  m_definedQueryGhosts &= ~(1<<1);
#endif
  m_physicsModels = a_physicsType;
  const int hasBoundary =
    (std::accumulate(a_periodic.begin(), a_periodic.end(), (int)0) == SpaceDim)
    ? 0 : 1;

  // FIXME: Might only be necessary for viscous
  if ((m_physicsModels & PhysicsViscous) ||
      (m_physicsModels & PhysicsThermPerf))  // Viscous
    {
      if (m_usePPMlimiter || (m_sgsModelType & SGSModelStretchedVortex) ||
          (m_faceInterpolationOrder == 5))   // If limited
        {
          m_numGhostList[NumGhostWcellExtTngBdry ] = hasBoundary*4;
                                                         // <W> in exterior
                                                         // ghost cells
          // Reduce cost of solving second-order face values for FCOR
          if (CRDparam::g_useFCOR)
            {
              m_numGhostList[NumGhostInertialFfaceAvg] = 0;
                                                         // <F> (inertial flux)
            }
          else
            {
              m_numGhostList[NumGhostInertialFfaceAvg] = 1;
                                                         // <F> (inertial flux)
            }                                            // on face
          if (m_faceDeconvolveLimit && (m_faceDeconvolveFlatten < 2))
          // Limit face deconvolutions
            {
              m_numGhostList[NumGhostWfaceAvgTan ] = 3;  // <W> on face
            }
          else
            {
              m_numGhostList[NumGhostWfaceAvgTan ] = 2;  // <W> on face
            }
                                                         // (tangent)
          m_numGhostList[NumGhostWfaceAvgNrm     ] = 2;  // <W> on face (normal)
          if (m_faceDeconvolveLimit && (m_faceDeconvolveFlatten < 2))
          // Limit face deconvolutions
            {
              m_numGhostList[NumGhostWpFaceAvgTan] = 4;  // <Wp> (native
                                                         // primitive state) on
                                                         // face in tangent Dir.
            }
          else
            {
              m_numGhostList[NumGhostWpFaceAvgTan] = 3;
            }
          m_numGhostList[NumGhostWcellAvg        ] = 6;  // <W> in cell
          // If SVS SGS LES consToPrim correction is used
          if (m_useConsToPrimCorrection)
            {
              m_numGhostList[NumGhostUcellAvg    ] = 8;  // <U> in cell
              m_numGhostList[NumGhostUcellPnt    ] = 6;  // (U) in cell
            }
          else
            {
              if (m_cellDeconvolveLimit && (m_cellDeconvolveFlatten < 2))
              // Limit cell deconvolutions
                {
                  m_numGhostList[NumGhostUcellAvg] = 8;  // <U> in cell
                }
              else
                {
                  m_numGhostList[NumGhostUcellAvg] = 7;  // <U> in cell
                }
              m_numGhostList[NumGhostUcellPnt    ] = 6;  // (U) in cell
            }
        }
      else                               // If not limited
        {
          m_numGhostList[NumGhostWcellExtTngBdry ] = hasBoundary*4;
                                                         // <W> in exterior
                                                         // ghost cells
          m_numGhostList[NumGhostInertialFfaceAvg] = 1;  // <F> (inertial flux)
                                                         // on face
          m_numGhostList[NumGhostWfaceAvgTan     ] = 2;  // <W> on face
                                                         // (tangent)
          m_numGhostList[NumGhostWfaceAvgNrm     ] = 2;  // <W> on face (normal)
          m_numGhostList[NumGhostWpFaceAvgTan    ] = 3;  // <Wp> (native
                                                         // primitive state) on
                                                         // face in tangent Dir.
          m_numGhostList[NumGhostWcellAvg        ] = 4;  // <W> in cell
          m_numGhostList[NumGhostUcellPnt        ] = 4;  // (U) in cell
          m_numGhostList[NumGhostUcellAvg        ] = 6;  // <U> in cell
        }
    }
  else                                   // Inertial only
    {
      m_numGhostList[NumGhostWcellExtTngBdry ] = 0;  // <W> in exterior
                                                     // ghost cells
      m_numGhostList[NumGhostInertialFfaceAvg] = 0;  // <F> (inertial flux) on
                                                     // face
      m_numGhostList[NumGhostWfaceAvgTan     ] = 0;  // <W> on face (tangent)
//    m_numGhostList[NumGhostInertialFfacePnt] = 0;
//    m_numGhostList[NumGhostWfacePntTan     ] = 0;
//    m_numGhostList[NumGhostWpFacePntTan    ] = 0;
      m_numGhostList[NumGhostWfaceAvgNrm     ] = 0;  // <W> on face (normal)
//    m_numGhostList[NumGhostWfacePntNrm     ] = 0;
//    m_numGhostList[NumGhostWpFacePntNrm    ] = 0;
//    m_numGhostList[NumGhostWpFaceAvgNrm    ] = 0;
      if (m_faceDeconvolveLimit && (m_faceDeconvolveFlatten < 2))
      // Limit face deconvolutions
        {
          m_numGhostList[NumGhostWpFaceAvgTan] = 2;  // <Wp> (native primitive
                                                     // state) on face in
                                                     // tangent directions
        }
      else
        {
          m_numGhostList[NumGhostWpFaceAvgTan] = 1;  // <Wp> (native primitive
                                                     // state) on face in
                                                     // tangent directions
        }
//    m_numGhostList[NumGhostFfromWpFaceAvg   ] = 1;
      m_numGhostList[NumGhostWcellAvg         ] = 4; // <W> in cell
//    m_numGhostList[NumGhostWcellPnt         ] = 4;
      m_numGhostList[NumGhostUcellPnt         ] = 4; // (U) in cell
      if (m_cellDeconvolveLimit && (m_cellDeconvolveFlatten < 2))
      // Limit cell deconvolutions
        {
          m_numGhostList[NumGhostUcellAvg     ] = 6; // <U> in cell
        }
      else
        {
          m_numGhostList[NumGhostUcellAvg     ] = 5; // <U> in cell
        }
    }
}

/*--------------------------------------------------------------------*/
//  Define turbulence model type
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineTurbulence(
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
  const Real                            a_startTimeAvgTime)
{
  m_turbModelType                  = a_turbModelType;
  m_sgsModelType                   = a_sgsModelType;
  m_globStreamwiseDir              = a_globStreamwiseDir;
  m_turbForce                      = a_turbForce;
  m_useConsToPrimCorrection        = a_useConsToPrimCorrection;
  m_useSGSCoarsening               = a_useSGSCoarsening;
  m_explicitFilterType             = a_explicitFilterType;
  m_spectralFilterProfile          = a_spectralFilterProfile;
  m_spectralFilterParam            = a_spectralFilterParam;
  m_turbForcingType                = a_turbForcingType;
  m_spectralForcingInterval        = a_spectralForcingInterval;
  m_spectralForcingEps             = a_spectralForcingEps;
  m_spectralForcingDt              = a_spectralForcingDt;
  m_sgskeCrsLevFilterRatio         = a_sgskeCrsLevFilterRatio;
  m_wallModelCrsRatio              = a_wallModelCrsRatio;
  m_enforceModeledWallShearStress  = a_enforceModeledWallShearStress;
  m_startTimeAvgTime               = a_startTimeAvgTime;
  for (int i = 0; i != SpaceDim; ++i)
    {
      m_spectralFilterDomainResolution[i] =
        a_spectralFilterDomainResolution[i];
    }
}

/*--------------------------------------------------------------------*/
//  Define initializing geometry
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineGeometry(const std::vector<int>&  a_domainBaseSize,
                                      const std::vector<Real>& a_domainLength,
                                      const std::vector<Real>& a_domainOrigin)
{
#ifndef NDEBUG
  m_definedComponents |= DefinedComponentGeometry;
#endif
  for (int i = 0; i != SpaceDim; ++i)
    {
      m_domainBaseSize[i] = a_domainBaseSize[i];
      m_domainLength[i]   = a_domainLength[i];
      m_domainOrigin[i]   = a_domainOrigin[i];
    }
  m_spectralFilterDomainResolution = m_domainBaseSize;
}

/*--------------------------------------------------------------------*/
//  Define initializing information about coordinate system type
/** Only used for some optimizations if domain is single-block
 *  Cartesian
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineCoordSysType(
  const CoordSysType a_coordSysType,
  const RealVect&    a_physicalLength)
{
  // Make sure defineGeometry has been called first
  CH_assert(m_definedComponents & DefinedComponentGeometry);

  m_coordSysType = a_coordSysType;

  bool cartBool = true;
  for (int i = 0; i != SpaceDim; ++i)
    {
      m_physicalLength[i] = a_physicalLength[i];
      if (m_physicalLength[i] != m_domainLength[i])
        {
          cartBool = false;
        }
    }
  if (a_coordSysType == CoordSysSingleBlockCartesian && cartBool)
    {
      m_cartesian = true;
    }
}

/*--------------------------------------------------------------------*/
//  Define initializing fluid parameters
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineFluid(const int  a_numStates,
                                   const Real a_rho,
                                   const Real a_T,
                                   const Real a_speed,
                                   const Real a_mu,
                                   const Real a_lambda,
                                   const Real a_gamma,
                                   const Real a_R,
                                   const Real a_grav,
                                   const Real a_Re,
                                   const Real a_K)
{
  m_numStates = a_numStates;
  m_rho       = a_rho;
  m_T         = a_T;
  m_speed     = a_speed;
  m_mu        = a_mu;
  m_lambda    = a_lambda;
  m_gamma     = a_gamma;
  m_R         = a_R;
  m_grav      = a_grav;
  m_Re        = a_Re;
  m_K         = a_K;
}

/*--------------------------------------------------------------------*/
//  Define initializing combustion parameters
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineCombustion(
  const int                       a_numSpecies,
  const Real                      a_reactionStartTime,
  const Real                      a_reactionEndTime,
  const bool                      a_useSpeciesCorrection,
  const int                       a_numReactions,
  const std::vector<std::string>& a_speciesNames)
{
  m_numSpecies = a_numSpecies;
  m_reactionStartTime = a_reactionStartTime;
  m_reactionEndTime = a_reactionEndTime;
  m_useSpeciesCorrection = a_useSpeciesCorrection;
  m_numReactions = a_numReactions;
  m_speciesNames = a_speciesNames;
}

/*--------------------------------------------------------------------*/
//  Define initializing small reference thermodynamics
/** \param[in]  a_smallr
 *                      Small reference density
 *  \param[in]  a_smallp
 *                      Small reference pressure
 *  These are usually calculated during IBC construction
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineSmall(const Real a_smallr, const Real a_smallp)
{
  m_smallr = a_smallr;
  m_smallp = a_smallp;
}

/*--------------------------------------------------------------------*/
//  Define initializing solver
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineSolver(const Real a_cfl,
                                    const Real a_initialDt,
                                    const Real a_maxDt,
                                    const Real a_initialCFL,
                                    const Real a_initialTime,
                                    CRDparam::TimeIntegrationMethod
                                    a_timeIntegrationMethod)
{
  m_cfl = a_cfl;
  m_initialDt = a_initialDt;
  m_maxDt = a_maxDt;
  m_initialCFL = a_initialCFL;
  m_initialTime = a_initialTime;
  m_timeIntegrationMethod = a_timeIntegrationMethod;
}

/*--------------------------------------------------------------------*/
//  Comments!!!
/*--------------------------------------------------------------------*/

void CRDparam::CRDparamVar::defineARK(const bool a_additiveRK,
                                      const int  a_ARKmaxLevel,
                                      const bool a_additiveRKUseAnalyticJac,
                                      const Real a_chemicalDtScale,
                                      const int  a_initERKSteps,
                                      const bool a_poutLinearSolverStats,
                                      const Real a_ARKNonlinearTol,
                                      const bool a_ARKUsePIDControl,
                                      const Real a_ARKPIDControlEps,
                                      const bool a_ARKExtrapInitGuess) {
  m_additiveRK = a_additiveRK;
  m_ARKmaxLevel = a_ARKmaxLevel;
  m_additiveRKUseAnalyticJac = a_additiveRKUseAnalyticJac;
  m_chemicalDtScale = a_chemicalDtScale;
  m_initERKSteps = a_initERKSteps;
  m_poutLinearSolverStats = a_poutLinearSolverStats;
  m_ARKNonlinearTol = a_ARKNonlinearTol;
  m_ARKUsePIDControl = a_ARKUsePIDControl;
  m_ARKPIDControlEps = a_ARKPIDControlEps;
  m_ARKExtrapInitGuess = a_ARKExtrapInitGuess;
}

/*--------------------------------------------------------------------*/
//  Define initializing limiting
/*--------------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineLimiter(const int  a_faceInterpolationOrder,
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
                                     const Real a_fifthOrderBlendingCoef)
{
#ifndef NDEBUG
  m_definedQueryGhosts &= ~(1<<0);
#endif
  m_faceInterpolationOrder      = a_faceInterpolationOrder;
  m_diffusiveDerivativeOrder    = a_diffusiveDerivativeOrder;
  m_reactionOrder               = a_reactionOrder;
  m_usePPMlimiter               = a_usePPMlimiter;
  m_limitFaceValues             = a_limitFaceValues;
  m_useFlattening               = a_useFlattening;
  m_cellConvolveFlatten         = a_cellConvolveFlatten;
  m_cellDeconvolveFlatten       = a_cellDeconvolveFlatten;
  m_faceConvolveFlatten         = a_faceConvolveFlatten;
  m_faceDeconvolveFlatten       = a_faceDeconvolveFlatten;
  m_cellConvolveLimit           = a_cellConvolveLimit;
  m_cellDeconvolveLimit         = a_cellDeconvolveLimit;
  m_faceConvolveLimit           = a_faceConvolveLimit;
  m_faceDeconvolveLimit         = a_faceDeconvolveLimit;
  m_useFCOR                     = a_useFCOR;
  m_useArtificialViscosity      = a_useArtificialViscosity;
  m_artificialViscosityCoef     = a_artificialViscosityCoef;
  m_useArtificialViscosity4thO  = a_useArtificialViscosity4thO;
  m_artificialViscosity4thOCoef = a_artificialViscosity4thOCoef;
  m_extraBoundLim               = a_extraBoundLim;
  m_noHOchecks                  = a_noHOchecks;
  m_wallCorrections             = a_wallCorrections;
  m_clipping                    = a_clipping;
  m_clippingHO                  = a_clippingHO;
  m_clippingPostSmooth          = a_clippingPostSmooth;
  m_fifthOrderBlendingCoef      = a_fifthOrderBlendingCoef;
}

/*--------------------------------------------------------------------*/
//  Define D/C flattening
/** \param[in] a_DCF    Base pointer to DCFlattening class. Allocate with
 *                      new and do not delete.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineDCF(const DCFlattening *const a_DCF)
{
  m_DCF = RefCountedPtr<const DCFlattening>(a_DCF);
}

/*--------------------------------------------------------------------*/
//  Define initializing AMR
/** \param[in]  a_minBoxSize
 *                      Minimum box size (usually set as block factor)
 *  \param[in]  a_maxBoxSize
 *                      Maximum box size
 *  \param[in]  a_refRatios
 *                      Refinement ratios between levels
 *  \param[in]  a_useSubcyling
 *                      T - use subcycling
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineAMR(const int               a_minBoxSize,
                                 const int               a_maxBoxSize,
                                 const std::vector<int>& a_refRatios,
                                 const bool              a_useSubcycling)
{
  m_minBoxSize = a_minBoxSize;
  m_maxBoxSize = a_maxBoxSize;
  const int numLevel = a_refRatios.size();
  m_refFromBase.resize(numLevel);
  int refTotal = 1;
  for (int i = 0; i != numLevel; ++i)
    {
      m_refFromBase[i] = refTotal;
      refTotal *= a_refRatios[i];
    }
  m_useSubcycling = a_useSubcycling;
}

/*--------------------------------------------------------------------*/
//  Define initializing physics class
/** \param[in] a_CRDPhysics
 *                      Base pointer to physics class.  Allocate with
 *                      new and do not delete.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineCRDPhysics(const CRDPhysics *const a_CRDPhysics)
{
  m_CRDPhysics[0] = RefCountedPtr<const CRDPhysics>(a_CRDPhysics);
}

/*--------------------------------------------------------------------*/
//  Define initializing a vector of physics classes
/** \param[in] a_numPhysics
 *                      Number of physics classes
 *  \param[in] a_CRDPhysics
 *                      C-style array of base pointers to physics
 *                      classes.  Allocate with new and do not delete.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineCRDPhysics(const int               a_numPhysics,
                                        const CRDPhysics *const a_CRDPhysics[])
{
  // If using this, remove g_CRDPhysics and fix all code to use
  // CRDparam::CRDPhysicsSet() instead.  g_CRDPhysics will segfault.
  m_CRDPhysics.assign(a_numPhysics, RefCountedPtr<const CRDPhysics>(nullptr));
  for (int i = 0; i != a_numPhysics; ++i)
    {
      m_CRDPhysics[i] = RefCountedPtr<const CRDPhysics>(a_CRDPhysics[i]);
    }
}

/*--------------------------------------------------------------------*/
//  Define initializing IBC
/** \param[in] a_CNSIBS
 *                      Base pointer to CNS initial and boundary
 *                      conditions class.  Allocate with new and do
 *                      not delete.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineIBC(const CNSIBC *const a_CNSIBC)
{
  m_CNSIBC = RefCountedPtr<const CNSIBC>(a_CNSIBC);
}

#ifdef CH_CTHR
/*--------------------------------------------------------------------*/
//  Define ThreadTeamArchitect
/** \param[in] a_CNSIBS
 *                      Base pointer to CNS initial and boundary
 *                      conditions class.  Allocate with new and do
 *                      not delete.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::defineThreads(
  ThreadTools::ThreadTeamArchitectChombo_t *const a_threads)
{
  m_threads = RefCountedPtr<ThreadTools::ThreadTeamArchitectChombo_t>(a_threads);
}
#endif

/*--------------------------------------------------------------------*/
//  Because the user can enter conflicting speed and Re, allow
//  adjustment
/** This should only be done during initialization.  Re and speed
 *  should still be considered constant.
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::adjustFluidSpeed(const Real a_speed, const Real a_length)
{
  m_speed = a_speed;
  m_Re = a_speed*a_length*m_rho/m_mu;
}

void
CRDparam::CRDparamVar::adjustFluidRe(const Real a_Re, const Real a_length)
{
  m_Re = a_Re;
  m_speed = a_Re*m_mu/(m_rho*a_length);
}

/*--------------------------------------------------------------------*/
//  Adjust refinement ratios to match restart file
/** This should only be done when restarting from a checkpoint file
 *//*-----------------------------------------------------------------*/

void
CRDparam::CRDparamVar::adjustRefFromBase(const int a_level, const int a_val)
{
  CH_assert(a_level < m_refFromBase.size());
  m_refFromBase[a_level] = a_val;
}


/*******************************************************************************
 *
 * Class CRDparamVar: static data member definitions
 *
 ******************************************************************************/

const int&      CRDparam::g_verbosity = CRDparam::CRDP.verbosity();
const bool&     CRDparam::g_verboseDt = CRDparam::CRDP.verboseDt();
const CRDparam::ProblemType& CRDparam::g_problemType =
  CRDparam::CRDP.problemType();
const CRDparam::CoordSysType& CRDparam::g_coordSysType =
  CRDparam::CRDP.coordSysType();
const RefCountedPtr<MultiBlockCoordSys>& CRDparam::g_coordSys =
  CRDparam::CRDP.coordSys();
const CRDparam::TurbModelType& CRDparam::g_turbModelType =
  CRDparam::CRDP.turbModelType();
const CRDparam::SGSModelType& CRDparam::g_sgsModelType =
  CRDparam::CRDP.sgsModelType();
const int&      CRDparam::g_physicsModels = CRDparam::CRDP.physicsModels();
const IntVect&  CRDparam::g_domainBaseSize = CRDparam::CRDP.domainBaseSize();
const RealVect& CRDparam::g_domainLength = CRDparam::CRDP.domainLength();
const RealVect& CRDparam::g_domainOrigin = CRDparam::CRDP.domainOrigin();
const RealVect& CRDparam::g_physicalLength = CRDparam::CRDP.physicalLength();
const IntVect&  CRDparam::g_periodicity = CRDparam::CRDP.periodicity();
const bool&     CRDparam::g_cartesian = CRDparam::CRDP.cartesian();
const int&      CRDparam::g_numStates = CRDparam::CRDP.numStates();
const int&      CRDparam::g_numTransportScalars =
  CRDparam::CRDP.numTransportScalars();
const Real&     CRDparam::g_rho = CRDparam::CRDP.rho();
const Real&     CRDparam::g_T = CRDparam::CRDP.T();
const Real&     CRDparam::g_K = CRDparam::CRDP.K();
const Real&     CRDparam::g_speed = CRDparam::CRDP.speed();
const Real&     CRDparam::g_small = CRDparam::CRDP.small();
const Real&     CRDparam::g_smallr = CRDparam::CRDP.smallr();
const Real&     CRDparam::g_smallp = CRDparam::CRDP.smallp();
const Real&     CRDparam::g_mu = CRDparam::CRDP.mu();
const Real&     CRDparam::g_lambda = CRDparam::CRDP.lambda();
const Real&     CRDparam::g_gamma = CRDparam::CRDP.gamma();
const Real&     CRDparam::g_R = CRDparam::CRDP.R();
const Real&     CRDparam::g_grav = CRDparam::CRDP.grav();
const Real&     CRDparam::g_Re = CRDparam::CRDP.Re();
const Real&     CRDparam::g_cfl = CRDparam::CRDP.cfl();
const Real&     CRDparam::g_initialDt = CRDparam::CRDP.initialDt();
const Real&     CRDparam::g_maxDt = CRDparam::CRDP.maxDt();
const Real&     CRDparam::g_initialCFL = CRDparam::CRDP.initialCFL();
const CRDparam::TimeIntegrationMethod& CRDparam::g_timeIntegrationMethod =
  CRDparam::CRDP.timeIntegrationMethod();
const bool&     CRDparam::g_additiveRK = CRDparam::CRDP.additiveRK();
const int&      CRDparam::g_ARKmaxLevel = CRDparam::CRDP.ARKmaxLevel();
const bool&     CRDparam::g_additiveRKUseAnalyticJac =
  CRDparam::CRDP.additiveRKUseAnalyticJac();
const Real&     CRDparam::g_chemicalDtScale = CRDparam::CRDP.chemicalDtScale();
const Real&     CRDparam::g_ARKNonlinearTol = CRDparam::CRDP.ARKNonlinearTol();
const int&      CRDparam::g_initERKSteps = CRDparam::CRDP.initERKSteps();
const bool&     CRDparam::g_poutLinearSolverStats =
  CRDparam::CRDP.poutLinearSolverStats();
const bool&     CRDparam::g_ARKUsePIDControl =
  CRDparam::CRDP.ARKUsePIDControl();
const Real&     CRDparam::g_ARKPIDControlEps =
  CRDparam::CRDP.ARKPIDControlEps();
const bool&     CRDparam::g_ARKExtrapInitGuess =
  CRDparam::CRDP.ARKExtrapInitGuess();
const Real&     CRDparam::g_initialTime = CRDparam::CRDP.initialTime();
const int&      CRDparam::g_minBoxSize = CRDparam::CRDP.minBoxSize();
const int&      CRDparam::g_maxBoxSize = CRDparam::CRDP.maxBoxSize();
const std::vector<int>& CRDparam::g_refFromBase =
  CRDparam::CRDP.refFromBase();
const bool&     CRDparam::g_useSubcycling = CRDparam::CRDP.useSubcycling();
const int&      CRDparam::g_faceInterpolationOrder =
  CRDparam::CRDP.faceInterpolationOrder();
const int&      CRDparam::g_diffusiveDerivativeOrder =
  CRDparam::CRDP.diffusiveDerivativeOrder();
const int&      CRDparam::g_reactionOrder = CRDparam::CRDP.reactionOrder();
const bool&     CRDparam::g_usePPMlimiter = CRDparam::CRDP.usePPMlimiter();
const bool&     CRDparam::g_limitFaceValues = CRDparam::CRDP.limitFaceValues();
const bool&     CRDparam::g_useFlattening = CRDparam::CRDP.useFlattening();
const int&      CRDparam::g_cellConvolveFlatten =
  CRDparam::CRDP.cellConvolveFlatten();
const int&      CRDparam::g_cellDeconvolveFlatten =
  CRDparam::CRDP.cellDeconvolveFlatten();
const int&      CRDparam::g_faceConvolveFlatten =
  CRDparam::CRDP.faceConvolveFlatten();
const int&      CRDparam::g_faceDeconvolveFlatten =
  CRDparam::CRDP.faceDeconvolveFlatten();
const int&      CRDparam::g_cellConvolveLimit =
  CRDparam::CRDP.cellConvolveLimit();
const int&      CRDparam::g_cellDeconvolveLimit =
  CRDparam::CRDP.cellDeconvolveLimit();
const int&      CRDparam::g_faceConvolveLimit =
  CRDparam::CRDP.faceConvolveLimit();
const int&      CRDparam::g_faceDeconvolveLimit =
  CRDparam::CRDP.faceDeconvolveLimit();
const bool&     CRDparam::g_useFCOR = CRDparam::CRDP.useFCOR();
const bool&     CRDparam::g_useArtificialViscosity =
  CRDparam::CRDP.useArtificialViscosity();
const Real&     CRDparam::g_artificialViscosityCoef =
  CRDparam::CRDP.artificialViscosityCoef();
const bool&     CRDparam::g_useArtificialViscosity4thO =
  CRDparam::CRDP.useArtificialViscosity4thO();
const Real&     CRDparam::g_artificialViscosity4thOCoef =
  CRDparam::CRDP.artificialViscosity4thOCoef();
const bool&     CRDparam::g_extraBoundLim = CRDparam::CRDP.extraBoundLim();
const bool&     CRDparam::g_noHOchecks = CRDparam::CRDP.noHOchecks();
const bool&     CRDparam::g_wallCorrections = CRDparam::CRDP.wallCorrections();
const bool&     CRDparam::g_clipping = CRDparam::CRDP.clipping();
const bool&     CRDparam::g_clippingHO = CRDparam::CRDP.clippingHO();
const bool&     CRDparam::g_clippingPostSmooth =
  CRDparam::CRDP.clippingPostSmooth();
const Real&     CRDparam::g_fifthOrderBlendingCoef =
  CRDparam::CRDP.fifthOrderBlendingCoef();
const int&      CRDparam::g_numSpecies = CRDparam::CRDP.numSpecies();
const Real&     CRDparam::g_reactionStartTime =
  CRDparam::CRDP.reactionStartTime();
const Real&     CRDparam::g_reactionEndTime = CRDparam::CRDP.reactionEndTime();
const bool&     CRDparam::g_useSpeciesCorrection =
  CRDparam::CRDP.useSpeciesCorrection();
const bool&     CRDparam::g_useTurbForce = CRDparam::CRDP.useTurbForce();
const bool&     CRDparam::g_useConsToPrimCorrection =
  CRDparam::CRDP.useConsToPrimCorrection();
const bool&     CRDparam::g_useSGSCoarsening =
  CRDparam::CRDP.useSGSCoarsening();
const int&      CRDparam::g_numReactions = CRDparam::CRDP.numReactions();
const std::vector<std::string>& CRDparam::g_speciesNames =
  CRDparam::CRDP.speciesNames();
const std::string& CRDparam::g_plotPrefix = CRDparam::CRDP.plotPrefix();
const int&      CRDparam::g_plotNumGhost = CRDparam::CRDP.plotNumGhost();
const bool&     CRDparam::g_plotJ = CRDparam::CRDP.plotJ();
const bool&     CRDparam::g_plotJU = CRDparam::CRDP.plotJU();
const bool&     CRDparam::g_plotDACFDCheck = CRDparam::CRDP.plotDACFDCheck();
const bool&     CRDparam::g_plotFlattening = CRDparam::CRDP.plotFlattening();
const bool&     CRDparam::g_plotExtraVars = CRDparam::CRDP.plotExtraVars();
const bool&     CRDparam::g_plotMappedDerivs =
  CRDparam::CRDP.plotMappedDerivs();
const bool&     CRDparam::g_plotWallDist = CRDparam::CRDP.plotWallDist();
const bool&     CRDparam::g_plotError = CRDparam::CRDP.plotError();
const bool&     CRDparam::g_plotTimeAvgTurb =
  CRDparam::CRDP.plotTimeAvgTurb();
const bool&     CRDparam::g_plotLoFaceCoordinates =
  CRDparam::CRDP.plotLoFaceCoordinates();
const bool&     CRDparam::g_plotHiFaceCoordinates =
  CRDparam::CRDP.plotHiFaceCoordinates();
const bool&     CRDparam::g_restartAddWallModel =
  CRDparam::CRDP.restartAddWallModel();
const int&     CRDparam::g_plotLoFaceAvgTimeAvgComps =
  CRDparam::CRDP.plotLoFaceAvgTimeAvgComps();
const int&     CRDparam::g_plotLoFaceAvgComps =
  CRDparam::CRDP.plotLoFaceAvgComps();
const CRDparam::ExplicitFilterType& CRDparam::g_explicitFilterType =
  CRDparam::CRDP.explicitFilterType();
const SpectralFilter::SpectralFilterProfile& CRDparam::g_spectralFilterProfile =
  CRDparam::CRDP.spectralFilterProfile();
const RealVect& CRDparam::g_globStreamwiseDir =
  CRDparam::CRDP.globStreamwiseDir();
const IntVect&  CRDparam::g_spectralFilterDomainResolution =
  CRDparam::CRDP.spectralFilterDomainResolution();
const Real&     CRDparam::g_spectralFilterParam =
  CRDparam::CRDP.spectralFilterParam();
const CRDparam::TurbForcingType& CRDparam::g_turbForcingType =
  CRDparam::CRDP.turbForcingType();
const Interval& CRDparam::g_spectralForcingInterval =
  CRDparam::CRDP.spectralForcingInterval();
const Real&     CRDparam::g_spectralForcingEps =
  CRDparam::CRDP.spectralForcingEps();
const Real&     CRDparam::g_spectralForcingDt =
  CRDparam::CRDP.spectralForcingDt();
const int&      CRDparam::g_sgskeCrsLevFilterRatio =
  CRDparam::CRDP.sgskeCrsLevFilterRatio();
const int&      CRDparam::g_wallModelCrsRatio =
  CRDparam::CRDP.wallModelCrsRatio();
const bool&     CRDparam::g_enforceModeledWallShearStress =
  CRDparam::CRDP.enforceModeledWallShearStress();
const Real&     CRDparam::g_startTimeAvgTime =
  CRDparam::CRDP.startTimeAvgTime();

//--Parameters that are variable for a run (these should only be used for
//--debugging or reporting)

int CRDparam::g_level = -1;
