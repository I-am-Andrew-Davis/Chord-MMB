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
 * \file ChordInput.cpp
 *
 * \brief Routines for input to Chord
 *
 *//*+*************************************************************************/

//----- System -----//

#ifdef CH_USE_CGNS
  #ifdef CH_MPI
    #include "pcgnslib.h"
  #else
    #include "cgnslib.h"
  #endif
#endif

//----- Chombo Library -----//

#include "ParmParse.H"
#include "CH_System.H"
#include "Misc.H"
#include "LevelGridMetrics.H"
// single-block mapping
#include "CartesianCS.H"
#include "WarpedCS.H"
#include "TwistedCS.H"
#include "LogCS.H"
#include "AnnulusCS.H"
#include "SkewCS.H"
#include "ExternalCS.H"
#include "SchwarzChristoffelRampCS.H"
#include "LogSchwarzChristoffelRampCS.H"
#include "SingleBlockCSAdaptor.H"
// multi-block mappings
#include "DoubleCartesianCS.H"
#include "ChannelSeparationCS.H"
#include "ExternalMultiBlockCS.H"
#include "CubedSphereShell3DCS.H"
#include "SmoothRampCS.H"
#include "JoukowskiAirfoilCS.H"

//----- Internal -----//

#include "ChordInput.H"
#include "CRDmsg.H"

// Precision to print floating point numbers
const int l_prec = 4;


/*==============================================================================
 * Helper functions
 *============================================================================*/

// Extract ramp parameters from input file
static void rampParameters(ParmParse& a_pp,
                           Real&      a_alpha,
                           Real&      a_leadLength,
                           Real&      a_rampLength)
{
  a_pp.get("ramp_angle", a_alpha);
  a_pp.get("lead_length", a_leadLength);
  a_pp.get("ramp_length", a_rampLength);
}

// Extract log shift parameters from input file
static void logshiftParameters(ParmParse& a_pp, RealVect& a_shift)
{
  std::vector<Real> COORDSYSshiftVec(SpaceDim, 0.);
  a_pp.queryarr("shift", COORDSYSshiftVec, 0, SpaceDim);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      a_shift[dir] = COORDSYSshiftVec[dir];
    }
}


/*============================================================================*/
//  Parse the input, check, and print
/**
 *
 *//*=========================================================================*/

void chordInput(ChomboParam& a_chparam)
{
  CRD::msg.setTerminateOnError(false);
  CRD::msg << "Parameters\n" << CRD::h1;

/*--------------------------------------------------------------------*
 * File parameters
 *--------------------------------------------------------------------*/

//--Read options about file and IO

  CRD::msg << "Parameters for 'file'" << CRD::h2;
  ParmParse ppFILE("file");

//--Read the restart file

  a_chparam.FILErestartFile = "";
  a_chparam.FILEdoRestart = false;
  if (ppFILE.contains("restart_file"))
    {
      ppFILE.get("restart_file", a_chparam.FILErestartFile);
      if (!CH_System::fileExists(a_chparam.FILErestartFile.c_str()))
        {
          CRD::msg << "Input: 'restart_file' " << a_chparam.FILErestartFile
                   << " does not exist!" << CRD::error;
        }
      CRD::msg << "Restart file\n" << a_chparam.FILErestartFile << CRD::var;
      a_chparam.FILEdoRestart = true;
    }
  else
    {
      CRD::msg << "Restart file\nnot provided" << CRD::var;
    }

//--Read the plot file prefix and interval

  a_chparam.FILEplotPrefix = "";
  a_chparam.FILEplotInterval = -1;
  a_chparam.FILEplotPeriod = -1.;
  // Local parameters
  int lclPlotNumGhost = 0;
  int lclPlotJ = 0;
  int lclPlotJU = 0;
  int lclPlotDACFDCheck = 0;
  int lclPlotFlattening = 0;
  int lclPlotExtraVars = 0;
  int lclPlotMappedDerivs = 0;
  int lclPlotWallDist = 0;
  int lclPlotError = 0;
  int lclPlotTimeAvgTurb = 0;
  int lclPlotLoFaceCoordinates = 0;
  int lclPlotHiFaceCoordinates = 0;
  int lclRestartAddWallModel = 0;
  int lclPlotLoFaceAvgTimeAvgComps = CRDparam::PlotUndefined;
  int lclPlotLoFaceAvgComps = CRDparam::PlotUndefined;
  if (ppFILE.contains("plot_prefix"))
    {
      // Query plot prefix
      ppFILE.get("plot_prefix", a_chparam.FILEplotPrefix);
      // Only test existance of path if '/' present
      size_t pos = a_chparam.FILEplotPrefix.rfind('/');
      if (pos != 0 && pos != std::string::npos)
        {
          std::string path = a_chparam.FILEplotPrefix.substr(0, pos);
          if (!CH_System::fileExists(path.c_str()))
            {
              CRD::msg << "Input: 'plot_prefix' path '" << path
                       << "' does not exist!" << CRD::error;
            }
        }
      CRD::msg << "Plot prefix\n" << a_chparam.FILEplotPrefix << CRD::var;
      // Query plot interval
      ppFILE.query("plot_interval", a_chparam.FILEplotInterval);
      ppFILE.query("plot_period", a_chparam.FILEplotPeriod);
      if (a_chparam.FILEplotInterval > 0 && a_chparam.FILEplotPeriod > 0.)
        {
          CRD::msg << "Input: 'plot_interval' and 'plot_period' cannot both"
                   << " be specified!" << CRD::error;
        }
      // Query to plot J
      ppFILE.query("plot_J", lclPlotJ);
      if (lclPlotJ != 0)
        {
          lclPlotJ = 1;
        }
      CRD::msg << "Plot J\n" << ((lclPlotJ) ? "true" : "false") << CRD::var;
      // Query to plot JU
      ppFILE.query("plot_JU", lclPlotJU);
      if (lclPlotJU != 0)
        {
          lclPlotJU = 1;
        }
      CRD::msg << "Plot JU\n" << ((lclPlotJU) ? "true" : "false") << CRD::var;
      // Query to plot UtoCheck
      ppFILE.query("plot_DACFDCheck", lclPlotDACFDCheck);
      if (lclPlotDACFDCheck != 0)
        {
          lclPlotDACFDCheck = 1;
        }
      CRD::msg << "Plot DACFDCheck\n"
               << ((lclPlotDACFDCheck) ? "true" : "false") << CRD::var;
      // Query to plot flattening coefficients
      ppFILE.query("plot_flattening", lclPlotFlattening);
      if (lclPlotFlattening != 0)
        {
          lclPlotFlattening = 1;
        }
      CRD::msg << "Plot flattening\n"
               << ((lclPlotFlattening) ? "true" : "false") << CRD::var;
      // Select number of ghosts to plot
      ppFILE.query("plot_ghost", lclPlotNumGhost);
      if (lclPlotNumGhost != 0 && lclPlotNumGhost != 1)
        {
          CRD::msg << "Input: 'plot_ghost' must be 0 or 1!" << CRD::error;
        }
      CRD::msg << "Plot number of ghosts\n" << lclPlotNumGhost << CRD::var;
      // Query to plot extra stuff
      ppFILE.query("plot_extra_vars",lclPlotExtraVars);
      CRD::msg << "Plot extra variables\n" <<
        ((lclPlotExtraVars) ? "true" : "false") << CRD::var;
      // Query to plot vorticity and divergence
      ppFILE.query("plot_mapped_derivs",lclPlotMappedDerivs);
      CRD::msg << "Plot mapped variables\n" <<
        ((lclPlotMappedDerivs) ? "true" : "false") << CRD::var;
      // Query to plot wall distance
      ppFILE.query("plot_wall_distance",lclPlotWallDist);
      CRD::msg << "Plot wall distance field\n" <<
        ((lclPlotWallDist) ? "true" : "false") << CRD::var;
      // Query to plot error from exact solution
      ppFILE.query("plot_error",lclPlotError);
      CRD::msg << "Plot error compared to exact solution\n" <<
        ((lclPlotError) ? "true" : "false") << CRD::var;
      // Query to plot time-averaged velocity field
      ppFILE.query("plot_time_averaged_turb_var",lclPlotTimeAvgTurb);
      CRD::msg << "Plot time-averaged velocity field\n" <<
        ((lclPlotTimeAvgTurb) ? "true" : "false") << CRD::var;
      // Query to plot low-face coordinates
      ppFILE.query("plot_lo_face_coordinates",lclPlotLoFaceCoordinates);
      CRD::msg << "Plot low-face coordinates\n" <<
        ((lclPlotLoFaceCoordinates) ? "true" : "false") << CRD::var;
      // Query to plot high-face coordinates
      ppFILE.query("plot_hi_face_coordinates",lclPlotHiFaceCoordinates);
      CRD::msg << "Plot high-face coordinates\n" <<
        ((lclPlotHiFaceCoordinates) ? "true" : "false") << CRD::var;
      // Query to restart case with added wall-model
      ppFILE.query("add_wall_model_on_restart",lclRestartAddWallModel);
      CRD::msg << "Restart case with wall-model added\n" <<
        ((lclRestartAddWallModel) ? "true" : "false") << CRD::var;
      // Setup plotting of face-centered, face-averaged, and time-averaged data
      if (ppFILE.contains("plot_low_face_avg") ||
          ppFILE.contains("plot_low_face_avg_time_avg"))
        {
          // First, check face-averaged time-averaged data
          if (ppFILE.contains("plot_low_face_avg_time_avg"))
            {
              const int numFacePlotComps =
                ppFILE.countval("plot_low_face_avg_time_avg");
              for (int comp = 0; comp != numFacePlotComps; ++comp)
                {
                  std::string compName;
                  int status = ppFILE.query(
                    "plot_low_face_avg_time_avg", compName, comp);
                  if (status == 0)
                    {
                      CRD::msg << "Input: Error reading output components for "
                               << "adding to plot_low_face_avg_time_avg"
                               << CRD::error;
                    }
                  if (compName.compare("primitive") == 0)
                    {
                      lclPlotLoFaceAvgTimeAvgComps |= CRDparam::PlotPrimitive;
                    }
                  else if (compName.compare("fluxes") == 0)
                    {
                      lclPlotLoFaceAvgTimeAvgComps |= CRDparam::PlotFluxes;
                    }
                  else if (compName.compare("turb") == 0)
                    {
                      lclPlotLoFaceAvgTimeAvgComps |=
                        CRDparam::PlotTurbulentComps;
                    }
                }
            }
          // Now, for the instantaneous face-averaged data
          if (ppFILE.contains("plot_low_face_avg"))
            {
              const int numFacePlotComps =
                ppFILE.countval("plot_low_face_avg");
              for (int comp = 0; comp != numFacePlotComps; ++comp)
                {
                  std::string compName;
                  int status = ppFILE.query(
                    "plot_low_face_avg", compName, comp);
                  if (status == 0)
                    {
                      CRD::msg << "Input: Error reading output components for "
                               << "adding to plot_low_face_avg" << CRD::error;
                    }
                  if (compName.compare("primitive") == 0)
                    {
                      lclPlotLoFaceAvgComps |= CRDparam::PlotPrimitive;
                    }
                  else if (compName.compare("fluxes") == 0)
                    {
                      lclPlotLoFaceAvgComps |= CRDparam::PlotFluxes;
                    }
                  else if (compName.compare("turb") == 0)
                    {
                      lclPlotLoFaceAvgComps |= CRDparam::PlotTurbulentComps;
                    }
                }
            }
        }
    }
  else
    {
      CRD::msg << "Plot prefix\nnot provided" << CRD::var;
      a_chparam.FILEplotInterval = -1;
    }
  CRD::msg << "Plot interval\n";
  if (a_chparam.FILEplotInterval < 0)
    {
      CRD::msg << "none";
    }
  else
    {
      CRD::msg << a_chparam.FILEplotInterval;
    }
  CRD::msg << CRD::var;
  CRD::msg << "Plot period\n";
  if (a_chparam.FILEplotPeriod < 0.)
    {
      CRD::msg << "none";
    }
  else
    {
      CRD::msg << a_chparam.FILEplotPeriod;
    }
  CRD::msg << CRD::var;
  CRDparam::CRDP.defineFiles(a_chparam.FILEplotPrefix,
                             lclPlotNumGhost,
                             lclPlotExtraVars,
                             lclPlotMappedDerivs,
                             lclPlotWallDist,
                             lclPlotJ,
                             lclPlotJU,
                             lclPlotDACFDCheck,
                             lclPlotFlattening,
                             lclPlotError,
                             lclPlotTimeAvgTurb,
                             lclPlotLoFaceCoordinates,
                             lclPlotHiFaceCoordinates,
                             lclRestartAddWallModel,
                             lclPlotLoFaceAvgTimeAvgComps,
                             lclPlotLoFaceAvgComps);

//--Read the checkpoint file prefix and interval

  a_chparam.FILEcheckpointPrefix = "";
  a_chparam.FILEcheckpointInterval = -1;
  if (ppFILE.contains("checkpoint_prefix"))
    {
      // Query checkpoint prefix
      ppFILE.get("checkpoint_prefix", a_chparam.FILEcheckpointPrefix);
      // Only test existance of path if '/' present
      size_t pos = a_chparam.FILEcheckpointPrefix.rfind('/');
      if (pos != 0 && pos != std::string::npos)
        {
          std::string path = a_chparam.FILEcheckpointPrefix.substr(0, pos);
          if (!CH_System::fileExists(path.c_str()))
            {
              CRD::msg << "Input: 'checkpoint_prefix' path '" << path
                       << "' does not exist!" << CRD::error;
            }
        }
      CRD::msg << "Checkpoint prefix\n" << a_chparam.FILEcheckpointPrefix
               << CRD::var;
      // Query checkpoint interval
      ppFILE.query("checkpoint_interval", a_chparam.FILEcheckpointInterval);
      if (a_chparam.FILEcheckpointInterval < -1)
        {
          CRD::msg << "Input: 'checkpoint_interval' must be >= -1!"
                   << CRD::error;
        }
    }
  else
    {
      CRD::msg << "Checkpoint prefix\nnot provided" << CRD::var;
      a_chparam.FILEcheckpointInterval = -1;
    }
  CRD::msg << "Checkpoint interval\n";
  if (a_chparam.FILEcheckpointInterval == -1)
    {
      CRD::msg << "none";
    }
  else
    {
      CRD::msg << a_chparam.FILEcheckpointInterval;
    }
  CRD::msg << CRD::var;


/*==============================================================================
 * Remaining variables
 *   - An internal simulation defines all remaining variables.  They are
 *     therefore defined here with default values.
 *   - All should be explicitly defined by a simulation (defaults should only
 *     apply to externally defined simulations).
 *   - All can be overridden in the input file.
 *============================================================================*/

  LocalParam lclparam;

//--Simulation paramters

  lclparam.SIMproblemType = CRDparam::ProblemUndefined;

//--Grid parameters

  lclparam.GRIDnumCells.assign(SpaceDim, 0);
  a_chparam.GRIDperiodic.assign(SpaceDim, 1);
  lclparam.GRIDdomainLength.assign(SpaceDim, 1.);
  lclparam.GRIDdomainOrigin.assign(SpaceDim, 0.);
  lclparam.GRIDphysicalLength = 0.;

//--Coordinate system parameters

  a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
  a_chparam.COORDSYSfactory = nullptr;

//--Solver parameters

  a_chparam.SOLfixedDt = 0.0;  // > 0 to fix
  a_chparam.SOLmaxDtGrowth = 1.1;
  a_chparam.SOLdtToleranceFactor = 1.1;
  lclparam.SOLcfl = 1.0;
  lclparam.SOLinitialDt = -1.0;
  lclparam.SOLmaxDt = -1.0;
  a_chparam.SOLinitialCFL = lclparam.SOLcfl;
  a_chparam.SOLinitialTime = 0.0;
  a_chparam.SOLmaxStep = 0;
  a_chparam.SOLmaxTime = 0.0;

//-- ARK parameters

  lclparam.ARKuse = false;
  lclparam.ARKuseAnalyticJac = false;
  lclparam.ARKchemDtScale = 1;
  lclparam.ARKmaxLevel = 0;
  lclparam.ARKinitERKSteps = 0;
  lclparam.ARKpoutLinearSolverStats = false;
  lclparam.ARKNonlinearTol = 1.e-8;
  lclparam.ARKUsePIDControl = false;
  lclparam.ARKPIDControlEps = 5.e-10;
  lclparam.ARKExtrapInitGuess = false;

  //--Limiter parameters

  lclparam.LIMITfaceInterpolationOrder = 4;
  lclparam.LIMITdiffusiveDerivativeOrder = 4;
  lclparam.LIMITreactionOrder = 4;
  lclparam.LIMITusePPMlimiter = 1;
  lclparam.LIMITlimitFaceValues = 1;
  lclparam.LIMITuseFlattening = 1;
  lclparam.LIMITcellConvolveFlatten   = 0;
  lclparam.LIMITcellDeconvolveFlatten = 0;
  lclparam.LIMITfaceConvolveFlatten   = 0;
  lclparam.LIMITfaceDeconvolveFlatten = 0;
  lclparam.LIMITcellConvolveLimit   = 0;
  lclparam.LIMITcellDeconvolveLimit = 0;
  lclparam.LIMITfaceConvolveLimit   = 0;
  lclparam.LIMITfaceDeconvolveLimit = 0;
  a_chparam.LIMITDCFlatTol = 0.2;
  a_chparam.LIMITFCORTol = 1.E-3;
  lclparam.LIMITuseFCOR = 0;
  lclparam.LIMITuseArtVisc = 1;
  lclparam.LIMITartViscCoef = 0.3;
  lclparam.LIMITuseArtVisc4thO = 1;
  lclparam.LIMITartViscCoef4thO = 0.3;
  lclparam.LIMITextraBoundLim = 0;
  lclparam.LIMITnoHOchecks = 0;
  lclparam.LIMITwallCorrections = 1;
  lclparam.LIMITclipping = false;
  lclparam.LIMITclippingHO = false;
  lclparam.LIMITclippingPostSmooth = false;
  lclparam.LIMITfifthOrderBlendingCoef = 0.;

//--AMR parameters

  a_chparam.AMRmaxLevel = 0;
  a_chparam.AMRtagBufferSize = 0;
  a_chparam.AMRbaseLevel = 0;
  a_chparam.AMRmaxGridSize = 32;
  a_chparam.AMRfillRatio = 0.75;
  a_chparam.AMRgridBufferSize = 0;
  a_chparam.AMRblockFactor = 8;
  a_chparam.AMRuseSubcycling = 1;
  a_chparam.AMRregridOnRestart = false;

  a_chparam.TAGfactory = NULL;

//--Fluid parameters (generally freestream ICAO STP)

  lclparam.FLUIDnumStates = SpaceDim + 2;
  lclparam.FLUIDrho    = 1.2250;     // (kg/m^3)
  lclparam.FLUIDT      = 288.15;     // (K)
  lclparam.FLUIDK      = 2.5326E-2;  // (W/m-K)
  lclparam.FLUIDspeed  = 0.;         // (m/s)
  // NACA A278141: nu = 1.4638E-5 (m^2/s)
  // P. Brimblecombe, "Air Composition and Chemistry": nu = 1.4607E-5 (m^2/s)
  lclparam.FLUIDmu     = 1.7894E-5;  // (kg/m s)
  lclparam.FLUIDlambda = 1.0/3.0;    //
  lclparam.FLUIDgamma  = 1.4;        //
  lclparam.FLUIDR      = 287.;       // (J/kg-K)
  lclparam.FLUIDgrav   = 9.80665;    // (m/s^2)
  lclparam.FLUIDRe     = 1.;

//--Combustion parameters

  lclparam.THERMnumSpecies = 0;
  lclparam.THERMreactionStartTime = -1.;
  lclparam.THERMreactionEndTime = 1.E6;
  lclparam.THERMuseSpeciesCorrection = true;
  lclparam.THERMnumReactions = 0;

//--Physics parameters

  lclparam.PHYSICSfluidModels =
    CRDparam::PhysicsInertial | CRDparam::PhysicsViscous;

//--Turbulence parameters

  lclparam.TURBmodelType = CRDparam::TurbModelUndefined;
  lclparam.TURBsgsModelType = CRDparam::SGSModelUndefined;
  lclparam.TURBglobStreamwiseDir = RealVect{1., 0., 0.};
  lclparam.TURBuseConsToPrimCorrection = false;
  lclparam.TURBuseSGSCoarsening = false;
  lclparam.TURBexplicitFilterType = CRDparam::ExplicitFilterNone;
  lclparam.TURBspectralFilterProfile = SpectralFilter::SpectralFilterSharp;
  lclparam.TURBspectralFilterDomainResolution.assign(SpaceDim, 0);
  lclparam.TURBspectralFilterParam = 1.;
  lclparam.TURBuseForcing = false;  // this turns on TaylorGreen forcing only
  lclparam.TURBforcingType = CRDparam::TurbForcingNone;
  // spectral forcing only for now
  lclparam.TURBsgskeCrsLevFilterRatio = 1.;
  lclparam.TURBwallModelCrsRatio = 1;
  lclparam.TURBenforceModeledWallShearStress = false;
  lclparam.TURBstartTimeAvgTime = 0;
  lclparam.TURBspectralForcingInterval = Interval(2, 4);

//--Thread parameters

  a_chparam.THREADuseThreads = false;
  a_chparam.THREADnumThreads = -1;
  a_chparam.THREADthreadsPerTeam = 1;
  a_chparam.THREADextraStackSize_MiB = 0;
  a_chparam.THREADuseHyperThreading = false;

/*--------------------------------------------------------------------*
 * Simulation
 *   - An internal simulation defines all remaining input variables
 *--------------------------------------------------------------------*/

//--Read options that define the simulation

  CRD::msg.newline();

  CRD::msg << "Parameters for 'simulation'" << CRD::h2;
  ParmParse ppSIM("sim");

  // Also might need the coordsys inputs
  ParmParse ppCOORDSYS("coordsys");

//--Type of problem
//**FIXME - these really need to set everything that is not default (e.g.,
//**        numCells)

  bool internalProblem = false;
  bool problemDefined = true;
  {
    std::string printProblemName("Undefined");
    if (ppSIM.contains("problem_type"))
      {
        internalProblem = true;

        std::string inputProblemName;
        ppSIM.get("problem_type", inputProblemName);
        if (inputProblemName == "External")
          {
            lclparam.SIMproblemType = CRDparam::ProblemExternal;
            printProblemName = "External";
            internalProblem = false;
          }
        else if (inputProblemName == "EulerAdvectionCube")
          {
            lclparam.SIMproblemType = CRDparam::ProblemEulerAdvectionCube;
            printProblemName = "Euler advection in cube";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
            lclparam.FLUIDrho = 1.4;
            lclparam.FLUIDT   = 1./(lclparam.FLUIDrho*lclparam.FLUIDR);
          }
        else if (inputProblemName == "EulerAdvectionCylinder")
          {
            lclparam.SIMproblemType = CRDparam::ProblemEulerAdvectionCylinder;
            printProblemName = "Euler advection in cylinder";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
          }
        else if (inputProblemName == "EulerShockTube")
          {
            lclparam.SIMproblemType = CRDparam::ProblemShockTube;
            printProblemName = "Euler Sod's shock tube";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
            //**FIXME periodicity should be read from the ibc files otherwise
            //**      there could be a mismatch between problem domain and ibc
            //**FIXME this should be changed to all non-periodic
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "NavierStokesShockTube")
          {
            lclparam.SIMproblemType = CRDparam::ProblemShockTube;
            printProblemName = "Navier-Stokes Sod's shock tube";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
            CRD::msg << "Input: 'problem_type' of NavierStokesShockTube not yet"
              "supported!" << CRD::error;
          }
        else if (inputProblemName == "EulerShockBox")
          {
            lclparam.SIMproblemType = CRDparam::ProblemEulerShockBox;
            printProblemName = "Euler shock box";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
            lclparam.FLUIDrho = 1.4;
            lclparam.FLUIDT   = 1./(lclparam.FLUIDrho*lclparam.FLUIDR);
          }
        else if (inputProblemName == "MachReflection")
          {
            lclparam.SIMproblemType = CRDparam::ProblemMachReflection;
            printProblemName = "Mach reflection on a ramp";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            // Check if stretching is used
            if (ppCOORDSYS.contains("shift"))
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockLogSchwarzChristoffelRamp;
              }
            else
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp;
              }
            lclparam.FLUIDrho = 1.4;
            lclparam.FLUIDT   = 1./(lclparam.FLUIDrho*lclparam.FLUIDR);
          }
        else if (inputProblemName == "SpecMachReflection")
          {
            lclparam.SIMproblemType = CRDparam::ProblemSpecMachReflection;
            printProblemName = "Mach reflection on a ramp with multi-species";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype =
              CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp;
            lclparam.FLUIDrho = 1.4;
            lclparam.FLUIDT   = 1./(lclparam.FLUIDrho*lclparam.FLUIDR);
          }
        else if (inputProblemName == "NavierStokesTransientCouette")
          {
            lclparam.SIMproblemType =
              CRDparam::ProblemNavierStokesTransientCouette;
            printProblemName = "Transient Couette flow";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[1] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "NavierStokesTransientFlatplate")
          {
            lclparam.SIMproblemType =
              CRDparam::ProblemNavierStokesTransientFlatplate;
            printProblemName = "Boundary layer growth on a flat plate";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            //**FIXME periodicity should be read from the ibc files otherwise
            //**      there could be a mismatch between problem domain and ibc
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            if (SpaceDim == 3)
              {
                a_chparam.GRIDperiodic[SpaceDim-1] = 1;
              }
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "NavierStokesTransientPoiseuille")
          {
            lclparam.SIMproblemType =
              CRDparam::ProblemNavierStokesTransientPoiseuille;
            printProblemName = "Transient Poiseuille flow";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsSource;
            //**FIXME periodicity should be read from the ibc files otherwise
            //**      there could be a mismatch between problem domain and ibc
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[SpaceDim - 1] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        //**FIXME This is only to test the inputs for combustion
        else if (inputProblemName == "CombustionTest")
          {
            lclparam.SIMproblemType = CRDparam::ProblemCombustion;
            printProblemName = "Combustion input test";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsSource |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.GRIDperiodic[0] = 1;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "LidDrivenCavity")
          {
            lclparam.SIMproblemType = CRDparam::ProblemLidDrivenCavity;
            printProblemName = "Lid driven cavity case";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemMixedCouette")
          {
            lclparam.SIMproblemType = CRDparam::ProblemMixedCouette;
            printProblemName = "Mixed Couette case";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.GRIDperiodic[0] = 1;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemJetInlet")
          {
            lclparam.SIMproblemType = CRDparam::ProblemJetInlet;
            printProblemName = "Jet inlet case";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemSpecShock")
          {
            lclparam.SIMproblemType = CRDparam::ProblemSpecShock;
            printProblemName = "Shock tube with species transport";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemSpecShockBox")
          {
            lclparam.SIMproblemType = CRDparam::ProblemSpecShockBox;
            printProblemName = "Shock box with species transport";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemShockBubble")
          {
            lclparam.SIMproblemType = CRDparam::ProblemShockBubble;
            printProblemName = "ProblemShockBubble";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemRMI")
          {
            lclparam.SIMproblemType = CRDparam::ProblemRMI;
            printProblemName = "ProblemRMI";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemBurnerTest")
          {
            lclparam.SIMproblemType = CRDparam::ProblemBurnerTest;
            printProblemName = "ProblemBurnerTest";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemReactionAdvection")
          {
            lclparam.SIMproblemType = CRDparam::ProblemReactionAdvection;
            printProblemName = "ProblemReactionAdvection";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemVortex")
          {
            lclparam.SIMproblemType = CRDparam::ProblemVortex;
            printProblemName = "Vortex test";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsThermPerf |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "WedgeDetonation")
          {
            lclparam.SIMproblemType = CRDparam::ProblemWedgeDetonation;
            printProblemName = "Wedge detonation";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsThermPerf |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.GRIDperiodic[1] = 0;
            // Check if stretching is used
            if (ppCOORDSYS.contains("shift"))
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockLogSchwarzChristoffelRamp;
              }
            else
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp;
              }
          }
        else if (inputProblemName == "ProblemChannelFlow")
          {
            lclparam.SIMproblemType = CRDparam::ProblemChannelFlow;
            printProblemName = "Channel flow problem";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            // a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ShuOsher")
          {
            lclparam.SIMproblemType = CRDparam::ProblemShuOsher;
            printProblemName = "Shu-Osher problem";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "CylinderFlow")
          {
            lclparam.SIMproblemType = CRDparam::ProblemCylinderFlow;
            printProblemName = "Flow over a cylinder";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemBluffBodyCombustion")
          {
            lclparam.SIMproblemType = CRDparam::ProblemBluffBodyCombustion;
            printProblemName = "ProblemBluffBodyCombustion";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous |
              CRDparam::PhysicsThermPerf;
            a_chparam.GRIDperiodic.assign(SpaceDim, 1);
            a_chparam.GRIDperiodic[0] = 0;
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemSphere")
          {
            lclparam.SIMproblemType = CRDparam::ProblemSphere;
            printProblemName = "Flow over a sphere";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        else if (inputProblemName == "ProblemShockShock")
          {
            lclparam.SIMproblemType = CRDparam::ProblemShockShock;
            printProblemName = "Shock-shock interaction";
            lclparam.PHYSICSfluidModels =
              CRDparam::PhysicsInertial |
              CRDparam::PhysicsViscous;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            a_chparam.COORDSYStype =
              CRDparam::CoordSysSingleBlockJoukowskiAirfoil;
          }
        else if (inputProblemName == "ObliqueWave")
          {
            lclparam.SIMproblemType = CRDparam::ProblemObliqueWave;
            printProblemName = "Oblique shock or expansion wave";
            lclparam.PHYSICSfluidModels = CRDparam::PhysicsInertial;
            a_chparam.GRIDperiodic.assign(SpaceDim, 0);
            // Check if stretching is used
            if (ppCOORDSYS.contains("shift"))
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockLogSchwarzChristoffelRamp;
              }
            else
              {
                a_chparam.COORDSYStype =
                  CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp;
              }
          }
        else // For any problem that inherits CNSIBCGeneralized
          {
            // Bool for whether periodicity must be specified from input file
            bool mustReadPeriod = false;
            if (inputProblemName == "ProblemGeneralized")
              {
                lclparam.SIMproblemType = CRDparam::ProblemGeneralized;
                printProblemName = "Generalized problem";
                mustReadPeriod = true;
              }
            else if (inputProblemName == "ProblemGeneralizedSingleBlock")
              {
                lclparam.SIMproblemType =
                  CRDparam::ProblemGeneralizedSingleBlock;
                printProblemName = "Generalized single block problem";
                mustReadPeriod = true;
              }
            else if (inputProblemName == "ProblemTurbulentCouette")
              {
                lclparam.SIMproblemType = CRDparam::ProblemTurbulentCouette;
                printProblemName = "Turbulent Couette problem";
              }
            else if (inputProblemName == "ProblemFlame")
              {
                lclparam.SIMproblemType = CRDparam::ProblemFlame;
                printProblemName = "General flame problem";
                mustReadPeriod = true;
              }
            else if (inputProblemName == "ProblemTaylorGreen")
              {
                lclparam.SIMproblemType = CRDparam::ProblemTaylorGreen;
                printProblemName = "Taylor-Green vortex problem";
              }
            else if (inputProblemName == "ProblemIsotropicTurbulence")
              {
                lclparam.SIMproblemType = CRDparam::ProblemIsotropicTurbulence;
                printProblemName = "Homogeneous isotropic turbulence problem";
              }
            else if (inputProblemName == "ProblemMixingLayer")
              {
                lclparam.SIMproblemType = CRDparam::ProblemMixingLayer;
                printProblemName = "Mixing layer problem";
              }
            else if (inputProblemName == "ProblemRiemannCube")
              {
                lclparam.SIMproblemType = CRDparam::ProblemRiemannCube;
                printProblemName = "3D Riemann problem";
              }
            else if (inputProblemName == "ProblemDetonation")
              {
                lclparam.SIMproblemType = CRDparam::ProblemDetonation;
                printProblemName = "ProblemDetonation";
                mustReadPeriod = true;
              }
            else if (inputProblemName == "ProblemGaussPulse")
              {
                lclparam.SIMproblemType = CRDparam::ProblemGaussPulse;
                printProblemName = "ProblemGaussPulse";
              }
            else if (inputProblemName == "ProblemShear")
              {
                lclparam.SIMproblemType = CRDparam::ProblemShear;
                printProblemName = "ProblemShear";
              }
            else if (inputProblemName == "ProblemVortex")
              {
                lclparam.SIMproblemType = CRDparam::ProblemVortex;
                printProblemName = "Vortex test";
              }
            else if (inputProblemName == "ProblemSpatiallyEvolvingShear")
              {
                lclparam.SIMproblemType =
                  CRDparam::ProblemSpatiallyEvolvingShear;
                printProblemName = "ProblemSpatiallyEvolvingShear";
              }
            else if (inputProblemName == "ProblemMMS")
              {
                lclparam.SIMproblemType = CRDparam::ProblemMMS;
                printProblemName = "Manufactured solution test";
              }
            else if (inputProblemName == "ProblemRecirculatingInletTFP")
              {
                lclparam.SIMproblemType =
                  CRDparam::ProblemRecirculatingInletTFP;
                printProblemName =
                  "Recirculating inlet flat-plate turbulent-boundary-layer";
              }
            else if (inputProblemName == "ProblemTemperatureDiffusion")
              {
                lclparam.SIMproblemType =
                  CRDparam::ProblemTemperatureDiffusion;
                printProblemName =
                  "Temperature diffusion test";
              }
            else if (inputProblemName == "ProblemChannelSeparation")
              {
                lclparam.SIMproblemType = CRDparam::ProblemChannelSeparation;
                printProblemName = "Smooth-ramp channel separation";
              }
            else
              {
                problemDefined = false;
              }
            ParmParse ppTHERM("therm");
            if (ppTHERM.contains("species"))
              {
                lclparam.PHYSICSfluidModels =
                  CRDparam::PhysicsInertial |
                  CRDparam::PhysicsThermPerf |
                  CRDparam::PhysicsViscous;
              }
            else
              {
                lclparam.PHYSICSfluidModels =
                  CRDparam::PhysicsInertial |
                  CRDparam::PhysicsViscous;
              }
            ParmParse ppIBC("ibc");
            // Problems are periodic by default
            std::vector<int> periodicity(SpaceDim, 1);
            if (mustReadPeriod)
              {
                ppIBC.getarr("periodicity", periodicity, 0, SpaceDim);
              }
            else
              {
                ppIBC.queryarr("periodicity", periodicity, 0, SpaceDim);
              }
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                a_chparam.GRIDperiodic[dir] = periodicity[dir];
              }
            a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockCartesian;
          }
        if (!problemDefined)
          {
            CRD::msg << "Input: 'problem_type' of " << inputProblemName
                     << " not recognized!" << CRD::error;
            printProblemName = "Unknown";
            internalProblem = false;
          }
      }
    CRD::msg << "Problem type\n" << printProblemName << CRD::var;
    // Define of simulation is deferred until grid parameters are read
  }

/*--------------------------------------------------------------------*
 * Grid parameters
 *--------------------------------------------------------------------*/

//--Read options that define the grid

  CRD::msg.newline();

  CRD::msg << "Parameters for 'grid'" << CRD::h2;
  ParmParse ppGRID("grid");

//--Check if reading from file, all other grid and CS reads will be ignored

  if (ppGRID.contains("file"))
    {
      a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockExternal;
      std::string grid_file;
      // Query file name
      ppGRID.query("file", grid_file);
      // Test existance of file
      if (!CH_System::fileExists(grid_file.c_str()))
        {
          CRD::msg << "Input: grid file '" << grid_file
                   << "' not found!" << CRD::error;
        }
      // Get the file extension
      std::string gridf_ext = grid_file.substr(grid_file.find_last_of(".") + 1);

      // Get the grid scaling
      std::vector<Real> scaleVec(SpaceDim, 0.);
      RealVect scale(RealVect::Unit);
      if (ppGRID.contains("scale"))
        {
          ppGRID.getarr("scale", scaleVec, 0, SpaceDim);
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              scale[dir] = scaleVec[dir];
            }
        }
      
      // Reading a CGNS file
      if (gridf_ext == "cgns" || gridf_ext == "CGNS" || gridf_ext == "Cgns")
        {
#ifdef CH_USE_CGNS
          // Open the CGNS file, read the number of zones, and close
          int cgerr;
          int idxFile;
#ifdef CH_MPI
          cgerr = cgp_open(grid_file.c_str(), CG_MODE_READ, &idxFile);
#else
          cgerr =  cg_open(grid_file.c_str(), CG_MODE_READ, &idxFile);
#endif
          CH_assert(cgerr == 0);
          int idxBase = 1;
          int numZones;
          cgerr = cg_nzones(idxFile, idxBase, &numZones);
          CH_assert(cgerr == 0);
#ifdef CH_MPI
          cgerr = cgp_close(idxFile);
#else
          cgerr =  cg_close(idxFile);
#endif
          CH_assert(cgerr == 0);
          std::string COORDSYSname;
          ppCOORDSYS.query("type", COORDSYSname);
          if (COORDSYSname == "ChannelSeparation_MMB")
          {
            RealVect scale(RealVect::Unit);
            // Dangerous, but provides very nice results if carefully done
            int blocksAreConformal = 0;
            ppGRID.query("is_conformal", blocksAreConformal);
            // Run only the inlet block
            int numBlocksChannel = 3;
            int onlyInlet = 0;
            ppGRID.query("only_simulate_channel_inlet", onlyInlet);
            if (onlyInlet)
              {
                numBlocksChannel = 1;
              }
            // these _dx are used to create the
            // spacing vector in the left, right, and z direction
            Real z_dx;
            Real l_dx;
            Real r_dx;
            ppGRID.query("z_dx", z_dx);
            ppGRID.query("l_dx", l_dx);
            ppGRID.query("r_dx", r_dx);
            // these _lengths are used to create the
            // total number of cells in the left, right, and z directions
            int z_length;
            int l_length;
            int r_length;
            ppGRID.query("z_length", z_length);
            ppGRID.query("l_length", l_length);
            ppGRID.query("r_length", r_length);

            // }
            CRD::msg << "Coordinate system\nChannel Separation" << CRD::var;
            //a_chparam.COORDSYStype = CRDparam::CoordSysMultiBlockCartesian;
            a_chparam.COORDSYStype = CRDparam::CoordSysExtrudedMultiBlock;
            a_chparam.COORDSYSfactory =

              new ChannelSeparationCSFactory(grid_file,
                                             scale,
                                             CRDparam::g_verbosity,
                                             z_dx,
                                             l_dx,
                                             r_dx,
                                             z_length,
                                             l_length,
                                             r_length,
                                             blocksAreConformal,
                                             numBlocksChannel);
          }

           else if (numZones == 1) // single block grid
            {
              a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockExternal;
              a_chparam.COORDSYSfactory =
                new SingleBlockCSAdaptorFactory(
                  new ExternalCSFactory(grid_file, scale));
            }
          else if (numZones >= 1) // multi-block grid
            {
              a_chparam.COORDSYStype = CRDparam::CoordSysMultiBlockExternal;
              a_chparam.COORDSYSfactory =
                new ExternalMultiBlockCSFactory(grid_file,
                                                scale,
                                                CRDparam::g_verbosity);
            }

          else
            {
              CRD::msg << "No zones in the first base of CGNS file"
                       << CRD::error;
            }
#else
          CRD::msg << "Input: CGNS file reading not enabled in this build! Compile with USE_CGNS=TRUE and the appropriate library\n" << CRD::error;
#endif  /* CH_USE_CGNS */
        }
      
      // Reading an Ogen file
      else if (gridf_ext == "hdf")
        {
#ifdef CH_USE_OGEN
          readOgen grid_data(grid_file);
          grid_data.getDim(lclparam.GRIDnumCells,
                           lclparam.GRIDdomainLength,
                           lclparam.GRIDdomainOrigin);
          
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new ExternalCSFactory(grid_file, scale));
#else
          CRD::msg << "Input: Ogen file reading not enabled in this build! Compile with USE_OGEN=TRUE and the appropriate libraries\n" << CRD::error;
#endif
        }
      
      // Did not recognize the grid file extension
      else
        {
          CRD::msg << "Input: Grid file extension not recognized!"
                   << CRD::error;
        }

      //--External Grid Scaling
      CRD::msg.setPrecFloatSN(l_prec);
      CRD::msg << "Physical Grid Scaling\n("<< scale[0];
      for (int dir = 1; dir != SpaceDim; ++dir)
        {
          CRD::msg << ',' << scale[dir];
        }
      CRD::msg << ')' << CRD::var;
      CRD::msg.setFloatDefault();
    }
  
//--Internal analytic mapped grid
  
//--Number of cells on the coarsest level

  if (ppGRID.contains("num_cells"))
    {
      ppGRID.getarr("num_cells", lclparam.GRIDnumCells, 0, SpaceDim);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (lclparam.GRIDnumCells[dir] < 5)  //**FIXME
            {
              CRD::msg << "Input: 'num_cells' must be >= 5 in direction " << dir
                       << '!' << CRD::error;
            }
        }
    }
  else if (!internalProblem)
    {
      CRD::msg << "Input: 'num_cells' must be specified if internal problem "
        "type not defined!" << CRD::error;
    }
  CRD::msg << "Number of cells\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << lclparam.GRIDnumCells[dir];
    }
  CRD::msg << ')' << CRD::var;

//--Domain length

  // This is the length in computational space!
  if (ppGRID.contains("domain_length"))
    {
      ppGRID.getarr("domain_length", lclparam.GRIDdomainLength, 0, SpaceDim);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          lclparam.GRIDphysicalLength[dir] = lclparam.GRIDdomainLength[dir];
          if (lclparam.GRIDdomainLength[dir] <= 0.)
            {
              CRD::msg << "Input: 'domain_length' must be > 0. in direction "
                       << dir << '!' << CRD::error;
            }
        }
    }
  else if (!internalProblem)
    {
      CRD::msg << "Input: 'domain_length' must be specified if internal "
        "problem type not defined!" << CRD::error;
    }
  CRD::msg.setPrecFloatSN(l_prec);
  CRD::msg << "Computational domain length\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << lclparam.GRIDdomainLength[dir];
    }
  CRD::msg << ')' << CRD::var;
  CRD::msg.setFloatDefault();

//--Domain origin (location of node (0,0,0))

  ppGRID.queryarr("domain_origin", lclparam.GRIDdomainOrigin, 0, SpaceDim);
  CRD::msg.setPrecFloatSN(l_prec);
  CRD::msg << "Computational domain origin\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << lclparam.GRIDdomainOrigin[dir];
    }
  CRD::msg << ')' << CRD::var;
  CRD::msg.setFloatDefault();

//--Periodicity

  if (ppGRID.contains("domain_periodic"))
    {
      std::vector<std::string> domainPeriodic(SpaceDim, "no");
      ppGRID.getarr("domain_periodic", domainPeriodic, 0, SpaceDim);
      stc::forEachElement<SpaceDim>([&]
                                    (const stc::array_size_type idx)
                                    {
                                      a_chparam.GRIDperiodic[idx] =
                                        (domainPeriodic[idx] == "yes");
                                    });
    }
  CRD::msg << "Domain periodicity\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << ((a_chparam.GRIDperiodic[dir] == 1) ? "periodic" : "none");
    }
  CRD::msg << ')' << CRD::var;

  Real ratioTest = 0.;
  D_TERM(,
         ratioTest += lclparam.GRIDdomainLength[0]/lclparam.GRIDdomainLength[1];
         ratioTest -= Real(lclparam.GRIDnumCells[0])/
         Real(lclparam.GRIDnumCells[1]);,
         ratioTest += lclparam.GRIDdomainLength[0]/lclparam.GRIDdomainLength[2];
         ratioTest -= Real(lclparam.GRIDnumCells[0])
         /Real(lclparam.GRIDnumCells[2]););
  if(std::abs(ratioTest) > 5.E-15)
    {
      CRD::msg << "Ratio of lengths to number of cells do not correspond! "
               << std::abs(ratioTest) << " ~ 0" << CRD::error;
    }

  CRDparam::CRDP.defineSimulation(lclparam.SIMproblemType,
                                  a_chparam.GRIDperiodic);
  CRDparam::CRDP.defineGeometry(lclparam.GRIDnumCells,
                                lclparam.GRIDdomainLength,
                                lclparam.GRIDdomainOrigin);

/*--------------------------------------------------------------------*
 * Coordinate system parameters
 * In some cases, the user can override the default selection of a
 * coordinate system for a given problem.  This is normally done
 * for experimental development work.
 *--------------------------------------------------------------------*/

//--Read options that define the coordinate system

  CRD::msg.newline();

  CRD::msg << "Parameters for 'coordsys'" << CRD::h2;
  // Coordsys already established above
  // ParmParse ppCOORDSYS("coordsys");

//--If the coordinate system is single-block Cartesian, allow the user to
//--modify it to some mapped variant in a hyper-cube.

  if (ppGRID.contains("file"))
    {
      // Output if CS comes from external file
      std::string mesh_file;
      ppGRID.query("file", mesh_file);
      
      CRD::msg << "Coordinate system\nfile '" 
               << mesh_file << "'" << CRD::var;
    }
  else if (a_chparam.COORDSYStype == CRDparam::CoordSysSingleBlockCartesian)
    {
      std::string COORDSYSname("cartesian");
      ppCOORDSYS.query("type", COORDSYSname);
      if (COORDSYSname == "cartesian")
        {
          std::vector<Real> COORDSYSoriginVec(SpaceDim, 0.);
          std::vector<Real> COORDSYSstretchVec(SpaceDim, 1.);
          ppCOORDSYS.queryarr("origin", COORDSYSoriginVec, 0, SpaceDim);
          RealVect COORDSYSorigin;     // Domain origin
          RealVect COORDSYSstretch;    // Stretching vector
          RealVect COORDSYSphysLength; // Length after stretching
          if (ppCOORDSYS.contains("physical_length"))
            {
              std::vector<Real> COORDSYSphysLengthVec(SpaceDim, 0.);
              CRD::msg.setPrecFloatSN(l_prec);
              ppCOORDSYS.getarr("physical_length",
                                COORDSYSphysLengthVec, 0, SpaceDim);
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  COORDSYSorigin[dir] = COORDSYSoriginVec[dir];
                  COORDSYSphysLength[dir] = COORDSYSphysLengthVec[dir];
                  COORDSYSstretch[dir] = COORDSYSphysLength[dir]/
                    lclparam.GRIDdomainLength[dir];
                }
              CRD::msg << ')' << CRD::var;
              CRD::msg.setFloatDefault();
            }
          else
            {
              ppCOORDSYS.queryarr("stretch", COORDSYSstretchVec, 0,
                                  SpaceDim);
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  COORDSYSorigin[dir] = COORDSYSoriginVec[dir];
                  COORDSYSstretch[dir] = COORDSYSstretchVec[dir];
                  COORDSYSphysLength[dir] =
                    lclparam.GRIDdomainLength[dir]*COORDSYSstretch[dir];
                }
            }
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new CartesianCSFactory(COORDSYSorigin, COORDSYSstretch),
              COORDSYSphysLength);

          CRD::msg << "Coordinate system\nCartesian" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Origin\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSorigin[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg << "  Stretch\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSstretch[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg << "  Physical length\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              lclparam.GRIDphysicalLength[dir] = COORDSYSphysLength[dir];
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSphysLength[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg.setFloatDefault();
        }
      else if (COORDSYSname == "warped")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockWarped;
          // Defaults for options
          std::vector<Real> COORDSYSscaleVec(SpaceDim, 0.1);
          Real COORDSYSrtol(RTOL);  // From WarpedCS.H
          Real COORDSYSatol(ATOL);  // From WarpedCS.H
          Real COORDSYSmaxIt(100);
          ppCOORDSYS.queryarr("scale", COORDSYSscaleVec, 0, SpaceDim);
          ppCOORDSYS.query("relative_tolerance", COORDSYSrtol);
          ppCOORDSYS.query("absolute_tolerance", COORDSYSatol);
          ppCOORDSYS.query("maximum_iterations", COORDSYSmaxIt);
          RealVect COORDSYSscale;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              COORDSYSscale[dir]  = COORDSYSscaleVec[dir];
            }
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new WarpedCSFactory(COORDSYSscale,
                                  COORDSYSrtol,
                                  COORDSYSatol,
                                  COORDSYSmaxIt),
              lclparam.GRIDphysicalLength);
          CRD::msg << "Coordinate system\nwarped" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Scale\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSscale[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg.setFloatDefault();
          CRD::msg << "  Relative tolerance\n" << COORDSYSrtol << CRD::var;
          CRD::msg << "  Absolute tolerance\n" << COORDSYSatol << CRD::var;
          CRD::msg << "  Maximum iterations\n" << COORDSYSmaxIt << CRD::var;
        }
      else if (COORDSYSname == "twisted")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockTwisted;
          // Defaults for options
          Real COORDSYSradius = 0.5;
          Real COORDSYStwist  = 0.5;
          ppCOORDSYS.query("radius", COORDSYSradius);
          ppCOORDSYS.query("twist", COORDSYStwist);
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new TwistedCSFactory(COORDSYSradius, COORDSYStwist),
              lclparam.GRIDphysicalLength);
          CRD::msg << "Coordinate system\ntwisted" << CRD::var;
          CRD::msg << "  Radius\n" << COORDSYSradius << CRD::var;
          CRD::msg << "  Twist\n" << COORDSYStwist << CRD::var;
        }
      else if (COORDSYSname == "logstretch")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockLogStretch;
          // Default for options
          std::vector<Real> COORDSYSshiftVec(SpaceDim, 0.);
          ppCOORDSYS.queryarr("shift", COORDSYSshiftVec, 0, SpaceDim);
          RealVect COORDSYSshift;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              COORDSYSshift[dir] = COORDSYSshiftVec[dir];
            }
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new LogCSFactory(COORDSYSshift),
              lclparam.GRIDphysicalLength);
          CRD::msg << "Coordinate system\nlog-stretch" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Shift\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSshift[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg.setFloatDefault();
        }
      else if (COORDSYSname == "annulus")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockAnnulus;
          // Default options
          Real COORDSYSinnerradius = 0.1;
          Real COORDSYSouterradius = 1.0;
          Real COORDSYSradstretch = 1.0;
          ppCOORDSYS.query("inner_radius", COORDSYSinnerradius);
          ppCOORDSYS.query("outer_radius", COORDSYSouterradius);
          ppCOORDSYS.query("radial_stretch", COORDSYSradstretch);
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new AnnulusCSFactory(COORDSYSinnerradius,
                                   COORDSYSouterradius,
                                   COORDSYSradstretch),
              lclparam.GRIDphysicalLength);
          CRD::msg << "Coordinate system\nannulus" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Inside radius\n" << COORDSYSinnerradius
                   << CRD::var;
          CRD::msg << "  Outside radius\n" << COORDSYSouterradius
                   << CRD::var;
          CRD::msg << "  Radial stretch\n" << COORDSYSradstretch
                   << CRD::var;
          CRD::msg.setFloatDefault();
        }
      else if (COORDSYSname == "skewed")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysSingleBlockSkewed;
          // Default options
          Real COORDSYSangle = 60;
          std::vector<Real> COORDSYSphysicalLengthVec(SpaceDim, 1.0);
          ppCOORDSYS.query("angle", COORDSYSangle);
          ppCOORDSYS.queryarr("physical_length", COORDSYSphysicalLengthVec, 0,
                              SpaceDim);
          RealVect COORDSYSphysicalLength;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              COORDSYSphysicalLength[dir] = COORDSYSphysicalLengthVec[dir];
            }
          a_chparam.COORDSYSfactory =
            new SingleBlockCSAdaptorFactory(
              new SkewCSFactory(COORDSYSangle,
                                COORDSYSphysicalLength),
              lclparam.GRIDphysicalLength);
          CRD::msg << "Coordinate system\nskew" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Angle (degrees)\n" << COORDSYSangle
                   << CRD::var;
          CRD::msg << "  Physical Length\n(";
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (dir != 0)
                {
                  CRD::msg << ',';
                }
              CRD::msg << COORDSYSphysicalLength[dir];
            }
          CRD::msg << ')' << CRD::var;
          CRD::msg.setFloatDefault();
        }

//**FIXME Some of these are not really Cartesian anymore.  They shouldn't be
//**      under a_chparam.COORDSYStype == CRDparam::CoordSysSingleBlockCartesian
//**      And when not under Cartesian, coordsys.type should be ignored from
//**      the input file.  See Joukowski, where the computational grid is still
//**      a single block, but the connectivity means it really should be treated
//**      differently.
//**      STOP PUTTING WEIRD GRIDS UNDER CARTESIAN!

      else if (COORDSYSname == "cartesian_MMB")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysMultiBlockCartesian;
          a_chparam.COORDSYSfactory = new DoubleCartesianCSFactory{true};
          CRD::msg << "Coordinate system\nmulti-block Cartesian" << CRD::var;
          // This requires a hyper-cube domain
          if (!stc::andEachElement<SpaceDim>(
                [size0 = CRDparam::g_domainBaseSize[0]]
                (const stc::array_size_type a_dir)
                { return CRDparam::g_domainBaseSize[a_dir] == size0; }))
            {
              CRD::msg << "Multiblock Cartesian domain must have equal "
                "dimensions" << CRD::error;
            }
          if ((CRDparam::g_domainBaseSize[0] % 2) != 0)
            {
              CRD::msg << "Multiblock Cartesian domain size must be divisible "
                "by 2" << CRD::error;
            }
        }
      else if (COORDSYSname == "sphere_outer")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysMultiBlockSphereShell;
          // Default options
          Real COORDSYSinnerradius = 0.1;
          Real COORDSYSouterradius = 1.0;
          Real COORDSYSradstretch = 1.0;
          Real COORDSYSalpha = -1.0;
          Real COORDSYSbeta = 1.0;
          ppCOORDSYS.query("inner_radius", COORDSYSinnerradius);
          ppCOORDSYS.query("outer_radius", COORDSYSouterradius);
          ppCOORDSYS.query("alpha", COORDSYSalpha);
          ppCOORDSYS.query("beta", COORDSYSbeta);
          a_chparam.COORDSYSfactory =
            new CubedSphereShell3DCSFactory(COORDSYSinnerradius,
                                            COORDSYSouterradius,
                                            COORDSYSalpha,
                                            COORDSYSbeta);
          CRD::msg << "Coordinate system\nmulti-block sphere" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Inside radius\n" << COORDSYSinnerradius
                   << CRD::var;
          CRD::msg << "  Outside radius\n" << COORDSYSouterradius
                   << CRD::var;
          CRD::msg << "  Radial stretch\n" << COORDSYSradstretch
                   << CRD::var;
          CRD::msg.setFloatDefault();
        }
      else if (COORDSYSname == "smooth_ramp_MMB")
        {
          a_chparam.COORDSYStype = CRDparam::CoordSysMultiBlockSmoothRamp;
          // Default options
          Real COORDSYSstretch = 0.;
          Real COORDSYSrampHeightRelax = 0.;
          Real COORDSYSzWidth = 0.128;
          int COORDSYSonlyInletBlock = 0;
          int COORDSYSnumXCellsInletBlock = 1;
          int COORDSYSnumXCellsRampBlockPreRamp = 1;
          int COORDSYSnumXCellsRampBlockRamp = 1;
          int COORDSYSnumXCellsRampBlockPostRamp = 1;
          int COORDSYSnumXCellsOutletBlock = 1;
          ppCOORDSYS.query("wall_normal_stretch", COORDSYSstretch);
          ppCOORDSYS.query("wall_normal_ramp_height_relaxation",
                           COORDSYSrampHeightRelax);
          ppCOORDSYS.query("spanwise_domain_width", COORDSYSzWidth);
          ppCOORDSYS.query("only_simulate_channel_inlet",
                           COORDSYSonlyInletBlock);
          ppCOORDSYS.query("num_x_cells_inlet_block",
                           COORDSYSnumXCellsInletBlock);
          ppCOORDSYS.query("num_x_cells_ramp_block_before_ramp",
                           COORDSYSnumXCellsRampBlockPreRamp);
          ppCOORDSYS.query("num_x_cells_total_ramp_block",
                           COORDSYSnumXCellsRampBlockRamp);
          ppCOORDSYS.query("num_x_cells_ramp_block_after_ramp",
                           COORDSYSnumXCellsRampBlockPostRamp);
          ppCOORDSYS.query("num_x_cells_outlet_block",
                           COORDSYSnumXCellsOutletBlock);
          a_chparam.COORDSYSfactory =
            new SmoothRampCSFactory(COORDSYSstretch,
                                    COORDSYSrampHeightRelax,
                                    COORDSYSzWidth,
                                    COORDSYSonlyInletBlock,
                                    COORDSYSnumXCellsInletBlock,
                                    COORDSYSnumXCellsRampBlockPreRamp,
                                    COORDSYSnumXCellsRampBlockRamp,
                                    COORDSYSnumXCellsRampBlockPostRamp,
                                    COORDSYSnumXCellsOutletBlock);
          CRD::msg << "Coordinate system\nmulti-block smooth-ramp" << CRD::var;
          CRD::msg.setPrecFloatSN(l_prec);
          CRD::msg << "  Wall-normal stretch\n" << COORDSYSstretch << CRD::var;
          CRD::msg << "  Wall-normal relaxation of ramp-height\n"
                   << COORDSYSrampHeightRelax << CRD::var;
          CRD::msg << "  Spanwise domain-width\n" << COORDSYSzWidth << CRD::var;
          CRD::msg << "  Only inlet block\n" << COORDSYSonlyInletBlock
                   << CRD::var;
          CRD::msg << "  Number of x-cells in inlet block\n"
                   << COORDSYSnumXCellsInletBlock << CRD::var;
          CRD::msg << "  Number of x-cells in ramp block before ramp\n"
                   << COORDSYSnumXCellsRampBlockPreRamp << CRD::var;
          CRD::msg << "  Number of x-cells in total ramp block\n"
                   << COORDSYSnumXCellsRampBlockRamp << CRD::var;
          CRD::msg << "  Number of x-cells in ramp block after ramp\n"
                   << COORDSYSnumXCellsRampBlockRamp << CRD::var;
          CRD::msg << "  Number of x-cells in outlet block\n"
                   << COORDSYSnumXCellsOutletBlock << CRD::var;
          CRD::msg.setFloatDefault();
        }
      else
        {
          CRD::msg << "Input: 'type' unrecognized!" << CRD::error;
        }
    }
  else if (a_chparam.COORDSYStype ==
           CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp)
    {
      Real alpha, leadLength, rampLength;
      rampParameters(ppCOORDSYS, alpha, leadLength, rampLength);
      Real length = leadLength + rampLength;
      Real height =
        lclparam.GRIDdomainLength[1]/lclparam.GRIDdomainLength[0]*length;
      lclparam.GRIDphysicalLength[0] = length;
      lclparam.GRIDphysicalLength[1] = height;
      CRD::msg << "Coordinate system\nSchwarz-Christoffel ramp" << CRD::var;
      CRD::msg.setPrecFloatSN(l_prec);
      CRD::msg << "  Angle of ramp (deg)\n" << alpha << CRD::var;
      CRD::msg << "  Lead length before ramp\n" << leadLength << CRD::var;
      CRD::msg << "  Ramp length (projected along x-axis)\n" << rampLength
               << CRD::var;
      CRD::msg << "  Physical Length\n(";
      for (const int dir : EachDir)
        {
          if (dir != 0) CRD::msg << ',';
          CRD::msg <<  lclparam.GRIDphysicalLength[dir];
        }
      CRD::msg << ')' << CRD::var;
      CRD::msg.setFloatDefault();
      a_chparam.COORDSYSfactory =
        new SingleBlockCSAdaptorFactory(
          new SchwarzChristoffelRampCSFactory(
            CRDparam::g_domainBaseSize[0],
            alpha*Pi/180.,
            leadLength,
            rampLength),
          lclparam.GRIDphysicalLength);
    }
  else if (a_chparam.COORDSYStype ==
           CRDparam::CoordSysSingleBlockLogSchwarzChristoffelRamp)
    {
      Real alpha, leadLength, rampLength;
      rampParameters(ppCOORDSYS, alpha, leadLength, rampLength);
      Real length = leadLength + rampLength;
      Real height =
        lclparam.GRIDdomainLength[1]/lclparam.GRIDdomainLength[0]*length;
      lclparam.GRIDphysicalLength[0] = length;
      lclparam.GRIDphysicalLength[1] = height;
      RealVect COORDSYSshift;
      logshiftParameters(ppCOORDSYS, COORDSYSshift);
      CRD::msg << "Coordinate system\nLog-stretch Schwarz-Christoffel ramp"
               << CRD::var;
      CRD::msg.setPrecFloatSN(l_prec);
      CRD::msg << "  Angle of ramp (deg)\n" << alpha << CRD::var;
      CRD::msg << "  Lead length before ramp\n" << leadLength << CRD::var;
      CRD::msg << "  Ramp length (projected along x-axis)\n" << rampLength
               << CRD::var;
      CRD::msg << "  Physical Length\n(";
      for (const int dir : EachDir)
        {
          if (dir != 0) CRD::msg << ',';
          CRD::msg <<  lclparam.GRIDphysicalLength[dir];
        }
      CRD::msg << ')' << CRD::var;
      CRD::msg << "  Shift\n" << COORDSYSshift << CRD::var;
      CRD::msg.setFloatDefault();
    }
  else if (a_chparam.COORDSYStype ==
           CRDparam::CoordSysSingleBlockJoukowskiAirfoil)
    {
      // Default options
      Real COORDSYSchord            = 1.;
      Real COORDSYSspan             = 1.;
      Real COORDSYSdomainRatio      = 4.;
      Real COORDSYScamberRatio      = 0.05;
      Real COORDSYSthicknessRatio   = 0.15;
      Real COORDSYSalpha            = 0.;
      Real COORDSYScellRatioRadial  = 0.1;
      Real COORDSYScellRatioAzimuth = 0.2;
      ppCOORDSYS.query("airfoil_chord", COORDSYSchord);
      ppCOORDSYS.query("airfoil_span", COORDSYSspan);
      ppCOORDSYS.query("domain_ratio", COORDSYSdomainRatio);
      ppCOORDSYS.query("camber_ratio", COORDSYScamberRatio);
      ppCOORDSYS.query("thickness_ratio", COORDSYSthicknessRatio);
      ppCOORDSYS.query("angle_of_attack", COORDSYSalpha);
      ppCOORDSYS.query("radial_cell_ratio", COORDSYScellRatioRadial);
      ppCOORDSYS.query("azimuthal_cell_ratio", COORDSYScellRatioAzimuth);
      a_chparam.COORDSYSfactory =
        new SingleBlockCSAdaptorFactory(
          new JoukowskiAirfoilCSFactory(COORDSYSchord,
                                        COORDSYSspan,
                                        COORDSYSdomainRatio,
                                        COORDSYScamberRatio,
                                        COORDSYSthicknessRatio,
                                        COORDSYSalpha,
                                        COORDSYScellRatioRadial,
                                        COORDSYScellRatioAzimuth),
          RealVect_unit,
          SingleBlockCSAdaptor::TopologyDonutX);
      CRD::msg << "Coordinate system\nJoukowski airfoil" << CRD::var;
      CRD::msg.setPrecFloatSN(l_prec);
      CRD::msg << "  Chord\n" << COORDSYSchord << CRD::var;
      CRD::msg << "  Span\n" << COORDSYSspan << CRD::var;
      CRD::msg << "  Domain ratio\n" << COORDSYSdomainRatio << CRD::var;
      CRD::msg << "  Camber ratio\n" << COORDSYScamberRatio << CRD::var;
      CRD::msg << "  Thickness ratio\n" << COORDSYSthicknessRatio << CRD::var;
      CRD::msg << "  Angle of attack\n" << COORDSYSalpha << CRD::var;
      CRD::msg << "  Radial cell ratio\n" << COORDSYScellRatioRadial
               << CRD::var;
      CRD::msg << "  Azimuthal cell ratio\n" << COORDSYScellRatioAzimuth
               << CRD::var;
      CRD::msg.setFloatDefault();
    }

  CRDparam::CRDP.defineCoordSysType(a_chparam.COORDSYStype,
                                    lclparam.GRIDphysicalLength);

/*--------------------------------------------------------------------*
 * Solver parameters
 *--------------------------------------------------------------------*/

//--Read options that define the solver

  CRD::msg.newline();

  CRD::msg << "Parameters for 'solver'" << CRD::h2;
  ParmParse ppSOL("sol");

//--Read options necessary to setup the amr class

//--Fixed dt
  std::string timeStepMethod("RK4");
  ppSOL.query("time_stepping_method", timeStepMethod);
  if (timeStepMethod == "RK2")
    {
      lclparam.SOLtimeIntegrationMethod = CRDparam::RK2;
    }
  else if (timeStepMethod == "RK4")
    {
      lclparam.SOLtimeIntegrationMethod = CRDparam::RK4;
    }
  else if (timeStepMethod == "ARK4")
    {
      lclparam.SOLtimeIntegrationMethod = CRDparam::ARK4;
    }
  else
    {
      CRD::msg << "Unknown time integration method: "
               << timeStepMethod
               << " needs to be one of: RK2, RK4, or ARK4."
               << CRD::error;
  }
  CRD::msg << "Time Integration Method\n" << timeStepMethod
           << CRD::var;

  ppSOL.query("fixed_dt", a_chparam.SOLfixedDt);
  if (a_chparam.SOLfixedDt < 0.)
    {
      CRD::msg << "Input: 'fixed_dt' must be >= 0.!" << CRD::error;
    }
  CRD::msg << "Fixed dt\n";
  if (a_chparam.SOLfixedDt == 0.)
    {
      CRD::msg << "variable";
    }
  else
    {
      CRD::msg << a_chparam.SOLfixedDt;
    }
  CRD::msg << CRD::var;

//--Limit on time step growth

  ppSOL.query("max_dt_growth", a_chparam.SOLmaxDtGrowth);
  if (a_chparam.SOLmaxDtGrowth < 1.0)
    {
      CRD::msg << "Input: 'max_dt_growth' must by >= 1.!" << CRD::error;
    }
  CRD::msg << "Limit on dt growth\n" << a_chparam.SOLmaxDtGrowth << CRD::var;

//--Let the time step grow by this factor above the "maximum" before reducing it

  ppSOL.query("dt_tolerance_factor", a_chparam.SOLdtToleranceFactor);
  if (a_chparam.SOLdtToleranceFactor < 1.0)
    {
      CRD::msg << "Input: 'dt_tolerance_factor' must by >= 1.!" << CRD::error;
    }
  CRD::msg << "dt tolerance factor\n" << a_chparam.SOLdtToleranceFactor
           << CRD::var;

//--CFL

  ppSOL.query("cfl", lclparam.SOLcfl);
  if (lclparam.SOLcfl <= 0.)
    {
      CRD::msg << "Input: 'cfl' must by > 0.!" << CRD::error;
    }
  CRD::msg << "CFL number\n" << lclparam.SOLcfl << CRD::var;

//--Initial time step

  ppSOL.query("initial_dt", lclparam.SOLinitialDt);
  if (lclparam.SOLinitialDt > 0.)
    {
      CRD::msg << "Initial dt\n" << lclparam.SOLinitialDt << CRD::var;
    }
  else if (lclparam.SOLinitialDt == 0.)
    {
      CRD::msg << "Input: 'initial_dt' cannot be equal to 0!" << CRD::error;
    }

//--Maximum time step

  ppSOL.query("max_dt", lclparam.SOLmaxDt);
  if (lclparam.SOLmaxDt > 0.)
    {
      if (lclparam.SOLinitialDt > 0. &&
          lclparam.SOLmaxDt < lclparam.SOLinitialDt)
        {
          CRD::msg << "Input: 'max_dt' (" << lclparam.SOLmaxDt
                   << ") cannot be less than 'initial_dt' ("
                   << lclparam.SOLinitialDt << ")!" << CRD::error;
        }
      else
        {
          CRD::msg << "Max dt\n" << lclparam.SOLmaxDt << CRD::var;
        }
    }
  else if (lclparam.SOLmaxDt == 0.)
    {
      CRD::msg << "Input: 'max_dt' cannot be equal to 0!" << CRD::error;
    }

//--Inital CFL

  if (!internalProblem)
    {
      a_chparam.SOLinitialCFL = lclparam.SOLcfl;
    }
  if (a_chparam.SOLinitialCFL > lclparam.SOLcfl)
    {
      a_chparam.SOLinitialCFL = lclparam.SOLcfl;
    }
  const int haveInitialCFL = ppSOL.query("initial_cfl", a_chparam.SOLinitialCFL);
  if (a_chparam.SOLinitialCFL <= 0. ||
      a_chparam.SOLinitialCFL > lclparam.SOLcfl)
    {
      CRD::msg << "Input: 'initial_cfl' must be > 0. and <= 'cfl = "
               << lclparam.SOLcfl << "'!" << CRD::error;
    }
  // If an SOLinitialCFL isn't specified and SOLcfl is bigger,
  // then we want the SOLinitialCFL to be the same as SOLcfl
  if (!haveInitialCFL && lclparam.SOLcfl > a_chparam.SOLinitialCFL)
    {
      a_chparam.SOLinitialCFL = lclparam.SOLcfl;
    }
  CRD::msg << "Initial CFL number\n" << a_chparam.SOLinitialCFL << CRD::var;

//--Initial solution time

  ppSOL.query("initial_time", a_chparam.SOLinitialTime);
  if (a_chparam.SOLinitialTime < 0)
    {
      CRD::msg << "Input: 'initial_time' must be >= 0.0!" << CRD::error;
    }
  CRD::msg << "Initial time\n" << a_chparam.SOLinitialTime << CRD::var;

//--Read options necessary to run the amr class

//--Stop after this number of steps

  ppSOL.query("max_step", a_chparam.SOLmaxStep);
  if (a_chparam.SOLmaxStep < 0)
    {
      CRD::msg << "Input: 'max_step' must be >= 0!" << CRD::error;
    }
  CRD::msg << "Maximum steps\n" << a_chparam.SOLmaxStep << CRD::var;

//--Maximum solution time

  if (ppSOL.contains("max_time"))
    {
      std::string SOLmaxTimeStr;
      ppSOL.get("max_time", SOLmaxTimeStr);
      size_t spos = 0;
      if (SOLmaxTimeStr[0] == '+')
        {
          spos = 1;
        }
      std::istringstream SOLmaxTimeIStrS(SOLmaxTimeStr.substr(spos));
      SOLmaxTimeIStrS >> a_chparam.SOLmaxTime;
      if (!SOLmaxTimeIStrS)
        {
          CRD::msg << "Input: Invalid formatting of value '" << SOLmaxTimeStr
                   << "' for keyword 'max_time'!" << CRD::error;
        }
      if (spos == 1)
        {
          a_chparam.SOLmaxTime += a_chparam.SOLinitialTime;
        }
    }
  if (a_chparam.SOLmaxTime < 0)
    {
      CRD::msg << "Input: 'max_time' must be >= 0.0!" << CRD::error;
    }
  CRD::msg << "Maximum time\n" << a_chparam.SOLmaxTime << CRD::var;

  CRDparam::CRDP.defineSolver(lclparam.SOLcfl,
                              lclparam.SOLinitialDt,
                              lclparam.SOLmaxDt,
                              a_chparam.SOLinitialCFL, //**FIXME make local?
                              a_chparam.SOLinitialTime,
                              lclparam.SOLtimeIntegrationMethod);

/*--------------------------------------------------------------------*
 * Additive Runge Kutta parameters
 *--------------------------------------------------------------------*/

  CRD::msg.newline();

  CRD::msg << "Parameters for 'ark'" << CRD::h2;
  ParmParse ppARK("ark");
  //--Whether to run with explicit RK time stepping or additive RK time stepping
  ppARK.query("use", lclparam.ARKuse);
  CRD::msg << "Time Step Method\n";
  if (lclparam.ARKuse) {
    CRD::msg << "Additive Runge Kutta";
  } else {
    CRD::msg << "Explicit Runge Kutta";
  }
  CRD::msg << CRD::var;

  if (lclparam.ARKuse) {
    // Query for the max AMR level
    ppARK.query("max_ark_level", lclparam.ARKmaxLevel);
    if (lclparam.ARKmaxLevel < 0) {
      CRD::msg << "Input: 'max_ark_level' must be >= 0!" << CRD::error;
    }
    CRD::msg << "Max level to apply ARK time stepping\n" << lclparam.ARKmaxLevel << CRD::var;

    ppARK.query("init_erk_steps", lclparam.ARKinitERKSteps);
    if (lclparam.ARKinitERKSteps < 0) {
      CRD::msg << "Input: 'init_erk_steps' must be >= 0!" << CRD::error;
    }
    CRD::msg << "Number of initial ERK steps\n" << lclparam.ARKinitERKSteps << CRD::var;

    //--If doing additive RK, whether to use the numerical Jacobian or
    //analytical Jacobian
    ppARK.query("use_analytic_jacobian",
                lclparam.ARKuseAnalyticJac);
    CRD::msg << "Using analytical jacobian?\n";
    if (lclparam.ARKuseAnalyticJac) {
      CRD::msg << "Yes";
    } else {
      CRD::msg << "No";
    }
    CRD::msg << CRD::var;

    // Check for time scales, used for additive RK scaling
    ppARK.query("chemical_dt_scale", lclparam.ARKchemDtScale);
    if (lclparam.ARKchemDtScale <= 0.) {
      CRD::msg << "Input: 'chemical_dt_scale' must be >= 0!" << CRD::error;
    }
    CRD::msg << "Chemical dt scale\n" << lclparam.ARKchemDtScale << CRD::var;

    ppARK.query("pout_ls_stats", lclparam.ARKpoutLinearSolverStats);
    CRD::msg << "Printing Linear Solver Stats\n";
    if(lclparam.ARKpoutLinearSolverStats) {
      CRD::msg << "Yes";
    } else {
      CRD::msg << "No";
    }
    CRD::msg << CRD::var;

    // Convergence tolerance of nonlinear solver
    ppARK.query("nonlinear_convergence_tolerance", lclparam.ARKNonlinearTol);
    if (lclparam.ARKNonlinearTol <= 0.) {
      CRD::msg << "Input: 'nonlinear_convergence_tolerance' must be > 0!" << CRD::error;
    }
    CRD::msg << "ARK nonlinear convergence tolerance\n" << lclparam.ARKNonlinearTol << CRD::var;

    ppARK.query("use_pid_control", lclparam.ARKUsePIDControl);
    CRD::msg << "Using PID Time Step Size Control\n";
    if (lclparam.ARKUsePIDControl) {
      CRD::msg << "Yes";
    } else {
      CRD::msg << "No";
    }
    CRD::msg << CRD::var;

    ppARK.query("pid_epsilon", lclparam.ARKPIDControlEps);
    if (lclparam.ARKPIDControlEps <= 0.) {
      CRD::msg << "Input: 'nonlinear_convergence_tolerance' must be > 0!" << CRD::error;
    }
    CRD::msg << "PID control epsilon\n" << lclparam.ARKPIDControlEps << CRD::var;

    ppARK.query("extrapolate_initial_guess", lclparam.ARKExtrapInitGuess);
    CRD::msg << "Extrapolating nonlinear solver's initial guess\n";
    if (lclparam.ARKExtrapInitGuess) {
      CRD::msg << "Yes";
    } else {
      CRD::msg << "No";
    }
    CRD::msg << CRD::var;

  }

  CRDparam::CRDP.defineARK(lclparam.ARKuse,
                           lclparam.ARKmaxLevel,
                           lclparam.ARKuseAnalyticJac,
                           lclparam.ARKchemDtScale,
                           lclparam.ARKinitERKSteps,
                           lclparam.ARKpoutLinearSolverStats,
                           lclparam.ARKNonlinearTol,
                           lclparam.ARKUsePIDControl,
                           lclparam.ARKPIDControlEps,
                           lclparam.ARKExtrapInitGuess);

/*--------------------------------------------------------------------*
 * Limiter parameters
 *--------------------------------------------------------------------*/

//--Read options that define the limiter

  CRD::msg.newline();

  CRD::msg << "Parameters for 'limiter'" << CRD::h2;
  ParmParse ppLIMIT("limit");

//--General statement on limiting (catch all that sets other values)

  if (ppLIMIT.contains("method"))
    {
      std::string printLimitMethod("unknown");
      std::string inputLimitMethod;
      ppLIMIT.get("method", inputLimitMethod);
      if (inputLimitMethod == "SecondOrderPPM")
        {
          printLimitMethod = "second-order PPM";
          lclparam.LIMITfaceInterpolationOrder = 2;
          lclparam.LIMITusePPMlimiter          = 1;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 1;
        }
      else if (inputLimitMethod == "FourthOrderPPM")
        {
          printLimitMethod = "fourth-order PPM";
          lclparam.LIMITfaceInterpolationOrder = 4;
          lclparam.LIMITusePPMlimiter          = 1;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 1;
        }
      else if (inputLimitMethod == "FifthOrderPPM")
        {
          printLimitMethod = "fifth-order PPM";
          CRD::msg << "Fifth-order PPM has demonstrated issues!" << CRD::warn;
          lclparam.LIMITfaceInterpolationOrder = 5;
          lclparam.LIMITusePPMlimiter          = 1;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 1;
        }
      else if (inputLimitMethod == "FifthOrderUpwind")
        {
          printLimitMethod = "fifth-order upwind";
          lclparam.LIMITfaceInterpolationOrder = 5;
          lclparam.LIMITusePPMlimiter          = 0;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 0;
        }
      else if (inputLimitMethod == "FourthOrderCentral")
        {
          printLimitMethod = "fourth-order central";
          lclparam.LIMITfaceInterpolationOrder = 4;
          lclparam.LIMITusePPMlimiter          = 0;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 0;
        }
      else if (inputLimitMethod == "FirstOrderGodunov")
        {
          printLimitMethod = "first-order Godunov";
          lclparam.LIMITfaceInterpolationOrder = 1;
          lclparam.LIMITusePPMlimiter          = 0;
          lclparam.LIMITlimitFaceValues        = 0;
          lclparam.LIMITuseFlattening          = 0;
        }
      else if (inputLimitMethod == "SecondOrderCentral")
        {
          printLimitMethod = "second-order central";
          lclparam.LIMITfaceInterpolationOrder   = 2;
          lclparam.LIMITdiffusiveDerivativeOrder = 2;
          lclparam.LIMITreactionOrder            = 2;
          lclparam.LIMITusePPMlimiter            = 0;
          lclparam.LIMITlimitFaceValues          = 0;
          lclparam.LIMITuseFlattening            = 0;
          lclparam.LIMITcellConvolveFlatten      = 2;
          lclparam.LIMITfaceConvolveFlatten      = 2;
          lclparam.LIMITcellDeconvolveFlatten    = 2;
          lclparam.LIMITfaceDeconvolveFlatten    = 2;
        }
      else
        {
            CRD::msg << "Input: 'limit.method' of " << inputLimitMethod
                     << " not recognized!" << CRD::error;
        }
      if (ppLIMIT.contains("PPM_limiter") ||
          ppLIMIT.contains("limit_face_values") ||
          ppLIMIT.contains("use_flattening") ||
          ppLIMIT.contains("face_order"))
        {
          printLimitMethod = "modified";
        }
      CRD::msg << "Limiting method\n" << printLimitMethod << CRD::var;
    }

//--Whether or not to use PPM cell limiting

  ppLIMIT.query("PPM_limiter", lclparam.LIMITusePPMlimiter);
  if (lclparam.LIMITusePPMlimiter != 0)
    {
      lclparam.LIMITusePPMlimiter = 1;
    }
  CRD::msg << "Use PPM limiter\n"
           << ((lclparam.LIMITusePPMlimiter) ? "true" : "false")
           << CRD::var;

//--Whether to use limit face values

  ppLIMIT.query("limit_face_values", lclparam.LIMITlimitFaceValues);
  if (lclparam.LIMITlimitFaceValues != 0)
    {
      lclparam.LIMITlimitFaceValues = 1;
    }
  CRD::msg << "Limit face values\n"
           << ((lclparam.LIMITlimitFaceValues) ? "true" : "false")
           << CRD::var;

//--Whether to use flattening

  ppLIMIT.query("use_flattening", lclparam.LIMITuseFlattening);
  if (lclparam.LIMITuseFlattening != 0)
    {
      lclparam.LIMITuseFlattening = 1;
    }
  CRD::msg << "Use flattening\n"
           << ((lclparam.LIMITuseFlattening) ? "true" : "false") << CRD::var;
  if (lclparam.LIMITuseFlattening == 1 && lclparam.LIMITusePPMlimiter == 0)
    {
      CRD::msg << "Input: PPM limiting must be selected to use flattening!"
               << CRD::error;
    }
  if (CRDparam::g_plotFlattening && lclparam.LIMITuseFlattening != 1)
    {
      CRD::msg << "Input: Option to plot flattening must be disabled if "
        "'use_flattening' disabled!" << CRD::error;
    }

//--Which order to construct faces at (1, 4, or 5)

  ppLIMIT.query("face_order", lclparam.LIMITfaceInterpolationOrder);
  switch (lclparam.LIMITfaceInterpolationOrder)
    {
    case 1:
    case 2:
    case 4:
    case 5:
      break;
    default:
      CRD::msg << "Input: 'face_order' must be 1, 2, 4, or 5!" << CRD::error;
    }
  CRD::msg << "Face interpolation order\n"
           << lclparam.LIMITfaceInterpolationOrder << CRD::var;

//--Determine order of diffusive derivatives

  ppLIMIT.query("diffusion_order", lclparam.LIMITdiffusiveDerivativeOrder);
  if (lclparam.LIMITdiffusiveDerivativeOrder != 4 &&
      lclparam.LIMITdiffusiveDerivativeOrder != 2)
    {
      CRD::msg << "Input: 'diffusive_order' must be 2 or 4!" << CRD::error;
    }
  CRD::msg << "Diffusive order\n" << lclparam.LIMITdiffusiveDerivativeOrder
           << CRD::var;

  //--Determine order of reaction source terms

  ppLIMIT.query("reaction_order", lclparam.LIMITreactionOrder);
  if (lclparam.LIMITreactionOrder != 4 &&
      lclparam.LIMITreactionOrder != 2)
    {
      CRD::msg << "Input: 'reaction_order' must be 2 or 4!" << CRD::error;
    }
  CRD::msg << "Reaction order\n" << lclparam.LIMITreactionOrder
           << CRD::var;

//--Determine if Convolution/Deconvolution limiting is used
  if (ppLIMIT.contains("deconvolution_limiting"))
    {
      std::string deconvLimitStr = "false";
      ppLIMIT.get("deconvolution_limiting", deconvLimitStr);
      if (deconvLimitStr == "true")
        {
          lclparam.LIMITcellDeconvolveLimit = 1;
          lclparam.LIMITfaceDeconvolveLimit = 1;
        }
    }
  if (ppLIMIT.contains("cell_deconvolution_limiting"))
    {
      std::string deconvLimitStr = "false";
      ppLIMIT.get("cell_deconvolution_limiting", deconvLimitStr);
      if (deconvLimitStr == "true")
        {
          lclparam.LIMITcellDeconvolveLimit = 1;
        }
    }
  if (ppLIMIT.contains("face_deconvolution_limiting"))
    {
      std::string deconvLimitStr = "false";
      ppLIMIT.get("face_deconvolution_limiting", deconvLimitStr);
      if (deconvLimitStr == "true")
        {
          lclparam.LIMITfaceDeconvolveLimit = 1;
        }
    }

  bool anyFlattenSet = false;
  {

//--Check to see if D/C flattening is used, only if slope flattening is used

    if (ppLIMIT.contains("total_dc_flattening"))
      {
        std::string dcFlattenStr = "false";
        ppLIMIT.get("total_dc_flattening", dcFlattenStr);
        if (dcFlattenStr == "true")
          {
            lclparam.LIMITcellConvolveFlatten = 1;
            lclparam.LIMITfaceConvolveFlatten = 1;
            lclparam.LIMITcellDeconvolveFlatten = 1;
            lclparam.LIMITfaceDeconvolveFlatten = 1;
          }
        else if (dcFlattenStr == "false")
          {
            lclparam.LIMITcellConvolveFlatten = 0;
            lclparam.LIMITfaceConvolveFlatten = 0;
            lclparam.LIMITcellDeconvolveFlatten = 0;
            lclparam.LIMITfaceDeconvolveFlatten = 0;
          }
        else if (dcFlattenStr == "2nd")
          {
            lclparam.LIMITcellConvolveFlatten = 2;
            lclparam.LIMITfaceConvolveFlatten = 2;
            lclparam.LIMITcellDeconvolveFlatten = 2;
            lclparam.LIMITfaceDeconvolveFlatten = 2;
          }
      }
    if (ppLIMIT.contains("convolution_flattening"))
      {
        std::string convFlattenStr = "false";
        ppLIMIT.get("convolution_flattening", convFlattenStr);
        if (convFlattenStr == "true")
          {
            lclparam.LIMITcellConvolveFlatten = 1;
            lclparam.LIMITfaceConvolveFlatten = 1;
          }
        else if (convFlattenStr == "false")
          {
            lclparam.LIMITcellConvolveFlatten = 0;
            lclparam.LIMITfaceConvolveFlatten = 0;
          }
        else if (convFlattenStr == "2nd")
          {
            lclparam.LIMITcellConvolveFlatten = 2;
            lclparam.LIMITfaceConvolveFlatten = 2;
            CRD::msg << "2nd-order convolution flattening is not recommended!"
                     << CRD::warn;
          }
      }
    if (ppLIMIT.contains("deconvolution_flattening"))
      {
        std::string deconvFlattenStr = "false";
        ppLIMIT.get("deconvolution_flattening", deconvFlattenStr);
        if (deconvFlattenStr == "true")
          {
            lclparam.LIMITcellDeconvolveFlatten = 1;
            lclparam.LIMITfaceDeconvolveFlatten = 1;
          }
        else if (deconvFlattenStr == "false")
          {
            lclparam.LIMITcellDeconvolveFlatten = 0;
            lclparam.LIMITfaceDeconvolveFlatten = 0;
          }
        else if (deconvFlattenStr == "2nd")
          {
            lclparam.LIMITcellDeconvolveFlatten = 2;
            lclparam.LIMITfaceDeconvolveFlatten = 2;
          }
      }
    std::string inputFlattenMethod = "none";
    ppLIMIT.query("flatten_cell_convolution", inputFlattenMethod);
    if (inputFlattenMethod == "false")
      {
        lclparam.LIMITcellConvolveFlatten = 0;
      }
    else if (inputFlattenMethod == "true")
      {
        lclparam.LIMITcellConvolveFlatten = 1;
      }
    else if (inputFlattenMethod == "2nd")
      {
        lclparam.LIMITcellConvolveFlatten = 2;
      }
    inputFlattenMethod = "none"; // Reset string
    ppLIMIT.query("flatten_cell_deconvolution", inputFlattenMethod);
    if (inputFlattenMethod == "false")
      {
        lclparam.LIMITcellDeconvolveFlatten = 0;
      }
    else if (inputFlattenMethod == "true")
      {
        lclparam.LIMITcellDeconvolveFlatten = 1;
      }
    else if (inputFlattenMethod == "2nd")
      {
        lclparam.LIMITcellDeconvolveFlatten = 2;
      }
    inputFlattenMethod = "none"; // Reset string
    ppLIMIT.query("flatten_face_convolution", inputFlattenMethod);
    if (inputFlattenMethod == "false")
      {
        lclparam.LIMITfaceConvolveFlatten = 0;
      }
    else if (inputFlattenMethod == "true")
      {
        lclparam.LIMITfaceConvolveFlatten = 1;
      }
    else if (inputFlattenMethod == "2nd")
      {
        lclparam.LIMITfaceConvolveFlatten = 2;
      }
    inputFlattenMethod = "none"; // Reset string
    ppLIMIT.query("flatten_face_deconvolution", inputFlattenMethod);
    if (inputFlattenMethod == "false")
      {
        lclparam.LIMITfaceDeconvolveFlatten = 0;
      }
    else if (inputFlattenMethod == "true")
      {
        lclparam.LIMITfaceDeconvolveFlatten = 1;
      }
    else if (inputFlattenMethod == "2nd")
      {
        lclparam.LIMITfaceDeconvolveFlatten = 2;
      }
    int sumFlattenTest = lclparam.LIMITcellDeconvolveFlatten +
      lclparam.LIMITcellConvolveFlatten +
      lclparam.LIMITfaceDeconvolveFlatten +
      lclparam.LIMITfaceConvolveFlatten;
    if (sumFlattenTest > 0)
      {
        anyFlattenSet = true;
      }
  }
  if (anyFlattenSet)
    {
      CRD::msg << "Flatten cell convolution\n";
      switch(lclparam.LIMITcellConvolveFlatten)
        {
        case 1: CRD::msg << "true" << CRD::var;
          break;
        case 2: CRD::msg << "2nd-order" << CRD::var;
          break;
        default: CRD::msg << "false" << CRD::var;
          break;
        }
      CRD::msg << "Flatten cell deconvolution\n";
      switch(lclparam.LIMITcellDeconvolveFlatten)
        {
        case 1: CRD::msg << "true" << CRD::var;
          break;
        case 2: CRD::msg << "2nd-order" << CRD::var;
          break;
        default: CRD::msg << "false" << CRD::var;
          break;
        }
      CRD::msg << "Flatten face convolution\n";
      switch(lclparam.LIMITfaceConvolveFlatten)
        {
        case 1: CRD::msg << "true" << CRD::var;
          break;
        case 2: CRD::msg << "2nd-order" << CRD::var;
          break;
        default: CRD::msg << "false" << CRD::var;
          break;
        }
      CRD::msg << "Flatten face deconvolution\n";
      switch(lclparam.LIMITfaceDeconvolveFlatten)
        {
        case 1: CRD::msg << "true" << CRD::var;
          break;
        case 2: CRD::msg << "2nd-order" << CRD::var;
          break;
        default: CRD::msg << "false" << CRD::var;
        }
      // DC flatten tolerance
      ppLIMIT.query("dc_flat_tol", a_chparam.LIMITDCFlatTol);
      CRD::msg << "DC flatten tolerance\n"
               << a_chparam.LIMITDCFlatTol << CRD::var;
    }
  else
    {
      CRD::msg << "DC flattening\nfalse" << CRD::var;
    }
  ppLIMIT.query("use_fcor", lclparam.LIMITuseFCOR);
  CRD::msg << "Use face construction order reduction\n"
           << ((lclparam.LIMITuseFCOR) ? "true" : "false") << CRD::var;
  if (lclparam.LIMITuseFCOR != 0)
    {
      lclparam.LIMITuseFCOR = 1;
      // FCOR tolerance
      ppLIMIT.query("fcor_tol", a_chparam.LIMITFCORTol);
      CRD::msg << "FCOR tolerance\n" << a_chparam.LIMITFCORTol << CRD::var;
    }

//--Whether to use artificial viscosity

  ppLIMIT.query("use_artificial_viscosity", lclparam.LIMITuseArtVisc);
  if (lclparam.LIMITuseArtVisc != 0)
    {
      lclparam.LIMITuseArtVisc = 1;
    }
  CRD::msg << "Use artificial viscosity\n"
           << ((lclparam.LIMITuseArtVisc) ? "true" : "false") << CRD::var;

//--Artificial viscosity coefficient

  ppLIMIT.query("artificial_viscosity_coef", lclparam.LIMITartViscCoef);
  CRD::msg << "Artificial viscosity coefficient\n" << lclparam.LIMITartViscCoef
           << CRD::var;

//--Whether to use fourth-order artificial viscosity

  ppLIMIT.query("use_fourth_order_artificial_viscosity",
                lclparam.LIMITuseArtVisc4thO);
  if (lclparam.LIMITuseArtVisc4thO != 0)
    {
      lclparam.LIMITuseArtVisc4thO = 1;
    }
  // Turning this message off since 4thO Art. Visc. is not used
  // CRD::msg << "Use fourth-order artificial viscosity\n"
  //          << ((lclparam.LIMITuseArtVisc4thO) ? "true" : "false") << CRD::var;

//--Artificial viscosity coefficient

  ppLIMIT.query("fourth_order_artificial_viscosity_coef",
                lclparam.LIMITartViscCoef4thO);
  CRD::msg << "Fourth-order artificial viscosity coefficient\n"
           << lclparam.LIMITartViscCoef4thO << CRD::var;

//--Forces reduced order interpolation on faces near boundaries

  ppLIMIT.query("extra_boundary_limiting", lclparam.LIMITextraBoundLim);
  if (lclparam.LIMITextraBoundLim != 0)
    {
      lclparam.LIMITextraBoundLim = 1;
      CRD::msg << "Extra boundary limiting\ntrue" << CRD::var;
    }

//--Whether to use HO checks used in face limiter

  ppLIMIT.query("no_ho_checks", lclparam.LIMITnoHOchecks);
  if (lclparam.LIMITnoHOchecks != 0)
    {
      lclparam.LIMITnoHOchecks = 1;
      if ((bool)lclparam.LIMITlimitFaceValues)
        {
          CRD::msg << "No HO (3rd derivative) checks\ntrue" << CRD::var;
        }
      else
        {
          CRD::msg << "Turning off 3rd derivative checks is only applicable"
                   << " with face limiting!" << CRD::warn;
        }
    }

//--Whether to use isentropic and acoustic corrections at walls (default: true)

  ppLIMIT.query("wall_acoustic_isentropic_corrections",
                lclparam.LIMITwallCorrections);
  if (lclparam.LIMITwallCorrections != 0)
    {
      lclparam.LIMITwallCorrections = 1;
    }
  CRD::msg << "Use wall isentropic and acoustic corrections\n"
           << ((lclparam.LIMITwallCorrections) ? "true" : "false") << CRD::var;

//--Whether to use interpolation clipping and redistribution
  ppLIMIT.query("use_clipping", lclparam.LIMITclipping);
  CRD::msg << "Use clipping and redistribution\n"
           << ((lclparam.LIMITclipping) ? "true" : "false") << CRD::var;

  if (lclparam.LIMITclipping)
    {
      ppLIMIT.query("use_HO_clipping", lclparam.LIMITclippingHO);
      CRD::msg << "High-order clipping and redistribution\n"
               << ((lclparam.LIMITclippingHO) ? "true" : "false") << CRD::var;

      ppLIMIT.query("use_clipping_post_smoothing", lclparam.LIMITclippingPostSmooth);
      CRD::msg << "Post-smoothing on clipping and redistribution\n"
               << ((lclparam.LIMITclippingPostSmooth) ? "true" : "false") << CRD::var;
    }

  if (lclparam.LIMITfaceInterpolationOrder == 5)
    {

//--Blending coefficient to smooth transition between fifth-order (coeff = 0)
//--and fourth-order (coeff = 1)
      ppLIMIT.query("fourth_order_fifth_order_blending",
                    lclparam.LIMITfifthOrderBlendingCoef);
      CRD::msg << "Fifth-order blending coefficient\n"
               << lclparam.LIMITfifthOrderBlendingCoef << CRD::var;
    }

  CRDparam::CRDP.defineLimiter(lclparam.LIMITfaceInterpolationOrder,
                               lclparam.LIMITdiffusiveDerivativeOrder,
                               lclparam.LIMITreactionOrder,
                               (bool)lclparam.LIMITusePPMlimiter,
                               (bool)lclparam.LIMITlimitFaceValues,
                               (bool)lclparam.LIMITuseFlattening,
                               lclparam.LIMITcellConvolveFlatten,
                               lclparam.LIMITcellDeconvolveFlatten,
                               lclparam.LIMITfaceConvolveFlatten,
                               lclparam.LIMITfaceDeconvolveFlatten,
                               lclparam.LIMITcellConvolveLimit,
                               lclparam.LIMITcellDeconvolveLimit,
                               lclparam.LIMITfaceConvolveLimit,
                               lclparam.LIMITfaceDeconvolveLimit,
                               (bool)lclparam.LIMITuseFCOR,
                               (bool)lclparam.LIMITuseArtVisc,
                               lclparam.LIMITartViscCoef,
                               (bool)lclparam.LIMITuseArtVisc4thO,
                               lclparam.LIMITartViscCoef4thO,
                               (bool)lclparam.LIMITextraBoundLim,
                               (bool)lclparam.LIMITnoHOchecks,
                               (bool)lclparam.LIMITwallCorrections,
                               lclparam.LIMITclipping,
                               lclparam.LIMITclippingHO,
                               lclparam.LIMITclippingPostSmooth,
                               lclparam.LIMITfifthOrderBlendingCoef);

/*--------------------------------------------------------------------*
 * Turbulence modeling parameters
 *--------------------------------------------------------------------*/

//--Read options that define the turbulence model

  CRD::msg.newline();

  ParmParse ppTURB("turb");
  if (ppTURB.contains("turb_model") &&
      (lclparam.PHYSICSfluidModels & CRDparam::PhysicsViscous))
    {
      CRD::msg.newline();
      CRD::msg << "Parameters for 'turb'" << CRD::h2;
      std::string inputTurbModel;
      ppTURB.get("turb_model", inputTurbModel);
      CRD::msg << "Model\n";
      if (inputTurbModel == "LES")
        {
          lclparam.TURBmodelType = CRDparam::TurbModelLES;
          CRD::msg << "LES" << CRD::var;
          CRD::msg << "SGS model\n";

//--Subgrid scale models for LES

          if (ppTURB.contains("sgs_model"))
            {
              std::string inputSGSModel;
              ppTURB.get("sgs_model", inputSGSModel);
              if (inputSGSModel == "Smagorinsky")
                {
                  lclparam.TURBsgsModelType = CRDparam::SGSModelSmagorinsky;
                  CRD::msg << "Smagorinsky" << CRD::var;
                }
              else if (inputSGSModel == "stretched-vortex" ||
                       inputSGSModel == "SVSGS" ||
                       inputSGSModel == "StretchedVortex")
                {
                  lclparam.TURBsgsModelType = CRDparam::SGSModelStretchedVortex;
                  CRD::msg << "stretched-vortex" << CRD::var;
                  ppTURB.query("use_sgs_ke_pressure_correction",
                               lclparam.TURBuseConsToPrimCorrection);
                  ppTURB.query("use_sgs_coarsening",
                               lclparam.TURBuseSGSCoarsening);
                  std::vector<Real> streamwiseDir(SpaceDim, 0.);
                  streamwiseDir[0] = 1.; // x-direction by default
                  ppTURB.queryarr("streamwise_direction",
                                  streamwiseDir, 0, SpaceDim);
                  RealVect streamVect = lclparam.TURBglobStreamwiseDir;
                  for (int dir = 0; dir != SpaceDim; ++dir)
                    {
                      streamVect[dir] = streamwiseDir[dir];
                    }
                  if (streamVect.vectorLength() < 1.e-20)
                    {
                      CRD::msg << "Streamwise direction too small. Setting "
                               << "to {1,0,0}" << CRD::warn;
                      streamVect = RealVect{1., 0., 0.};
                    }
                  lclparam.TURBglobStreamwiseDir =
                    streamVect/(streamVect.vectorLength());
                  ppTURB.query("sgs_ke_coarsest_coarsening_ratio",
                               lclparam.TURBsgskeCrsLevFilterRatio);
                  if ((lclparam.TURBsgskeCrsLevFilterRatio > 1) &&
                      !CRDparam::CoordSysSingleBlockCartesian)
                    {
                      CRD::msg << "Input: coarsest level coarsening only "
                               << "works with Cartesian single-block cases!"
                               << CRD::abort;
                    }
                  ppTURB.query("wall_model_coarsening_ratio",
                               lclparam.TURBwallModelCrsRatio);
                  if (lclparam.TURBwallModelCrsRatio <= 0)
                    {
                      CRD::msg << "Input: wall_model_coarsening_ratio must "
                               << "be >= 1" << CRD::abort;
                    }
                }
              else
                {
                  CRD::msg << "Input: 'sgs_model' unrecognized!"
                           << CRD::error;
                }
            }
          else
            {
              lclparam.TURBsgsModelType = CRDparam::SGSModelSmagorinsky;
              CRD::msg << "default (Smagorinsky)" << CRD::var;
            }
        }
      else if (inputTurbModel == "SA")
        {
          lclparam.TURBmodelType = CRDparam::TurbModelSA;
          CRD::msg << "SA" << CRD::var;
        }
      else if (inputTurbModel == "none" || inputTurbModel == "None")
        {
          CRD::msg << "none" << CRD::var;
        }
      else
        {
          CRD::msg << "Input: 'turb_model' unrecognized!" << CRD::error;
        }
    }
  else if (ppTURB.contains("turb_model"))
    {
      CRD::msg << "Must use viscous physics to apply turbulence model!"
               << CRD::error;
    }
  if (ppTURB.contains("turb_forcing") &&
      (lclparam.PHYSICSfluidModels & CRDparam::PhysicsViscous))
    {
      lclparam.TURBuseForcing = true;
    }
  else if (ppTURB.contains("turb_forcing"))
    {
      CRD::msg << "Must use viscous physics to apply turbulent forcing!"
               << CRD::error;
    }
  // Law-of-the-wall shear-stress can be enforced without turbulence-model
  ppTURB.query("enforce_modeled_wall_shear_stress",
               lclparam.TURBenforceModeledWallShearStress);

//--Turbulence forcing

  if (ppTURB.contains("turb_forcing") &&
      (lclparam.PHYSICSfluidModels & CRDparam::PhysicsViscous))
    {
      std::string inputTurbForcing;
      ppTURB.get("turb_forcing", inputTurbForcing);
      CRD::msg << "Forcing type\n";

      if (inputTurbForcing=="spectral" || inputTurbForcing=="PetersenLivescu")
        {
          lclparam.TURBforcingType = CRDparam::TurbForcingSpectral;
          CRD::msg << "spectral (Petersen & Livescu 2010)" << CRD::var;

          if (ppTURB.contains("spectral_forcing_interval"))
            {
              std::vector<int> input_arr;
              ppTURB.getarr("spectral_forcing_interval", input_arr, 0, 2);
              if (input_arr[0] < 1 || input_arr[0] >= input_arr[1])
                {
                  CRD::msg << "Input: 'spectral_forcing_interval' must be a "
                           "valid interval within the range [1, infty).  "
                           "Given [" << input_arr[0] << "," << input_arr[1]
                           << ")!" << CRD::error;
                }
              lclparam.TURBspectralForcingInterval = Interval(input_arr[0],
                                                              input_arr[1]);
              CRD::msg << "Forcing interval\n["
                       << lclparam.TURBspectralForcingInterval.begin() << ','
                       << lclparam.TURBspectralForcingInterval.end() << ')'
                       << CRD::var;
            }
          if (!ppTURB.contains("spectral_forcing_eps"))
            {
              CRD::msg << "Input: 'spectral_forcing_eps' must be "
                "specified for spectral forcing!" << CRD::error;
            }
          else
            {
              ppTURB.get("spectral_forcing_eps",
                         lclparam.TURBspectralForcingEps);
              CRD::msg << "Forcing KE injection rate\n"
                       << lclparam.TURBspectralForcingEps
                       << " J/m^3-s"
                       << CRD::var;
            }
          if (!ppTURB.contains("spectral_forcing_dt"))
            {
              CRD::msg << "Input: 'spectral_forcing_dt' must be "
                "specified for spectral forcing!" << CRD::error;
            }
          else
            {
              ppTURB.get("spectral_forcing_dt",
                         lclparam.TURBspectralForcingDt);
              CRD::msg << "Forcing source refresh rate\n"
                       << lclparam.TURBspectralForcingDt
                       << " s"
                       << CRD::var;
            }
        }
      else if (inputTurbForcing == "linear" || inputTurbForcing == "Rosales" )
        {
          lclparam.TURBforcingType = CRDparam::TurbForcingLinear;
          CRD::msg << "linear (Rosales 2005)" << CRD::var;
        }
      else
        {
          lclparam.TURBforcingType = CRDparam::TurbForcingNone;
          CRD::msg << "None" << CRD::var;
        }
    }
  else if (ppTURB.contains("turb_forcing"))
    {
      CRD::msg << "Must use viscous physics to apply turbulent forcing!"
               << CRD::error;
    }

//--Explicit filtering

  if (ppTURB.contains("explicit_filter"))
    {
      CRD::msg << "Explicit filter\n";
      std::string explicitFilterName;
      ppTURB.get("explicit_filter", explicitFilterName);
      if (explicitFilterName == "spectral")
        {
          CRD::msg << "spectral" << CRD::var;
          lclparam.TURBexplicitFilterType = CRDparam::ExplicitFilterSpectral;
          if (!ppTURB.contains("spectral_domain_resolution"))
            {
              CRD::msg << "Input: 'spectral_domain_resolution' must be "
                "specified for spectral 'explicit_filter'!" << CRD::error;
            }
          else
            {
              ppTURB.getarr("spectral_domain_resolution",
                            lclparam.TURBspectralFilterDomainResolution,
                            0, SpaceDim);
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if (lclparam.TURBspectralFilterDomainResolution[dir] < 1 ||
                      lclparam.TURBspectralFilterDomainResolution[dir] >
                      lclparam.GRIDnumCells[dir])
                    {
                      CRD::msg << "Input: 'spectral_domain_resolution' must"
                        " be in range [1:num_cells].  Given " << lclparam
                        .TURBspectralFilterDomainResolution[dir] << " for "
                        "direction " << dir << '!' << CRD::error;
                    }
                }
              CRD::msg << "Spectral domain resolution\n(";
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if (dir != 0)
                    {
                      CRD::msg << ',';
                    }
                  CRD::msg << lclparam.TURBspectralFilterDomainResolution[dir];
                }
              CRD::msg << ')' << CRD::var;
            }
          CRD::msg << "Spectral filter profile\n";
          if (ppTURB.contains("spectral_filter_profile"))
            {
              std::string spectralFilterProfile;
              ppTURB.get("spectral_filter_profile", spectralFilterProfile);
              if (spectralFilterProfile == "sharp")
                {
                  CRD::msg << "sharp" << CRD::var;
                  lclparam.TURBspectralFilterProfile =
                    SpectralFilter::SpectralFilterSharp;
                }
              else if (spectralFilterProfile == "sharp_isotropic")
                {
                  CRD::msg << "sharp_isotropic" << CRD::var;
                  lclparam.TURBspectralFilterProfile =
                    SpectralFilter::SpectralFilterSharpIsotropic;
                }
              else if (spectralFilterProfile == "tanh")
                {
                  CRD::msg << "tanh" << CRD::var;
                  lclparam.TURBspectralFilterProfile =
                    SpectralFilter::SpectralFilterTanh;
                  ppTURB.query("spectral_filter_width",
                               lclparam.TURBspectralFilterParam);
                }
              else if (spectralFilterProfile == "gaussian")
                {
                  CRD::msg << "gaussian" << CRD::var;
                  lclparam.TURBspectralFilterProfile =
                    SpectralFilter::SpectralFilterGaussian;
                  ppTURB.query("spectral_filter_width",
                               lclparam.TURBspectralFilterParam);
                  if (lclparam.TURBspectralFilterParam >
                      lclparam.TURBspectralFilterDomainResolution[0])
                    {
                      CRD::msg << "Input: 'gaussian_filter_rolloff' must be "
                        << "<= 'spectral_domain_resolution[0]'. Given "
                        << lclparam.TURBspectralFilterParam
                        << " for 'spectral_filter_width' and "
                        << lclparam.TURBspectralFilterDomainResolution[0]
                        << " for 'spectral_domain_resolution[0]'. Setting "
                        << "'gaussian_filter_rolloff to "
                        << "'spectral_domain_resolution[0]'" << CRD::var;
                      lclparam.TURBspectralFilterParam =
                        lclparam.TURBspectralFilterDomainResolution[0];
                    }
                }
              else
                {
                  CRD::msg << "Input: 'spectral_filter_profile' "
                           << spectralFilterProfile << " not supported!"
                           << CRD::error;
                }
            }
          else
            {
              lclparam.TURBspectralFilterProfile =
                SpectralFilter::SpectralFilterSharp;
              CRD::msg << "default (sharp)" << CRD::var;
            }
          // Must be 3D, MPI, with a periodic rectangular Cartesian grid.  No
          // AMR.
          if (SpaceDim != 3)
            {
              CRD::msg << "Input: Spectral 'explicit_filter' can only be used "
                "in 3-D!" << CRD::error;
            }
#ifndef CH_MPI
          CRD::msg << "Input: Spectral 'explicit_filter' can only be used with "
            " MPI!" << CRD::error;
#endif
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (a_chparam.GRIDperiodic[dir] == 0)
                {
                  CRD::msg << "Input: Spectral 'explicit_filter' can only be "
                    "applied to periodic grids.  Direction " << dir << " is "
                    "not periodic!" << CRD::error;
                }
            }
          if (a_chparam.COORDSYStype != CRDparam::CoordSysSingleBlockCartesian)
            {
              CRD::msg << "Input: Spectral 'explicit_filter' can only be "
                "applied to rectangular Cartesian grids!" << CRD::error;
            }
          if (a_chparam.AMRmaxLevel > 0)
            {
              CRD::msg << "Input: Spectral 'explicit_filter' is not supported "
                "on AMR grids!" << CRD::error;
            }
        }
      else
        {
          CRD::msg << "none" << CRD::var;
        }
    }

//--Data processing
  // Query for averaging length
  ppTURB.query("start_averaging_time",
               lclparam.TURBstartTimeAvgTime);
  CRD::msg << "Start time for time-avg\n" << lclparam.TURBstartTimeAvgTime
           << CRD::var;

  CRDparam::CRDP.defineTurbulence(lclparam.TURBmodelType,
                                  lclparam.TURBsgsModelType,
                                  lclparam.TURBglobStreamwiseDir,
                                  lclparam.TURBuseForcing,
                                  lclparam.TURBuseConsToPrimCorrection,
                                  lclparam.TURBuseSGSCoarsening,
                                  lclparam.TURBexplicitFilterType,
                                  lclparam.TURBspectralFilterProfile,
                                  lclparam.TURBspectralFilterDomainResolution,
                                  lclparam.TURBspectralFilterParam,
                                  lclparam.TURBforcingType,
                                  lclparam.TURBspectralForcingInterval,
                                  lclparam.TURBspectralForcingEps,
                                  lclparam.TURBspectralForcingDt,
                                  lclparam.TURBsgskeCrsLevFilterRatio,
                                  lclparam.TURBwallModelCrsRatio,
                                  lclparam.TURBenforceModeledWallShearStress,
                                  lclparam.TURBstartTimeAvgTime);

/*--------------------------------------------------------------------*
 * Physics parameters
 *--------------------------------------------------------------------*/

//--Read options that define physics models

  CRD::msg.newline();

  CRD::msg << "Parameters for 'physics'" << CRD::h2;
  ParmParse ppPHYSICS("physics");

//--Fluid models

  if (ppPHYSICS.contains("fluid_models"))
    {
      lclparam.PHYSICSfluidModels = 0;  // Reset
      int numModel = ppPHYSICS.countval("fluid_models");
      std::vector<std::string> inputFluidModels(numModel);
      ppPHYSICS.getarr("fluid_models", inputFluidModels, 0, numModel);
      for (int idxModel = 0; idxModel != numModel; ++idxModel)
        {
          if (inputFluidModels[idxModel] == "inertial")
            {
              lclparam.PHYSICSfluidModels |= CRDparam::PhysicsInertial;
            }
          else if (inputFluidModels[idxModel] == "viscous")
            {
              lclparam.PHYSICSfluidModels |= CRDparam::PhysicsViscous;
            }
          else if (inputFluidModels[idxModel] == "multispecies")
            {
              lclparam.PHYSICSfluidModels |= CRDparam::PhysicsThermPerf;
            }
          else if (inputFluidModels[idxModel] == "source")
            {
              lclparam.PHYSICSfluidModels |= CRDparam::PhysicsSource;
            }
          else
            {
              CRD::msg << "Input: 'fluid_models' named "
                       << inputFluidModels[idxModel] << " not recognized!"
                       << CRD::error;
            }
        }
    }

  // Errors
  if (lclparam.PHYSICSfluidModels == 0)
    {
      CRD::msg << "Input: no fluid models defined!" << CRD::error;
    }

  // Output
  CRD::msg << "Fluid models:" << CRD::body;
  if (lclparam.PHYSICSfluidModels & CRDparam::PhysicsInertial)
    {
      CRD::msg << "Inertial (hyperbolic terms in Navier-Stokes)" << CRD::body;
    }
  if (lclparam.PHYSICSfluidModels & CRDparam::PhysicsViscous)
    {
      CRD::msg << "Viscous (elliptic terms in Navier-Stokes)" << CRD::body;
    }
  if (lclparam.PHYSICSfluidModels & CRDparam::PhysicsThermPerf)
    {
      CRD::msg << "Multispecies (species transport and reactions)" << CRD::body;
    }
  if (lclparam.PHYSICSfluidModels & CRDparam::PhysicsSource)
    {
      CRD::msg << "Body force source term" << CRD::body;
    }

  // IMPORTANT: defineLimiter must be called before definePhysics!
  // IMPORTANT: defineTurbulence must be called before definePhysics!
  //**FIXME IBC can still change domain boundaries.  Also, multiblock?
  CRDparam::CRDP.definePhysics(lclparam.PHYSICSfluidModels,
                               a_chparam.GRIDperiodic);

/*--------------------------------------------------------------------*
 * Fluid parameters
 *--------------------------------------------------------------------*/

//--Read options that define fluid parameters

  CRD::msg.newline();

  CRD::msg << "Parameters for 'fluid'" << CRD::h2;
  ParmParse ppFLUID("fluid");

//--Universal gas constant

  ppFLUID.query("R", lclparam.FLUIDR);
  if (lclparam.FLUIDR <= 0.)
    {
      CRD::msg << "Input: 'R' must by > 0.!" << CRD::error;
    }

//--Gamma

  ppFLUID.query("gamma", lclparam.FLUIDgamma);
  if (lclparam.FLUIDgamma <= 0.)
    {
      CRD::msg << "Input: 'gamma' must by > 0.!" << CRD::error;
    }

//--Density

  ppFLUID.query("rho", lclparam.FLUIDrho);
  if (lclparam.FLUIDrho <= 0.)
    {
      CRD::msg << "Input: 'rho' must by > 0.!" << CRD::error;
    }

//--Temperature

  ppFLUID.query("T", lclparam.FLUIDT);
  if (lclparam.FLUIDT <= 0.)
    {
      CRD::msg << "Input: 'T' must by > 0.!" << CRD::error;
    }
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      CRD::msg << "Gas constant\nfunction of species" << CRD::var;
      CRD::msg << "Specific heat ratio\nfunction of temperature" << CRD::var;
    }
  else
    {
      CRD::msg << "Gas constant\n" << lclparam.FLUIDR << CRD::var;
      CRD::msg << "Specific heat ratio\n" << lclparam.FLUIDgamma << CRD::var;
      CRD::msg << "Density\n" << lclparam.FLUIDrho << CRD::var;
      CRD::msg << "Temperature\n" << lclparam.FLUIDT << CRD::var;
    }

//--Thermal conductivity and dynamic viscosity

  ppFLUID.query("K", lclparam.FLUIDK);
  ppFLUID.query("mu", lclparam.FLUIDmu);
  if (lclparam.FLUIDK < 0. && lclparam.FLUIDmu < 0.)
    {
      CRD::msg << "Thermal conductivity\nfunction of temperature" << CRD::var;
      CRD::msg << "Dynamic viscosity\nfunction of temperature" << CRD::var;
    }
  else if (lclparam.FLUIDK >= 0. && lclparam.FLUIDmu >= 0.)
    {
      CRD::msg << "Thermal conductivity\nconstant " <<
        lclparam.FLUIDK << CRD::var;
      CRD::msg << "Dynamic viscosity\nconstant " <<
        lclparam.FLUIDmu << CRD::var;
    }
  else
    {
      CRD::msg << "Input: 'K' and 'mu', both must be either constant or "
               << "solved for (= -1). Cannot set one and solve the other"
               << CRD::error;
    }

//--Reynolds number

  ppFLUID.query("Re", lclparam.FLUIDRe);
  if (lclparam.FLUIDRe < 0.)
    {
      CRD::msg << "Input: 'Re' must by >= 0.!" << CRD::error;
    }
  CRD::msg << "Reynolds number\n" << lclparam.FLUIDRe << CRD::var;

  CRDparam::CRDP.defineFluid(lclparam.FLUIDnumStates,
                             lclparam.FLUIDrho,
                             lclparam.FLUIDT,
                             lclparam.FLUIDspeed,
                             lclparam.FLUIDmu,
                             lclparam.FLUIDlambda,
                             lclparam.FLUIDgamma,
                             lclparam.FLUIDR,
                             lclparam.FLUIDgrav,
                             lclparam.FLUIDRe,
                             lclparam.FLUIDK);

/*--------------------------------------------------------------------*
 * AMR parameters
 *--------------------------------------------------------------------*/

//--Read options that define AMR parameters

  CRD::msg.newline();

  CRD::msg << "Parameters for 'AMR'" << CRD::h2;
  ParmParse ppAMR("amr");

//--Update grids on restart

  ppAMR.query("regrid_on_restart", a_chparam.AMRregridOnRestart);
  if (a_chparam.AMRregridOnRestart)
    {
      CRD::msg << "Doing regrid on restart\n" << CRD::var;
    }

//--Maximum number of AMR levels

  if (ppAMR.contains("max_level"))
    {
      ppAMR.get("max_level", a_chparam.AMRmaxLevel);
      if (a_chparam.AMRmaxLevel < 0)
        {
          CRD::msg << "Input: 'max_level' must be >= 0!" << CRD::error;
        }
    }
  else if (!internalProblem)
    {
      CRD::msg << "Input: 'max_level' must be specified if internal problem "
        "type not defined!" << CRD::error;
    }
  CRD::msg << "Maximum number of AMR levels\n" << a_chparam.AMRmaxLevel
           << CRD::var;
  // Size of some vector arrays related to number of levels.
  int numReadLevels = std::max(a_chparam.AMRmaxLevel, 1);

//--Refinement ratios

  if (ppAMR.contains("ref_ratio"))
    {
      a_chparam.AMRrefRatios.clear();
      a_chparam.AMRrefRatios.assign(a_chparam.AMRmaxLevel+1, 1);
      if (a_chparam.AMRmaxLevel > 0)
        {
          ppAMR.getarr("ref_ratio", a_chparam.AMRrefRatios, 0,
                       a_chparam.AMRmaxLevel);
          for (int iLev = 0; iLev != a_chparam.AMRmaxLevel; ++iLev)
            {
              if (a_chparam.AMRrefRatios[iLev] < 1 ||
                  !Misc::isPower2(a_chparam.AMRrefRatios[iLev]))
                {
                  CRD::msg << "Input: 'ref_ratio' must be >= 1 and power of 2 "
                    "for level " << iLev << '!' << CRD::error;
                }
            }
        }
    }
  else if (!internalProblem)
    {
      CRD::msg << "Input: 'ref_ratio' must be specified if internal problem "
        "type not defined!" << CRD::error;
    }
  CRD::msg << "Refinement ratios\n";
  for (int iLev = 0; iLev <= a_chparam.AMRmaxLevel; ++iLev)
    {
      if (iLev != 0)
        {
          CRD::msg << ' ';
        }
      CRD::msg << a_chparam.AMRrefRatios[iLev];
    }
  CRD::msg << CRD::var;

//--Regrid interval

  if (ppAMR.contains("regrid_interval"))
    {
      a_chparam.AMRregridIntervals.clear();
      a_chparam.AMRregridIntervals.assign(numReadLevels, 1);
      if (a_chparam.AMRmaxLevel > 0)
        {
          ppAMR.getarr("regrid_interval", a_chparam.AMRregridIntervals, 0,
                       numReadLevels);
        }
      for (int iLev = 0; iLev != numReadLevels; ++iLev)
        {
          if ((a_chparam.AMRregridIntervals[iLev] < 1) &&
              (a_chparam.AMRregridIntervals[iLev] != -1))
            {
              CRD::msg << "Input: 'regrid_interval' must be >= 1 or set to -1 "
                "to disable for grid level " << iLev << '!' << CRD::error;
            }
        }
    }
  else if (!internalProblem)
    {
      CRD::msg << "Input: 'regrid_interval' must be specified if internal "
        "problem type not defined!" << CRD::error;
    }
  CRD::msg << "Regrid intervals\n";
  for (int iLev = 0; iLev != numReadLevels; ++iLev)
    {
      if (iLev != 0)
        {
          CRD::msg << ' ';
        }
      CRD::msg << a_chparam.AMRregridIntervals[iLev];
    }
  CRD::msg << CRD::var;

//--Tag buffer size

  ppAMR.query("tag_buffer_size", a_chparam.AMRtagBufferSize);
  if (a_chparam.AMRtagBufferSize < 0)
    {
      CRD::msg << "Input: 'tag_buffer_size' must be > 0!" << CRD::error;
    }
  CRD::msg << "Tag buffer size\n" << a_chparam.AMRtagBufferSize << CRD::var;

//--Base level (this and lower levels are fixed)

  ppAMR.query("base_level", a_chparam.AMRbaseLevel);
  if (a_chparam.AMRbaseLevel < 0)
    {
      CRD::msg << "Input: 'base_level' must be > 0!" << CRD::error;
    }
  CRD::msg << "Base level\n" << a_chparam.AMRbaseLevel << CRD::var;
  if (a_chparam.AMRbaseLevel > 0)
    {
      CRD::msg << "Notice: When setting base level, it is advised that tags "
        "for base and lower levels not be solution dependent." << CRD::warn;
      CRD::msg << "Notice: Revising regrid intervals to -1 for levels < base "
        "level to -1." << CRD::warn;
      for (int iLev = 0; iLev != a_chparam.AMRbaseLevel; ++iLev)
        {
          a_chparam.AMRregridIntervals[iLev] = -1;
        }
      CRD::msg << "Revised regrid intervals\n";
      for (int iLev = 0; iLev != numReadLevels; ++iLev)
        {
          if (iLev != 0)
            {
              CRD::msg << ' ';
            }
          CRD::msg << a_chparam.AMRregridIntervals[iLev];
        }
      CRD::msg << CRD::var;
    }

//--Max grid size

  ppAMR.query("max_grid_size", a_chparam.AMRmaxGridSize);
  if (a_chparam.AMRmaxGridSize < 0)
    {
      CRD::msg << "Input: 'max_grid_size' must be > 0!" << CRD::error;
    }
  CRD::msg << "Maximum box size\n" << a_chparam.AMRmaxGridSize << CRD::var;

//--AMR fill ratio

  ppAMR.query("fill_ratio", a_chparam.AMRfillRatio);
  if (a_chparam.AMRfillRatio < 0. || a_chparam.AMRfillRatio > 1.)
    {
      CRD::msg << "Input: 'fill_ratio' must be >= 0. and <= 1.!" << CRD::error;
    }
  CRD::msg << "Fill ratio\n" << a_chparam.AMRfillRatio << CRD::var;

//--Grid buffer size

//  We don't know the number of ghosts until defineLimiter and definePhysics are
//  called
  const int codeReqBufferSize = LevelGridMetrics::bufferSize4thO(
    a_chparam.AMRrefRatios,
    a_chparam.AMRmaxLevel,
    CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg));
  if (ppAMR.contains("grid_buffer_size"))
    {
      CRD::msg << "Input: Do not manually set 'grid_buffer_size' as it is best "
        "determined automatically by the software." << CRD::warn;
      ppAMR.get("grid_buffer_size", a_chparam.AMRgridBufferSize);
      if (a_chparam.AMRgridBufferSize < 0)
        {
          CRD::msg << "Input: 'grid_buffer_size' must be > 0!" << CRD::error;
        }
      CRD::msg << "Input: 'grid_buffer_size' of " << a_chparam.AMRgridBufferSize
               << " is less than minimum buffer size, " << codeReqBufferSize
               << CRD::warn;
    }
  else
    {
      a_chparam.AMRgridBufferSize = codeReqBufferSize;
    }
  CRD::msg << "Grid buffer size\n" << a_chparam.AMRgridBufferSize << CRD::var;

//--Block factor

  ppAMR.query("block_factor", a_chparam.AMRblockFactor);
  if (a_chparam.AMRblockFactor < 1 || !Misc::isPower2(a_chparam.AMRblockFactor))
    {
      CRD::msg << "Input: 'block_factor' must be >= 1 and a power of 2!"
               << CRD::error;
    }
  CRD::msg << "Block factor (minimum box size)\n" << a_chparam.AMRblockFactor
           << CRD::var;
  if (a_chparam.AMRblockFactor <
      CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg))
    {
      CRD::msg << "'block_factor' must be >= "
               << CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg)
               << ", which is the number of interior ghost cells!"
               << CRD::error;
    }

//--Use subcycling

  ppAMR.query("use_subcycling", a_chparam.AMRuseSubcycling);
  bool useSubcyclingBool = false;
  if (a_chparam.AMRuseSubcycling != 0)
    {
      a_chparam.AMRuseSubcycling = 1;
      useSubcyclingBool = true;
    }
  CRD::msg << "Use subcycling\n"
           << ((a_chparam.AMRuseSubcycling) ? "true" : "false") << CRD::var;

  CRDparam::CRDP.defineAMR(a_chparam.AMRblockFactor,
                           a_chparam.AMRmaxGridSize,
                           a_chparam.AMRrefRatios,
                           useSubcyclingBool);

/*--------------------------------------------------------------------*
 * Combustion parameters
 *--------------------------------------------------------------------*/

//--Read options that define the combustion

  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      CRD::msg.newline();

      CRD::msg << "Parameters for 'therm'" << CRD::h2;
      ParmParse ppTHERM("therm");

      lclparam.THERMnumSpecies = ppTHERM.countval("species");

      if (lclparam.THERMnumSpecies <= 0)
        {
          CRD::msg << "Input: 'num_species' must be > 0!" << CRD::error;
        }
      CRD::msg << "Number of species\n" << lclparam.THERMnumSpecies << CRD::var;

      const int numSpecies = lclparam.THERMnumSpecies;
      lclparam.THERMspeciesNames.assign(numSpecies,"U");
      for (int row = 0; row != numSpecies; ++row)
        {
          ppTHERM.get("species", lclparam.THERMspeciesNames[row], row);
          CRD::msg << "Species " << row << "\n"
                   << lclparam.THERMspeciesNames[row] << CRD::var;
        }
      // Check to ensure the same species is not present multiple times
      for (int i = 0; i != numSpecies; ++i)
        {
          std::string type1 = lclparam.THERMspeciesNames[i];
          for (int j = 0; j != numSpecies; ++j)
            {
              std::string type2 = lclparam.THERMspeciesNames[j];
              if (i != j && type1 == type2)
                {
                  CRD::msg << "Input: 'species' the species "
                           << lclparam.THERMspeciesNames[i]
                           << " is repeated in species list" << CRD::error;
                }
            }
        }
      ppTHERM.query("use_species_correction",
                    lclparam.THERMuseSpeciesCorrection);

//--Read the reaction file

      ppTHERM.query("num_of_reactions",lclparam.THERMnumReactions);
      if (lclparam.THERMnumReactions < 0)
        {
          CRD::msg << "Input: 'num_of_reactions' must be >= 0!" << CRD::error;
        }
      CRD::msg << "Number of reactions\n" << lclparam.THERMnumReactions
               << CRD::var;

      if (lclparam.THERMnumReactions > 0)
        {
          // Check for reaction start and end times
          ppTHERM.query("reaction_start_time", lclparam.THERMreactionStartTime);
          lclparam.THERMreactionEndTime = a_chparam.SOLmaxTime;
          ppTHERM.query("reaction_end_time", lclparam.THERMreactionEndTime);
          if (ppTHERM.contains("reaction_start_time"))
            {
              if (lclparam.THERMreactionStartTime >
                  lclparam.THERMreactionEndTime)
                {
                  CRD::msg << "Input: 'reaction_start_time' is less than "
                           << "'reaction_end_time'!" << CRD::warn;
                }
              CRD::msg << "Reaction start time\n"
                       << lclparam.THERMreactionStartTime << CRD::var;
              CRD::msg << "Reaction end time\n"
                       << lclparam.THERMreactionEndTime << CRD::var;
            }
        }

      CRDparam::CRDP.defineCombustion(lclparam.THERMnumSpecies,
                                      lclparam.THERMreactionStartTime,
                                      lclparam.THERMreactionEndTime,
                                      lclparam.THERMuseSpeciesCorrection,
                                      lclparam.THERMnumReactions,
                                      lclparam.THERMspeciesNames);
    }
  else if ((CRDparam::g_mu < 0. || CRDparam::g_K < 0.) &&
           (lclparam.PHYSICSfluidModels & CRDparam::PhysicsViscous))
    {
      CRD::msg << "Cannot solve for 'mu' or 'K' unless using species "
               << "transport physics" << CRD::error;
    }


/*--------------------------------------------------------------------*
 * Threading parameters
 *--------------------------------------------------------------------*/

//--Read options that define the threading

  CRD::msg.newline();

  CRD::msg << "Parameters for 'threads'" << CRD::h2;
  {
    ParmParse ppThreads("threads");
    a_chparam.THREADuseThreads = ppThreads.countval("useThreads");
    ppThreads.query("useThreads", a_chparam.THREADuseThreads);

    CRD::msg << "Use threads\n"
             << (a_chparam.THREADuseThreads ? "true" : "false") << CRD::var;

    if (a_chparam.THREADuseThreads)
      {
        ppThreads.query("numThreads", a_chparam.THREADnumThreads);
        ppThreads.query("threadsPerTeam", a_chparam.THREADthreadsPerTeam);
        ppThreads.query("extraStackSize_MiB",
                        a_chparam.THREADextraStackSize_MiB);
        ppThreads.query("useHyperThreading", a_chparam.THREADuseHyperThreading);

        CRD::msg << "Number of threads\n" << a_chparam.THREADnumThreads
                 << CRD::var;
        CRD::msg << "Number of threads per team\n" << a_chparam.THREADthreadsPerTeam
                 << CRD::var;
        CRD::msg << "extra stacksize (MiB)\n"
                 << a_chparam.THREADextraStackSize_MiB << CRD::var;
        CRD::msg << "Use hyper threading\n"
                 << (a_chparam.THREADuseHyperThreading  ? "true" : "false")
                 << CRD::var;
      }
  }

  // Set to terminate on error
  CRD::msg.setTerminateOnError(true);
}
