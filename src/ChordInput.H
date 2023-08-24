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
 * \file ChordInput.H
 *
 * \brief Argument processing and checking
 *
 *//*+*************************************************************************/

#ifndef _CHORDINPUT_H_
#define _CHORDINPUT_H_

//----- Standard Library -----//

#include <string>
#include <vector>

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class MultiBlockCoordSysFactory;
#include "NamespaceFooter.H"

class TagLevelFactory;

//--Options used to initialize Chombo AMR classes.  These are destroyed once
//--the classes are constructed.  Persistent options are instead placed in
//--class CRDparam

struct ChomboParam
{
  // FILE
  std::string FILErestartFile;
  bool FILEdoRestart;
  std::string FILEplotPrefix;
  int FILEplotInterval;
  Real FILEplotPeriod;
  std::string FILEcheckpointPrefix;
  int FILEcheckpointInterval;
  // GRID
  std::vector<int> GRIDperiodic;
  // COORDSYS
  CRDparam::CoordSysType COORDSYStype;
  // Instead of just giving the type (above), the CoordSys factory can be
  // allocated and defined.  This allows for reading "custom" input for a
  // specific coordinate system.
  MultiBlockCoordSysFactory* COORDSYSfactory;
  // SOLUTION
  Real SOLfixedDt;
  Real SOLmaxDtGrowth;
  Real SOLdtToleranceFactor;
  Real SOLinitialCFL;
  Real SOLinitialTime;
  int SOLmaxStep;
  Real SOLmaxTime;
  // AMR
  int AMRmaxLevel;
  std::vector<int> AMRrefRatios;
  std::vector<int> AMRregridIntervals;
  int AMRtagBufferSize;
  int AMRbaseLevel;
  int AMRmaxGridSize;
  Real AMRfillRatio;
  int AMRgridBufferSize;
  int AMRblockFactor;
  int AMRuseSubcycling;
  bool AMRregridOnRestart;
  // TAG
  TagLevelFactory* TAGfactory;
  // LIMITER
  Real LIMITDCFlatTol;
  Real LIMITFCORTol;
  // THREADS
  bool THREADuseThreads;
  int THREADnumThreads;
  int THREADthreadsPerTeam;
  size_t THREADextraStackSize_MiB;
  bool THREADuseHyperThreading;
};

//--Options used internally in ChordInput.  These are immediately placed in
//--class CRDparam

struct LocalParam
{
  // SIMULATION
  CRDparam::ProblemType SIMproblemType;
  // GRID
  std::vector<int> GRIDnumCells;
  std::vector<Real> GRIDdomainLength;
  std::vector<Real> GRIDdomainOrigin;
  RealVect GRIDphysicalLength;
  // SOLUTION
  Real SOLcfl;
  Real SOLinitialDt;
  Real SOLmaxDt;
  CRDparam::TimeIntegrationMethod SOLtimeIntegrationMethod;
  // ARK
  Real ARKchemDtScale;
  bool ARKuse;
  bool ARKuseAnalyticJac;
  int  ARKmaxLevel;
  int  ARKinitERKSteps;
  bool ARKpoutLinearSolverStats;
  bool ARKpoutDtInfo;
  Real ARKNonlinearTol;
  bool ARKUsePIDControl;
  Real ARKPIDControlEps;
  bool ARKExtrapInitGuess;
  // LIMITER
  int LIMITfaceInterpolationOrder;
  int LIMITdiffusiveDerivativeOrder;
  int LIMITreactionOrder;
  int LIMITusePPMlimiter;
  int LIMITlimitFaceValues;
  int LIMITuseFlattening;
  int LIMITcellConvolveFlatten;
  int LIMITcellDeconvolveFlatten;
  int LIMITfaceConvolveFlatten;
  int LIMITfaceDeconvolveFlatten;
  int LIMITcellConvolveLimit;
  int LIMITcellDeconvolveLimit;
  int LIMITfaceConvolveLimit;
  int LIMITfaceDeconvolveLimit;
  int LIMITuseFCOR;
  int LIMITuseArtVisc;
  Real LIMITartViscCoef;
  int LIMITuseArtVisc4thO;
  Real LIMITartViscCoef4thO;
  int LIMITextraBoundLim;
  int LIMITnoHOchecks;
  int LIMITuseIRS;
  Real LIMITbetaIRS;
  int LIMITwallCorrections;
  bool LIMITclipping;
  bool LIMITclippingHO;
  bool LIMITclippingPostSmooth;
  Real LIMITfifthOrderBlendingCoef;
  // FLUID
  int  FLUIDnumStates;
  Real FLUIDrho;
  Real FLUIDT;
  Real FLUIDK;
  Real FLUIDspeed;
  Real FLUIDmu;
  Real FLUIDlambda;
  Real FLUIDgamma;
  Real FLUIDR;
  Real FLUIDgrav;
  Real FLUIDRe;
  // PHYSICS
  int PHYSICSfluidModels;
  // THERM
  int THERMnumSpecies;
  Real THERMreactionStartTime;
  Real THERMreactionEndTime;
  bool THERMuseSpeciesCorrection;
  int THERMnumReactions;
  std::vector<std::string> THERMspeciesNames;
  std::string THERMreactionFile;
  std::string THERMthermoFile;
  std::string THERMtransFile;
  // TURB
  CRDparam::TurbModelType TURBmodelType;
  CRDparam::SGSModelType TURBsgsModelType;
  RealVect TURBglobStreamwiseDir;
  bool TURBuseForcing;
  bool TURBuseConsToPrimCorrection;
  bool TURBuseSGSCoarsening;
  CRDparam::ExplicitFilterType TURBexplicitFilterType;
  SpectralFilter::SpectralFilterProfile TURBspectralFilterProfile;
  std::vector<int> TURBspectralFilterDomainResolution;
  Real TURBspectralFilterParam;
  CRDparam::TurbForcingType TURBforcingType;
  Interval TURBspectralForcingInterval;
  Real TURBspectralForcingEps;
  Real TURBspectralForcingDt;
  int  TURBsgskeCrsLevFilterRatio;
  int  TURBwallModelCrsRatio;
  bool TURBenforceModeledWallShearStress;
  Real TURBstartTimeAvgTime;
};

/// Parse the input, check, and print
void chordInput(ChomboParam& a_chparam);

//FIXME: This should be removed when CNSIBC files read in coordsys
/// Return parameters for Mach reflection case
void rampParameters(Real& a_alpha,
                    Real& a_leadLength,
                    Real& a_rampLength);

void logshiftParameters(RealVect& a_shift);

template <typename T>
struct ModDirect
{
  T operator()(const T& a_val) const
    {
      return a_val;
    }
};

template <typename T>
struct ModConvertBool
{
  ModConvertBool(const T& a_valFalse)
  :
  m_valFalse(a_valFalse)
    { }
  bool operator()(const T& a_val) const
    {
      return (a_val != m_valFalse);
    }
  T m_valFalse;
};

template <typename T>
struct ModAdd
{
  ModAdd(const T& a_valAdd)
  :
  m_valAdd(a_valAdd)
    { }
  T operator()(const T& a_val) const
    {
      return (a_val + m_valAdd);
    }
  T m_valAdd;
};

/// Load a SpaceDim bool array
template <typename T, typename S, typename Op = ModDirect<S> >
struct SpaceDimArray
{
  static void loadFromScalar(T *const a_p, const S& a_val, Op op = Op())
    {
      D_TERM6(a_p[0] = op(a_val);,
              a_p[1] = op(a_val);,
              a_p[2] = op(a_val);,
              a_p[3] = op(a_val);,
              a_p[4] = op(a_val);,
              a_p[5] = op(a_val);)
    }
  static void loadFromArray(T *const a_p, const S *const a_q, Op op = Op())
    {
      D_TERM6(a_p[0] = op(a_q[0]);,
              a_p[1] = op(a_q[1]);,
              a_p[2] = op(a_q[2]);,
              a_p[3] = op(a_q[3]);,
              a_p[4] = op(a_q[4]);,
              a_p[5] = op(a_q[5]);)
    }
};

/// Case insensitive string equality
inline bool iequals(const std::string& a, const std::string& b)
{
    return std::equal(a.begin(), a.end(),
                      b.begin(), b.end(),
                      [](char a, char b) {
                          return tolower(a) == tolower(b);
                      });
}

#endif /* ! defined _CHORDINPUT_H_ */
