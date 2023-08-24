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
 * \file LevelCNSOp.H
 *
 * \brief AMRLevel operator for compressible Navier-Stokes
 *
 *//*+*************************************************************************/

#ifndef _LEVELCNSOP_H_
#define _LEVELCNSOP_H_

#include <atomic>
#include <memory>

//----- Internal -----//

#include "PatchCNSOp.H"
#include "PatchMappedFunc.H"
#include "MMBSingleLevel.H"
#include "SpectralUtil.H"
#include "SpectralForcing.H"
#include "NonLinearSolver.H"
#include "ARKUtil.H"

//----- Chombo Library -----//

#include "CoarseAverage.H"
#include "FourthOrderFineInterp.H"

#include "UsingNamespace.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class LevelFluxRegister;
class LevelGridMetrics;
class NewFourthOrderCoordSys;
class FArrayBox;
template <typename T> class LevelData;
struct SubcycleParams;
#include "NamespaceFooter.H"
typedef NewFourthOrderCoordSys BlockCoordSys;
typedef FourthOrderFineInterp Interp;


/*******************************************************************************
 */
///  Level operator for compressible Navier-Stokes
/**
 *//*+*************************************************************************/

/* DESIGN CONSIDERATIONS
 * - LevelGridMetrics should really be AMRLevelGridMetrics.  We are going to
 *   start building this class with full LevelGridMetrics but really we should
 *   only have access to the data m_N and m_J.  So try to avoid using LGM as
 *   much as is possible.
 * - Multiblock operators should be removed from LevelGridMetrics since they
 *   are mostly single-level operators.
 *   This is mostly a convenience class.  As a set of operators, data should
 *   not be defined here.  Data should arrive through function arguments or
 *   set as references during construction.
 * - We are trying to define data that we cannot get through function arguments
 *   as references.  So once AMRLevel defines them, they are defined here as
 *   well.  Examples include LGM and the DBL
 */

class LevelCNSOp
{


/*==============================================================================
 * Flags to pass for computing terms of the NSE
 *============================================================================*/

public:
  enum Terms
    {
     None = 0, // Don't compute any of the terms (???)
     NonStiff = 1, // Compute the nonstiff terms (inviscid, viscous, linear source)
     Stiff = 2 // Compute the stiff terms (nonlinear reacting source)
    };

public:

  typedef LevelData<FArrayBox> RHS;
  typedef LevelData<FArrayBox> SOLN;

  /// Cache giving a partitioned box for a particular box size and Num. threads
  using BoxPartCache = ComputeCacheUMap<
    stc::IVec<SpaceDim+1>,                                     // Key
    std::pair<IntVect, IntVect>,                               // Cached value
    std::function<std::pair<IntVect, IntVect>(const IntVect&,  // Box size
                                              const int,       // Num threads
                                              const int)>>;    // Min part size

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  LevelCNSOp(const DisjointBoxLayout& a_boxes,
             LevelGridMetrics&        a_levelGridMetrics,
             LevelData<FArrayBox>&    a_U,
             const Copier&            a_UExchangeCopier,
             const Copier&            a_JUExchangeCopier,
             const bool&              a_hasCoarserGrid,
             const bool&              a_hasFinerGrid,
             BoxPartCache&            a_boxPartCache);

  // Use synthesized destructor


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  LevelCNSOp(const LevelCNSOp&);

//--Assignment not permitted

  LevelCNSOp &operator=(const LevelCNSOp&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Define the level operator (weak construction)
  void define(const int       a_level,
              const RealVect& a_dx);

  /// Define the level operator (weak construction) in reverse traversal
  void reverseDefine();

  /// Sets the time interpolator for this level operator.
  //  Note:  Unlike most const references herein, this cannot be set in the
  //  constructor because it is difficult to have it not be an rvalue which
  //  leads to taking a reference to a temporary when m_levelOp is constructed.
  void setTimeInterpolatorPtr(const TimeInterpolatorRK4* a_timeInterpolator);

  /// Find unit normals for applying the Riemann problem on mapped grids
  void defineUnitNormals(LevelData<FluxBox>& a_NLev);

  ///  Define data structures for computing SGS kinetic energy estimate
  void defineSGSKE();

  /// Compute new timestep from latest estimate of the solution
  Real computeNewDt(const LevelData<FArrayBox>& a_U,
                    const Real&                 a_time,
                    const Real&                 a_currentDt,
                    const Real&                 a_prevStepDt) const;

  /// Add artificial viscosity to a_Unew
  void addArtificialViscosity(
    LevelData<FArrayBox>&       a_JUnew,
    LevelFluxRegister&          a_finerFluxRegister,
    LevelFluxRegister&          a_coarserFluxRegister,
    const LevelData<FArrayBox>& a_Uold,
    const LevelData<FArrayBox>& a_WOld,
    const Real                  a_weight,
    const Real                  a_time) const;

//--LevelRK4 required methods

  /// Evaluate a_RHS (which is d(JU)/dt) at the current time based on a_JU
  void evalRHS(RHS&               a_RHS,
               SOLN&              a_JU,
               const int          a_stage,
               const Real         a_stageTime,
               const Real         a_stageWeight,
               LevelFluxRegister& a_fnFluxRegister,
               LevelFluxRegister& a_crFluxRegister,
               const Real         a_dt,
               const Real         a_timeOld,
               const Real         a_crTimeOld,
               const Real         a_crTimeNew,
               SubcycleParams     a_subcycleParams,
               const int          a_termFlags,
               SOLN&              a_WOld);

  /// Update the solution (a_JU += dt*a_RHS)
  void updateODE(SOLN&       a_JU,
                 const RHS&  a_RHS,
                 Real        a_dt,
                 const int   a_stage = -1,
                 bool        a_toNewSoln = true,
                 const SOLN* a_oldSoln = nullptr) const;

  /// Define a_newJU to match a_JU, including ghost cells
  void defineSolnData(SOLN&       a_newJU,
                      const SOLN& a_JU) const;

  /// Define type a_newRHS based on a_JU, including required ghost cells
  void defineRHSData(RHS&        a_newRHS,
                     const SOLN& a_JU) const;

  /// Copy a_srcJU to a_dstJU
  void copySolnData(SOLN&       a_dstJU,
                    const SOLN& a_srcJU) const;

  /// Species correction
  void speciesCorrection(LevelData<FArrayBox>& a_JUnew) const;

  /// Exchange operation for restart
  void restartExchange(LevelData<FArrayBox>& a_U) const;

  // Increment operation, does a_lhs += a_scale*a_rhs
  void increment(SOLN&       a_lhs,
                 const RHS&  a_rhs,
                 const Real& a_scale);

  // Solve the Ax=b system in ARK time stepping
  void solve(SOLN&              a_newSoln,
             const SOLN&        a_prevStageSoln,
             const SOLN&        a_prevTimeSoln,
             const RHS&         a_rhs,
             int                a_stage,
             Real               a_time,
             Real               a_stageweight,
             LevelFluxRegister& a_finerFluxRegister,
             LevelFluxRegister& a_coarserFluxRegister,
             const Real         a_dt,
             const Real         a_timeOld,
             const Real         a_timeCoarseOld,
             const Real         a_timeCoarseNew,
             SubcycleParams     a_subcycleParams,
             SOLN&              a_WOld);

  // The terms that will be treated as non-stiff by the ARK time stepping method
  void evalNonStiff(RHS&               a_RHS,
                    SOLN&              a_JU,
                    const int          a_stage,
                    const Real         a_stageTime,
                    const Real         a_stageWeight,
                    LevelFluxRegister& a_fnFluxRegister,
                    LevelFluxRegister& a_crFluxRegister,
                    const Real         a_dt,
                    const Real         a_timeOld,
                    const Real         a_crTimeOld,
                    const Real         a_crTimeNew,
                    SubcycleParams     a_subcycleParams,
                    SOLN&              a_WOld);

  // The terms that will be treated as stiff by the the ARK time stepping method
  void evalStiff(FArrayBox&         a_RHS,
                 const FArrayBox&   a_JU,
                 const Box&         a_box,
                 const Real         a_stageTime,
                 const DataIterator a_dit,
                 FArrayBox&         a_WOld);

  // The terms that will be treated as stiff by the ARK time stepping method
  void evalStiff(RHS&               a_RHS,
                 SOLN&              a_JU,
                 const int          a_stage,
                 const Real         a_stageTime,
                 const Real         a_stageWeight,
                 LevelFluxRegister& a_fnFluxRegister,
                 LevelFluxRegister& a_crFluxRegister,
                 const Real         a_dt,
                 const Real         a_timeOld,
                 const Real         a_crTimeOld,
                 const Real         a_crTimeNew,
                 SubcycleParams     a_subcycleParams,
                 SOLN&              a_WOld);

  // Calculate the source term Jacobian for a single cell
  void
  calcRxnJacobian(FArrayBox&          a_reactionJacobian,
                  const FArrayBox&    a_JU,
                  const Box&          a_box,
                  const Real          a_stageTime,
                  const Real          a_dt,
                  const DataIterator& a_dit,
                  FArrayBox&          a_WOld);

  /// Get number of components in flattening array (for plotting)
  int numDebugPlotVar() const
    { return c_numDebugPlotVar; }

  /// Get flattening coefficients/debugging variable for plotting
  const LevelData<FArrayBox>& getDebugVar() const
    { return m_flattening; }

  /// Get SGS KE estimate (for plotting up-to-date instantaneous value)
  const LevelData<FArrayBox>& getSGSKineticEnergyVar() const
    { return m_sgske; }

  /// Get J*SGS KE estimate (for plotting up-to-date instantaneous value)
  const LevelData<FArrayBox>& getJSGSKineticEnergyVar() const
    { return m_Jsgske; }

  /// Get time-averaged variables for plotting
  const LevelData<FArrayBox>& getTimeAvgVar() const
    { return m_timeAvgData; }

  /// Get face-averaged variables for plotting
  const LevelData<FluxBox>& getFaceAvgVar() const
    { return m_faceAvgPlotData; }

  /// Get face-averaged time-averaged variables for plotting
  const LevelData<FluxBox>& getFaceAvgTimeAvgVar() const
    { return m_faceAvgTimeAvgPlotData; }

  /// Get face unit-normals outside LevelCNSOp for initializing turb variables
  const LayoutData<FluxBox>& getUnitNormals() const
    {
      CH_assert(m_unitNormalsDefined);
      return m_unitNormals;
    }

//--Functions related to modeling

  /// Compute global sums for turbulence forcing
  void computeTurbForcingSums(
    const LevelData<FArrayBox>& a_WcellAvgLvl,
    Real&                       a_globalKE,
    Real&                       a_globalHelicity) const;

  /// Set wall-model data from finer level
  void setWallModelFromFiner(LevelData<FArrayBox>& a_JU);

  /// Compute coarsest level LES SGS KE for single-block problems
  void estimateCrsLevSingleBlockSGSKE(LevelData<FArrayBox>& a_JU);

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
  /// Add spectral forcing perturbations to solution field
  Real addSpectralForce(LevelData<FArrayBox>&    a_U,
                        LevelData<FArrayBox>&    a_Source,
                        Real                     a_dt);
#else
  void addSpectralForce(LevelData<FArrayBox>&    a_U,
                        LevelData<FArrayBox>&    a_Source,
                        Real                     a_dt)
    { }
#endif

  //-- Functions related to ARK

  ///  Compute the max norm of JU-\hat{JU}
  Real computeErrorEstimate(const LevelData<FArrayBox>& a_JU,
                            const LevelData<FArrayBox>& a_JUhat) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Returns cached guidelines on how to partition a box
  std::pair<IntVect, IntVect> cachedBoxPart(const Box& a_box,
                                            const int  a_numPart);

  /// Fill in ghost cells of m_U to RK4 intermediates
  void fillGhostsRK4AndComputeU(LevelData<FArrayBox>& a_JU,
                                const int             a_stage,
                                const Real            a_dt,
                                const Real            a_timeOld,
                                const Real            a_crTimeOld,
                                const Real            a_crTimeNew,
                                SubcycleParams        a_subcycleParams);


/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  Real computeFirstDt(const LevelData<FArrayBox>& a_U,
                      const Real&                 a_time) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

//--References to data defined in AMRLevelCNS.  These are all set during
//--construction and are fully defined by AMRLevelCNS.

  const DisjointBoxLayout& m_boxes;   ///< Boxes used on the level
  LevelGridMetrics& m_levelGridMetrics;
                                      ///< Grid metrics for the level.  Non-
                                      ///< const for timeIntermediate.
  //**This is not available through LevelRK4 but a switch to DenseRK4Output may
  //**fix it so we can pass m_U as an argument to evalRHS.
  LevelData<FArrayBox>& m_U;          ///< Conservative state in physical space
  const Copier& m_UExchangeCopier;    ///< Copier for \<U\>
  const Copier& m_JUExchangeCopier;   ///< Copier for \<JU\>
  const TimeInterpolatorRK4* m_timeInterpolator;
                                      ///< RK4 time interpolator

  // Intermediate data on a stage
  LevelData<FArrayBox> m_WcellPntLvl; ///< Point values of the full primitive
                                      ///< state in all cells
  LevelData<FArrayBox> m_WcellAvgLvl; ///< Average values of the full primitive
                                      ///< state in all cells
  LevelData<FluxBox> m_WfaceAvgLvl;   ///< Average values of the full primitive
                                      ///< state on all faces

  // These describe whether or not a level is actually used (i.e., has a
  // non-zero number of cells) rather than just if those levels were defined.
  const bool& m_hasCoarserGrid;       ///< T - A coarser level exists and
                                      ///<     is in use (has grids)
  const bool& m_hasFinerGrid;         ///< T - A finer level exists and
                                      ///<     is in use (has grids)
  // Cache for using threads
  BoxPartCache& m_boxPartCache;       ///< A cache of guidelines for
                                      ///< partitioning a box into a specific
                                      ///< number of partitions.

//--New data for this class

  PatchCNSOp m_patchOp;               ///< CNS patch operator

  LayoutData<FluxBox> m_unitNormals;  ///< Unit normal basis for applying
                                      ///< Riemann solvers on mapped grids
  LevelData<FArrayBox> m_flattening;  ///< Flattening coefficients, if used.
  LevelData<FArrayBox> m_timeAvgData; ///< Time-averaged data, if used.
  LevelData<FluxBox> m_faceAvgPlotData;
                                      ///< If used, gathers instantaneous
                                      ///< face-averaged data for plotting
  LevelData<FluxBox> m_faceAvgTimeAvgPlotData;
                                      ///< If used, gathers time-averaged
                                      ///< face-averaged data for plotting
  MMBSingleLevel m_MMBSingleLevelOp;  ///< Encapsulated mapped-multiblock
                                      ///< operations for a single level

  RealVect m_dx;                      ///< Mesh spacing.  Somewhat ambiguous for
                                      ///< multiblock.
  int m_stage;                        ///< Index of the RK stage
  int m_level;                        ///< Index of the level
  mutable Real m_minConvDt;           ///< Local minimum convective time step
  mutable Real m_minDiffDt;           ///< Local minimum diffusive time step
  mutable Real m_minChemDt;           ///< Local minimum chemical time step
  mutable Real m_minLocalDt;          ///< Local minimum time step
  mutable FixedQueue<Real, 3> m_UminusUhat;
                                      ///< Stores the last three U - \hat{U}
                                      /// values used for PID step-size control
  mutable IntVect m_minConvDtCell;    ///< Cell with local minimum convective
                                      ///< time step
  mutable IntVect m_minDiffDtCell;    ///< Cell with local minimum diffusive
                                      ///< time step
  mutable IntVect m_minChemDtCell;    ///< Cell with local minimum chemical time
                                      ///< step
  mutable IntVect m_minPidDtCell;     ///< Cell with local minimum PID-controller
                                      ///  time step
  mutable int m_cnlim;                ///< Species that limits the time step
                                      ///< size when doing chemical reactions
  std::atomic_flag m_lockMinLocalDt;  ///< For threadsafe access to
                                      ///< m_minLocalDt
  mutable bool m_firstDt;             ///< Bool to tell if this is the first
                                      ///< calculation of dt
  Real m_prevDt;                      ///< Time-step size from previous step
  bool m_unitNormalsDefined;          ///< True if unitNormals has been defined
  mutable bool m_firstTau;            ///< Bool to tell if it is first
                                      ///< chemical time step

  unique_ptr<NonLinearSolverI> m_nlsolver;
                                      ///< The nonlinear solver for
                                      ///< semi-implicit time marching methods

//--Data and operators for coarsened LES SGS model calculations

  LevelData<FArrayBox> m_currModelJU; ///< JU on this level for model calc
  LevelData<FArrayBox> m_partJWcellAvgLvl;
                                      ///< Cell-avg J*density + J*velocity on
                                      ///< this level
  LevelData<FArrayBox> m_cellAvgJWGradNtJLvl;
                                      ///< Cell-avg J*grad W N^t/J on this level
  LevelData<FArrayBox> m_cellPntJ;    ///< Cell-centered J on the current block
  LevelData<FArrayBox> m_crsCellPntJ; ///< Cell-centered J on coarser mesh
  LevelData<FArrayBox> m_crsXLvl;     ///< Cell-centered X on coarser mesh
  LevelData<FArrayBox> m_crsDeltaC;   ///< Cell-centered delta-c on crsr mesh
  LevelData<FArrayBox> m_sgske;       ///< Plot sgs kinetic energy
  LayoutData<FluxBox>  m_faceCoord;   ///< Face-centered coordinates
                                      ///< This is a temporary speedup
  LevelData<FArrayBox> m_Jsgske;      ///< Plot sgs kinetic energy
  LayoutData<FluxBox>  m_facePntDeltaC;
                                      ///< Face-centered J - temporary speedup
  LevelData<FArrayBox> m_crsCellAvgJW;///< Coarsened <JU> from the base mesh
  LevelData<FArrayBox> m_crsCellAvgJPhysGrad;
                                      ///< Coarsened base mesh <J grad W N^t/J>
  LevelData<FArrayBox> m_crsCellAvgJSGSKE;
                                      ///< Coarsened <J-SGSKE> from base mesh
  LevelData<FArrayBox> m_fnCellAvgSGSKE;
                                      ///< The <SGSKE> interpolated to base mesh
  LevelData<FArrayBox> m_fnCellAvgJSGSKE;
                                      ///< The <J-SGSKE> interp to base mesh

  CoarseAverage m_cellAvgJWAvgOp;     ///< Average down <JW> for SGS KE estimate
  CoarseAverage m_cellAvgJPhysGradWAvgOp;
                                      ///< Average down <J grad W N^t/J> for
                                      ///< SGS KE estimate
  Copier m_crsCellAvgJWExchange;      ///< Exchange copier for <JU> and <J>
  Copier m_crsCellAvgSGSKEExchange;   ///< Exchange copier for coarse <SGSKE>
  Copier m_fnCellAvgSGSKEExchange;    ///< Exchange copier for fine <SGSKE>
  Copier m_fnCellAvgJSGSKEExchange;   ///< Exchange copier for fine <J-SGSKE>
  Copier m_crsCellAvgJGradWExchange;  ///< Exchange copier for crs
                                      ///< <J grad W N^t/J>
  Interp m_crsToFnSGSKEInterpolator;  ///< Coarse-to-fine <SGSKE> interpolator
  Real   m_totalTimeAfterFilterInit;  ///< Total time after initialization of
                                      ///< sharp time-averaging filter

//--Data for temporally-evolving characteristic boundary conditions
  LayoutData<stc::Vector<FArrayBox, 2*SpaceDim>> m_bndryCellData;
                                      ///< Boundary cell-averaged data at
                                      ///< previous time-step

//--NtJ for wall-model calculations only along wall-boundaries
  LayoutData<stc::Vector<FArrayBox, 2*SpaceDim>> m_bndryNtJ;
                                      ///< NtJ at interior face for wall-model
                                      ///< Necessary because NtJ doesn't exist
                                      ///< out to 4 ghost faces or cells and
                                      ///< LGM.m_N tangential components only
                                      ///< exist out to 3 ghost cells

//--Special operators for specific cases

public:

  SpectralFilter m_spectralFilter;    ///< Spectrally filter the state data.
                                      ///< Used for research into LES.  Requires
                                      ///< single-level periodic rectangular
                                      ///< domain.
  SpectralForcing m_spectralForcing;  ///< Used for spectrally-truncated
                                      ///< deterministic forcing.
                                      ///< Requires periodic rectangular
                                      ///< domain.
  Real m_spectralRefreshTime;
  Real m_globalKE;                    ///< Domain sum of kinetic energy for
                                      ///< turbulence forcing
  Real m_globalHelicity;              ///< Domain sum of helicity for
                                      ///< turbulence forcing

protected:

  // Extra plot variables used for debugging
  enum : int
  {
    c_numDebugPlotVar = 0
  };
};


/*******************************************************************************
 *
 * Class LevelCNSOp: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Returns cached guidelines on how to partition a box
/** \param[in]  a_box   Box to partition
 *  \param[in]  a_numPart
 *                      Number of partitions
 *  \return             A pair of IntVects.  The first is the size
 *                      of each sub-box.  The second is the number of
 *                      partitions in each direction that were made to
 *                      a_box.
 *//*-----------------------------------------------------------------*/

inline std::pair<IntVect, IntVect>
LevelCNSOp::cachedBoxPart(const Box& a_box, const int  a_numPart)
{
  // Make a key { a_box.size(), a_numPart }
  const stc::IVec<SpaceDim+1> key = stc::IVec<SpaceDim+1>(a_box.size()) +
    a_numPart*stc::make_IVec<SpaceDim+1>::basis(SpaceDim);
  return m_boxPartCache(key, a_box.size(), a_numPart, CRDparam::g_minBoxSize);
}

#endif  /* ! defined _LEVELCNSOP_H_ */
