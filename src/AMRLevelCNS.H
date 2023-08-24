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
 * \file AMRLevelCNS.H
 *
 * \brief AMRLevel operator for compressible Navier-Stokes
 *
 *//*+*************************************************************************/

#ifndef _AMRLEVELCNS_H_
#define _AMRLEVELCNS_H_

//----- Chombo Library -----//

#include "CoarseAverage.H"
#include "AMRLevel.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "LevelARK4.H"
#include "UsingNamespace.H"
#include "NamespaceVar.H"

//----- Internal -----//

#include "LevelCNSOp.H"

//-- Forward declarations from Chombo

#include "NamespaceHeader.H"
class LevelGridMetrics;
class MultiBlockCoordSys;
class MultiBlockCoordSysFactory;
class LevelFluxRegister;
class TimeInterpolatorARK4;
class TimeInterpolatorRK4;
class TimeInterpolatorRK2;
#include "NamespaceFooter.H"

//-- Forward declarations

class AMRLevelCNSFactory; 
class TagLevel;


/*******************************************************************************
 */
///  AMRLevel operator for compressible Navier-Stokes
/**
 *   This class derives from AMRLevel and performs the necessary AMR operations
 *   to allow for operations on a level.
 *
 *//*+*************************************************************************/

class AMRLevelCNS : public CH_XDIR::AMRLevel
{

//--Friend classes

  friend class AMRLevelCNSFactory;
  friend class CH_XDIR::LevelGridMetrics;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default/factory constructor
  AMRLevelCNS(const AMRLevelCNSFactory *const a_factory);

  /// Destructor
  virtual ~AMRLevelCNS();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  AMRLevelCNS(const AMRLevelCNS&);

//--Assignment not permitted

  AMRLevelCNS &operator=(const AMRLevelCNS&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Define the new AMR level (weak construction)
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio) override;

  /// Provide a vector of block domains for the level
  virtual std::pair<int, const Box*> blockDomainVector() const override;

  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newBoxes) override;

  /// Compute the grid metrics for mapped grids
  virtual void postInitialGrid(const bool a_restart) override;

  /// Initialize data
  virtual void initialData() override;

  /// Things to do after initialization
  virtual void postInitialize() override;

  /// Advance by one timestep and return new timestep.
  virtual Real advance() override;

//   /// add dissipation to m_Unew:  Kreiss-Oliger artificial viscosity
// //   virtual void addDissipation();

  /// Things to do after a timestep -- reflux
  virtual void postTimeStep() override;

  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) override;

  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) override;

  /// Pre-regrid operations necessary for updating metric terms
  virtual void preRegrid(int                         a_baseLevel,
                         const Vector<Vector<Box> >& a_newBoxes) override;

  /// Fine-to-coarse traversal after a regrid
  virtual void postRegrid(int a_baseLevel) override;

  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newBoxes) override;

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const override;

  /// Write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const override;

  /// Read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle) override;

  /// Read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle) override;

  /// Write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const override;

  /// Write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const override;

  /// Write mapped-grid info
  virtual void writeMappedPlotFile() const;
#endif

  /// Conclude by testing conservation and reporting error if necessary
  virtual void conclude(int a_step) const override;

  /// Compute dt
  virtual Real computeDt() override;

  /// Compute dt using initial data
  virtual Real computeInitialDt() override;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Setup menagerie of data structures and operators for the level
  void levelOpSetup();

  /// Reverse-order setup of data structures and operators for the level
  void reverseLevelOpSetup();

  /// Do we have a pointer to a coarser level?
  bool hasCoarserPtr() const;

  //  After levelOpSetup is called, you can use the flag m_hasCoarserGrid to
  //  see if the coarser level is defined and has grids

  /// Get the next coarser level
  //  Returns level if ptr exists, level may not be used
  AMRLevelCNS* getCoarserLevel() const;

  /// Do we have a pointer to a finer level?
  bool hasFinerPtr() const;

  //  After levelOpSetup is called, you can use the flag m_hasFinerGrid to
  //  see if the finer level is defined and has grids

  /// Get the next finer level
  //  Returns level if ptr exists, level may not be used
  AMRLevelCNS* getFinerLevel() const;

  /// Compute a sum of the components over all levels
  void computeSum(Vector<Real>& a_sum) const;

  // /// Compute and write out norm of conserved quantities
  // void reportNorm() const;

  /// Report norms of errors (requires exact solution)
  void reportError() const;

  /// Initialize argument <U> and <JU> at time m_time using IBC class
  void loadIBCData(LevelData<FArrayBox>& a_U,
                   LevelData<FArrayBox>& a_JU);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

//--Parameters initialized by factory constructor

  const AMRLevelCNSFactory *const m_factory;
                                      ///< We hold our own factory object
                                      ///< so we can use it during weak
                                      ///< construction (ideally this would be
                                      ///< held in AMR and provided to the
                                      ///< define function).
                                      ///< Note: the factory object gives itself
                                      ///< to the new_amrlevel routine and we
                                      ///< then store it here.  This is weird
                                      ///< but there is not much choice.  The
                                      ///< pointer to the factory object should
                                      ///< be maintained in the driver code.
  int m_numGhost;                     ///< Number of ghost cells (in each
                                      ///< direction).  Currently fixed to 5
                                      ///< for hyperbolic algorithm.
  IntVect m_ghostVect;                ///< Ghost cell vector.
  RealVect m_dx;                      ///< Mesh spacing.  Somewhat ambiguous for
                                      ///< multiblock.
  int m_levelStep;                    ///< Time step index for the level
  int m_begFinerLevelStep;            ///< Time step index for the finer level
                                      ///< at beginning of the coarser time
                                      ///< step.  Used to detect reduction in
                                      ///< time step on finer levels which could
                                      ///< indicate an instability
  LevelCNSOp::BoxPartCache m_boxPartCache;
                                      ///< For a given number of threads, this
                                      ///< cache stores the breakdown of any
                                      ///< box.  This cache persists for the
                                      ///< duration of the run

//--Parameters initialized by define

  RefCountedPtr<MultiBlockCoordSys> m_multiBlockCoordSys;
                                      ///< A pointer to the instance of the
                                      ///< coordinate system for this level.
  RefCountedPtr<TagLevel> m_tagLevel; ///< Method for tagging cells on a level
  LevelGridMetrics m_levelGridMetrics;
                                      ///< Grid metrics for the level
  MultiBlockUtil m_mbUtil;            ///< Multiblock utilities

//--Parameters initialized by initialGrid

  DisjointBoxLayout m_boxes;          ///< Boxes used on the level
  LevelData<FArrayBox> m_U;           ///< \<U\> is placed here so that we
                                      ///< can take a reference to it early in
                                      ///< construction and pass it to the level
                                      ///< operator.
  LevelData<FArrayBox> m_WOld;        ///< Old primitive data, used to stabilize
                                      ///  consToPrim
  MappedLevelData m_data;             ///< Level data (conserved variables
                                      ///< \<JU\> and \<U\>)

//--Parameters initialized by levelOpSetup

  // These describe whether or not a level is actually used (i.e., has a
  // non-zero number of cells) rather than just if those levels were defined.
  bool m_hasCoarserGrid;              ///< T - A coarser level exists and
                                      ///<     is in use (has grids)
  bool m_hasFinerGrid;                ///< T - A finer level exists and
                                      ///<     is in use (has grids)
  // Averaging from coarse to fine level
  CoarseAverage m_averageOp;          ///< Averages from a finer level to
                                      ///< this (coarser) level
  LevelFluxRegister* m_fluxRegisterPtr;
                                      ///< Multiblock-capable flux register
                                      ///< between this (coarser) level and a
                                      ///< finer level
  // Interpolation in time
  TimeInterpolatorRK4* m_timeInterpolator;
                                      ///< The time interpolator for
                                      ///< interpolating coarse data to a fine
                                      ///< level

  // spectral forcing for HIT and spectral-nudging eDNS research simulations
  LevelData<FArrayBox> m_spectralSource;
                                      ///< Persistent memory for spectral
                                      ///< forcing

//--May reference previous values

  LevelCNSOp m_levelOp;               ///< Level operator (ideally no AMR
                                      ///< details)

  LevelARK4<LevelData<FArrayBox>,     ///< The Additive Runge Kutta time
            LevelData<FArrayBox>,     ///< integrator
            TimeInterpolatorARK4,
            LevelCNSOp> m_ark_time_integrator;

//--Miscellaneous

  Real m_dtNew;                       ///< Time step estimate from latest
                                      ///< solution state
  Real m_prevStepDt;                  ///< Time step size from previous step
  bool m_exchangeForTagging;          ///< Necessary for when a restart occurs
                                      ///< on a step where tagging occurs for
                                      ///< finer levels
  bool m_justRegridded;               ///< True when a regrid was just performed on
                                      ///  this level. Reset to false after a step

//--Testing/validation

  static Vector<Real> s_JUConsvRef;   ///< Sum of JU at start of solution
                                      ///< used for checking conservation
};


/*******************************************************************************
 *
 * Class AMRLevelCNS: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Do we have a pointer to a coarser level?
/** This only finds the bounds on levels.  The coarser level may not
 *  be used.  If you need to ensure the level is defined and has
 *  grids, use AMRLevel::hasCoarserLevel() instead.
 *  \return             T - coarser level exists
 *//*-----------------------------------------------------------------*/

inline bool
AMRLevelCNS::hasCoarserPtr() const
{
  return (m_coarser_level_ptr != NULL);
}

/*--------------------------------------------------------------------*/
//  Do we have a pointer to a finer level?
/** This only finds the bounds on levels.  The finer level may not
 *  be used.  If you need to ensure the level is defined and has
 *  grids, use AMRLevel::hasFinerLevel() instead.
 *  \return             T - finer level exists
 *//*-----------------------------------------------------------------*/

inline bool
AMRLevelCNS::hasFinerPtr() const
{
  return (m_finer_level_ptr != NULL);
}

#endif  /* ! defined _AMRLEVELCNS_H_ */
