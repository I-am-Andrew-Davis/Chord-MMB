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
 * \file AMRLevelCNS.cpp
 *
 * \brief Member functions for AMRLevelCNS
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "MultiBlockCoordSys.H"
#include "MultiBlockFluxRegister.H"
#include "LoadBalance.H"
#include "computeSum.H"
#include "computeNorm.H"
#include "FourthOrderUtil.H"
#include "TimeInterpolatorRK4.H"
#include "TimeInterpolatorARK4.H"
#include "TimeInterpolatorRK2.H"
#include "LevelRK4_v2.H"
#include "LevelRK2.H"
#include "MOLUtilities.H"
#include "NodeAMRIO.H"
#include "CONSTANTS.H"  /* DBG REMOVE WITH FFTW INIT */

//----- Internal -----//

#include "AMRLevelCNS.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "AMRLevelCNSFactory.H"
#include "CRDPhysics.H"
#include "CNSIBC.H"
#include "TagLevel.H"
#include "LevelMappedFunc.H"
#include "DataTemp.H"
#include "SetCentersF_F.H"


/*******************************************************************************
 *
 * Class AMRLevelCNS: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Function prototypes for functions that are bound
 *============================================================================*/

/*--------------------------------------------------------------------*/
/// Partition a box into sub-boxes
/** \param[in]  a_boxSize
 *                      Size of the box to partition
 *  \param[in]  a_numPart
 *                      Desired number of partitions
 *  \param[in]  a_minPartSize
 *                      Minimum size of a sub-box in a direction
 *  \return             A pair of IntVects.  The first is the size
 *                      of each sub-box.  The second is the number of
 *                      partitions in each direction that were made to
 *                      a_boxSize.
 *  \note
 *  <ul>
 *    <li> a_boxSize must be evenly divisible by a_minPartSize
 *  </ul>
 *//*-----------------------------------------------------------------*/

std::pair<IntVect, IntVect> partitionBox(const IntVect& a_boxSize,
                                         const int      a_numPart,
                                         const int      a_minPartSize)
{
  if (a_numPart == 1 || a_boxSize.product() == 0)
    {
      return std::make_pair(a_boxSize, IntVect_unit);
    }
  else
    {
      IntVect splitDim = a_boxSize/a_minPartSize;
      // Enforce sanity for now...
      CH_assert(splitDim*a_minPartSize == a_boxSize);
      while (splitDim.product() > a_numPart)
        {
          IntVect partSize = a_boxSize/splitDim;
          // Neglect any direction where partSize == a_boxSize
          for (const int d : EachDir)
            {
              if (splitDim[d] == 1)
                {
                  partSize[d] = std::numeric_limits<int>::max();
                }
            }
          const int dir = stc::indexMinElem(partSize);
          splitDim[dir] /= 2;
        }
      const IntVect partSize = a_boxSize/splitDim;
      CH_assert(splitDim*partSize == a_boxSize);
      CH_assert(a_numPart >= splitDim.product());
      return std::make_pair(partSize, splitDim);
    }
}

/*==============================================================================
 * Definition of static variables
 *============================================================================*/

Vector<Real> AMRLevelCNS::s_JUConsvRef;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default/factory constructor
/** \param[in]  a_factory
 *                      Pointer to the factory object which is being
 *                      used to construct this class (!)
 *  This constructor is called by the AMRLevelCNSFactory object.
 *//*-----------------------------------------------------------------*/

AMRLevelCNS::AMRLevelCNS(const AMRLevelCNSFactory *const a_factory)
  :
  m_factory(a_factory),
  m_numGhost(CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg)),
  m_ghostVect(CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg)*IntVect::Unit),
  // Hard-code 4th order spatial accuracy
  m_dx(RealVect::Unit),
  m_levelStep(0),
  m_boxPartCache(std::function<
                 std::pair<IntVect, IntVect>(const IntVect&,
                                             const int,
                                             const int)>(partitionBox)),
  m_levelGridMetrics(CRDparam::g_CRDPhysics->numConservative(), 4),
  m_averageOp(),
  m_fluxRegisterPtr(nullptr),
  m_timeInterpolator(nullptr),
  m_levelOp(m_boxes,             // All arguments are just a reference grab
            m_levelGridMetrics,  // and do not necessarily pass information yet.
            m_U,
            m_data.refUExchangeCopier(),
            m_data.refJUExchangeCopier(),
            m_hasCoarserGrid,
            m_hasFinerGrid,
            m_boxPartCache),
  m_dtNew(-1.),
  m_prevStepDt(-1.),
  m_exchangeForTagging(true),
  m_justRegridded(false)
{ }

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

AMRLevelCNS::~AMRLevelCNS()
{
  if (m_fluxRegisterPtr != nullptr)
    {
      delete m_fluxRegisterPtr;
    }
  if(m_timeInterpolator != nullptr)
    {
      delete m_timeInterpolator;
    }
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the new AMR level (weak construction)
/** \param[in]  a_coarserLevelPtr
 *                      Pointer to next coarser level object.
 *  \param[in]  a_problemDomain
 *                      Problem domain of this level
 *  \param[in]  a_level Index of this level.  The base level is zero
 *  \param[in]  a_refRatio
 *                      The refinement ratio between this level and
 *                      the next finer level
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::define(AMRLevel*            a_coarserLevelPtr,
                    const ProblemDomain& a_problemDomain,  // Empty
                    int                  a_level,
                    int                  a_refRatio)
{
  CH_TIME("AMRLevelCNS::define");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::define (level: " << a_level
           << ")" << CRD::end;

  // Define the base class
  // Defines m_coarser_level_ptr, m_problem_domain, m_level, m_ref_ratio,
  // m_finer_level_ptr(NULL), m_isDefined(true).  Pointers are null if
  // at extents of AMR hierarchy.
  // WARNING: m_problem_domain must never be used!  Use
  //   m_multiBlockCoordSys->problemDomain(Box) or
  //   m_multiBlockCoordSys->levelDomain() instead.
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,  // Empty
                   a_level,
                   a_refRatio);

  AMRLevel::verbosity(CRDparam::g_verbosity);

  m_factory->define_amrlevel(a_level,               // I
                             a_refRatio,            // I
                             m_dx,                  // O
                             m_multiBlockCoordSys,  // O
                             m_tagLevel);           // O

  // Get setup information from the next coarser level
  LevelGridMetrics* coarserLGMptr = NULL;
  if (hasCoarserPtr())
    {
      coarserLGMptr = &(getCoarserLevel()->m_levelGridMetrics);
    }

  // Define the flux register class
  //**FIXME Can this have a define so as to not be a pointer?
  // m_fluxRegisterPtr = new MultiBlockFluxRegister(m_multiBlockCoordSys,
  //                                                m_coordSysFine);

  //**FIXME Can these share a coordinate system? The only trick is LGM deletes
  //        its CS on destruction

  // Define the grid metrics class
  //**FIXME do we want access to cell NtJ all the time?
  //**FIXME why do we need multiblock vector data? Exchanges of U all happen
  //        in physical space and velocity is independent of the grid
  const int numCellCoordGhosts = 1;
  m_levelGridMetrics.define(
    this,
    m_multiBlockCoordSys.operator->(),
    coarserLGMptr,
    m_dx,
    m_ghostVect,
    LevelGridMetrics::TransverseNOpAverage,
    false, // (CRDparam::g_CRDPhysics->velocityInterval().size() > 0),
    true,  // solve for m_cellNtJ for cell centered mappings
    numCellCoordGhosts // number of ghost cells for cell-centered coordinates
    );
  // set clipping options
  m_levelGridMetrics.defineInterpOptions(CRDparam::g_clipping, CRDparam::g_clippingHO, CRDparam::g_clippingPostSmooth);

  // Define the multiblock utilities class
  m_mbUtil.define(m_multiBlockCoordSys.operator->());
}

/*--------------------------------------------------------------------*/
//  Provide a vector of block domains for the level
/** \return             The number of blocks and a box representing
 *                      each domain.
 *//*-----------------------------------------------------------------*/

std::pair<int, const Box*>
AMRLevelCNS::blockDomainVector() const
{
  return std::make_pair(
    m_multiBlockCoordSys->numBlocks(),
    m_multiBlockCoordSys->mappingBlocks().constStdVector().data());
}

/*--------------------------------------------------------------------*/
//  Initialize grids
/** \param[in]  a_newBoxes
 *                      A vector of boxes defining the new grid level
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::initialGrid(const Vector<Box>& a_newBoxes)
{
  CH_TIME("AMRLevelCNS::initialGrid");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::initialGrid (level: " << m_level
           << ")" << CRD::end;

//--Define the grid

  int numBoxes = a_newBoxes.size();

  // Save the original grids in the base AMRLevel class (required)
  m_level_grids = a_newBoxes;

  // Load balance
  Vector<int> procIDs;
  LoadBalance(procIDs, a_newBoxes);
  if (CRDparam::g_verbosity >= 4)
    {
      CRD::msg << CRD::fv4 << "AMRLevelCNS::initialGrid: processor map:"
               << CRD::end;
      for (int iBox = 0; iBox != numBoxes; ++iBox)
        {
          CRD::msg << CRD::fv4 << "  Box " << iBox << ", proc: "
                   << procIDs[iBox] << ", cells: " << a_newBoxes[iBox].volume()
                   << CRD::end;
        }
    }

  {
    DisjointBoxLayout dbl(a_newBoxes,
                          procIDs,
                          m_multiBlockCoordSys->levelDomain());
    dbl.close();
    m_boxes = dbl;
    // This will abort if a box is not contained in a block
    m_boxes.defineLocalBlocks(m_multiBlockCoordSys->mappingBlocks());
  }

//--More set up for LevelGridMetrics

  // Pass to grid metrics
  m_levelGridMetrics.initialGrid(&m_boxes);

  const int numPrimitive = CRDparam::g_CRDPhysics->numPrimitive();
  const int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostWfromUcellAvg);
  m_WOld.define(m_boxes, numPrimitive, IntVect::Unit*numGhosts);
  // Define a distance field along with grid metrics
  // only need this for turbulence models
  if ((CRDparam::g_turbModelType) &&
      (!(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)))
    {
      // get the boundary cell boxes, and direction of boundary
      Vector<Box> boundaryBoxes;
      Vector<int> faceDir;
      CRDparam::g_CNSIBC->getAllWallBoundaryFaces(boundaryBoxes,
                                                  faceDir,
                                                  0.0,
                                                  m_level,
                                                  m_levelGridMetrics);
      // setup the wall distance function as part of LGM
      m_levelGridMetrics.defineDistances(boundaryBoxes,
                                         faceDir);
    }

//--Define MappedLevelData

  // Once LevelGridMetrics knows the layout, we can define the data
  const LevelData<FArrayBox>* coarseUPtr = 0;
  if (hasCoarserPtr())
    {
      coarseUPtr = &(getCoarserLevel()->m_data.rawU());
    }
  m_data.define(&m_levelGridMetrics,
                &m_U,
                0,  // 0 means m_data will allocate the LevelData
                0,  // 0 means m_data will allocate the LevelData
                coarseUPtr,
                CRDparam::g_CRDPhysics->numConservative(),      // Num states
                CRDparam::g_CRDPhysics->velocityInterval(),  // Velocity loc
                m_numGhost);

//--Allocate the flux register

  //**FIXME Can this have a define so as to not be a pointer?
  if (hasFinerPtr())
    {
      if (m_fluxRegisterPtr != nullptr)
        {
          delete m_fluxRegisterPtr;
        }
      m_fluxRegisterPtr = new MultiBlockFluxRegister(
        m_multiBlockCoordSys.operator->(),
        getFinerLevel()->m_multiBlockCoordSys.operator->());
    }

//--Setup the remaining operators on the level

  levelOpSetup();
}

/*--------------------------------------------------------------------*/
//  Compute the grid metrics for mapped grids
/** postInitialGrid is call after initialGrid and after
 *  readCheckpointLevel if this is a restart.  Consequently grids
 *  should be defined.  postInitialGrid traverses from fine to coarse.
 *  \param[in]  a_restart
 *                      T - This solution is a restart
 *                      F - Initialized solution
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::postInitialGrid(const bool a_restart)
{
  CH_TIME("AMRLevelCNS::postInitialGrid");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::postInitialGrid (level: " << m_level
           << ")" << CRD::end;

  // Define the metric terms
  m_levelGridMetrics.postInitialGrid(NULL);

//--Level setup for fine-to-coarse traversal

  reverseLevelOpSetup();

//--Things to do if this is a restart

  if (a_restart)
    {
      // Normally performed in initialData()

      // If restarting and adding the LES wall-model, initialize the data
      if (CRDparam::g_restartAddWallModel)
        {
          LevelData<FArrayBox>& JU = m_data.rawJU();
          CRDparam::g_CNSIBC->restartWallModel(
            JU, m_levelGridMetrics, m_levelOp.getUnitNormals());
        }

      // Set this level as coarsened-fine data for the interpolator to the next
      // finer level.
      if (m_hasFinerGrid)
        {
          // Find the number of ghost cells on the coarsened-fine mesh required
          // to fill the invalid ghost cells in the fine mesh.
          const IntVect& interpolatorCrFnGhostVect =
            m_levelGridMetrics.interpolatorCrFnNumGhost(true);
          // Assumes ghosts in all directions are equal
          const int numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
          m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
            m_data.getU(0, -1, numInterpolatorCrFnGhost),
            CRDparam::g_CRDPhysics->velocityInterval());
        }

      // Normally performed in postInitialize()

      // This is our conservation reference
      if (m_level == 0)
        {
          s_JUConsvRef.resize(CRDparam::g_CRDPhysics->numConservative());
          computeSum(s_JUConsvRef);
          m_levelOp.restartExchange(m_data.getU(1,1));
          m_exchangeForTagging = false;
        }
      else
        {
          // Exchange must occur before tagging for a restart on finer levels
          m_exchangeForTagging = true;
        }

      // Compute the wall distance if needed
      if ((CRDparam::g_turbModelType) &&
          (!(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)))
        {
          // get the boundary cell boxes, and direction of boundary
          Vector<Box> boundaryBoxes;
          Vector<int> faceDir;
          CRDparam::g_CNSIBC->getAllWallBoundaryFaces(boundaryBoxes,
                                                      faceDir,
                                                      0.0,
                                                      m_level,
                                                      m_levelGridMetrics);
          // setup the wall distance function as part of LGM
          m_levelGridMetrics.defineDistances(boundaryBoxes,
                                             faceDir);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Initialize data
/** This is only called for levels with grids
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::initialData()
{
  CH_TIME("AMRLevelCNS::initialData");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::initialData (level: " << m_level
           << ")" << CRD::end;

  LevelData<FArrayBox>& U = m_data.rawU();
  LevelData<FArrayBox>& JU = m_data.rawJU();

  // Sets <U> and <JU> from IBC class
  loadIBCData(U, JU);

  // Invalidate all ghost cells
  m_data.invalidate();

  // Set coarsened-fine data for the interpolator from the next coarser level.
  // We have to do this from the finer level because the space-interpolator is
  // not defined until now.  Of course, don't bother if there is actually no
  // grid on this level.
  if (m_hasCoarserGrid && m_boxes.size() > 0)
    {
      // Find the number of ghost cells on the coarsened-fine mesh required
      // to fill the invalid ghost cells in the fine mesh.
      const IntVect& interpolatorCrFnGhostVect =
        getCoarserLevel()->m_levelGridMetrics.interpolatorCrFnNumGhost(true);
      // Assumes ghosts in all directions are equal
      const int numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
      // All extra-block ghosts are filled but use this value for clarity
      getCoarserLevel()->m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
        getCoarserLevel()->m_data.getU(0, -1, numInterpolatorCrFnGhost),
        CRDparam::g_CRDPhysics->velocityInterval());
    }
}

/*--------------------------------------------------------------------*/
//  Things to do after initialization
/**
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::postInitialize()
{
  CH_TIME("AMRLevelCNS::postInitialize");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::postInitialize (level: " << m_level
           << ")" << CRD::end;

//--Average from finer level data (as in postTimeStep)

  if (m_hasFinerGrid)
    {
      // Volume weighted average from finer level data
      AMRLevelCNS* finerAMRLevel = getFinerLevel();
      m_averageOp.averageToCoarse(
        m_data.getJU(),
        finerAMRLevel->m_data.getJU());
      // Set wall-model data from finer level
      m_levelOp.setWallModelFromFiner(m_data.getJU());
      m_data.invalidate();
      // Set coarsened-fine data for the interpolator for the next finer level.
      // We rely on proper-nesting to keep this independent from the next
      // coarser level
      // Find the number of ghost cells on the coarsened-fine mesh required
      // to fill the invalid ghost cells in the fine mesh.
      const IntVect& interpolatorCrFnGhostVect =
        m_levelGridMetrics.interpolatorCrFnNumGhost(true);
      // Assumes ghosts in all directions are equal
      const int numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
      // All extra-block ghosts are filled but use this value for clarity
      m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
        m_data.getU(0, -1, numInterpolatorCrFnGhost),
        CRDparam::g_CRDPhysics->velocityInterval());
    }

  // This is our conservation reference
  if (m_level == 0)
    {
      s_JUConsvRef.resize(CRDparam::g_CRDPhysics->numConservative());
      computeSum(s_JUConsvRef);
    }

  // only do FFTs on coarsest level where we are certain a uniform grid exists
  // (also to keep the FFTs cheap)
  if (m_level == 0 &&
      CRDparam::g_turbForcingType == CRDparam::TurbForcingSpectral)
    {
      m_levelOp.m_spectralRefreshTime = CRDparam::g_spectralForcingDt;
      m_levelOp.m_spectralForcing.calcSpectralForce(m_boxes,
                                                    m_data.getJU(),
                                                    m_spectralSource);

      // Now interpolate coarse perturb to all finer AMR levels here.
    }
}

/*--------------------------------------------------------------------*/
//  Advance by one timestep and return new timestep
/**
 *//*-----------------------------------------------------------------*/

Real
AMRLevelCNS::advance()
{
  CH_TIME("AMRLevelCNS::advance");
  CRD::msg << CRD::fv2 << "AMRLevelCNS::advance (level: " << m_level
           << ", time: " << m_time << ", dt: " << m_dt << ")" << CRD::end;
  CRDparam::g_level = m_level;

  if (m_hasFinerGrid)
    {
      // At beginning of this level time step, store the finer level step for
      // comparison in postTimeStep
      m_begFinerLevelStep = getFinerLevel()->m_levelStep;
    }

  // Subcycle parameters (only meaningful for AMR levels)
  const SubcycleParams& subcycleParams = AMRLevel::getSubcycleParams();

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;  //**FIXME see if we can get rid of this and use
                              //**FIXME nullptr instead

  // Set arguments to dummy values and then fix if real values are available
  TimeInterpolatorRK4* finerTI = nullptr;

  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarserGrid)
    {
      AMRLevelCNS *const coarserAMRLevel = getCoarserLevel();

      // Recall that my flux register goes between my level and the next
      // finer level
      coarserFR = coarserAMRLevel->m_fluxRegisterPtr;

      tCoarserNew = coarserAMRLevel->m_time;
      tCoarserOld = tCoarserNew - coarserAMRLevel->m_dt;
    }

  // A finer level exists
  if (m_hasFinerGrid)
    {
      // We set intermediate time data from this level so the time interpolator
      // on the fine level can do the interpolations
      finerTI = getFinerLevel()->m_timeInterpolator;
      // Recall that my flux register goes between my level and the next
      // finer level
      finerFR = m_fluxRegisterPtr;
      finerFR->setToZero();
    }

//--Fill the invalid ghost of <U> (and 1 layer of <JU>) if there is a coarser
//--level.  This is akin to calling fillGhostsRK4AndComputeU() for the first
//--stage in LevelMappedConsOperator.cpp.  <U> will be required right now if
//--artificial viscosity is added
//**FIXME Perhaps we can also do the computation of <U> for the first stage and
//**avoid it in LevelMappedConsOperator.
//**FIXME If this is the first subcycle, we should be able to just use the
//**coarsened-fine data previously set by the call to postTimeStep() for the
//**next coarser level
//**FIXME If this is not the first subcycle, we should still have the correct
//**coarsened-fine data, set at the end of the previous advance() on this level

  if (m_hasCoarserGrid)
    {
      // Start of time-step as time interpolation coefficient in range [0:1] of
      // coarse time-step
      Real dtCoarse = tCoarserNew - tCoarserOld;
      Real alpha = (m_time - tCoarserOld)/dtCoarse;
      if (subcycleParams.isFirstSubcycle)
        {
          // We are at the beginning time boundary of the coarser level.  Reset
          // coarsened-fine data from the coarser level.  Note that at this
          // point, U^(l-1) has data at the end time boundary of the coarser
          // level.  However, the time interpolator has separately cached this
          // information.  It will detect alpha (is exactly) = 0. and simply
          // copy.
          alpha = 0.;
        }
      CH_assert((alpha >= 0.) && (alpha < 1.));

      // This will internally store the coarsened-fine data used to fill
      // ghosts cells for this level.  We interpolate the coarsened-fine data
      // to the *begin* time for the full time-step on *this* level.
      const Real dtRatio = m_dt/dtCoarse;
      m_levelGridMetrics.timeIntermediate(
        *m_timeInterpolator,
        alpha,
        dtRatio,
        0,  // Always stage zero since we are always outside of stages
        CRDparam::g_CRDPhysics->velocityInterval(),
        m_data.getJU(),
        true);
    }

  // Save <JU> in old location.  This means <JU> new is already initialized to
  // the current and this operation can be avoided in RK4LevelAdvance
  m_data.copyJUNewToOld();

//--Add in artificial viscosity

  if (CRDparam::g_useArtificialViscosity)
    {
      LevelData<FArrayBox>& UOld  = m_data.getU(m_numGhost, m_numGhost);
      LevelData<FArrayBox>& JUNew = m_data.getJU();
      // The updates the flux registers and JUNew
      m_levelOp.addArtificialViscosity(JUNew,
                                       *finerFR,
                                       *coarserFR,
                                       UOld,
                                       m_WOld,
                                       m_dt,
                                       m_time);
    }

//--Advance conservation law by one time step using 4th-order Runge-Kutta.

  if (CRDparam::g_additiveRK
      && (m_level <= CRDparam::g_ARKmaxLevel)
      && (m_levelStep >= CRDparam::g_initERKSteps)
      && (CRDparam::g_timeIntegrationMethod == CRDparam::ARK4))
    {
      CRD::msg << CRD::fv2 << "AMRLevelCNS::advance: level " << m_level
               << ": Doing ARK time step" << CRD::end;
      // Either ARK4::Extrapolate2ndOrder or ARK4::PreviousStageValues
      ARK4::StageValuePredictorMethod initGuessMethod =
        ARK4::PreviousStageValues;
      if (CRDparam::g_ARKExtrapInitGuess && !m_justRegridded)
        {
          initGuessMethod = ARK4::Extrapolate2ndOrder;
        }
      TimeInterpolatorARK4* ARKfinerTI =
        dynamic_cast<TimeInterpolatorARK4*>(finerTI);
      m_ark_time_integrator.advance
                (m_data.getJU(),
                 m_data.getJUOld(),
                 m_time,
                 m_dt,
                 m_prevStepDt,
                 ARKfinerTI,
                 m_levelOp,
                 initGuessMethod,
                 CRDparam::g_ARKUsePIDControl,
                 *finerFR,
                 *coarserFR,
                 m_dt,
                 m_time,
                 tCoarserOld,
                 tCoarserNew,
                 subcycleParams,
                 m_WOld);
    }
  else if (CRDparam::g_timeIntegrationMethod == CRDparam::RK4)
    {
      CRD::msg << CRD::fv2 << "AMRLevelCNS::advance: level " << m_level
               << ": Doing explicit RK4 time step" << CRD::end;
      // The time interpolator from the finer level is filled with information
      // from the stages taken on this level.  We'll use it later when we
      // calculate on the finer level.
      // Construct the flags for which terms to calculate
      int termFlags =
        LevelCNSOp::Terms::NonStiff | LevelCNSOp::Terms::Stiff;
      RK4LevelAdvance<LevelData<FArrayBox>,
                      LevelData<FArrayBox>,
                      TimeInterpolatorRK4,
                      LevelCNSOp>
        (m_data.getJU(),
         m_data.getJUOld(),
         finerTI,
         m_time,
         m_dt,
         false,
         m_levelOp,
         // Additional arguments forwarded to m_levelOp
         *finerFR,
         *coarserFR,
         m_dt,
         m_time,
         tCoarserOld,
         tCoarserNew,
         subcycleParams,
         termFlags,
         m_WOld);
    }
  else if(CRDparam::g_timeIntegrationMethod == CRDparam::RK2)
    {
      CRD::msg << CRD::fv2 << "AMRLevelCNS::advance: level " << m_level
               << ": Doing explicit RK2 time step" << CRD::end;
      // The time interpolator from the finer level is filled with information
      // from the stages taken on this level.  We'll use it later when we
      // calculate on the finer level.
      // Construct the flags for which terms to calculate
      int termFlags =
        LevelCNSOp::Terms::NonStiff | LevelCNSOp::Terms::Stiff;
      TimeInterpolatorRK2* RK2TimeInterp =
        dynamic_cast<TimeInterpolatorRK2*>(finerTI);
      RK2LevelAdvance<LevelData<FArrayBox>,
                      LevelData<FArrayBox>,
                      TimeInterpolatorRK2,
                      LevelCNSOp>
        (m_data.getJU(),
         m_data.getJUOld(),
         RK2TimeInterp,
         m_time,
         m_dt,
         false,
         m_levelOp,
         // Additional arguments forwarded to m_levelOp
         *finerFR,
         *coarserFR,
         m_dt,
         m_time,
         tCoarserOld,
         tCoarserNew,
         subcycleParams,
         termFlags,
         m_WOld);

    }

  {
    LevelData<FArrayBox>& JUNew = m_data.getJU();
    m_levelOp.speciesCorrection(JUNew);
  }

  // Update the time and store the new timestep
  m_time += m_dt;

//--So that we can compute <U> properly in the invalid ghost cells on this
//--level, we have to set the CrFnU to the current time in the interpolator
//**FIXME If this is the last subcycle, then postTimeStep will also update
//**CrFnU.  Can we avoid doing this twice?  ComputeNewDt gets in the way.

  // if (m_hasCoarserGrid)
  //   {
  //     //**FIXME Is this necessary or can we just use the results from the last
  //     //**      stage?
  //     // End of time-step as time interpolation coefficient in range [0:1] of
  //     // coarse time-step
  //     if (subcycleParams.isLastSubcycle)
  //       {
  //         // We have reached the end time boundary of the coarser level.  Set
  //         // coarsened-fine data from U (which is still new U) from the coarser
  //         // level
  //         AMRLevelCNS *const coarserAMRLevel = getCoarserLevel();
  //         coarserAMRLevel->m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
  //           coarserAMRLevel->m_data.getU(0, -1),
  //           CRDparam::g_CRDPhysics->velocityInterval());
  //       }
  //     else
  //       {
  //         //**FIXME -- not better to use stage 3 of previous time step?
  //         // We need to interpolate the state at the end of the time step
  //         Real dtCoarse = tCoarserNew - tCoarserOld;
  //         // End of time-step as time interpolation coefficient in range [0:1]
  //         // of coarse time-step
  //         Real alpha = (m_time - tCoarserOld)/dtCoarse;
  //         CH_assert((alpha > 0.) && (alpha < 1.));
  //         const Real dtRatio = m_dt/dtCoarse;
  //         m_levelGridMetrics.timeIntermediate(
  //           m_timeInterpolator,
  //           alpha,
  //           dtRatio,
  //           0,  // Always stage zero since we are always outside of stages.
  //               // Note: we are at the first stage of the next subcycle
  //           CRDparam::g_CRDPhysics->velocityInterval());
  //       }
  //   }

//--Store current level SGS kinetic energy estimate in m_U and m_JU
  if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
    {
      LevelData<FArrayBox>& JUNew = m_data.getJU();
      LevelData<FArrayBox>& UNew = m_data.getU();
      const LevelData<FArrayBox>& SGSKE = m_levelOp.getSGSKineticEnergyVar();
      const LevelData<FArrayBox>& JSGSKE = m_levelOp.getJSGSKineticEnergyVar();
      const int sgsKEcomp = CRDparam::g_CRDPhysics->turbConsInterval().begin();
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& keFab = SGSKE[dit];
          const FArrayBox& JkeFab = JSGSKE[dit];
          FArrayBox& JUFab = JUNew[dit];
          FArrayBox& UFab = UNew[dit];
          MD_BOXLOOP(m_boxes[dit], i)
            {
              JUFab[MD_IX(i, sgsKEcomp)] = JkeFab[MD_IX(i, 0)];
              UFab[MD_IX(i, sgsKEcomp)] = keFab[MD_IX(i, 0)];
            }
        }
    }

//--Special processing of new state

  if (CRDparam::g_explicitFilterType == CRDparam::ExplicitFilterSpectral)
    {
      m_levelOp.m_spectralFilter.apply(
        CRDparam::g_spectralFilterDomainResolution,
        m_data.getJU(),
        Interval(0, CRDparam::g_CRDPhysics->numConservative() - 1));
    }

  // All we really know is <JU> on valid cells.
  // The next finer level may perform a regrid as the next operation:
  //   1) If it is on the first sub-cycle, it needs coarse data from this level
  //      from the old time.  This was already stored in a previous
  //      postTimeStep() on this level.  The 'false' argument below indicates
  //      that this data is still valid.
  //   2) If it is on a later sub-cycle, it will have updated the coarsened-
  //      fine information it requires in the block above.
  m_data.invalidate(false);

  ++m_levelStep;

  // Return an estimate of the next time step.
  //**FIXME We probably need to think carefully about where this is done -- here
  //**it forces lots of message passing to get <U>.  It would be nice to do this
  //**in post-time step.
  m_dtNew = m_levelOp.computeNewDt(m_data.getU(1, 1),
                                   m_time,
                                   m_dt,
                                   m_prevStepDt);
  m_prevStepDt = m_dt;
  return m_dtNew;  // This is actually discarded in AMR
}

/*--------------------------------------------------------------------*/
//  Things to do after a timestep -- reflux
/**
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::postTimeStep()
{
  CH_TIME("AMRLevelCNS::postTimeStep");
  CRD::msg << CRD::fv2 << "AMRLevelCNS::postTimeStep (level: " << m_level
           << ")" << CRD::end;
  if (m_hasFinerGrid)
    {
      AMRLevelCNS* finerAMRLevel = getFinerLevel();

      // Test if finer level took extra time steps than refRatio
      if (finerAMRLevel->m_levelStep - m_begFinerLevelStep != refRatio())
        {
          CRD::msg << "Level " << m_level+1 << " required "
                   << finerAMRLevel->m_levelStep - m_begFinerLevelStep
                   << " time-steps instead of " << refRatio() << CRD::warn;
        }

      // Reflux
      Real scale = -1.0/m_dx[0];  //**FIXME dxVect
      m_fluxRegisterPtr->reflux(m_data.getJU(), scale);

      // Set m_JUnew at cell 'i' to the mean of amrConsFinerPtr->m_JUnew
      // over all finer cells within 'i'.
      m_averageOp.averageToCoarse(
        m_data.getJU(),
        finerAMRLevel->m_data.getJU());

      // Set wall-model data from finer level
      m_levelOp.setWallModelFromFiner(m_data.getJU());

      // Average metric terms (these only average if required)
      m_levelGridMetrics.postTimeStep();
      m_data.invalidate();  // Since <JU> was just modified

      // Set coarsened-fine data for the interpolator for the next finer level.
      // We rely on proper-nesting to keep this independent from the next
      // coarser level.  This is required for tagging (maybe move there?)
      // Find the number of ghost cells on the coarsened-fine mesh required
      // to fill the invalid ghost cells in the fine mesh.
      const IntVect& interpolatorCrFnGhostVect =
        m_levelGridMetrics.interpolatorCrFnNumGhost(true);
      // Assumes ghosts in all directions are equal
      const int numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
      // All extra-block ghosts are filled but use this value for clarity
      m_levelGridMetrics.presetThis2FnInterpolatorCrFnLevU(
        m_data.getU(0, -1, numInterpolatorCrFnGhost),
        CRDparam::g_CRDPhysics->velocityInterval());
    }

  // if (m_level == 0)
  //   {
  //     reportError();
  //   }

  // Only do FFTs on coarsest level where we are certain a uniform grid exists
  // (also to keep the FFTs cheap)
  if (m_level == 0 &&
      CRDparam::g_turbForcingType == CRDparam::TurbForcingSpectral)
    {
      m_levelOp.addSpectralForce(m_data.getJU(), m_spectralSource, m_dt);

      if (m_levelOp.m_spectralRefreshTime < (m_time + 0.5*m_dt))
        {
          m_levelOp.m_spectralForcing.calcSpectralForce(m_boxes,
                                                        m_data.getU(),
                                                        m_spectralSource);

          // Now interpolate coarse perturb to all finer AMR levels here.
          // Also don't forget to add interpolation in regrid() or something.

          while (m_levelOp.m_spectralRefreshTime < (m_time + 0.5*m_dt))
            {
              m_levelOp.m_spectralRefreshTime += CRDparam::g_spectralForcingDt;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Create tags for regridding
/**
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::tagCells(IntVectSet& a_tags)
{
  CH_TIME("AMRLevelCNS::tagCells");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::tagCells (level: " << m_level << ")"
           << CRD::end;

  IntVectSet localTags;
  // If we are restarting and tagging occurs right away, the ghost cells
  // on the finer levels have not been filled
  // FIXME: This should be moved to another location so that it is not
  // constantly being checked
  if (m_exchangeForTagging)
    {
      m_levelOp.restartExchange(m_data.getU(1,1));
      m_exchangeForTagging = false;
    }
  m_tagLevel->tagCells(localTags,
                       m_data,
                       m_levelGridMetrics,
                       *m_multiBlockCoordSys,
                       m_time,
                       m_level);

  a_tags = localTags;
}

/*--------------------------------------------------------------------*/
//  Create tags at initialization
/**
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::tagCellsInit(IntVectSet& a_tags)
{
  CH_TIME("AMRLevelCNS::tagCellsInit");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::tagCellsInit (level: " << m_level
           << ")" << CRD::end;

  tagCells(a_tags);
}


/*--------------------------------------------------------------------*/
//  Pre-regrid operations necessary for updating metric terms
/** \param[in]  a_baseLevel
 *                      Index of the base of the regrid (boxes do not
 *                      change on this level).
 *  \param[in]  a_newBoxes
 *                      Vectors of boxes defining the new grid levels
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::preRegrid(int                         a_baseLevel,
                       const Vector<Vector<Box> >& a_newBoxes)
{
  CH_TIME("AMRLevelCNS::preRegrid");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::preRegrid (level: " << m_level << ")"
           << CRD::end;

  if (m_level == finestLevel())
    {
      // Set the DBL on all levels that have boxes
      AMRLevelCNS* walkAMRLevelCNS = this;
      Vector<int> procIDs;
      for (int level = m_level; level > a_baseLevel; --level)
        {
          const int numBoxes = a_newBoxes[level].size();
          if (numBoxes > 0)
            {
              LoadBalance(procIDs, a_newBoxes[level]);
              if (CRDparam::g_verbosity >= 4)
                {
                  CRD::msg << CRD::fv4 << "AMRLevelCNS::preRegrid: processor "
                    "map for level " << level << ':' << CRD::end;
                  for (int iBox = 0; iBox != numBoxes; ++iBox)
                    {
                      CRD::msg << CRD::fv4 << "  Box " << iBox << ", proc: "
                               << procIDs[iBox] << ", cells: "
                               << a_newBoxes[level][iBox].volume() << CRD::end;
                    }
                }
            }
          DisjointBoxLayout dbl(
            a_newBoxes[level],
            procIDs,
            walkAMRLevelCNS->m_multiBlockCoordSys->levelDomain());
          dbl.close();

          walkAMRLevelCNS->m_boxes = dbl;
          walkAMRLevelCNS->m_boxes.defineLocalBlocks(
            walkAMRLevelCNS->m_multiBlockCoordSys->mappingBlocks());
          walkAMRLevelCNS = walkAMRLevelCNS->getCoarserLevel();
        }
    }

  // <U> needs to be up-to-date here (on the old mesh), in the valid cells of
  // the coarse level and in the valid cells + 1 layer of valid ghost cells for
  // this level.  However, if this is the finest level on the new hierarchy,
  // then we don't need any <U>s.  Since this routine traverses fine-to-coarse,
  // we retrieve <U> on the coarser level, if it exists.

  // We do not require <U> and coarser <U> if this is the finest level on
  // the new hierarchy (note that new levels can only be added one at a
  // time).  UPtr will be dereferenced but may not be used.  Start by pointing
  // it to rawU.
  LevelData<FArrayBox>* UPtr             = &(m_data.rawU());
  const LevelData<FArrayBox>* coarseUPtr = nullptr;
  DisjointBoxLayout* coarseBoxes = nullptr;
  if (m_level < AMRLevel::finestLevel())
    {
      if (m_hasCoarserGrid)
        {
//**FIXME As written, this is almost certainly not correct
//**      We probably want to exchange const coarse U using MMB flux correction
//**      code
          //**FIXME We need (1,1) for multiblock until the BlockRegister is used
          //**      to do snapback
          // coarseUPtr = &(getCoarserLevel()->m_data.getU());
          coarseUPtr = &(getCoarserLevel()->m_data.getU(1, 1));
        }
      UPtr = &(m_data.getU(1, 0));
    }
  if (m_hasCoarserGrid)
    {
      coarseBoxes = &(getCoarserLevel()->m_boxes);
    }
  // Use pointers/references here to clearly avoid side-effects
  LevelData<FArrayBox>& JUNew = m_data.getJU();
  m_levelGridMetrics.preRegrid(a_baseLevel,
                               m_boxes,
                               coarseBoxes,
                               coarseUPtr,
                               *UPtr,
                               JUNew);
  const int numPrimitive = CRDparam::g_CRDPhysics->numPrimitive();
  const int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostWfromUcellAvg);
  m_WOld.define(m_boxes, numPrimitive, IntVect::Unit*numGhosts);

  // And now everything is invalid because <JU> may have changed
  m_data.invalidate();
}

/*--------------------------------------------------------------------*/
//  Set up data on this level after regridding
/** \param[in]  a_newBoxes
 *                      A vector of boxes defining the new grid level
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::regrid(const Vector<Box>& a_newBoxes)
{
  CH_TIME("AMRLevelCNS::regrid");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::regrid (level: " << m_level << ")"
           << CRD::end;

//--If we don't receive any boxes, then this level is going away.  This has been
//--separated out to clearly illustrate what should happen

  if (a_newBoxes.size() == 0)
    {
      CRD::msg << CRD::fv3 << "AMRLevelCNS::regrid: level " << m_level
               << " is being created with 0 cells. " << CRD::end;

      // Sanity check: LGM should have also detected that this level is going
      // away and the metric terms should be undefined
      CH_assert(!m_levelGridMetrics.metricsDefined());
      CH_assert(m_levelGridMetrics.getBoxes().size() == 0);

      // Save the original grids in the base AMRLevel class (required)
      m_level_grids = a_newBoxes;

      // Load balance and construction of DBL still must be done (if this level
      // persists, normally it would have been done in preRegrid).
      {
        Vector<int> procIDs;
        LoadBalance(procIDs, a_newBoxes);
        DisjointBoxLayout dbl(a_newBoxes,
                              procIDs,
                              m_multiBlockCoordSys->levelDomain());
        dbl.close();
        m_boxes = dbl;
      }

      // Everything is built but with a zero-sized box layout
      m_data.reshapeJUNew();
      m_data.reshapeJUOld();
      m_data.reshapeU();
      m_data.invalidate();

      // levelOpSetup() will correctly set m_hasFiner and m_hasCoarser for this
      // level and fix up m_hasFiner for the coarser level.
      levelOpSetup();
      reverseLevelOpSetup();
      return;
    }

  int numBoxes = a_newBoxes.size();

  // Given a properly construct base, all boxes should be in the domain

  // First, intersect the given boxes with each of the blocks in our
  // coordinate system.
  // const Vector<Box>& blocks = m_multiBlockCoordSys->mappingBlocks();
  // Vector<Box> validBoxes;
  for (int iBox = 0; iBox < numBoxes; ++iBox)
    {
      const Box& box = a_newBoxes[iBox];
      CH_assert(m_multiBlockCoordSys->whichBlock(box) >= 0);
      // for (int iBlock = 0; iBlock < blocks.size(); ++iBlock)
      //   {
      //     const Box& blockBox = blocks[iBlock];
      //     Box intersect = box & blockBox;
      //     if (!intersect.isEmpty())
      //       {
      //         validBoxes.push_back(intersect);
      //       }
      //   }
    }
  // numBoxes = validBoxes.size();

  // Save the original grids in the base AMRLevel class (required)
  //**FIXME can we do this in preRegrid?
  m_justRegridded = true;
  m_level_grids = a_newBoxes;

  // Save data for later
  m_data.copyJUNewToOld();

  // Reshape is simple because MappedLevelData already knows about the grid
  // layout from LevelGridMetrics.  Reshape state with new grids
  m_data.reshapeJUNew();

  // Note that U is now on a different layout from JU.  Any attempts to fill
  // <U> (beyond a simple exchange) will therefore fail.

//--Fill <JU>

  // Interpolate from coarser level
  if (AMRLevel::hasCoarserLevel())  // Use from base class since levelOpSetup()
    {                               // has not yet been called.  Of coarse, why
                                    // would the base level ever be regridded?
      AMRLevelCNS *const coarserAMRLevel = getCoarserLevel();
      // Use references here to clearly avoid side-effects.  With proper
      // nesting, we only need <U> computed from valid <JU>.  I.e., 1 layer
      // away from coarse-fine interfaces is sufficient.
      // We also need extra-block ghosts to be filled.
      //**FIXME suggest a getU with (valid ghosts, invalid ghosts,
      //**                           extra-block ghosts)
      LevelData<FArrayBox>& crU     = coarserAMRLevel->m_data.getU(0, -1);
      coarserAMRLevel->m_levelGridMetrics.multiblockExchangeU(
        crU, CRDparam::g_CRDPhysics->velocityInterval());
      LevelData<FArrayBox>& crJUNew = coarserAMRLevel->m_data.getJU();
      // We don't required any ghosts in JU here
      m_levelGridMetrics.regrid(m_data.getJU(),
                                crU,
                                crJUNew,
                                CRDparam::g_CRDPhysics->velocityInterval());
    }

  // Copy from old state
  m_data.copyJUOldToNew();

  // Reshape JUold to be consistent with JUnew
  m_data.reshapeJUOld();

  // And <U> needs to be reshaped too
  m_data.reshapeU();

  // And now everything is invalid because <JU> may have changed
  m_data.invalidate();

  // Setup the remaining operators on the level
  levelOpSetup();
}

/*--------------------------------------------------------------------*/
//  Fine-to-coarse traversal after a regrid
/** \param[in]  a_newBoxes
 *                      A vector of boxes defining the new grid level
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::postRegrid(int a_baseLevel)
{
  CH_TIME("AMRLevelCNS::postRegrid");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::postRegrid (level: " << m_level << ")"
           << CRD::end;

  reverseLevelOpSetup();
}

#ifdef CH_USE_HDF5

/*--------------------------------------------------------------------*/
//  Write checkpoint header
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCNS::writeCheckpointHeader");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::writeCheckpointHeader (level: "
           << m_level << ")" << CRD::end;

  //--We write out all components of m_JU or m_U

  HDF5HeaderData header;
  if (CRDparam::g_plotDACFDCheck)
    {
      // This part only for DA-CFD system

      // Set up the number of components
      const int numStates = CRDparam::g_CRDPhysics->numOutputVar();
      int numStatesOutput = numStates + 1;  // U + T + J

      header.m_int["num_components"] = numStatesOutput;

      // Set up the component names
      char compStr[32];

      int comploc = 0;
      // U
      for (int comp = 0; comp != numStates; ++comp)
        {
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] =
            CRDparam::g_CRDPhysics->consvStateName(comp);
          ++comploc;
        }

      // J
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "J";
      ++comploc;

      // Sanity check
      CH_assert(comploc == numStatesOutput);

    }
  else
    {
      // DEFAULT

      // Set up the number of components
      const int numStates = CRDparam::g_CRDPhysics->numConservative();
      int numStatesOutput = numStates;  // JU

      header.m_int["num_components"] = numStatesOutput;

      // Set up the component names
      char compStr[32];
      char nameStr[64];

      for (int comp = 0; comp < numStates; ++comp)
        {
          sprintf(compStr, "component_%d", comp);
          sprintf(nameStr, "J-%s",
                  CRDparam::g_CRDPhysics->consvStateName(comp));
          header.m_string[compStr] = nameStr;
        }
    }

  // Write the header
  header.writeToFile(a_handle);
  CRD::msg << CRD::fv3 << header << CRD::end;
}

/*--------------------------------------------------------------------*/
//  Write checkpoint data for this level
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCNS::writeCheckpointLevel");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::writeCheckpointLevel (level: "
           << m_level << ")" << CRD::end;

  // Setup the level string
  char levelStr[32];
  sprintf(levelStr, "%d", m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  const ProblemDomain& domain = m_multiBlockCoordSys->levelDomain();
  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = 0;  // Not stored in Chord
  header.m_real["dx"]              = m_dx[0];  //**FIXME
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = domain.domainBox();

  // Setup the periodicity info
  D_TERM6(
    header.m_int ["is_periodic_0"] = (int)domain.isPeriodic(0);,
    header.m_int ["is_periodic_1"] = (int)domain.isPeriodic(1);,
    header.m_int ["is_periodic_2"] = (int)domain.isPeriodic(2);,
    header.m_int ["is_periodic_3"] = (int)domain.isPeriodic(3);,
    header.m_int ["is_periodic_4"] = (int)domain.isPeriodic(4);,
    header.m_int ["is_periodic_5"] = (int)domain.isPeriodic(5););

  // Write the header for this level
  header.writeToFile(a_handle);
  CRD::msg << CRD::fv3 << header << CRD::end;


  // Write the data for this level, current <U> and then <J>
  if (CRDparam::g_plotDACFDCheck)
    {
      // The number of ghosts we output
      const IntVect outputGhost = CRDparam::g_plotNumGhost*IntVect::Unit;
      // Number of outputs for the plot files
      const int numOutStates = CRDparam::g_CRDPhysics->numOutputVar();
      // Write the boxes and construct the dataset for multiple levels
      LevelData<FArrayBox> U(m_data.getU(1, 1).getBoxes(),
                             numOutStates,
                             outputGhost);

      pout() << "For DACFD checkpoint, numOutStates = " << numOutStates << std::endl;

      CRDparam::g_CRDPhysics->outputLevelData(U,
                                              m_data.getU(1, 1),
                                              m_WOld,
                                              m_levelGridMetrics);
      int numCompWrite = U.nComp();
      CH_assert(numCompWrite == numOutStates);

      // // Component for temperature
      // numCompWrite += 1;

      // Components for J
      CH_assert(m_levelGridMetrics.m_J.ghostVect() >= outputGhost);
      CH_assert(m_levelGridMetrics.m_J.nComp() == 1);
      numCompWrite += 1;

      WriteMultiData<FArrayBox> writer(a_handle,
                                       U.getBoxes(),
                                       numCompWrite,
                                       "data",
                                       CH_HDF5::IOPolicyDefault,
                                       outputGhost);


      m_levelOp.speciesCorrection(U);

      int iCompWrite = 0;
      writer.writeData(U,
                       U.interval(),
                       U.interval());
      iCompWrite += U.nComp();

      //FIXME (temperature)

      writer.writeData(m_levelGridMetrics.m_J,
                       m_levelGridMetrics.m_J.interval(),
                       Interval(iCompWrite, iCompWrite));
      iCompWrite += 1;

    }
  else
    {
      // Write the data for this level, current <JU>
      write(a_handle, m_data.rawJU().boxLayout());
      write(a_handle, m_data.rawJU(), "data");
    }

}

/*--------------------------------------------------------------------*/
//  Read checkpoint header
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::readCheckpointHeader(HDF5Handle& a_handle)
{
  CH_TIME("AMRLevelCNS::readCheckpointHeader");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::readCheckpointHeader (level: "
           << m_level << ")" << CRD::end;

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  CRD::msg << CRD::fv4 << "HDF5 header data: " << header << CRD::end;

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointHeader: checkpoint file does "
        "not have num_components" << CRD::error;
    }

  const int numStatesPhys = CRDparam::g_CRDPhysics->numConservative();
  int numStatesInput = numStatesPhys;
  int numStatesRead = header.m_int["num_components"];
  if (CRDparam::g_plotDACFDCheck)
    {
      // numStatesInput = numStatesInput + 2; // T + J
      numStatesInput = CRDparam::g_CRDPhysics->numOutputVar() + 1;
    }
  if (CRDparam::g_restartAddWallModel)
    {
      numStatesInput = CRDparam::g_CRDPhysics->numConservative() - 1 - SpaceDim;
    }

  if (numStatesInput != numStatesRead)
    {
      CRD::msg << "AMRLevelCNS::readCheckpointHeader: num_components in "
        "checkpoint file does not match solver physics" << CRD::error;
    }

  // Get the component names
  std::string stateNameRead;
  char compStr[32];
  char stateNamePhys[64];
  if (CRDparam::g_plotDACFDCheck)
    {
      // U + T
      int numStatesDA = CRDparam::g_CRDPhysics->numOutputVar();
      for (int comp = 0; comp < numStatesDA; ++comp)
        {
          sprintf(compStr, "component_%d", comp);
          if (header.m_string.find(compStr) == header.m_string.end())
            {
              CRD::msg << "AMRLevelCNS::readCheckpointHeader: checkpoint file "
                "does not have enough component names" << CRD::error;
            }

          stateNameRead = header.m_string[compStr];
          sprintf(stateNamePhys,
                  CRDparam::g_CRDPhysics->consvStateName(comp));
          if (std::strcmp(stateNameRead.c_str(), stateNamePhys) != 0)
            {
              CRD::msg << "AMRLevelCNS::readCheckpointHeader: state_name in "
                "checkpoint (" << stateNameRead << ") does not match name in "
                "solver (" << stateNamePhys << CRD::error;
            }
        }

      // J
      for (int comp = numStatesDA; comp < numStatesInput; ++comp)
        {
          sprintf(compStr, "component_%d", comp);
          header.m_string[compStr] = "J";

        }
    }
  else if (CRDparam::g_restartAddWallModel)
    {
      // Don't worry about checking states for now
    }
  else
    {
      // DEFAULT (JU)
      for (int comp = 0; comp < numStatesPhys; ++comp)
        {
          sprintf(compStr, "component_%d", comp);
          if (header.m_string.find(compStr) == header.m_string.end())
            {
              CRD::msg << "AMRLevelCNS::readCheckpointHeader: checkpoint file "
                "does not have enough component names" << CRD::error;
            }

          stateNameRead = header.m_string[compStr];
          sprintf(stateNamePhys,
                  "J-%s",
                  CRDparam::g_CRDPhysics->consvStateName(comp));
          if (std::strcmp(stateNameRead.c_str(), stateNamePhys) != 0)
            {
              CRD::msg << "AMRLevelCNS::readCheckpointHeader: state_name in "
                "checkpoint (" << stateNameRead << ") does not match name in "
                "solver (" << stateNamePhys << CRD::error;
            }
        }
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointHeader: checkpoint file does "
        "not contain iteration" << CRD::error;
    }
  m_levelStep = header.m_int ["iteration"];
}

/*--------------------------------------------------------------------*/
//  Read checkpoint data for this level
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_TIME("AMRLevelCNS::readCheckpointLevel");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::readCheckpointLevel (level: "
           << m_level << ")" << CRD::end;

  // Read the header for this level
  char levelStr[32];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  CRD::msg << CRD::fv4 << "HDF5 header data: " << header << CRD::end;


//--Read preliminary data

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain "
        "ref_ratio at level " << m_level << '!' << CRD::error;
    }
  // the checkpoint refinement ratio stored in the checkpoint is meaningless
  // unless a finer level was actually defined. The ref_ratio is instead taken
  // from the input file
  if (hasFinerLevel())
    {
      m_ref_ratio = header.m_int["ref_ratio"];
    }
  CH_assert(m_ref_ratio > 0);
  CRD::msg << CRD::fv2 << "Refinement ratio\n: " << m_ref_ratio << CRD::varend;
  if (m_level > 0)
    {
      // The ref ratio we are interested in is from the coarser level.  As
      // this routine reads levels from coarse->fine, we know the coarser level
      // is defined
      const int readRefRatio = getCoarserLevel()->m_ref_ratio;
      const int inputRefRatio = (CRDparam::g_refFromBase[m_level]/
                                 CRDparam::g_refFromBase[m_level-1]);
      if (readRefRatio != inputRefRatio)
        {
          const int refRatioChange = readRefRatio/inputRefRatio;
          CRD::msg << "Modifying reference ratios from input file to match "
            "checkpoint file at level " << m_level  << '!' << CRD::warn;
          for (int i = m_level; i != CRDparam::numAMRLevel(); ++i)
            {
              CRDparam::CRDP.adjustRefFromBase(
                i,
                CRDparam::g_refFromBase[m_level]*refRatioChange);
            }
        }
    }

  // In Chord, we neglect the tag buffer size

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain dx "
        "at level " << m_level << '!' << CRD::error;
    }
  m_dx = header.m_real["dx"]*RealVect::Unit;  //**FIXME dxVect
  CRD::msg << CRD::fv2 << "Cell spacing\n: " << m_dx << CRD::varend;

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain dt "
        "at level " << m_level << '!' << CRD::error;
    }
  m_dt = header.m_real["dt"];
  CRD::msg << CRD::fv2 << "Time step\n: " << m_dt << CRD::varend;

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain "
        "time at level " << m_level << '!' << CRD::error;
    }
  m_time = header.m_real["time"];
  CRD::msg << CRD::fv2 << "Time\n: " << m_time << CRD::varend;

  // // Get the problem domain
  // if (header.m_box.find("prob_domain") == header.m_box.end())
  //   {
  //     CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain "
  //       "prob_domain at level " << m_level << '!' << CRD::error;
  //   }
  // Box domainBox = header.m_box["prob_domain"];  //**FIXME multiblock

  // // Get the periodicity info -- this is more complicated than it really
  // // needs to be in order to preserve backward compatibility.
  // bool isPeriodic[SpaceDim];
  // D_TERM6(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
  //           isPeriodic[0] = (header.m_int["is_periodic_0"] == 1);
  //         else
  //           isPeriodic[0] = false; ,

  //         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
  //           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
  //         else
  //           isPeriodic[1] = false; ,

  //         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
  //           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
  //         else
  //           isPeriodic[2] = false; ,

  //         if (!(header.m_int.find("is_periodic_3") == header.m_int.end()))
  //           isPeriodic[3] =  (header.m_int["is_periodic_3"] == 1);
  //         else
  //           isPeriodic[3] = false; ,

  //         if (!(header.m_int.find("is_periodic_4") == header.m_int.end()))
  //           isPeriodic[4] =  (header.m_int["is_periodic_4"] == 1);
  //         else
  //           isPeriodic[4] = false; ,

  //         if (!(header.m_int.find("is_periodic_5") == header.m_int.end()))
  //           isPeriodic[5] =  (header.m_int["is_periodic_5"] == 1);
  //         else
  //           isPeriodic[5] = false;);

  // Set empty problem domain.  m_problem_domain must not be used!
  m_problem_domain = ProblemDomain{};

  // Get the grids
  Vector<Box> readBoxes;
  const int gridStatus = read(a_handle, readBoxes);
  if (gridStatus != 0)
    {
      CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain "
        "Vector<Box> at level " << m_level << '!' << CRD::error;
    }
  const int numBoxes = readBoxes.size();

  // Save the original grids in the base AMRLevel class (required)
  m_level_grids = readBoxes;

  // Load balance
  Vector<int> procIDs;
  LoadBalance(procIDs, readBoxes);
  if (CRDparam::g_verbosity >= 4)
    {
      CRD::msg << CRD::fv4 << "AMRLevelCNS::readCheckpointLevel: processor map:"
               << CRD::end;
      for (int iBox = 0; iBox != numBoxes; ++iBox)
        {
          CRD::msg << CRD::fv4 << "  Box " << iBox << ", proc: "
                   << procIDs[iBox] << ", cells: " << readBoxes[iBox].volume()
                   << CRD::end;
        }
    }

  {
    DisjointBoxLayout dbl(readBoxes,
                          procIDs,
                          m_multiBlockCoordSys->levelDomain());
    dbl.close();
    m_boxes = dbl;
    m_boxes.defineLocalBlocks(m_multiBlockCoordSys->mappingBlocks());
  }

//--Set up for LevelGridMetrics

  // Pass to grid metrics
  m_levelGridMetrics.initialGrid(&m_boxes);
  const int numPrimitive = CRDparam::g_CRDPhysics->numPrimitive();
  const int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostWfromUcellAvg);
  m_WOld.define(m_boxes, numPrimitive, IntVect::Unit*numGhosts);

//--Define MappedLevelData

  // Once LevelGridMetrics knows the layout, we can define the data
  const LevelData<FArrayBox>* coarseUPtr = 0;
  if (hasCoarserPtr())
    {
      coarseUPtr = &(getCoarserLevel()->m_data.rawU());
    }
  m_data.define(&m_levelGridMetrics,
                &m_U,
                0,  // 0 means m_data will allocate the LevelData
                0,  // 0 means m_data will allocate the LevelData
                coarseUPtr,
                CRDparam::g_CRDPhysics->numConservative(),   // Num states
                CRDparam::g_CRDPhysics->velocityInterval(),  // Velocity loc
                m_numGhost);

//--Allocate the flux register

  //**FIXME Can this have a define so as to not be a pointer?
  if (hasFinerPtr())
    {
      if (m_fluxRegisterPtr != nullptr)
        {
          delete m_fluxRegisterPtr;
        }
      m_fluxRegisterPtr = new MultiBlockFluxRegister(
        m_multiBlockCoordSys.operator->(),
        getFinerLevel()->m_multiBlockCoordSys.operator->());
    }

//--Setup the remaining operators on the level

  levelOpSetup();

//--Read the data
  if (CRDparam::g_plotDACFDCheck)
    {
      // U and J (MODIFIED by Yijun 11/22/2019)
      LevelData<FArrayBox>& a_U = m_data.rawU();
      LevelData<FArrayBox>& a_JU = m_data.rawJU();

      // The number of ghosts we output
      const IntVect inputGhost = CRDparam::g_plotNumGhost*IntVect::Unit;
      // Number of J-component values
      const int numStates = CRDparam::g_CRDPhysics->numConservative();
      // New LevelData for reading all the data (U+J) at this level
      LevelData<FArrayBox> a_checkUJ (a_U.getBoxes(), numStates+3, inputGhost+1);
      // NEW LevelData for J (FIXME)
      LevelData<FArrayBox> a_J (a_U.getBoxes(), 1, inputGhost+1);

      // Interval for copyTo
      const int compShift =  CRDparam::g_CRDPhysics->numOutputVar() - numStates;
      const Interval JIntvSrc(numStates+compShift, numStates+compShift);

      const int dataStatus = read<FArrayBox>(a_handle,
                                             a_checkUJ,
                                             "data",
                                             m_boxes,
                                             a_checkUJ.interval(),
                                             false);
      // Copy information to U
      a_checkUJ.copyTo(a_U.interval(),
                       a_U,
                       a_U.interval());
      // Copy information to J
      a_checkUJ.copyTo(JIntvSrc,
                       a_J,
                       a_J.interval());

      // Call the exchange function
      a_U.exchange();
      a_J.exchange();

      //--Compute <JU> from <U> and <J>. One-sided derivates are used for <U> at
      //--the physical boundaries
      for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok(); ++dit)
        {
          const Box interiorBox = m_levelGridMetrics.getBoxes()[dit];
          FArrayBox& JUFab = a_JU[dit];
          const FArrayBox& UFab = a_U[dit];
          // const FArrayBox& JFab = m_levelGridMetrics.m_J[dit];
          const FArrayBox& JFab = a_J[dit];
          fourthOrderCellProd(JUFab, UFab, JFab, interiorBox,
                              m_multiBlockCoordSys->problemDomain(interiorBox),
                              true);
        }
      if (dataStatus != 0)
        {
          CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain all "
            "state data at level " << m_level << '!' << CRD::error;
        }
    }
  else if (CRDparam::g_restartAddWallModel)
    {
      // JU
      const int numReadStates =
        CRDparam::g_CRDPhysics->numConservative() - 1 - SpaceDim;
      LevelData<FArrayBox>& JU = m_data.rawJU();
      // Set the new values to zero for now
      for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok();
           ++dit)
        {
          FArrayBox& JUFab = JU[dit];
          const int startNewComp = numReadStates;
          const int numNewComp = 1 + SpaceDim;
          for (int newComp = 0; newComp != numNewComp; ++newComp)
            {
              JUFab.setVal(0., newComp+startNewComp);
            }
        }
      Interval readIntv(0, numReadStates - 1);
      const int dataStatus = read<FArrayBox>(a_handle,
                                             JU,
                                             "data",
                                             m_boxes,
                                             readIntv,
                                             false);  // Do not redefine JU
      if (dataStatus != 0)
        {
          CRD::msg << "AMRLevelCNS::readCheckpointLevel:"
                   << " file does not contain all state data at level "
                   << m_level << '!' << CRD::error;
        }
    }
  else
    {
      // JU
      LevelData<FArrayBox>& JU = m_data.rawJU();
      const int dataStatus = read<FArrayBox>(a_handle,
                                             JU,
                                             "data",
                                             m_boxes,
                                             JU.interval(),
                                             false);  // Do not redefine JU
      if (dataStatus != 0)
        {
          CRD::msg << "AMRLevelCNS::readCheckpointLevel: file does not contain all "
            "state data at level " << m_level << '!' << CRD::error;
        }
    }

  m_data.invalidate();  // Only know <JU> on valid real cells
}

/*--------------------------------------------------------------------*/
//  Write plotfile header
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCNS::writePlotHeader");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::writePlotHeader (level: " << m_level
           << ")" << CRD::end;

  // Setup the number of components
  const int numOutStates = CRDparam::g_CRDPhysics->numOutputVar();
  // Setup the number of components for the J-components
  const int numStates = CRDparam::g_CRDPhysics->numConservative();
  HDF5HeaderData header;
  int numStatesOutput = numOutStates;  // U
  if (CRDparam::g_plotJU)
    {
      numStatesOutput += numStates;    // JU
    }
  if (CRDparam::g_plotJ)
    {
      ++numStatesOutput;               // J
    }
  if (CRDparam::g_plotMappedDerivs)
    {
      ++numStatesOutput;               // Magnitude of gradient of density
      ++numStatesOutput;               // divergence of velocity
      numStatesOutput += PatchMappedFunc::m_numCurlComps;
                                       // Vorticity (curl of velocity)
    }
  if (CRDparam::g_plotWallDist && m_levelGridMetrics.distanceDefined())
    {
      ++numStatesOutput;               // distance field
    }
  if (CRDparam::g_plotFlattening || m_levelOp.numDebugPlotVar() > 0)
    {
      numStatesOutput +=               // Flattening/Debug
        (int)(CRDparam::g_plotFlattening) + m_levelOp.numDebugPlotVar();
    }
  if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
    {
      ++numStatesOutput;               // Subgrid-scale kinetic energy estimate
    }
  if (CRDparam::g_plotError)
    {
      numStatesOutput += numStates;    // Error of conservative state with
                                       // respect to exact solution
    }
  if (CRDparam::g_plotTimeAvgTurb)
    {
      numStatesOutput += 1;            // Wall shear-stress from extrapolation
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          numStatesOutput += 1;        // Wall shear-stress from LES model
        }
      numStatesOutput += SpaceDim*SpaceDim;
                                       // Time-averaged velocity on the low-face
                                       // of each cell in all directions
                                       // These are face-averaged values
      numStatesOutput += SpaceDim*SpaceDim;
                                       // Time-averaged velocity on the low-face
                                       // of each cell in all directions
                                       // These are face-centered values
      numStatesOutput += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
                                       // Time-avg velocity stress tensor on
                                       // low-face of each cell in all dirs
                                       // These are face-averaged values
      numStatesOutput += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
                                       // Time-avg velocity stress tensor on
                                       // low-face of each cell in all dirs
                                       // These are face-centered values

      // Time-avg density and temperature in the cells
      numStatesOutput += 1 + 1;

      // New components
      numStatesOutput += 2; // rho*u and rho*u*u
      numStatesOutput += 1; // <tau_{wall-model}^2>
      // Time-avg <SGS> stress tensor on x-faces
      numStatesOutput += (SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
      numStatesOutput += (SpaceDim-1); // Actual wall shear-stress
                                       // Momentum flux on wall -- only use the
                                       // wall-tangential components
    }
  if (CRDparam::g_plotLoFaceCoordinates)
    {
      numStatesOutput += SpaceDim*SpaceDim;
                                       // SpaceDim number of coordinates times
                                       // SpaceDim number of faces
    }
  if (CRDparam::g_plotHiFaceCoordinates)
    {
      numStatesOutput += SpaceDim*SpaceDim;
                                       // SpaceDim number of coordinates times
                                       // SpaceDim number of faces
    }
  if (CRDparam::g_plotLoFaceAvgComps)
    {
      int numAvgComps = 0;
      int plotVariables = CRDparam::g_plotLoFaceAvgComps;
      plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          numAvgComps += CRDparam::g_CRDPhysics->numPrimitive();
        }
      if (plotVariables & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
        }
      if (plotVariables & CRDparam::PlotTurbulentComps)
        {
          numAvgComps += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        }
      numStatesOutput += SpaceDim*numAvgComps;
    }
  if (CRDparam::g_plotLoFaceAvgTimeAvgComps)
    {
      int numAvgComps = 0;
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotPrimitive)
        {
          numAvgComps += CRDparam::g_CRDPhysics->numPrimitive();
        }
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
        }
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotTurbulentComps)
        {
          numAvgComps += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        }
      numStatesOutput += SpaceDim*numAvgComps;
    }
  header.m_int["num_components"] = numStatesOutput;

//--Set up the component names

  char compStr[64];
  int comploc = 0;
  // U
  for (int comp = 0; comp != numOutStates; ++comp)
    {
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] =
        CRDparam::g_CRDPhysics->consvStateName(comp);
      ++comploc;
    }
  // JU
  if (CRDparam::g_plotJU)
    {
      char nameStr[64];
      for (int comp = 0; comp != numStates; ++comp)
        {
          sprintf(compStr, "component_%d", comploc);
          sprintf(nameStr,
                  "J-%s",
                  CRDparam::g_CRDPhysics->consvStateName(comp));
          header.m_string[compStr] = nameStr;
          ++comploc;
        }
    }
  // J
  if (CRDparam::g_plotJ)
    {
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "J";
      ++comploc;
    }
  // Derivatives
  if (CRDparam::g_plotMappedDerivs)
    {
      // magnitude of gradient of density
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "grad-density_magnitude";
      ++comploc;

      // divergence
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "div-velocity";
      ++comploc;

      // vorticity
      switch (SpaceDim)
        {
        case 2:
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "vorticity";
          ++comploc;
          break;
        case 3:
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "x-vorticity";
          ++comploc;
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "y-vorticity";
          ++comploc;
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "z-vorticity";
          ++comploc;
          break;
        }
    }
  // Wall distance
  if (CRDparam::g_plotWallDist && m_levelGridMetrics.distanceDefined())
    {
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "wall-distance";
      ++comploc;
    }
  // Flattening/Debug
  if (CRDparam::g_plotFlattening || m_levelOp.numDebugPlotVar() > 0)
    {
      if (CRDparam::g_plotFlattening)
        {
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "flattening-coef";
          ++comploc;
        }
      char nameStr[64];
      for (int c = 0, cend = m_levelOp.numDebugPlotVar(); c != cend; ++c)
        {
          sprintf(compStr, "component_%d", comploc);
          sprintf(nameStr, "debug%d", c);
          header.m_string[compStr] = nameStr;
          ++comploc;
        }
    }
  // Subgrid-scale kinetic energy estimate
  if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
    {
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "Subgrid-scale-kinetic-energy";
      ++comploc;
    }
  // Error of conservative state with respect to exact solution
  if (CRDparam::g_plotError)
    {
      char nameStr[64];
      for (int comp = 0; comp != numStates; ++comp)
        {
          sprintf(compStr, "component_%d", comploc);
          sprintf(nameStr,
                  "Error-%s",
                  CRDparam::g_CRDPhysics->consvStateName(comp));
          header.m_string[compStr] = nameStr;
          ++comploc;
        }
    }
  // Time-averaged state
  if (CRDparam::g_plotTimeAvgTurb)
    {
      char nameStr[64];
      // Wall shear-stress from extrapolation
      sprintf(compStr, "component_%d", comploc);
      header.m_string[compStr] = "TimeAvg_extrap_wall_shear-stress";
      ++comploc;
      // Wall shear-stress from LES model
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          sprintf(compStr, "component_%d", comploc);
          header.m_string[compStr] = "TimeAvg_model_wall_shear-stress";
          ++comploc;
        }
      // velocity components
      const char* velName[] = {"u","v","w"};
      const char* dirName[] = {"x","y","z"};
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "TimeAvg-FaceAvg-%sVel-%sFace",
                      dirName[dir], dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
      // velocity components
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "TimeAvg-FacePnt-%sVel-%sFace",
                      dirName[dir], dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
      // Time-avg velocity stress tensor on low-face of each cell in all dirs
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int stressComp = dir; stressComp != SpaceDim; ++stressComp)
                {
                  sprintf(compStr, "component_%d", comploc);
                  sprintf(nameStr, "TimeAvg-FaceAvg-%s%s-%sFace",
                          velName[dir], velName[stressComp], dirName[faceDir]);
                  header.m_string[compStr] = nameStr;
                  ++comploc;
                }
            }
        }
      // Time-avg velocity stress tensor on low-face of each cell in all dirs
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int stressComp = dir; stressComp != SpaceDim; ++stressComp)
                {
                  sprintf(compStr, "component_%d", comploc);
                  sprintf(nameStr, "TimeAvg-FacePnt-%s%s-%sFace",
                          velName[dir], velName[stressComp], dirName[faceDir]);
                  header.m_string[compStr] = nameStr;
                  ++comploc;
                }
            }
        }

      // Time-avg density and temperature in the cells
      sprintf(compStr, "component_%d", comploc);
      sprintf(nameStr, "TimeAvg-CellAvg-density");
      header.m_string[compStr] = nameStr;
      ++comploc;
      sprintf(compStr, "component_%d", comploc);
      sprintf(nameStr, "TimeAvg-CellAvg-pressure");
      header.m_string[compStr] = nameStr;
      ++comploc;

      // Time-avg rho*u and rho*u*u
      sprintf(compStr, "component_%d", comploc);
      sprintf(nameStr, "TimeAvg-rhoU-xFace");
      header.m_string[compStr] = nameStr;
      ++comploc;
      sprintf(compStr, "component_%d", comploc);
      sprintf(nameStr, "TimeAvg-rhoUU-xFace");
      header.m_string[compStr] = nameStr;
      ++comploc;

      // Time-avg eta^2
      sprintf(compStr, "component_%d", comploc);
      sprintf(nameStr, "TimeAvg-eta-Sqrd");
      header.m_string[compStr] = nameStr;
      ++comploc;

      // Time-avg sgs-stresses
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          for (int stressComp = dir; stressComp != SpaceDim; ++stressComp)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "TimeAvg-FaceAvg-SGS-%s%s-xFace",
                      velName[dir], velName[stressComp]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }

      // Time-avg actual wall shear-stress
      for (int dir = 0; dir != (SpaceDim - 1); ++dir)
        {
          sprintf(compStr, "component_%d", comploc);
          sprintf(nameStr, "TimeAvg-actualWallTau-dir-%d", dir);
          header.m_string[compStr] = nameStr;
          ++comploc;
        }
    }
  // Low-side face coordinates
  if (CRDparam::g_plotLoFaceCoordinates)
    {
      char nameStr[64];
      const char* dirName[] = {"x","y","z"};
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int coordComp = 0; coordComp != SpaceDim; ++coordComp)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "Lo-%s-Face-%s-Coord",
                      dirName[faceDir], dirName[coordComp]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
    }
  // High-side face coordinates
  if (CRDparam::g_plotHiFaceCoordinates)
    {
      char nameStr[64];
      const char* dirName[] = {"x","y","z"};
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          for (int coordComp = 0; coordComp != SpaceDim; ++coordComp)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "Hi-%s-Face-%s-Coord",
                      dirName[faceDir], dirName[coordComp]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
    }
  // Low-side face-averaged variables
  if (CRDparam::g_plotLoFaceAvgComps)
    {
      char nameStr[64];
      const char* dirName[] = {"x","y","z"};
      const char* velName[] = {"u","v","w"};
      int plotVariables = CRDparam::g_plotLoFaceAvgComps;
      plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          const int numPrim = CRDparam::g_CRDPhysics->numPrimitive();
          for (int primComp = 0; primComp != numPrim; ++primComp)
            {
              for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                {
                  sprintf(compStr, "component_%d", comploc);
                  sprintf(nameStr, "low-%s-face-avg-%s", dirName[faceDir],
                          CRDparam::g_CRDPhysics->primStateName(primComp));
                  header.m_string[compStr] = nameStr;
                  ++comploc;
                }
            }
        }
      if (plotVariables & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              const int numFlux = CRDparam::g_CRDPhysics->numFluxes();
              for (int fluxComp = 0; fluxComp != numFlux; ++fluxComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-%s-Inertialflux",
                              dirName[faceDir],
                              CRDparam::g_CRDPhysics->consvStateName(fluxComp));
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              const int numFlux = CRDparam::g_CRDPhysics->numFluxes();
              for (int fluxComp = 0; fluxComp != numFlux; ++fluxComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-%s-Viscousflux",
                              dirName[faceDir],
                              CRDparam::g_CRDPhysics->consvStateName(fluxComp));
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
        }
      if (plotVariables & CRDparam::PlotTurbulentComps)
        {
          // Low-side face-averaged velocity stress tensor in all dirs
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-%s%s", dirName[faceDir],
                              velName[dir], velName[sComp]);
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          // Low-side face-averaged p^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-pressSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged temperature^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tempSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged rho^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-rhoSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged SGS stress tensor in all dirs
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-SGS-%s%s",
                              dirName[faceDir],
                              velName[dir], velName[sComp]);
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          // Low-side face-averaged magnitude of wall shear-stress
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tau_w_mag", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged streamwise wall shear-stress
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tau_w_streamwise",
                      dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged magnitude of wall shear-stress squared
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tau_w_sqrd", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
    }
  // Low-side face-averaged time-averaged variables
  if (CRDparam::g_plotLoFaceAvgTimeAvgComps)
    {
      char nameStr[64];
      const char* dirName[] = {"x","y","z"};
      const char* velName[] = {"u","v","w"};
      const int plotVariables = CRDparam::g_plotLoFaceAvgTimeAvgComps;
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          const int numPrim = CRDparam::g_CRDPhysics->numPrimitive();
          for (int primComp = 0; primComp != numPrim; ++primComp)
            {
              for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                {
                  sprintf(compStr, "component_%d", comploc);
                  sprintf(nameStr, "low-%s-face-avg-tavg-%s", dirName[faceDir],
                          CRDparam::g_CRDPhysics->primStateName(primComp));
                  header.m_string[compStr] = nameStr;
                  ++comploc;
                }
            }
        }
      if (plotVariables & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              const int numFlux = CRDparam::g_CRDPhysics->numFluxes();
              for (int fluxComp = 0; fluxComp != numFlux; ++fluxComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-tavg-%s-Inertialflux",
                              dirName[faceDir],
                              CRDparam::g_CRDPhysics->consvStateName(fluxComp));
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              const int numFlux = CRDparam::g_CRDPhysics->numFluxes();
              for (int fluxComp = 0; fluxComp != numFlux; ++fluxComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-tavg-%s-Viscousflux",
                              dirName[faceDir],
                              CRDparam::g_CRDPhysics->consvStateName(fluxComp));
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
        }
      if (plotVariables & CRDparam::PlotTurbulentComps)
        {
          // Low-side face-averaged velocity stress tensor in all dirs
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-tavg-%s%s",
                              dirName[faceDir],
                              velName[dir], velName[sComp]);
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          // Low-side face-averaged p^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-pressSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged temperature^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tempSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged rho^2 in all dirs
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-rhoSq", dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged SGS stress tensor in all dirs
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                    {
                      sprintf(compStr, "component_%d", comploc);
                      sprintf(nameStr, "low-%s-face-avg-tavg-SGS-%s%s",
                              dirName[faceDir],
                              velName[dir], velName[sComp]);
                      header.m_string[compStr] = nameStr;
                      ++comploc;
                    }
                }
            }
          // Low-side face-averaged magnitude of wall shear-stress
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tavg-tau_w_mag",
                      dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged streamwise wall shear-stress
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tau_w_streamwise",
                      dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
          // Low-side face-averaged magnitude of wall shear-stress squared
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              sprintf(compStr, "component_%d", comploc);
              sprintf(nameStr, "low-%s-face-avg-tavg-tau_w_sqrd",
                      dirName[faceDir]);
              header.m_string[compStr] = nameStr;
              ++comploc;
            }
        }
    }

  // Sanity check
  CH_assert(comploc == numStatesOutput);

  // Write the header
  header.writeToFile(a_handle);
  CRD::msg << CRD::fv3 << header << CRD::end;

  // Write expressions
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  CRDparam::g_CRDPhysics->expressions(expressions);
  expressions.writeToFile(a_handle);
}

/*--------------------------------------------------------------------*/
//  Write plotfile data for this level
/** \param[in]  a_handle
 *                      Chombo class that holds a number of hid_t
 *                      objects.  Essentially the dataset to write.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_TIME("AMRLevelCNS::writePlotLevel");
  CRD::msg << CRD::fv2 << "AMRLevelCNS::writePlotLevel (level: " << m_level
           << ")" << CRD::end;

  // Write out the mapped-grid geometry info (only valid at level 0)
  if (m_level == 0)
    {
      writeMappedPlotFile();
    }

  // Setup the level string
  char levelStr[32];
  sprintf(levelStr, "%d", m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  // Internal data
  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];  //**FIXME
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] =
    m_multiBlockCoordSys->levelDomain().domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);
  CRD::msg << CRD::fv3 << header << CRD::end;

  // The number of ghosts we output
  const IntVect outputGhost = CRDparam::g_plotNumGhost*IntVect::Unit;
  // Number of outputs for the plot files
  const int numOutStates = CRDparam::g_CRDPhysics->numOutputVar();
  // Number of J-component values
  const int numStates = CRDparam::g_CRDPhysics->numConservative();

  // Components for U
  LevelData<FArrayBox> U(m_data.getU(1, 1).getBoxes(),
                         numOutStates,
                         outputGhost);
  CRDparam::g_CRDPhysics->outputLevelData(U,
                                          m_data.getU(1, 1),
                                          m_WOld,
                                          m_levelGridMetrics);
  int numCompWrite = U.nComp();
  CH_assert(numCompWrite == numOutStates);

  // Do species correction which modifies U
  m_levelOp.speciesCorrection(U);

  // Components for JU
  LevelData<FArrayBox> *JU = 0;
  if (CRDparam::g_plotJU)
    {
      // Apply species correction to JU as well
      JU->define(m_data.getJU(1,1).getBoxes(),
                 m_data.getJU(1,1).nComp(),
                 m_data.getJU(1,1).ghostVect());
      m_data.getJU(1,1).copyTo(*JU);
      m_levelOp.speciesCorrection(*JU);
      const int numJUcomp = JU->nComp();
      CH_assert(numJUcomp == numStates);
      numCompWrite += numJUcomp;
    }

  // Components for J
  if (CRDparam::g_plotJ)
    {
      CH_assert(m_levelGridMetrics.m_J.ghostVect() >= outputGhost);
      CH_assert(m_levelGridMetrics.m_J.nComp() == 1);
      numCompWrite += 1;
    }

  // Derivatives
  if (CRDparam::g_plotMappedDerivs)
    {
      // magnitude of gradient of density
      ++numCompWrite;
      // divergence
      ++numCompWrite;
      // vorticity
      numCompWrite += PatchMappedFunc::m_numCurlComps;
    }

  // Wall distance
  if (CRDparam::g_plotWallDist && m_levelGridMetrics.distanceDefined())
    {
      CH_assert(m_levelGridMetrics.m_distance.ghostVect() >= outputGhost);
      CH_assert(m_levelGridMetrics.m_distance.nComp() == 1);
      numCompWrite += 1;
    }

  // Flattening coefficients and extra debug variables
  if (CRDparam::g_plotFlattening || m_levelOp.numDebugPlotVar() > 0)
    {
      numCompWrite +=
        (int)(CRDparam::g_plotFlattening) + m_levelOp.numDebugPlotVar();
    }

  // Subgrid-scale kinetic energy estimate
  if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
    {
      ++numCompWrite;
    }

  // Error of conservative state with respect to exact solution
  if (CRDparam::g_plotError)
    {
      numCompWrite += numStates;
    }

  // Time-averaged solution data
  if (CRDparam::g_plotTimeAvgTurb)
    {
      numCompWrite += 1;               // Wall shear-stress from extrapolation
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          numCompWrite += 1;           // Wall shear-stress from LES model
        }
      numCompWrite += SpaceDim*SpaceDim;
                                       // Time-averaged velocity on the low-face
                                       // of each cell in all directions
                                       // These are face-averaged values
      numCompWrite += SpaceDim*SpaceDim;
                                       // Time-averaged velocity on the low-face
                                       // of each cell in all directions
                                       // These are face-centered values
      numCompWrite += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
                                       // Time-avg velocity stress tensor on
                                       // low-face of each cell in all dirs
                                       // These are face-averaged values
      numCompWrite += SpaceDim*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
                                       // Time-avg velocity stress tensor on
                                       // low-face of each cell in all dirs
                                       // These are face-centered values

      // Time-avg density and pressure in the cells
      numCompWrite += 1 + 1;

      // New components
      numCompWrite += 2; // rho*u and rho*u*u
      numCompWrite += 1; // <tau_{wall-model}^2>
      // Time-avg <SGS> stress tensor on x-faces
      numCompWrite += (SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
      numCompWrite += (SpaceDim-1); // Actual wall shear-stress
                                    // Momentum flux on wall -- only use the
                                    // wall-tangential components
    }

  // Low-side face coordinates
  if (CRDparam::g_plotLoFaceCoordinates)
    {
      numCompWrite += SpaceDim*SpaceDim;
                                       // SpaceDim number of coordinates times
                                       // SpaceDim number of faces
    }

  // High-side face coordinates
  if (CRDparam::g_plotHiFaceCoordinates)
    {
      numCompWrite += SpaceDim*SpaceDim;
                                       // SpaceDim number of coordinates times
                                       // SpaceDim number of faces
    }

  // Low-side face-averaged variables
  if (CRDparam::g_plotLoFaceAvgComps)
    {
      int numAvgComps = 0;
      int plotVariables = CRDparam::g_plotLoFaceAvgComps;
      plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          numAvgComps += CRDparam::g_CRDPhysics->numPrimitive();
        }
      if (plotVariables & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
        }
      if (plotVariables & CRDparam::PlotTurbulentComps)
        {
          numAvgComps += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        }
      numCompWrite += SpaceDim*numAvgComps;
    }

  // Low-side face-averaged time-averaged variables
  if (CRDparam::g_plotLoFaceAvgTimeAvgComps)
    {
      int numAvgComps = 0;
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotPrimitive)
        {
          numAvgComps += CRDparam::g_CRDPhysics->numPrimitive();
        }
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotFluxes)
        {
          if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
          if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
            {
              numAvgComps += CRDparam::g_CRDPhysics->numFluxes();
            }
        }
      if (CRDparam::g_plotLoFaceAvgTimeAvgComps & CRDparam::PlotTurbulentComps)
        {
          numAvgComps += 6 + 2*(SpaceDim*(SpaceDim - 1.)/2. + SpaceDim);
        }
      numCompWrite += SpaceDim*numAvgComps;
    }

//--Begin writing data

  // Write the boxes and construct the dataset for multiple levels
  WriteMultiData<FArrayBox> writer(a_handle,
                                   U.getBoxes(),
                                   numCompWrite,
                                   "data",
                                   CH_HDF5::IOPolicyDefault,
                                   outputGhost);


  int iCompWrite = 0;
  // Write U
  {
    writer.writeData(U,
                     U.interval(),
                     U.interval());
    iCompWrite += U.nComp();
  }

  // Write JU
  if (CRDparam::g_plotJU)
    {
      writer.writeData(*JU,
                       JU->interval(),
                       Interval(iCompWrite, iCompWrite + JU->nComp() - 1));
      iCompWrite += JU->nComp();
    }

  // Write J
  if (CRDparam::g_plotJ)
    {
      writer.writeData(m_levelGridMetrics.m_J,
                       m_levelGridMetrics.m_J.interval(),
                       Interval(iCompWrite, iCompWrite));
      iCompWrite += 1;
    }

  // Write derivatives
  if (CRDparam::g_plotMappedDerivs)
    {
      const int vortComp = PatchMappedFunc::m_numCurlComps;
      // magGradDensity (1), divergence(1), vorticity(1 or 3)
      const int numMapDerivs = 1 + 1 + vortComp;
      LevelData<FArrayBox> mappedDerivs(m_boxes, numMapDerivs, outputGhost);
      int comp = 0;
      //FIXME: should consider using consToPrim to operate on velocity rather
      // than manually going from momentum to velocity

      // Write Magnitude of gradient of density
      LevelMappedFunc::magGradient2OPS(
        mappedDerivs,
        comp,
        outputGhost,
        m_data.getU(CRDparam::g_plotNumGhost + 1,
                    CRDparam::g_plotNumGhost + 1),
        Interval(CRDparam::g_CRDPhysics->densityIndex(),
                 CRDparam::g_CRDPhysics->densityIndex()),
        m_levelGridMetrics);
      ++comp;

      // Write divergence of velocity
      LevelMappedFunc::divergence2OPS(
        mappedDerivs,
        comp,
        outputGhost,
        m_data.getU(CRDparam::g_plotNumGhost + 1,
                    CRDparam::g_plotNumGhost + 1),
        CRDparam::g_CRDPhysics->velocityInterval(),
        m_levelGridMetrics,
        CRDparam::g_CRDPhysics->densityIndex());
      ++comp;

      // Write Vorticity
      LevelMappedFunc::curl2OPS(
        mappedDerivs,
        comp,
        outputGhost,
        m_data.getU(CRDparam::g_plotNumGhost + 1,
                    CRDparam::g_plotNumGhost + 1),
        CRDparam::g_CRDPhysics->velocityInterval(),
        m_levelGridMetrics,
        CRDparam::g_CRDPhysics->densityIndex());
      comp += vortComp;

      CH_assert(comp == numMapDerivs);
      writer.writeData(mappedDerivs,
                       mappedDerivs.interval(),
                       Interval(iCompWrite, iCompWrite + numMapDerivs - 1));
      iCompWrite += numMapDerivs;
    }

  // Write wall distance
  if (CRDparam::g_plotWallDist && m_levelGridMetrics.distanceDefined())
    {
      writer.writeData(m_levelGridMetrics.m_distance,
                       m_levelGridMetrics.m_distance.interval(),
                       Interval(iCompWrite, iCompWrite));
      iCompWrite += 1;
    }

  // Write flattening and extra debug variables
  if (CRDparam::g_plotFlattening || m_levelOp.numDebugPlotVar() > 0)
    {
      int intv0 = (int)(!CRDparam::g_plotFlattening !=
                        !CRDparam::g_useFlattening);
      Interval intv(intv0, intv0 + (int)CRDparam::g_plotFlattening +
                    m_levelOp.numDebugPlotVar() - 1);
      writer.writeData(m_levelOp.getDebugVar(),
                       intv,
                       Interval(iCompWrite, iCompWrite + intv.size() - 1));
      iCompWrite += intv.size();
    }

  // Write subgrid-scale kinetic energy estimate
  if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
    {
      writer.writeData(m_levelOp.getSGSKineticEnergyVar(),
                       Interval(0,0),
                       Interval(iCompWrite, iCompWrite));
      ++iCompWrite;
    }

  // Error of conservative state with respect to exact solution
  if (CRDparam::g_plotError)
    {
      LevelData<FArrayBox> consError(m_boxes, numStates, outputGhost);
      const LevelData<FArrayBox>& UData = m_data.getU(1,1);
      UData.copyTo(consError);
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = m_boxes[dit];
          const BlockDomain& domain =
            m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
          // Cell-averaged exact solution only required on disjointBox
          Box UavgBox = disjointBox;
          // Cell-centered exact solution required on disjointBox + 1
          Box UpntBox = grow(disjointBox, 1);
          UpntBox &= domain;

          // Get <Ux>
          FABSTACKTEMP(Ux, UpntBox, numStates);

          if (CRDparam::g_CNSIBC->haveExactSol())
            {
              const int stat = CRDparam::g_CNSIBC->exactSol(
                Ux,
                UavgBox,
                disjointBox,
                m_levelGridMetrics,
                m_levelOp.getUnitNormals()[dit],
                dit(),
                m_time,
                m_level);
              if (stat)
                {
                  CRD::msg << "Error obtaining exact solution" << CRD::error;
                }
            }
          else
            {
              Ux.setVal(0.);
            }

          // Compute the error
          consError[dit] -= Ux;
       }
      // Print the error
      Interval intv(0,numStates-1);
      writer.writeData(consError,
                       intv,
                       Interval(iCompWrite, iCompWrite + intv.size() - 1));
      iCompWrite += intv.size();
    }

  // Write time-averaged data
  if (CRDparam::g_plotTimeAvgTurb)
    {
      int numTimeAvgStates = m_levelOp.getTimeAvgVar().nComp();
      writer.writeData(m_levelOp.getTimeAvgVar(),
                       Interval(0,numTimeAvgStates-1),
                       Interval(iCompWrite, iCompWrite+numTimeAvgStates-1));
      iCompWrite += numTimeAvgStates;
    }

  // Low-side face coordinates
  if (CRDparam::g_plotLoFaceCoordinates)
    {
      const int numFaceCoordinates = SpaceDim*SpaceDim;
      LevelData<FArrayBox> faceCoords(m_boxes, numFaceCoordinates, outputGhost);
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          // Compute the low-side face-coordinates
          for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
            {
              const Box disjointBox = m_boxes[dit];
              const BlockDomain& domain =
                m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
              const BlockCoordSys& blockCoordSys =
                *(m_levelGridMetrics.getCoordSys(disjointBox));
              // Get the cells on which the face coordinates will be computed
              Box coordBox = grow(disjointBox, outputGhost);
              coordBox &= domain;
              // Get the face box on which to compute the coordinates
              Box coordFaceBox = coordBox;
              coordFaceBox.surroundingNodes(faceDir);
              FABSTACKTEMP(faceXFab, coordFaceBox, SpaceDim);
              FABSTACKTEMP(faceXiFab, coordFaceBox, SpaceDim);
              CRDparam::g_CNSIBC->getFaceCoordinates(
                coordFaceBox, faceXiFab, faceXFab, faceDir, blockCoordSys);
              // Now, shift the low faces into the cells
              faceXFab.shiftHalf(faceDir, 1);
              int destComp = faceDir*SpaceDim;
              faceCoords[dit].copy(
                faceXFab, coordBox, 0, coordBox, destComp, SpaceDim);
            }
        }
      writer.writeData(faceCoords,
                       Interval(0,numFaceCoordinates-1),
                       Interval(iCompWrite, iCompWrite+numFaceCoordinates-1));
      iCompWrite += numFaceCoordinates;
    }

  // High-side face coordinates
  if (CRDparam::g_plotHiFaceCoordinates)
    {
      const int numFaceCoordinates = SpaceDim*SpaceDim;
      LevelData<FArrayBox> faceCoords(m_boxes, numFaceCoordinates, outputGhost);
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          // Compute the low-side face-coordinates
          for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
            {
              const Box disjointBox = m_boxes[dit];
              const BlockDomain& domain =
                m_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
              const BlockCoordSys& blockCoordSys =
                *(m_levelGridMetrics.getCoordSys(disjointBox));
              // Get the cells on which the face coordinates will be computed
              Box coordBox = grow(disjointBox, outputGhost);
              coordBox &= domain;
              // Get the face box on which to compute the coordinates
              Box coordFaceBox = coordBox;
              coordFaceBox.surroundingNodes(faceDir);
              FABSTACKTEMP(faceXFab, coordFaceBox, SpaceDim);
              FABSTACKTEMP(faceXiFab, coordFaceBox, SpaceDim);
              CRDparam::g_CNSIBC->getFaceCoordinates(
                coordFaceBox, faceXiFab, faceXFab, faceDir, blockCoordSys);
              // Now, shift the high faces into the cells
              faceXFab.shiftHalf(faceDir, -1);
              int destComp = faceDir*SpaceDim;
              faceCoords[dit].copy(
                faceXFab, coordBox, 0, coordBox, destComp, SpaceDim);
            }
        }
      writer.writeData(faceCoords,
                       Interval(0,numFaceCoordinates-1),
                       Interval(iCompWrite, iCompWrite+numFaceCoordinates-1));
      iCompWrite += numFaceCoordinates;
    }

  // Write low-side face-averaged variables
  if (CRDparam::g_plotLoFaceAvgComps)
    {
      int numAvgComp = m_levelOp.getFaceAvgVar().nComp();
      int numPrintComp = SpaceDim*numAvgComp;
      LevelData<FArrayBox> loFaceAvgData(m_boxes, numPrintComp, outputGhost);
      const LevelData<FluxBox>& avgData = m_levelOp.getFaceAvgVar();
      // Compute the low-side face-coordinates
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = m_boxes[dit];

          // Fill loFaceAvgData with avgData (low-face data into cells)
          FArrayBox& loFaceAvgDataFab = loFaceAvgData[dit];
          for (int avgComp = 0; avgComp != numAvgComp; ++avgComp)
            {
              for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                {
                  const int printComp = faceDir + SpaceDim*avgComp;
                  const FArrayBox& avgDataFab = avgData[dit][faceDir];
                  MD_BOXLOOP(disjointBox, i)
                    {
                      loFaceAvgDataFab[MD_IX(i, printComp)] =
                        avgDataFab[MD_IX(i, avgComp)];
                    }
                }
            }
        }
      writer.writeData(loFaceAvgData,
                       Interval(0,numPrintComp-1),
                       Interval(iCompWrite, iCompWrite + numPrintComp - 1));
      iCompWrite += numPrintComp;
    }

  // Write low-side face-averaged time-averaged variables
  if (CRDparam::g_plotLoFaceAvgTimeAvgComps)
    {
      int numAvgComp = m_levelOp.getFaceAvgTimeAvgVar().nComp();
      int numPrintComp = SpaceDim*numAvgComp;
      LevelData<FArrayBox> loFaceAvgData(m_boxes, numPrintComp, outputGhost);
      const LevelData<FluxBox>& avgData = m_levelOp.getFaceAvgTimeAvgVar();
      // Compute the low-side face-coordinates
      for (DataIterator dit = m_boxes.dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = m_boxes[dit];

          // Fill loFaceAvgData with avgData (low-face data into cells)
          FArrayBox& loFaceAvgDataFab = loFaceAvgData[dit];
          for (int avgComp = 0; avgComp != numAvgComp; ++avgComp)
            {
              for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
                {
                  const int printComp = faceDir + SpaceDim*avgComp;
                  const FArrayBox& avgDataFab = avgData[dit][faceDir];
                  MD_BOXLOOP(disjointBox, i)
                    {
                      loFaceAvgDataFab[MD_IX(i, printComp)] =
                        avgDataFab[MD_IX(i, avgComp)];
                    }
                }
            }
        }
      writer.writeData(loFaceAvgData,
                       Interval(0,numPrintComp-1),
                       Interval(iCompWrite, iCompWrite + numPrintComp - 1));
      iCompWrite += numPrintComp;
    }

  // Sanity check
  CH_assert(iCompWrite == numCompWrite);
}

/*--------------------------------------------------------------------*/
//  Write mapped-grid info
/** This routine writes the mapping between computational and physical
 *  node locations.  It is called for level 0 and writes the mapping
 *  for all levels
 *  \note
 *  <ul>
 *    <li> FIXME - This should integrated into the regular plot files
 *    <li> FIXME - This should be done per level
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::writeMappedPlotFile() const
{
  CH_TIME("AMRLevelCNS::writeMappedPlotFile");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::writeMappedPlotFile" << CRD::end;

  // Only done on level 0
  CH_assert(m_level == 0);

  // Gather AMR Levels and create node-centered dataset of node locations of
  // mapped grids

  Vector<AMRLevel*> vAMRLevels;
  {
    // Cast away const to call this function
    vAMRLevels = const_cast<AMRLevelCNS*>(this)->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();

  Vector<int> vRefRatios(numLevels, 0);
  Vector<DisjointBoxLayout> vGrids(numLevels);
  Vector<LevelData<NodeFArrayBox>*> vNodeLoc(numLevels, NULL);

//--Prepare the data

  for (int iLev = 0; iLev < numLevels; ++iLev)
    {
      AMRLevelCNS* AMRlevelPtr = static_cast<AMRLevelCNS*>(vAMRLevels[iLev]);
      vRefRatios[iLev] = AMRlevelPtr->m_ref_ratio;
      vGrids[iLev] = AMRlevelPtr->m_boxes;
      vNodeLoc[iLev] = new LevelData<NodeFArrayBox>(
        vGrids[iLev], SpaceDim, CRDparam::g_plotNumGhost*IntVect::Unit);

      const RealVect& levelDx = AMRlevelPtr->m_dx;

      for (DataIterator dit = vGrids[iLev].dataIterator(); dit.ok(); ++dit)
        {
          const Box& box = vGrids[iLev][dit];
          const BlockCoordSys& blockCoordSys =
            *(AMRlevelPtr->m_levelGridMetrics.getCoordSys(box));

          NodeFArrayBox& nodeFab = (*vNodeLoc[iLev])[dit];
          FArrayBox& XFab = nodeFab.getFab();
          const Box& nodeBox = XFab.box();
          FABSTACKTEMP(XiFab, nodeBox, SpaceDim);
          // Set location in computational space
          FORT_SETCORNERSVEC(CHF_FRA(XiFab),
                             CHF_CONST_REALVECT(levelDx),
                             CHF_BOX(nodeBox));
          // Get location in physical space
          blockCoordSys.realCoord(XFab, XiFab, nodeBox);
        }
    }

//--Write the data

  // Create names
  Vector<string> varNames(SpaceDim);
  D_TERM(varNames[0] = "x";,
         varNames[1] = "y";,
         varNames[2] = "z";)

  // Create filename
  char suffix[64];
  sprintf(suffix, "%06d.%dd.map.hdf5", m_levelStep, SpaceDim);
  string fileName(CRDparam::g_plotPrefix);
  fileName += suffix;

  // Now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(fileName,
                        vGrids,
                        vNodeLoc,
                        varNames,
                        m_multiBlockCoordSys->levelDomain().domainBox(),
                        m_dx[0],  //**FIXME dx
                        m_dt,
                        m_time,
                        vRefRatios,
                        numLevels);

  // Clean up memory
  for (int iLev = 0; iLev < numLevels; ++iLev)
    {
      delete vNodeLoc[iLev];
    }
}
#endif

/*--------------------------------------------------------------------*/
//  Conclude by testing conservation and reporting error if necessary
/** \param[in]  a_step  The last step in the simulation
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::conclude(const int a_step) const
{
  CH_TIME("AMRLevelCNS::conclude");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::conclude (level: " << m_level
           << ")" << CRD::end;
  if (m_level == 0)
    {
      computeSum(s_JUConsvRef);
      if (CRDparam::g_CNSIBC->haveExactSol())
        {
          reportError();
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute dt
/** Computed previously.  Retrieve stored value
 *//*-----------------------------------------------------------------*/

Real
AMRLevelCNS::computeDt()
{
  CRD::msg << CRD::fv3 << "AMRLevelCNS::computeDt (level: " << m_level
           << ")" << CRD::end;
  return m_dtNew;
}

/*--------------------------------------------------------------------*/
//  Compute dt using initial data
/**
 *//*-----------------------------------------------------------------*/

Real
AMRLevelCNS::computeInitialDt()
{
  CRD::msg << CRD::fv3 << "AMRLevelCNS::computeInitialDt (level: " << m_level
           << ")" << CRD::end;

  // Need 1 layer of ghost cells in U
  if (CRDparam::g_initialDt == -1.)
    {
      m_dtNew = (CRDparam::g_initialCFL/CRDparam::g_cfl)*
        m_levelOp.computeNewDt(m_data.getU(1, 1),
                               m_time,
                               m_dt,
                               m_prevStepDt);
    }

  else
    {
      m_dtNew = CRDparam::g_initialDt;
    }
  return m_dtNew;
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Setup menagerie of data structures and operators for the level
/** Defines if coarser or finer grids exist.  Defines the
 *  averager, flux register, and level operator.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::levelOpSetup()
{
  CH_TIME("AMRLevelCNS::levelOpSetup");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::levelOpSetup (level: " << m_level
           << ")" << CRD::end;

  // In general, this can only look at coarser levels since during
  // initialization, the levels are traversed in order of coarse to fine.
  AMRLevelCNS* coarserAMRLevel = getCoarserLevel();

  // The AMRLevel::has* routines test for pointer, definition, and grids.
  m_hasCoarserGrid = AMRLevel::hasCoarserLevel();
  // During initialization, traversing coarse to fine, this will most likely be
  // set to false even if a finer level will be set up next.  So we fix up
  // what a coarser level thinks from the finer level (see below)
  m_hasFinerGrid   = AMRLevel::hasFinerLevel();
  // Fix up what a coarser level thinks from the finer level
  if (m_hasCoarserGrid)
    {
      coarserAMRLevel->m_hasFinerGrid = (m_boxes.size() > 0);
    }

  // The flux register is allocated if a finer level is possible (but not
  // necessarily used)
  if (hasFinerPtr())
    {
      m_fluxRegisterPtr->undefine();
    }

  // If there is no grid on this level, don't set up anything else
  if (m_boxes.size() == 0)
    {
      return;
    }

//--Set up some interlevel operators on the coarser level now that we know
//--a finer level (this one) exists.

  if (m_hasCoarserGrid)
    {
      // Maintain averager
      coarserAMRLevel->m_averageOp.define(
        m_boxes,
        coarserAMRLevel->m_boxes,
        CRDparam::g_CRDPhysics->numConservative(),
        coarserAMRLevel->refRatio());

      // Maintain flux registers
      coarserAMRLevel->m_fluxRegisterPtr->define(
        m_boxes,
        coarserAMRLevel->m_boxes,
        m_multiBlockCoordSys->levelDomain(),
        coarserAMRLevel->refRatio(),
        CRDparam::g_CRDPhysics->numConservative());
      coarserAMRLevel->m_fluxRegisterPtr->setToZero();

      // Maintain time interpolator
      if (m_timeInterpolator != nullptr)
        {
          delete m_timeInterpolator;
        }
      if (CRDparam::g_timeIntegrationMethod == CRDparam::ARK4)
        {
          m_timeInterpolator = new TimeInterpolatorARK4;
        }
      else if (CRDparam::g_timeIntegrationMethod == CRDparam::RK4)
        {
          m_timeInterpolator = new TimeInterpolatorRK4;
        }
      else if (CRDparam::g_timeIntegrationMethod == CRDparam::RK2)
        {
          m_timeInterpolator = new TimeInterpolatorRK2;
        }
      else
        {
          CRD::msg << "Error finding correct time interpolator to use."
                   << " I am looking for: " << CRDparam::g_timeIntegrationMethod
                   << " but did not find it." << CRD::error;
        }
      m_levelOp.setTimeInterpolatorPtr(m_timeInterpolator);

      if (m_levelGridMetrics.isMultiBlock())
        {
          /* This triggers a hack where the "fine" grid is faked to completely
             cover the coarse grid.  In other words, we are tricking
             m_timeInterpolator into not really doing a copyTo so that we can
             instead to a custom one later in FourthOrderMappedFineInterp.
             The custom one handles extra-block ghosts.
          */
          DisjointBoxLayout dummyFineGrids;
          refine(dummyFineGrids,
                 coarserAMRLevel->m_boxes,
                 coarserAMRLevel->refRatio());
          m_timeInterpolator->define(dummyFineGrids,
                                     coarserAMRLevel->m_boxes,
                                     m_multiBlockCoordSys->levelDomain(),
                                     coarserAMRLevel->refRatio(),
                                     CRDparam::g_CRDPhysics->numConservative(),
                                     0);  // Proper nesting means we shouldn't
                                          // need any ghosts
        }
      else
        {
          // Find the number of ghost cells on the coarsened-fine mesh required
          // to fill the invalid ghost cells in the fine mesh.
          const IntVect& interpolatorCrFnGhostVect =
            coarserAMRLevel->m_levelGridMetrics.interpolatorCrFnNumGhost(true);
          // Assumes ghosts in all directions are equal
          const int numInterpolatorCrFnGhost = interpolatorCrFnGhostVect[0];
          //**FIXME dxVect
          m_timeInterpolator->define(m_boxes,
                                     coarserAMRLevel->m_boxes,
                                     m_multiBlockCoordSys->levelDomain(),
                                     coarserAMRLevel->refRatio(),
                                     CRDparam::g_CRDPhysics->numConservative(),
                                     numInterpolatorCrFnGhost);
        }
    }

//--Define the level operator

  if(CRDparam::g_additiveRK)
    {
      m_ark_time_integrator.define(m_boxes,
                                   CRDparam::g_CRDPhysics->numConservative(),
                                   m_data.getJU().ghostVect());
    }

  m_levelOp.define(m_level, m_dx);

//--Define the persistent data for spectral forcing

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
  if (CRDparam::g_turbForcingType == CRDparam::TurbForcingSpectral)
    {
      m_spectralSource.define(m_boxes, 3, IntVect::Zero);
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Reverse-order setup of data structures and operators for the level
/** Allows defines in a fine-to-coarse traversal.  These traversals
 *  happen in postInitialGrid and postRegrid.
 *
 *  Example: MMBSingleLevel::define depends on MBRegions but the
 *  latter is defined during inter-level operations.  It is not
 *  defined on the coarser level until the finer level is accessed.
 *  In reverse, MMBRegions has been defined everywhere.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::reverseLevelOpSetup()
{
  CH_TIME("AMRLevelCNS::reverseLevelOpSetup");
  CRD::msg << CRD::fv3 << "AMRLevelCNS::reverseLevelOpSetup (level: " << m_level
           << ")" << CRD::end;

  m_levelOp.reverseDefine();
}

/*--------------------------------------------------------------------*/
//  Get the next coarser level
/** The base class stores AMRLevel objects so they must be cast to
 *  the derived type.
 *  \return             Pointer of type AMRLevelCNS to coarser level
 *                      or NULL if this is the coarsest level.
 *//*-----------------------------------------------------------------*/

AMRLevelCNS*
AMRLevelCNS::getCoarserLevel() const
{
  AMRLevelCNS* coarserPtr = NULL;
  if (hasCoarserPtr())
    {
      coarserPtr = dynamic_cast<AMRLevelCNS*>(m_coarser_level_ptr);
      if (coarserPtr == NULL)
        {
          CRD::msg << "AMRLevelCNS::getCoarserLevel(): dynamic cast failed"
                   << CRD::error;
        }
    }
  return coarserPtr;
}

/*--------------------------------------------------------------------*/
//  Get the next finer level
/** The base class stores AMRLevel objects so they must be cast to
 *  the derived type.
 *  \return             Pointer of type AMRLevelCNS to finer level
 *                      or NULL if this is the finest level.
 *//*-----------------------------------------------------------------*/

AMRLevelCNS*
AMRLevelCNS::getFinerLevel() const
{
  AMRLevelCNS* finerPtr = NULL;
  if (hasFinerPtr())
    {
      finerPtr = dynamic_cast<AMRLevelCNS*>(m_finer_level_ptr);
      if (finerPtr == NULL)
        {
          CRD::msg << "AMRLevelCNS::getFinerLevel(): dynamic cast failed"
                   << CRD::error;
        }
    }
  return finerPtr;
}

/*--------------------------------------------------------------------*/
//  Compute a sum of the components over all levels
/** \param[out] a_sum   Sum for each component
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::computeSum(Vector<Real>& a_sum) const
{
  CH_TIME("AMRLevelCNS::computeSum");
  CH_assert(m_level == 0);

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // Cast away const to call this function
    vAMRLevels = const_cast<AMRLevelCNS*>(this)->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();
  // Gather ref ratios
  Vector<int> vRefRatios(numLevels - 1);
  for (int iLev = 0; iLev < numLevels - 1; ++iLev)
    {
      vRefRatios[iLev] =
        static_cast<AMRLevelCNS*>(vAMRLevels[iLev])->refRatio();
    }

//--Sum of all components for testing conservation

  Vector<LevelData<FArrayBox>*> vPhi(numLevels);
  for (int iLev = 0; iLev != numLevels; ++iLev)
    {
      vPhi[iLev] = &(static_cast<AMRLevelCNS*>(vAMRLevels[iLev])
                     ->m_data.getJU());
    }
  const int numComp = CRDparam::g_CRDPhysics->numConservative();
  for (int iComp = 0; iComp != numComp; ++iComp)
    {
      const Interval compIntv(iComp, iComp);
      a_sum[iComp] = ::computeSum(vPhi,
                                  vRefRatios,
                                  1.,
                                  compIntv,
                                  0);
    }
  CRD::msg << CRD::fv2 << "AMRLevelCNS::computeSum: " << CRD::end;
  CRD::msg.setPrecFloatSN(15);
  for (int i = 0; i != a_sum.size(); ++i)
    {
      CRD::msg << CRD::fv2 << a_sum[i] << " ";
    }
  CRD::msg.setFloatDefault();
  CRD::msg << CRD::fv2 << CRD::end;
  return;
}

/*--------------------------------------------------------------------*/
//  Report norms of errors (requires exact solution)
/** \note
 *  <ul>
 *    <li> Error is printed to output.
 *    <li> The solution \<JU\> is replaced by the error in the
 *         process.
 *    <li> There is a capability for printing the location of the
 *         maximum error on each level.  Set 'maxLocComp' to a valid
 *         component index.  This is intended only for debugging.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::reportError() const
{
  CH_TIME("AMRLevelCNS::reportError");
  CH_assert(m_level == 0);
  const int maxLocComp = -1;  // Set this to a valid component index to receive
                              // detailed information on the location of the
                              // maximum error on each level

  // Gather AMR Levels
  Vector<AMRLevel*> vAMRLevels;
  {
    // Cast away const to call this function
    vAMRLevels = const_cast<AMRLevelCNS*>(this)->getAMRLevelHierarchy();
  }
  int numLevels = vAMRLevels.size();
  // Gather ref ratios
  Vector<int> vRefRatios(numLevels - 1);
  for (int iLev = 0; iLev < numLevels - 1; ++iLev)
    {
      vRefRatios[iLev] =
        static_cast<AMRLevelCNS*>(vAMRLevels[iLev])->refRatio();
    }

  if (!CRDparam::g_CNSIBC->haveExactSol())
    {
      CRD::msg << "No exact solution known.  Computing norm of solution "
        "instead of error." << CRD::warn;
    }

  Vector<LevelData<FArrayBox>*> vPhi(numLevels);
  for (int iLev = 0; iLev != numLevels; ++iLev)
    {
      AMRLevelCNS* AMRlevelPtr = static_cast<AMRLevelCNS*>(vAMRLevels[iLev]);
      vPhi[iLev] = &(AMRlevelPtr->m_data.getJU());

      IntVect maxLoc = IntVect::Zero;
      Real maxErr = 0.;
      for (DataIterator dit = vPhi[iLev]->dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = vPhi[iLev]->disjointBoxLayout()[dit];
          const BlockDomain& domain =
            AMRlevelPtr->m_multiBlockCoordSys->problemDomain(disjointBox);
          Box UavgBox = grow(disjointBox, 1);
          UavgBox &= domain;
          Box UpntBox = grow(disjointBox, 2);
          UpntBox &= domain;

          // Get <Ux>
          FABSTACKTEMP(Ux,
                       UpntBox,
                       CRDparam::g_CRDPhysics->numConservative());

          if (CRDparam::g_CNSIBC->haveExactSol())
            {
              const int stat = CRDparam::g_CNSIBC->exactSol(
                Ux,
                UavgBox,
                disjointBox,
                AMRlevelPtr->m_levelGridMetrics,
                AMRlevelPtr->m_levelOp.getUnitNormals()[dit],
                dit(),
                m_time,
                m_level);
              if (stat)
                {
                  CRD::msg << "Error obtaining exact solution" << CRD::error;
                }
            }
          else
            {
              Ux.setVal(0.);
            }

          // <JUx>
          FABSTACKTEMP(JUx,
                       disjointBox,
                       CRDparam::g_CRDPhysics->numConservative());
          const FArrayBox& JFab = AMRlevelPtr->m_levelGridMetrics.m_J[dit];
          fourthOrderCellProd(JUx,
                              Ux,
                              JFab,
                              disjointBox,
                              domain,
                              true);

          // Compute the error
          (*vPhi[iLev])[dit] -= JUx;

          // {
//             IntVect loU = disjointBox.smallEnd();
//             IntVect hiU = disjointBox.bigEnd();
//             dumpFAB2DSlicePretty(&Ux, 1, loU, hiU, 6, pout());
// // dumpFAB2DSlicePretty(&(*vPhi[iLev])[dit], 0, loU, hiU, 6, pout());
//           }

          if (maxLocComp >= 0 &&
              maxLocComp < CRDparam::g_CRDPhysics->numConservative())
            {
              for (BoxIterator bit(disjointBox); bit.ok(); ++bit)
                {
                  const Real pntErr = abs((*vPhi[iLev])[dit](bit(), 0));
                  if (pntErr > maxErr)
                    {
                      maxErr = pntErr;
                      maxLoc = bit();
                    }
                }
            }
        }
      if (maxLocComp >= 0 &&
          maxLocComp < CRDparam::g_CRDPhysics->numConservative())
        {
          CRD::msg.setMaxPrecFloatSN();
          CRD::msg << "Max error on level " << iLev << "\n" << maxErr
                   << CRD::var;
          CRD::msg << "Max error location " << iLev << "\n" << maxLoc
                   << CRD::var;
          CRD::msg.setFloatDefault();
        }
    }

  // Compute and print the error norms
  CRD::msg.setMaxPrecFloatSN();
  CRD::msg << "Error norms at time " << m_time << CRD::h2;
  Real errorNorm[3];
  for (int comp = 0;
       comp != CRDparam::g_CRDPhysics->numConservative();
       ++comp)
    {
      const Interval compIntv(comp, comp);
      errorNorm[0] = ::computeNorm(vPhi,
                                   vRefRatios,
                                   m_dx[0],  //** FIXME dx
                                   compIntv,
                                   0);
      errorNorm[1] = ::computeNorm(vPhi,
                                   vRefRatios,
                                   m_dx[0],  //** FIXME dx
                                   compIntv,
                                   1);
      errorNorm[2] = ::computeNorm(vPhi,
                                   vRefRatios,
                                   m_dx[0],  //** FIXME dx
                                   compIntv,
                                   2);
      CRD::msg<< "Component " << comp << ": "
              << CRDparam::g_CRDPhysics->consvStateName(comp) << CRD::body;
      CRD::msg<< "L1  (error)\n" << errorNorm[1] << CRD::var;
      CRD::msg<< "L2  (error)\n" << errorNorm[2] << CRD::var;
      CRD::msg<< "max (error)\n" << errorNorm[0] << CRD::var;
    }
  CRD::msg.setFloatDefault();
}

/*--------------------------------------------------------------------*/
//  Initialize argument <U> and <JU> at time m_time using IBC class
/** \param[out] a_U     Initialized to average <U>
 *  \param[out] a_JU    Initialized using product rule
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNS::loadIBCData(LevelData<FArrayBox>& a_U,
                         LevelData<FArrayBox>& a_JU)
{
  CH_TIME("AMRLevelCNS::loadIBCData");

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
//-- Adding random HIT IC here for now. SpectralForcing could get added
//-- to CRDPhysics, same as TurbModeling, so that CNSIBC
//-- can call the random spectral IC generator and spectral force generator
//-- itself in CNSIBC::initialize() and CNSIBC::addSourceTerm()
  if (m_level == 0 &&
      CRDparam::g_turbForcingType == CRDparam::TurbForcingSpectral)
    {
      // get an unnormalized random velocity field
      // use same grid size and MPI task count with this rseed for testing
      const long int rseed = CHprocID();
      // // use this rseed for a new random IC every run
      // const long int rseed = std::time(NULL);
      m_levelOp.m_spectralForcing.calcRandomIC(a_U, rseed);

      // exchange ghost cells (is this the right way to do this here?)
      m_levelOp.restartExchange(a_U);
    }
#endif /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */

  // // Assumed constant but we only know the time now
  // const_cast<CNSIBC* const>(CRDparam::g_CNSIBC)->setTime(m_time);
  // Sets <U> on valid+1 except a physical boundaries so we can calculate <JU>
  // next
  CRDparam::g_CNSIBC->initialize(
    a_U, m_levelGridMetrics, m_levelOp.getUnitNormals(), m_time, m_level);

//--Compute <JU> from <U> and <J>.  One-sided derivates are used for <U> at
//--the physical boundaries

  for (DataIterator dit = m_levelGridMetrics.getDataIterator(); dit.ok(); ++dit)
    {
      const Box interiorBox = m_levelGridMetrics.getBoxes()[dit];
      FArrayBox& JUFab = a_JU[dit];
      const FArrayBox& UFab = a_U[dit];
      const FArrayBox& JFab = m_levelGridMetrics.m_J[dit];
      fourthOrderCellProd(JUFab, UFab, JFab, interiorBox,
                          m_multiBlockCoordSys->problemDomain(interiorBox),
                          true);

//**DBG-begin
// FFTW testing -- leaving here for now.  Will move soon to test.
#if 0
      // Component 0
      JUFab.setVal(1., 0);
      // Component 1
      constexpr double L = 2*Pi;
      const int n = interiorBox.size(0);
      const double dL = L/n;
      for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
            {
              const double x = iv[0]*dL;
              const double y = iv[1]*dL;
              double u = 1.0;
              double val = u;
              for (int k = 1; k < (n/2+1); ++k)
                {
                  // c decreases from ~1 to 1/(n/2)
                  const double c = u*((double)(n - k - n/2 + 1))/(n/2);
                  val += c*(std::cos(k*x) +
                            std::sin(k*y) +
                            std::cos(k*x)*std::sin(k*y));
                }
              JUFab(iv, 1) = val;
            }
        }
      // Component 2
      for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
            {
              const double y = iv[1]*dL;
              const double z = iv[2]*dL;
              double u = 1.0;
              double val = u;
              for (int k = 1; k < (n/2+1); ++k)
                {
                  // c decreases from ~1 to 1/(n/2)
                  const double c = u*((double)(n - k - n/2 + 1))/(n/2);
                  val += c*(std::cos(k*y) +
                            std::sin(k*z) +
                            std::cos(k*y)*std::sin(k*z));
                }
              JUFab(iv, 2) = val;
            }
        }
      // Component 3
      for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
            {
              const double z = iv[2]*dL;
              const double x = iv[0]*dL;
              double u = 1.0;
              double val = u;
              for (int k = 1; k < (n/2+1); ++k)
                {
                  // c decreases from ~1 to 1/(n/2)
                  const double c = u*((double)(n - k - n/2 + 1))/(n/2);
                  val += c*(std::cos(k*z) +
                            std::sin(k*x) +
                            std::cos(k*z)*std::sin(k*x));
                }
              JUFab(iv, 3) = val;
            }
        }
      // Component 4
      for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
            {
              const double x = iv[0]*dL;
              const double y = iv[1]*dL;
              const double z = iv[2]*dL;
              double u = 1.0;
              double val = u;
              for (int k = 1; k < (n/2+1); ++k)
                {
                  // c decreases from ~1 to 1/(n/2)
                  const double c = u*((double)(n - k - n/2 + 1))/(n/2);
                  val += c*(std::cos(k*x) +
                            std::sin(k*y) +
                            std::cos(k*z) +
                            std::cos(k*x)*std::sin(k*y));
                }
              JUFab(iv, 4) = val;
            }
        }
#endif
//**DBG-end
    }
}
