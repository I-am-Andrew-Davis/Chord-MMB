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
 * \file CNSIBC.H
 *
 * \brief Abstract base class for initial and boundary conditions
 *
 *//*+*************************************************************************/

#ifndef _CNSIBC_H_
#define _CNSIBC_H_

//----- Standard Library -----//

#include <unordered_map>

//----- Chombo Library -----//

#include "REAL.H"
#include "PhysMappedIBC.H"
#include "Vector.H"
#include "CH_Hash.H"
#include "ParmParse.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class FArrayBox;
class FluxBox;
class ProblemDomain;
template <typename T> class LevelData;
class LevelGridMetrics;
class NewFourthOrderCoordSys;
#include "NamespaceFooter.H"

//--Forward declarations

typedef NewFourthOrderCoordSys BlockCoordSys;
class ViscousTensor4thOrderOp;
class TagLevelFactory;
class CRDPhysics;

// Forward declaration of BoundaryIndex for defining isHashable which is then
// runtime-checked in the definition of struct BoundaryIndex.
struct BoundaryIndex;
#include "BaseNamespaceHeader.H"
namespace CH_Hash
{

/// Remark that BoundaryIndex can be byte-hashed and give bytes to hash
template <>
struct isHashable<BoundaryIndex> : public std::true_type
{
  static constexpr int c_hashSize = 2*sizeof(int) +
    isHashable<Side::LoHiSide>::c_hashSize;
  // It is also runtime-asserted that BoundaryIndex is not padded in
  // BoundaryIndex::assertValid()
};

}
#include "BaseNamespaceFooter.H"


/*******************************************************************************
 */
/// Index for a boundary in a multi-block scheme
/**
 *//*+*************************************************************************/

struct BoundaryIndex
{
  /// Default constructor
  BoundaryIndex()
    { }

  /// Constructor
  BoundaryIndex(int a_block, int a_dir, Side::LoHiSide a_side)
    :
    m_block(a_block),
    m_dir(a_dir),
    m_side(a_side)
    {
      assertValid();
    }

  /// Weak construction
  void define(int a_block, int a_dir, Side::LoHiSide a_side)
    {
      m_block = a_block;
      m_dir = a_dir;
      m_side = a_side;
      assertValid();
    }

  /// Equilvalence
  bool operator==(const BoundaryIndex &other) const
    {
      return (this->m_block == other.m_block &&
              this->m_dir == other.m_dir &&
              this->m_side == other.m_side);
    }

  /// Checks
  void assertValid() const
    {
      CH_assert(m_block >= 0 && m_block <= CRDparam::g_coordSys->numBlocks());
      CH_assert(m_dir >= 0 && m_dir < SpaceDim);
      CH_assert(m_side != Side::Invalid);
      // Ensure no padding in BoundaryIndex
      CH_assert(reinterpret_cast<const char*>(this) +
                CH_Hash::isHashable<BoundaryIndex>::c_hashSize ==
                reinterpret_cast<const char*>(&m_side) +
                CH_Hash::isHashable<Side::LoHiSide>::c_hashSize);
    }

  // Member data initially invalid 
  int m_block = -1;
  int m_dir = -1;
  Side::LoHiSide m_side = Side::Invalid;
};

/// Output of a boundary index
std::ostream& operator<<(std::ostream&       a_os,
                         const BoundaryIndex a_idxB);


/*******************************************************************************
 */
/// Information about a particular BC
/**
 *//*+*************************************************************************/

struct BCInfo
{
  // Defaults are invalid
  CRDparam::DomainBCType m_type = CRDparam::DomainBCTypeUndefined;
                                      ///< BC type
  int m_order = -1;                   ///< Order of accuracy at this BC
  int m_idxState = -1;                ///< New state index
  Real m_relaxCBCStateParam = 0.;     ///< Relaxation parameter for CBC state
  Real m_relaxCBCWaveParam = 0.;      ///< Relaxation parameter for CBC waves
};


/*******************************************************************************
 */
/// BC and initial conditions
/** This class serves as abstract base for setting initial and boundary
 *  conditions. All IBC derive from this base class.
 *
 *  \note
 *  <ul>
 *    <li> Constructors should not require too much information as this class
 *         may itself read input to define parameters.  Reference states are
 *         required for setting the Fortran common block.
 *    <li> Expect domainBC to eventually be superseded by block domain BC.
 *  </ul>
 *
 *//*+*************************************************************************/

class CNSIBC
{
public:


/*==============================================================================
 * Types
 *============================================================================*/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor sets reference values and BC at domain extents
  CNSIBC();

  // Copy not permitted
  CNSIBC(const CNSIBC&) = delete;

  // Assignment not permitted
  CNSIBC &operator=(const CNSIBC&) = delete;

  /// Destructor
  virtual ~CNSIBC();


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Allocations of new physics states.  Customize for each derivative class.
  static std::vector<CRDPhysics*> allocatePhysics();

  /// Define common reference state
  void defineCommonRef(const Real a_rhoRef, const Real a_uRef);

  /// Set domainBC on all boundaries
  void setAllDomainBC(const BCInfo& a_domainBC);

  /// Set domain BC and order on one boundary
  void setDomainBC(const BoundaryIndex& a_bcIndex,
                   const BCInfo& a_domainBC);

  // /// Set domain BC state values, referenced by
  // void setDomainBCstate(const int a_stateIndex,
  //                       const BCstate& a_domainBC);
  
  /// Get BC info from a boundary
  const BCInfo& getDomainBC(const BoundaryIndex& a_bcIndex) const;

  // /// Get BC values from a boundary
  // const BCstate& getDomainBCstate(const BoundaryIndex& a_bcIndex) const;

  // /// Get BC values from a boundary
  // const BCstate& getDomainBCstate(const BCInfo& a_bc) const;
  
  /// Add components in an interval to a vector
  void addIntervaltoVec(Vector<int>&    a_compVec,
                        const Interval& a_interval) const;

  /// Return a name describing the IBC
  virtual const char *const IBCName() const = 0;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Write out each specified boundary condition
  virtual void writeBoundaryConditions() const;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);

//--Level operators

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;

  /// Initialize wall-model state after restart
  virtual void restartWallModel(LevelData<FArrayBox>&      a_JU,
                                LevelGridMetrics&          a_gridMetrics,
                                const LayoutData<FluxBox>& a_unitNormals) const;

  /// Initialize inlet region data structures
  virtual void initializeInletDataStructures(
    const LevelData<FArrayBox>& a_U,
    const int                   a_level,
    const DisjointBoxLayout&    a_disjointBoxLayout,
    const LevelGridMetrics&     a_gridMetrics,
    const bool                  a_hasFinerGrid) const;

  /// Copy the interior plane to the inlet plane and rescale it
  virtual void copyInteriorToInletAndRescale(
    const LevelData<FArrayBox>& a_U,
    const int                   a_level,
    const Real                  a_t,
    const Real                  a_dt,
    const int                   a_stage,
    const DisjointBoxLayout&    a_disjointBoxLayout,
    const LevelGridMetrics&     a_gridMetrics) const;

//--Patch operators

  /// Solve for the cell averaged gradients in the exterior ghost cells
  virtual void cellAvgGradBC(const FArrayBox&     a_PhicellAvgFab,
                             FArrayBox&           a_GradPhicellAvgFab,
                             const Box&           a_box,
                             const BlockDomain&   a_domain,
                             const RealVect&      a_dxVect,
                             const bool           a_solveTemp,
                             const int            a_level) const;
 
  /// Set boundary slopes
  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time);

  /// Adjust boundary fluxes to account for artificial viscosity
  virtual void artViscBC(
    FArrayBox&                                   a_NtFdir,
    const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
    const FArrayBox&                             a_U,
    const FArrayBox&                             a_unitNormalBasis,
    const FArrayBox&                             a_divVel,
    const FArrayBox&                             a_csq,
    const FArrayBox&                             a_dxFace,
    const Interval&                              a_momIntv,
    const Real                                   a_alpha,
    const Real                                   a_beta,
    const Box&                                   a_loFaceBox,
    const int                                    a_hasLo,
    const Box&                                   a_hiFaceBox,
    const int                                    a_hasHi,
    const int                                    a_dir,
    const Box&                                   a_disjointBox,
    LevelGridMetrics&                            a_gridMetrics,
    const Real                                   a_time,
    const int                                    a_level) const;

  /// This function is called to get boundary faces of a_box.
  virtual void getBoundaryFaces(Box&                  a_boundaryBox,
                                const Box&            a_box,
                                const BlockDomain&    a_domain,
                                const int&            a_dir,
                                const Side::LoHiSide& a_side) const;

  /// This function is called to get interior faces of a_box.
  virtual void getInteriorBoundaryFaces(Box&                  a_boundaryBox,
                                        const Box&            a_box,
                                        const BlockDomain&    a_domain,
                                        const int&            a_dir,
                                        const Side::LoHiSide& a_side) const;

  /// This function retrieves boxes that are at 2 boundary intersections
  virtual void getBoundaryEdgeBox(Box&                  a_boundaryBox,
                                  const Box&            a_box,
                                  const BlockDomain&    a_domain,
                                  const int&            a_dir1,
                                  const Side::LoHiSide& a_side1,
                                  const int&            a_dir2,
                                  const Side::LoHiSide& a_side2) const;

  /// Return boundary adjacent cells for a specific boundary type
  virtual void getBoundaryAdjCells(
    Vector<Box>&            a_boundAdjCells,
    const Box&              a_box,
    const BlockDomain&      a_domain,
    const int&              a_dir,
    const Side::LoHiSide&   a_side,
    const int&              a_level,
    const Real&             a_time,
    const int               a_numCells,
    LevelGridMetrics&       a_gridMetrics,
    CRDparam::DomainBCType& a_boundaryType) const;

  /// This function is called to get all boundary faces in domain
  void getAllWallBoundaryFaces(Vector<Box>&      a_boundaryBox,
                               Vector<int>&      a_faceDir,
                               const Real        a_time,
                               const Real        a_level,
                               LevelGridMetrics& a_gridMetrics) const;

  /// Set the primitive state at a domain boundary
  virtual void setBCprimState(FArrayBox&            a_WfaceAvgBndryDirFab,
                              FArrayBox&            a_WfaceAvgDirFab,
                              FArrayBox&            a_WcellAvgFab,
                              const FArrayBox&      a_bndryFaceAvgSlipVelDirFab,
                              FArrayBox&            a_bndryCellFab,
                              const Box&            a_boundaryFaceBox,
                              const Box&            a_boundaryFaceGhostBox,
                              const Box&            a_disjointBox,
                              const FArrayBox&      a_unitNormalBasisFab,
                              BCInfo&               a_bcInfo,
                              Box&                  a_bcBox,
                              const int             a_dir,
                              const Side::LoHiSide& a_side,
                              LevelGridMetrics&     a_gridMetrics,
                              const Real            a_time,
                              const Real            a_prevDt,
                              const int             a_level) const;

  /// Set the primitive state in domain ghost cells
  virtual void extrapolateDomainGhostCells(
    const FArrayBox&      a_WfaceAvgBndryDirFab,
    const FArrayBox&      a_WfaceAvgDirFab,
    FArrayBox&            a_WcellAvgFab,
    const Box&            a_boundaryFaceBox,
    const Box&            a_boundaryFaceGhostBox,
    const Box&            a_disjointBox,
    const FArrayBox&      a_unitNormalBasisFab,
    const int             a_dir,
    const Side::LoHiSide& a_side,
    const BCInfo&         a_bcInfo,
    const Box&            a_bcBox,
    LevelGridMetrics&     a_gridMetrics,
    const Real            a_time,
    const int             a_level) const;
  
  /// Set the primitive state at a domain boundary -- DEPRECATED
  virtual void setBCprimState(
    FArrayBox&            a_WavgFace,
    const FArrayBox&      a_WavgCell,
    const FArrayBox&      a_boundarySlipVelocity,
    const Box&            a_boundaryFaceBox,
    const Box&            a_boundaryFaceGhostBox,
    const Box&            a_disjointBox,
    const FArrayBox&      a_unitNormalBasisFab,
    Vector<BCInfo>&       a_bcTypes,
    Vector<Box>&          a_bcBoxes,
    const int             a_dir,
    const Side::LoHiSide& a_side,
    LevelGridMetrics&     a_gridMetrics,
    const Real            a_time,
    const int             a_level) const;

  /// Set the primitive state at a domain boundary -- DEPRECATED
  virtual void applyBCStates(
    FArrayBox&            a_WfaceBdryFab,
    FArrayBox&            a_WavgFace,
    FArrayBox&            a_WavgCell,
    const FArrayBox&      a_boundarySlipVelocity,
    const Box&            a_boundaryFaceBox,
    const Box&            a_boundaryFaceGhostBox,
    const Box&            a_disjointBox,
    const FArrayBox&      a_unitNormalBasisFab,
    const int             a_dir,
    const Side::LoHiSide& a_side,
    const Vector<BCInfo>& a_bcTypes,
    const Vector<Box>&    a_bcBoxes,
    LevelGridMetrics&     a_gridMetrics,
    const Real            a_time,
    const int             a_level) const;

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells (defaults to a no-op)
  virtual int exactSol(FArrayBox&              a_Ux,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       const LevelGridMetrics& a_gridMetrics,
                       const FluxBox&          a_unitNormals,
                       const DataIndex&        a_didx,
                       const Real              a_time,
                       const int               a_level) const;

  /// Add a source term to the flow
  virtual void addSourceTerm(FArrayBox&           a_sourceFab,
                             FArrayBox&           a_invDtFab,
                             const FArrayBox&     a_Wcell,
                             const FArrayBox&     a_UcellAvg,
                             const FluxBox&       a_WfaceAvgFxb,
                             const BlockDomain&   a_domain,
                             LevelGridMetrics&    a_gridMetrics,
                             const Real           a_time,
                             const Real           a_stageWeight,
                             const int            a_level,
                             const Box&           a_disjointBox,
                             const Box&           a_solveBox,
                             const DataIndex&     a_dataIndx,
                             const Real           a_globalKE,
                             const Real           a_globalHelicity) const;

  /// Assign mass fractions for various problems
  int assignMassFractions(std::vector<Real>& a_mfVals,
                          const std::string& a_cnInput,
                          const std::string& a_mfInput);

  /// Get the computational coordinates on cells
  virtual void getCellCompCoordinates(
    const Box&           a_box,
    FArrayBox&           a_XiFab,
    const BlockCoordSys& a_blockCoordSys) const;
  
  /// Get the computational and physical space coordinates on cells
  virtual void getCellCoordinates(const Box&           a_box,
                                  FArrayBox&           a_XiFab,
                                  FArrayBox&           a_XFab,
                                  const BlockCoordSys& a_blockCoordSys) const;

  /// Get the computational coordinates on faces
  virtual void getFaceCompCoordinates(
    const Box&           a_box,
    FArrayBox&           a_XiFab,
    const int            a_dir,
    const BlockCoordSys& a_blockCoordSys) const;

  /// Get the computational and physical space coordinates on faces
  virtual void getFaceCoordinates(const Box&           a_box,
                                  FArrayBox&           a_XiFab,
                                  FArrayBox&           a_XFab,
                                  const int            a_dir,
                                  const BlockCoordSys& a_blockCoordSys) const;

  /// Get the computational and physical space coordinates at nodes
  virtual void getNodeCoordinates(const Box&           a_box,
                                  FArrayBox&           a_XiFab,
                                  FArrayBox&           a_XFab,
                                  const BlockCoordSys& a_blockCoordSys) const;

  /// Set wall velocity on entire box
  virtual void setBoxWallVelocity(FArrayBox&           a_wallVelocity,
                                  const FArrayBox&     a_unitNormalBasisFab,
                                  const Box&           a_boundaryFaceBox,
                                  const Box&           a_disjointBox,
                                  LevelGridMetrics&    a_gridMetrics,
                                  const ProblemDomain& a_domain,
                                  const Real           a_time,
                                  const int            a_level,
                                  const int            a_dir,
                                  const Side::LoHiSide a_side) const;

  /// Set wall velocity on specific boundary -- user can specialize this
  virtual void setBndryWallVelocity(FArrayBox&           a_wallVelocity,
                                    const FArrayBox&     a_unitNormalBasisFab,
                                    const Box&           a_boundaryFaceBox,
                                    const Box&           a_disjointBox,
                                    LevelGridMetrics&    a_gridMetrics,
                                    const ProblemDomain& a_domain,
                                    const Real           a_time,
                                    const int            a_level,
                                    const int            a_dir,
                                    const Side::LoHiSide a_side,
                                    const BCInfo&        a_domT) const;

  /// Spatially average data if necessary
  virtual void spatiallyAverageData(
    LevelData<FArrayBox>& a_avgData,
    const int             a_stage,
    const Real            a_t) const;

  /// Output pnt-values and domain-sums specific to simulation
  virtual void inSituSumPntProcessing(
    const LevelData<FluxBox>&   a_faceAvgPlotData,
    const LevelData<FluxBox>&   a_WfaceAvgFxb,
    const LevelData<FArrayBox>& a_WcellAvgFab,
    const LevelData<FArrayBox>& a_WcellPntFab,
    const LayoutData<FluxBox>&  a_unitNormalsFxb,
    const LevelGridMetrics&     a_gridMetrics,
    const Real                  a_time,
    const int                   a_stage,
    const int                   a_level) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const;

  /// Set the primitive state at wall BC
  virtual void setWallBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_boundarySlipVelocity,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const;

  /// Set inlet profiles for velocity, temperature, species
  virtual void setCNSCBCProfiles(
    FArrayBox&                    a_BCProfile,
    const Box&                    a_boundaryBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const BoundaryIndex&          a_bcIdx,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const BCInfo&                 a_domT) const;

  /// Set exterior face-state for relaxed characteristic boundary conditions
  virtual void setRelaxedCBCPrimState(
    FArrayBox&              a_WfaceAvgExterior,
    FArrayBox&              a_bndryCellFab,
    FArrayBox&              a_WcellAvgFab,
    const FArrayBox&        a_WfaceAvgDirFab,
    const FArrayBox&        a_unitNormalBasisFab,
    const BoundaryIndex&    a_bcIdx,
    const BCInfo&           a_bcInfo,
    const Box&              a_disjointBox,
    const Box&              a_totalFaceBox,
    const LevelGridMetrics& a_gridMetrics,
    const int               a_dir,
    const Side::LoHiSide    a_side,
    const Real              a_time,
    const Real              a_prevDt,
    const int               a_level) const;

  /// Compute characteristic wave amplitudes
  virtual void computeCharWaveAmp(FArrayBox&       a_waveAmp,
                                  const FArrayBox& a_refState,
                                  const FArrayBox& a_derivState,
                                  const int        a_dir,
                                  const Box&       a_cellBox) const;

  /// Determine which components to use a Neumann condition for
  /// and which variable to solve for (rho, p, or T)
  virtual int setExtraBCs(Vector<int>&                  a_NeumannComp,
                          Vector<int>&                  a_interiorExtrap,
                          Vector<int>&                  a_faceValExtrap,
                          FArrayBox&                    a_WgradFace,
                          const BoundaryIndex&          a_bcIdx,
                          const Box&                    a_disjointBox,
                          LevelGridMetrics&             a_gridMetrics,
                          const Real                    a_time,
                          const int                     a_level,
                          const BCInfo&                 a_domT) const;

  /// Characteristic BC's using Motheau et al. 2017
  void applyCNSCBC(FArrayBox&                    a_WavgFace,
                   const Box&                    a_boundaryFaceBox,
                   FArrayBox&                    a_WavgCell,
                   const FArrayBox&              a_unitNormalBasisFab,
                   const BoundaryIndex&          a_bcIdx,
                   const Box&                    a_disjointBox,
                   LevelGridMetrics&             a_gridMetrics,
                   const Real                    a_time,
                   const int                     a_level,
                   const bool                    a_edgeBool,
                   const BCInfo&                 a_domT) const;

  /// Set boxes for inflow or outflow conditions in mixed boundaries
  virtual int setMixedBC(const Box&                    a_boundaryFaceBox,
                         const BoundaryIndex&          a_bcIdx,
                         const Box&                    a_disjointBox,
                         LevelGridMetrics&             a_gridMetrics,
                         const Real                    a_time,
                         const int                     a_level,
                         Vector<Box>&                  a_boxVect,
                         Vector<BCInfo>&               a_domainBC) const;

  /// Compute the cell-averaged gradients in the ghost cells
  void computeCellAvgGradBdry(FArrayBox&           a_GradPhicellAvgFab,
                              const FArrayBox&     a_PhicellAvgFab,
                              const Box&           a_combinedBox,
                              const int            a_dir,
                              const Side::LoHiSide a_side,
                              const int            a_numComp,
                              const RealVect&      a_dxVect) const;

  /// Compute the cell-averaged gradients in the ghost cells adjacent
  /// to the corners
  void computeCornerGrad(FArrayBox&           a_GradPhiCellAvgFab,
                         const FArrayBox&     a_PhicellAvgFab,
                         const Box&           a_cornerBox,
                         const Box&           a_biasedBox,
                         const int            a_gradDir,
                         const Side::LoHiSide a_side,
                         const int            a_numComp,
                         const Real           a_dx) const;

  // Solve for the wall state
  void computeWallPrimState(FArrayBox&       a_Wface,
                            const Box        a_boundaryFaceBox,
                            const FArrayBox& a_Wcell,
                            const FArrayBox& a_boundarySlipVelocity,
                            const FArrayBox& a_gammaCellFab,
                            const FArrayBox& a_unitNormalBasisFab,
                            const RealVect&  a_wallVel,
                            const int        a_viscousSlip,
                            const int        a_dir,
                            const int        a_lohiSign) const;

  /// Extrapolate ghost cell values based on imposed face values
  void extrapGhostFromFaceVals(const Box&            a_ghostBox1,
                               FArrayBox&            a_Wcell,
                               const FArrayBox&      a_Wface,
                               const Side::LoHiSide& a_side,
                               const int             a_dir,
                               const Vector<int>&    a_compList) const;

  /// Extrapolate ghost cell values and reset face values based on
  /// Neumann conditions
  void extrapGhostNeumann(const Box&            a_ghostBox1,
                          FArrayBox&            a_Wcell,
                          FArrayBox&            a_Wface,
                          const FArrayBox&      a_WgradFace,
                          const Side::LoHiSide& a_side,
                          const int             a_dir,
                          const Real            a_dx,
                          const Vector<int>&    a_compList) const;

  /// Extrapolate ghost cell values using interior values
  void extrapGhostFromInterior(const Box&            a_ghostBox1,
                               FArrayBox&            a_Wcell,
                               const Side::LoHiSide& a_side,
                               const int             a_dir,
                               const Vector<int>&    a_compList) const;


/*==============================================================================
 * Deprecated member functions
 * These override some base class functions and force an error.  Do not use
 * them.
 *============================================================================*/

public:

  /// Initialize a level (DEPRECATED)
  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// Set domainBC on all boundaries to a single type
  void setAllDomainBC(const int a_domainBC);

  /// Set domainBC on all boundaries
  void setAllDomainBC(const CRDparam::DomainBCType* a_domainBC);

  /// Set domainBC on one boundary
  void setDomainBC(const int  a_dir,
                   const int  a_side,
                   const int  a_domainBC);
  
  /// Set BC order on all boundaries
  void setAllDomainBCOrder(const int a_order);

  /// Set BC order on one boundary
  void setDomainBCOrder(const int  a_dir,
                        const int  a_side,
                        const int  a_order);
  
  /// (DEPRECATED) use updated interface instead
  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// (DEPRECATED) use updated interface instead
  /// Set the primitive state at wall BC
  virtual void setWallBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_boundarySlipVelocity,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// (DEPRECATED) use updated interface instead
  /// Set boxes for inflow or outflow conditions in mixed boundaries
  virtual int setMixedBC(const Box&                      a_boundaryFaceBox,
                         const int                       a_dir,
                         const Side::LoHiSide&           a_side,
                         const Box&                      a_disjointBox,
                         LevelGridMetrics&               a_gridMetrics,
                         const Real                      a_time,
                         const int                       a_level,
                         Vector<Box>&                    a_boxVect,
                         Vector<int>&                    a_orderBC,
                         Vector<CRDparam::DomainBCType>& a_domainType) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::unordered_map<BoundaryIndex, BCInfo,
                     CH_Hash::google_CityHash64<BoundaryIndex> > m_blockBCInfo;
                                      ///< Domain boundary conditions
  // std::unordered_map<int, BCstate,
  //                    CH_Hash::google_CityHash64<int> >           m_blockBCstate;
                                      ///< Domain boundary values
                                      ///< these are separated from the BC
                                      ///< conditions because it is possible
                                      ///< for several block boundaries to
                                      ///< share boundary states
  bool m_readInput;                   ///< Input file should only be read once.
                                      ///< Set to false in constructor and
                                      ///< test/set to true when reading input
                                      ///< file.
  int m_idxStateRef;                  ///< A reference or default state set to
                                      ///< the state named "reference".  If
                                      ///< "reference" does not exist, it is
                                      ///< created and set to default normalized
                                      ///< conditions.

  // !!! FIXME !!! These should part of BCvalues only, requires reworking
  // CombustionReference
  // Parameters for CNSBC smoothing  
  Real m_etaMax = 0.287;   ///< Smoothing coefficient for CNSCBC inflows
  Real m_etaCN  = 0.287;   ///< Smoothing coefficient for CNSCBC
  ///< inflows for mass fractions
  Real m_CBCsigma = 1./4.; ///< Smoothing coefficient for CNSCBC outflows
  Real m_CBCbeta  = 0.5;   ///< Smoothing coefficient for CNSCBC
  ///< outflows transverse conditions

public:

  static bool s_firstRKStage;         ///< If RK first stage
  static bool s_lastRKStage;          ///< If RK last stage
};

/*******************************************************************************
 *
 * Class CNSIBC: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get BC info from a boundary
/** \param[in] a_bcIndex  Index for a physical block boundary
 *  \return               Domain BC info on this block boundary
 *//*-----------------------------------------------------------------*/

inline const BCInfo&
CNSIBC::getDomainBC(const BoundaryIndex& a_bcIndex) const
{
  a_bcIndex.assertValid();
  auto boundary = CRDparam::g_coordSys->boundary(a_bcIndex.m_block,
                                                 a_bcIndex.m_dir,
                                                 a_bcIndex.m_side);
  if (!boundary.isDomainBoundary())
    {
      CRD::msg << "getDomainBC not valid for " << a_bcIndex
               << " because it is a multi-block interface." << CRD::error;
    }
  return m_blockBCInfo.at(a_bcIndex);
}

/*--------------------------------------------------------------------*/
//  Get BC values from a boundary
/** \param[in] a_bcIndex  Index for a physical block boundary
 *  \return               Domain BC info on this block boundary
 *//*-----------------------------------------------------------------*/

// inline const BCstate&
// CNSIBC::getDomainBCstate(const BoundaryIndex& a_bcIndex) const
// {
//   auto bc = getDomainBC(a_bcIndex);
//   return getDomainBCstate(bc);
// }

/*--------------------------------------------------------------------*/
//  Get BC values from a boundary
/** \param[in] a_bcIndex  Index for a physical block boundary
 *  \return               Domain BC info on this block boundary
 *//*-----------------------------------------------------------------*/

// inline const BCstate&
// CNSIBC::getDomainBCstate(const BCInfo& a_bc) const
// {
//   // this check should be here, but requires fully committing to
//   // the generalzied IBC approach
//   CH_assert(a_bc.m_stateIdx >= 0);
//   return m_blockBCstate.at(a_bc.m_stateIdx);
// }

/*--------------------------------------------------------------------*/
//  Add the components of an interval to a vector
/** 
 *  \param[out] a_compVec
 *                      The vector of components to add to
 *  \param[in]  a_interval
 *                      The interval to add to the vector
 *//*-----------------------------------------------------------------*/

inline void
CNSIBC::addIntervaltoVec(Vector<int>&    a_compVec,
                         const Interval& a_interval) const
{
  const int numIntComp = a_interval.size();
  if (a_interval.end() == -1)
    {
      return;
    }
  if (a_compVec[0] == -1)
    {
      a_compVec[0] = a_interval.begin();
    }
  else
    {
      a_compVec.push_back(a_interval.begin());
    }
  for(int comp = 1; comp != numIntComp; ++comp)
    {
      int val = a_interval.begin() + comp;
      a_compVec.push_back(val);
    }
}
  
/*--------------------------------------------------------------------*/
//  Provided the names of the inputs for species names and mass fractions,
//  this function will assign those values to the member data and ensure they
//  sum to 1.
/** 
 *  \param[out] a_mfVals
 *                      Member data vector of the mass fractions
 *  \param[in]  a_cnInput
 *                      Input string for the mass fraction names
 *                      Example: If our input file read
 *                      ibc.left_specs = O2 N2
 *                      Then a_cnInput = "left_specs"
 *  \param[in]  a_mfInput
 *                      Input string for the mass fraction values
 *                      Example: If our input file read
 *                      ibc.left_mfs = 0.233 0.767
 *                      Then a_mfInput = "left_mfs"
 *  \return     mCheck  Returns 0 if masses sum to 1 and 1 if not
 *//*-----------------------------------------------------------------*/

inline int
CNSIBC::assignMassFractions(std::vector<Real>& a_mfVals,
                            const std::string& a_cnInput,
                            const std::string& a_mfInput)
{
  // Return massCheck = 1 if the mass fractions do not sum to 1.
  int mCheck = 0;
  ParmParse ppIBC("ibc");
  const int numSpecies = CRDparam::g_numSpecies;
  // Assign vector to 0
  a_mfVals.assign(numSpecies, 0.);
  CH_assert(a_mfVals.size() == numSpecies);
  const char* cnInp = a_cnInput.c_str();
  const char* mfInp = a_mfInput.c_str();
  const int numMassSpec = ppIBC.countval(cnInp);
  std::vector<std::string> specNames(numMassSpec);
  std::vector<Real> specMF(numMassSpec);
  ppIBC.getarr(cnInp, specNames, 0, numMassSpec);
  ppIBC.getarr(mfInp, specMF, 0, numMassSpec);
  Real massTest = 0.;
  for(int spec = 0; spec != numSpecies; ++spec)
    {
      bool specCheck = false;
      for(int speci = 0; speci != numMassSpec; ++speci)
        {
          if (specNames[speci] == CRDparam::g_speciesNames[spec])
            {
              a_mfVals[spec] = specMF[speci];
              specCheck = true;
              break;
            }
        }
      if (!specCheck)
        {
          CRD::msg << "CNSIBC::assignMassFractions:The species "
                   << specNames[spec] << " was not specified as a species"
                   << CRD::error;
        }
      massTest += a_mfVals[spec];
    }
  // If error is greater than 0.01, return an error
  if (std::abs(1. - massTest) > 0.01)
    {
      mCheck = 1;
    }
  // Otherwise, renormalize the values
  else if (std::abs(1. - massTest) > 1.E-15)
    {
      for(int speci = 0; speci != numSpecies; ++speci)
        {
          a_mfVals[speci] = a_mfVals[speci]/massTest;
        }
      CRD::msg << "Initial mass fractions " << a_cnInput << " had to be "
               << "normalized so sum = 1 to compensate difference of "
               << massTest << "!" << CRD::warn;
      mCheck = 0;
    }
  return mCheck;
}
#endif /* ! defined _CNSIBC_H_ */
