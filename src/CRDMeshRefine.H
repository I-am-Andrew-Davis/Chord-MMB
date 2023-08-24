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
 * \file CRDMeshRefine.H
 *
 * \brief Version of BRMeshRefine used for multiblock in Chord
 *
 *//*+*************************************************************************/

#ifndef _CRDMESHREFINE_H_
#define _CRDMESHREFINE_H_

//----- Standard Library -----//

//----- Chombo Library -----//

#include "BRMeshRefine.H"
#include "MultiBlockCoordSys.H"

#include "UsingNamespace.H"

//----- Internal -----//



/*******************************************************************************
 */
/// CRDMeshRefine
/**
 *  Adapts BRMeshRefine for multiblock.
 *  \tparam MBCSFactory The function object that provides a MultiBlockCoordSys*
 *          for a level (often a lambda)
 *
 *//*+*************************************************************************/

template <typename MBCSFactory>
class CRDMeshRefine : public BRMeshRefine
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  CRDMeshRefine(const ProblemDomain& a_baseDomain,
                const Vector<int>&   a_refRatios,
                const Real           a_fillRatio,
                const int            a_blockFactor,
                const int            a_bufferSize,
                const int            a_maxSize,
                const int            a_baseLevel,
                MBCSFactory&&        a_MBCSfactory);

  /// Destructor
  virtual ~CRDMeshRefine()
    { }


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Is the argument box properly nested?
  virtual bool properlyNested(const Box&           a_boxBF,
                              const ProblemDomain& a_domainBF,
                              const IntVectSet&    a_pndBF,
                              int                  a_totalBufferSizeBF,
                              const int            a_level) const override;

  /// Add tags from a finer level to this level
  virtual void buildSupport(const ProblemDomain& a_lvldomain,
                            Vector<Box>&         a_lvlboxes,
                            IntVectSet&          a_modifiedTags,
                            const int            a_level) override;

  /// Clip grown mesh from a finer level to the problem domain (null Op.)
  virtual void clipBox(Box&                 a_box,
                       const ProblemDomain& a_domain) const override;

  /// Hijack makePNDs to allow use of a higher-level mesh as the PNDvoid 
  void makePNDs(Vector<IntVectSet>&          a_pnds,
                Vector<int>&                 a_totalBufferSize,
                const int                    a_baseLevel,
                const int                    a_topLevel,
                const Vector<ProblemDomain>& a_domains,
                const Vector<Vector<Box>>&   a_oldMeshes,
                const Vector<int>&           a_bufferSize) const override;

  /// Get the multiblock coordinate system for the level
  const MultiBlockCoordSys* getLevelMBCoordSys(const int a_level) const;


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  MBCSFactory&& m_MBCSfactory;        ///< Function object acting as factory to
                                      ///< create pointers to multiblock
                                      ///< coordinate systems
  mutable int m_prevIdxBlk;           ///< Last used block index in
                                      ///< properlyNested.  Recheck for next
                                      ///< use.
  const int m_baseLevel;              ///< Desired level for building PND
};


/*******************************************************************************
 *
 * Class CRDMeshRefine: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_baseDomain
 *                      The problem domain for the base level (0)
 *  \param[in]  a_refRatios
 *                      Refinement ratios betweeen all levels
 *  \param[in]  a_fillRatio
 *                      Amount of a box that should be filled with
 *                      tagged cells
 *  \param[in]  a_blockFactor
 *                      All boxes can be coarsened by this amount.
 *                      Also, the minimum box size in any direction.
 *  \param[in]  a_bufferSize
 *                      The grid buffer size needed so that operations
 *                      on level l+1 are independent from those on
 *                      level l-1.
 *  \param[in]  a_maxSize
 *                      The maximum box size in any direction
 *  \param[in]  a_baseLevel
 *                      If non-zero, the PND is set from this level
 *                      The level 0:a_baseLevel should not be
 *                      regridded
 *  \param[in]  a_MBCSfactory
 *                      Factory object that creates a multiblock
 *                      coordinate system for any level
 *//*-----------------------------------------------------------------*/

template <typename MBCSFactory>
CRDMeshRefine<MBCSFactory>::CRDMeshRefine(
  const ProblemDomain& a_baseDomain,
  const Vector<int>&   a_refRatios,
  const Real           a_fillRatio,
  const int            a_blockFactor,
  const int            a_bufferSize,
  const int            a_maxSize,
  const int            a_baseLevel,
  MBCSFactory&&        a_MBCSfactory)
:
  BRMeshRefine(a_baseDomain,
               a_refRatios,
               a_fillRatio,
               a_blockFactor,
               a_bufferSize,
               a_maxSize),
  m_MBCSfactory(std::forward<MBCSFactory>(a_MBCSfactory)),
  m_prevIdxBlk(-1),
  m_baseLevel(a_baseLevel)
{ }


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Is the argument box properly nested?
/**
 *  With sufficient partitioning, the boxes should eventually fit
 *  inside the proper nesting domain.  The role of this function is to
 *  simply determine whether the candidate box, given by a_boxBF, is
 *  properly nested or not.
 *
 *  All arguments to this function have been coarsened by the blocking
 *  factor
 *  \param[in]  a_boxBF The box to test
 *  \param[in]  a_domainBF
 *                      The problem domain
 *  \param[in]  a_pndBF An IntVectset representing the proper nesting
 *                      domain, essentially a union of all boxes on
 *                      the base level.
 *  \param[in]  a_totalBufferSizeBF
 *                      The total buffer size from this level to the
 *                      base level
 *  \param[in]  a_level Grid level
 *  \return             T - the box is properly nested.
 *//*-----------------------------------------------------------------*/

template <typename MBCSFactory>
bool
CRDMeshRefine<MBCSFactory>::properlyNested(
  const Box&           a_boxBF,
  const ProblemDomain& a_domainBF,
  const IntVectSet&    a_pndBF,
  int                  a_totalBufferSizeBF,
  const int            a_level) const
{
  if (m_PNDMode == 0)
    {
      MayDay::Error("CRDMeshRefine only implemented for PNDMode != 0");
    }
  else
    {
      Box growBox(a_boxBF);
      growBox.grow(a_totalBufferSizeBF*m_refineDirs);
      // pout() << "PN: tb: " << a_totalBufferSizeBF << ": " << growBox
      //        << std::endl;

      // Perform all tests at standard resolution
      growBox.refine(m_level_blockfactors[a_level]);
      Box box = refine(a_boxBF, m_level_blockfactors[a_level]);
      // pout() << "PN on: " << a_boxBF << ' ' << box << std::endl;
      // pout() << "PN grow on: " << growBox << std::endl;

      const MultiBlockCoordSys* coordSys = getLevelMBCoordSys(a_level);
      // There's no option here but calling whichBlock since there are no rules
      // about a_box.  However, see if a previous value can be used.
      int idxBlk = -1;
      if (m_prevIdxBlk >= 0)
        {
          if (coordSys->mappingBlocks()[m_prevIdxBlk].contains(box))
            {
              idxBlk = m_prevIdxBlk;
            }
        }
      if (idxBlk == -1)
        {
          idxBlk = coordSys->whichBlock(box);
        }
      if (idxBlk == -1)
        {
          // pout() << "PN done" << std::endl;
          return false;  // Not contained within a single block
        }
      m_prevIdxBlk = idxBlk;
      Box block = coordSys->mappingBlocks()[idxBlk];
      // Crop within the block
      Box growBoxBF = growBox & block;
      // Coarsen by the block factor
      growBoxBF.coarsen(m_level_blockfactors[a_level]);
      if (!a_pndBF.contains(growBoxBF))
        {
          // Nevermind block boundaries, not even in the blocks PND
          // pout() << "PN done" << std::endl;
          return false;
        }
  
//**The pnd should really be an interpolation stencil into ghosts and then the
//**MBEx stencil used to fill those ghosts
//**1) Query each ghost cell for a valid cell and stencil
//**2) Is each stencil cell in the PND?
//**3) Are the stencil cells supported
//**4) But stencil cells are not part of the meta-data

//**A) Assume that to fill ghosts, only info from adjacent boxes is required.
//**B) Find necessary ghosts, find stencil cells required to fill those ghosts,
//**   and make sure those cells are properly nested in their blocks.

      if (!block.contains(growBox))
        {
          for (const int dir : EachDir)
            {
              for (const auto side : EachSide)
                {
                  if (Side::sign(side)*(block.sideEnd(side)[dir] -
                                        growBox.sideEnd(side)[dir]) < 0)
                    {
                      const BlockBoundary& boundary =
                        coordSys->boundary(idxBlk, dir, side);
                      // Corners are neglected but should have been properly
                      // considered in buildSupport.
                      if (boundary.isInterface())
                        {
                          Box trfmGrowBox =
                            boundary.getTransformation().transformFwd(growBox);
                          // Grow a bit more because of the stencil for MBEx
                          trfmGrowBox.grow(2);
                          // Crop within the neighbor block
                          trfmGrowBox &=
                            coordSys->mappingBlocks()[boundary.neighbor()];
                          // Coarsen by the block factor
                          Box trfmGrowBoxBF = trfmGrowBox;
                          trfmGrowBoxBF.coarsen(m_level_blockfactors[a_level]);
                          if (!a_pndBF.contains(trfmGrowBoxBF))
                            {
                              // pout() << "PN done" << std::endl;
                              return false;
                            }
                        }
                    }
                }
            }
        }
    }

  // pout() << "PN done" << std::endl;
  return true;
}

/*--------------------------------------------------------------------*/
//  Add tags from a finer level to this level
/**
 *  a_lvlboxes are the (new) boxes from a finer level that have been
 *  coarsened and grown by the buffer size.  These need to be added to
 *  a_modifiedTags.  The addition of the buffer can extend the boxes
 *  across multiblock and periodic domains.  We have to take that into
 *  consideration here when adding to a_modifiedTags.
 *
 *  All arguments to this function have been coarsened by the blocking
 *  factor
 *
 *  \param[in]  a_domainBF
 *                      The domain
 *  \param[in]  a_boxesBF
 *                      Boxes from a finer level that have been
 *                      coarsened and grown by the buffer size.
 *  \param[in]  a_modifiedTagsBF
 *                      User defined tags for the level
 *  \param[out] a_modifiedTagsBF
 *                      Add boxes from finer level
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

template <typename MBCSFactory>
void
CRDMeshRefine<MBCSFactory>::buildSupport(const ProblemDomain& a_domainBF, 
                                         Vector<Box>&         a_boxesBF, 
                                         IntVectSet&          a_modifiedTagsBF,
                                         const int            a_level)
{
  // pout() << "BS entry" << std::endl;
  struct BoxData
  {
    BoxData(const Box& a_box, const int a_idxBlk) noexcept
      :
      m_box(a_box),
      m_idxBlk(a_idxBlk)
      { }
    Box m_box;
    int m_idxBlk;
  };

  // UMap_Box_int = boxesBacking1;
  Vector<BoxData> boxesBacking1;
  boxesBacking1.reserve(a_boxesBF.size());
  Vector<BoxData> boxesBacking2;
  boxesBacking2.reserve(64);

  // Get the CS
  const MultiBlockCoordSys* coordSys = getLevelMBCoordSys(a_level);

  // Let's get the block lookup out of the way
  for (const Box& boxBF : a_boxesBF)
    {
      // Perform all tests at standard resolution
      Box box = boxBF;
      box.refine(m_level_blockfactors[a_level]);
      // std::cout << "BS box: " << boxBF << " --level " << a_level << "-- "  << box << std::endl;
      
      const int idxBlk = coordSys->whichBlockOverlap(box);
      // If this fails, the blocks need to be moved further apart on the lattice
      if (idxBlk == -1)
        {
          MayDay::Error("Box + buffer reaches into more than one block.  Move "
                        "blocks further apart on the lattice.  This is an "
                        "internal error.");
        }
      boxesBacking1.emplace_back(box, idxBlk);
    }

//**This works fine for Cartesian boxes, but is not MMB general

//**For each cell that reaches into the ghosts, you need to find it's stencil
//**and all these are tags.

  Vector<BoxData>* boxes = &boxesBacking1;
  Vector<BoxData>* nextBoxes = &boxesBacking2;
  // std::cout << "Level: " << a_level << "  " << m_level_blockfactors[a_level]
  //           << "  " << m_blockFactor << std::endl;

  int idxBB = 0;  // Number of block boundaries crossed
  while (boxes->size() > 0)
    {
      // std::cout << "In boxes" << std::endl;
      // for (BoxData& boxData : *boxes)
      //   {
      //     std::cout << "  " << boxData.m_box << std::endl;
      //   }
      for (BoxData& boxData : *boxes)
        {
          Box box = boxData.m_box;
          const int idxBlk = boxData.m_idxBlk;

          const Box& block = coordSys->mappingBlocks()[idxBlk];
          if (block.contains(box))
            {
              // Whole box is in one block
              a_modifiedTagsBF |= coarsen(box, m_level_blockfactors[a_level]);
              continue;
            }
          else
            {
              a_modifiedTagsBF |= coarsen(box & block,
                                          m_level_blockfactors[a_level]);
// std::cout << "On box " << box << std::endl;
              // Must deal with block extending across block boundaries.  It's
              // okay if the nextBox's produced from box partly overlap.
              for (const int dir : EachDir)
                {
                  for (const auto side : EachSide)
                    {
                      if (Side::sign(side)*(block.sideEnd(side)[dir] -
                                            box.sideEnd(side)[dir]) < 0)
                        {
                          Box nextBox = box;
                          if (idxBB == 0)
                            {
                              nextBox.grow(2);  // MBEx stencil
                            }
                          // Crop
                          switch (side)
                            {
                            case Side::Lo:
                              nextBox.setBig(dir, block.smallEnd(dir) - 1);
                              break;
                            case Side::Hi:
                              nextBox.setSmall(dir, block.bigEnd(dir) + 1);
                              break;
                            default:
                              CH_assert(false);
                            }
                          // std::cout << "  Dir " << dir << ", side " << side
                          //           << " crop: " << nextBox << std::endl;
                          CH_assert(!nextBox.isEmpty());
                          // Now transform nextBox across the boundary
                          const BlockBoundary& boundary =
                            coordSys->boundary(idxBlk, dir, side);
                          // Physical boundaries may not have been clipped
                          if (boundary.isInterface())
                            {
                          // // Physical boundaries should have already been
                          // // clipped
                          // CH_assert(boundary.isInterface());
                              Box trfmNextBox =
                                boundary.getTransformation().transformFwd(
                                  nextBox);
                              // std::cout << "  Trfm " << trfmNextBox << std::endl;
                              nextBoxes->emplace_back(trfmNextBox,
                                                      boundary.neighbor());
                            }
                        }
                    }
                }
            }
        }
      Vector<BoxData>* tmp = boxes;
      boxes = nextBoxes;
      nextBoxes = tmp;
      nextBoxes->clear();
      ++idxBB;
    }
  // if (dbgc >= 4) abort();
  // pout() << "BS done" << std::endl;
}

/*--------------------------------------------------------------------*/
//  Clip grown mesh from a finer level to the problem domain (null
//  Op.)
/**
 *  Boxes from the fine level have been coarsend to this level and
 *  grown by the buffer size.  We want the boxes to extend across
 *  periodic boundaries and multiblock interfaces.  This
 *  characteristic is used in splitBox above to extend tags across
 *  those interfaces.  However, the grown boxes must be clipped at
 *  physical boundaries.
 *
 *  However, physical boundaries are detected and ignored in
 *  buildSupport.  Hence, this routine does not need to do anything.
 *
 *  All arguments to this function have been coarsened by the blocking
 *  factor
 *
 *  \param[in]  a_box   Box to clip
 *  \param[in]  a_domain
 *                      Problem domain
 *//*-----------------------------------------------------------------*/

template <typename MBCSFactory>
inline void
CRDMeshRefine<MBCSFactory>::clipBox(Box&                 a_box,
                                    const ProblemDomain& a_domain) const
{
  // Physical boundaries are detected and ignored in buildSupport.
  // Hence, this routine does not need to do anything.
}

/*--------------------------------------------------------------------*/
//  Hijack makePNDs to allow use of a higher-level mesh as the PND
/** \param[out] a_pnds  IntVectSet of the pnd on the base level
 *                      (blocked)
 *  \param[out] a_totalBufferSize
 *                      The grid buffer to the base pnd (blocked)
 *  \param[in]  a_baseLevel
 *                      Base of the regrid (the mesh on this level
 *                      does not change)
 *  \param[in]  a_topLevel
 *                      Highest level with tags.  A mesh can be
 *                      created for a_topLevel+1
 *  \param[in]  a_domains
 *                      Problem domains
 *  \param[in]  a_oldMeshes
 *                      Meshes before regrid (source for the PND)
 *  \param[in]  a_bufferSize
 *                      Grid buffer for each level (blocked)
 *                      
 *  \note
 *  Current expectations:
 *  <ol>
 *    <li> Meshes at levels <= m_baseLevel are not dynamic.  They
 *         should have their regrid level set to -1.  I.e. for a
 *         static PND on level 1, 'amr.regrid_interval = -1 2'
 *    <li> m_baseLevel > 1 has not been tested.
 *  </ol>
 *  Other notes:
 *  <ul>
 *    <li> m_baseLevel > 1 has not been tested.
 *    <li> Dynamic appearance and disappearance of mesh levels has not
 *         been explored but should work.
 *  </ul>
*//*------------------------------------------------------------------*/

template <typename MBCSFactory>
void 
CRDMeshRefine<MBCSFactory>::makePNDs(
  Vector<IntVectSet>&          a_pnds,
  Vector<int>&                 a_totalBufferSize,
  const int                    a_baseLevel,
  const int                    a_topLevel,
  const Vector<ProblemDomain>& a_domains,
  const Vector<Vector<Box>>&   a_oldMeshes,
  const Vector<int>&           a_bufferSize) const
{
  bool forceStandardPND = false;
  if (m_baseLevel > 0 && a_topLevel >= m_baseLevel)
    {
      if (a_oldMeshes.size() <= m_baseLevel)
        {
          forceStandardPND = true;
        }
      else
        {
          if (a_oldMeshes[m_baseLevel].size() == 0)
            {
              forceStandardPND = true;
            }
        }
      if (forceStandardPND)
        {
          CRD::msg << "Reverting to standard mesh for base PND!" << CRD::warn;
        }
    }
  if (m_baseLevel == 0 ||
      a_topLevel < m_baseLevel ||
      a_baseLevel >= m_baseLevel ||
      forceStandardPND)
    {
      // Standard PND
      MeshRefine::makePNDs(a_pnds,
                           a_totalBufferSize,
                           a_baseLevel,
                           a_topLevel,
                           a_domains,
                           a_oldMeshes,
                           a_bufferSize);
    }
  else
    {
      // Build the PND from the mesh on level m_baseLevel.  Get the boxes from
      // level m_baseLevel and coarsen them as if they were on level a_baseLevel
      IntVectSet mesh;
      int refToBase = 1;
      for (int lvl = m_baseLevel - 1; lvl >= a_baseLevel; --lvl)
        {
          refToBase *= m_nRefVect[lvl];
        }
      const IntVect crFactor =
        this->inRefineDirs(m_level_blockfactors[a_baseLevel]*refToBase);
      const Vector<Box>& pndMesh = a_oldMeshes[m_baseLevel];
      CH_assert(pndMesh.size() > 0);
      for (int idxBox = 0, idxBox_end = pndMesh.size(); idxBox != idxBox_end;
           ++idxBox)
        {
          Box box = pndMesh[idxBox];
          CRD::msg << CRD::fv4 << "PND" << m_baseLevel << " box: " << box
                   << CRD::end;
          box.coarsen(crFactor);
          mesh |= box;
        }
      MeshRefine::makePNDs(a_pnds,
                           a_totalBufferSize,
                           a_baseLevel,
                           a_topLevel,
                           a_domains,
                           mesh,
                           a_bufferSize);
    }
}

/*--------------------------------------------------------------------*/
//  Get the multiblock coordinate system for the level
/** \param[in]  a_level Grid level
 *  \return     Pointer to multiblock coordinate system for the level
 *//*-----------------------------------------------------------------*/

template <typename MBCSFactory>
inline const MultiBlockCoordSys*
CRDMeshRefine<MBCSFactory>::getLevelMBCoordSys(const int a_level) const
{
  return m_MBCSfactory(a_level);
}

#endif  /* ! defined _CRDMESHREFINE_H_ */
