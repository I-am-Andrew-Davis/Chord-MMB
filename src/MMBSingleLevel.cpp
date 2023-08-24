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
 * \file MMBSingleLevel.cpp
 *
 * \brief Member functions for MMBSingleLevel
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

//----- Chombo Library -----//

#include "LevelGridMetrics.H"

//----- Internal -----//

#include "MMBSingleLevel.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Class MMBAddFaceOp: member definitions
 *
 * Used in the copier to transform data as it is read in from buffers
 *
 ******************************************************************************/

class MBAssignFaceOp: public LDOperator<FArrayBox>
{

public:

//--Constructors/destructors

  /// Construct
  MBAssignFaceOp(
    const MultiBlockCoordSys* const               a_coordSysPtr,
    const BoxLayout&                              a_lvlStoLayout,
    const LayoutData<MultiBlockRegions::FaceTag>& a_lvlSto2Origin,
    const Real                                    a_scale = 1.0,
    const bool                                    a_fluxSign = true)
    :
    m_coordSysPtr(a_coordSysPtr),
    m_lvlStoLayout(a_lvlStoLayout),
    m_lvlSto2Origin(a_lvlSto2Origin),
    m_scale(a_scale),
    m_fluxSign(a_fluxSign)
    { }

  // Use synthesized destructor, copy, and assignment

//--Member functions

  /// Use API 1 (which gives the dataIndex for the destination)
  virtual int API() const override
    { return 1; }

  /// API 0 not supported
  virtual void linearIn(FArrayBox& arg,
                        void* buf,
                        const Box& regionTo,
                        const Interval& comps) const override
    {
      MayDay::Error("API 0 not supported for MBAssignFaceOp");
    }

  /// Operation for distributed memory
  virtual void linearIn1(FArrayBox&       a_dstFab,
                         void*            a_srcBuf,
                         const Box&       a_dstRegion,
                         const Interval&  a_comps,
                         const DataIndex& a_dstDidx) const override
    {
      Real* buffer = static_cast<Real*>(a_srcBuf);

      MultiBlockRegions::FaceTag faceTag;
      IndicesTransformation trfm = getTransform(a_dstDidx, faceTag);
      Box srcRegion = trfm.transformFwd(a_dstRegion);

      FArrayBox tmp(srcRegion, a_comps.size(), buffer);
      trfmOp(a_dstFab, a_comps, a_dstRegion, tmp, a_comps, faceTag.dir(), trfm);
    }

  /// API 0 not supported
  void op(FArrayBox&       a_dest,
          const Box&       a_regionFrom,
          const Interval&  a_Cdest,
          const Box&       a_regionTo,
          const FArrayBox& a_src,
          const Interval&  a_Csrc) const override
    {
      MayDay::Error("API 0 not supported for MBAssignFaceOp");
    }

  /// Operation for shared memory
  virtual void op1(FArrayBox&       a_dstFab,
                   const Box&       a_srcRegion,
                   const Interval&  a_dstComps,
                   const Box&       a_dstRegion,
                   const FArrayBox& a_srcFab,
                   const Interval&  a_srcComps,
                   const DataIndex& a_dstDidx) const override
    {
      CH_assert(a_srcRegion.numPts() == a_dstRegion.numPts());
      CH_assert(a_srcComps.size()    == a_dstComps.size());

      MultiBlockRegions::FaceTag faceTag;
      IndicesTransformation trfm = getTransform(a_dstDidx, faceTag);
      CH_assert(trfm.transformFwd(a_dstRegion) == a_srcRegion);
      trfmOp(a_dstFab, a_dstComps, a_dstRegion, a_srcFab, a_srcComps,
             faceTag.dir(), trfm);
    }

private:

  /// The op for copying from a src to dst
  void trfmOp(FArrayBox&                   a_dstFab,
              const Interval&              a_dstComps,
              const Box&                   a_dstRegion,
              const FArrayBox&             a_srcFab,
              const Interval&              a_srcComps,
              const int                    a_dir,
              const IndicesTransformation& a_trfm) const
    {
      const IndicesTransformation invTr = a_trfm.inverse();
      const int fluxSign = m_fluxSign ? invTr.getSign()[a_dir] : 1;

      const Real scale = fluxSign*m_scale;

      for (int cDst = a_dstComps.begin(), cDst_end = a_dstComps.end() + 1,
             cSrc = a_srcComps.begin(); cDst != cDst_end; ++cDst, ++cSrc)
        {
          MD_BOXLOOP(a_dstRegion, i)
            {
              const IntVect ivDst = MD_GETIV(i);
              const IntVect ivSrc = a_trfm.transformFwd(ivDst);
              a_dstFab[MD_IX(i, cDst)] = scale*a_srcFab[MD_IV(ivSrc, cSrc)];
            }
        }
    }

  /// Get the transformation for this copy
  /** \param[in]  a_dstDidx
   *                      The dataIndex for data we are copying into
   *  \param[out] a_faceTag
   *                      The face we are working on
   */
  IndicesTransformation getTransform(
    const DataIndex&            a_dstDidx,
    MultiBlockRegions::FaceTag& a_faceTag) const
    {
      int dstIdxBlk = m_lvlStoLayout.blockIndex(a_dstDidx);
      // When lvlSto was created, we retained relations to the original layout,
      // which designates the source data index and face.  Both of these are in
      // the face tag.
      a_faceTag = m_lvlSto2Origin[a_dstDidx];
      const BlockBoundary& bb =
        m_coordSysPtr->boundary(dstIdxBlk, a_faceTag.dir(), a_faceTag.side());
      return bb.getTransformation();
    }

//--Data

  const MultiBlockCoordSys* m_coordSysPtr;
                                      ///< Multiblock coordinate system
  const BoxLayout& m_lvlStoLayout;    ///< Layout of the compact storage
  const LayoutData<MultiBlockRegions::FaceTag>& m_lvlSto2Origin;
                                      ///< The FaceTag for the origination
                                      ///< DataIndex and face
  const Real m_scale;                 ///< Scaling factor (default 1)
  const bool m_fluxSign;              ///< Modify scale by sign (default true)
                                      ///< Set to true for fluxes and false for
                                      ///< state
};


/*******************************************************************************
 *
 * Class MMBSingleLevel: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor
/** \param[in]  a_levelGridMetrics
 *                      Level grid metrics class used to retrieve CS
 *//*-----------------------------------------------------------------*/

MMBSingleLevel::MMBSingleLevel(const LevelGridMetrics& a_levelGridMetrics)
  :
  m_levelGridMetrics(a_levelGridMetrics)
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Weak construction
/** Define the storage used for exchanging fluxes for averaging on a
 *  single level across multiblock boundaries
 *  \param[in]  a_maxComp
 *                      A maximum number of components to average.
 *                      This amount is always exchanged.
 *//*-----------------------------------------------------------------*/

void
MMBSingleLevel::define(const int a_maxComp)
{
  m_lvlSto.define(m_levelGridMetrics.getMBRegions().lvlStoLayout(), a_maxComp);
}

/*--------------------------------------------------------------------*/
//  Set the flux on a connected block boundary
/** \param[in]  a_disjointBox
 *                      The disjoint box (or a subset) on which a_flux
 *                      is based
 *  \param[in]  a_flux  Flux to add to the register.  These are fluxes
 *                      on the faces of connected blocks that need to
 *                      be averaged.
 *  \param[in]  a_didx  Data index for the box on the level
 *  \param[in]  a_intv  The interval to average in a_flux
 *//*-----------------------------------------------------------------*/

void
MMBSingleLevel::setFlux(const Box&       a_disjointBox,
                        const FluxBox&   a_flux,
                        const DataIndex& a_didx,
                        const Interval&  a_intv)
{
  CH_assert(a_flux.box().contains(a_disjointBox));
  CH_assert(a_intv.size() <= m_lvlSto.nComp());

  const stc::Vector<Vector<Box>, 2*SpaceDim>& lvlLocations =
    getMBRegions().getLvlLocations(a_didx);
  for (const int dir : EachDir)
    {
      for (const auto side : EachSide)
        {
          const int idxFace = MultiBlockRegions::FaceTag::indexFace(dir, side);
          if (lvlLocations[idxFace].size() > 0)
            {
              MultiBlockRegions::FaceTag faceTag(dir, side, a_didx);
              const DataIndex& didxSto = getMBRegions().getLvlStoDidx(faceTag);
              FArrayBox& fluxSto = m_lvlSto[didxSto];
              // Create an alias so shifting is thread-safe
              FArrayBox fluxRcl(a_flux[dir].interval(), const_cast<FArrayBox&>(a_flux[dir]));
              // Shift so faces overlay interior cells
              fluxRcl.shiftHalf(dir, -Side::sign(side));
              for (Box box : lvlLocations[idxFace])
                {
                  CH_assert(fluxSto.contains(box));
                  // 'box' may be larger than fluxRcl if the latter was
                  // partitioned for threading
                  box &= a_disjointBox;
                  if (!box.isEmpty())
                    {
                      CH_assert(fluxRcl.contains(box));
                      for (int cR = a_intv.begin(), cR_end = a_intv.end() + 1,
                             cS = 0; cR != cR_end; ++cR, ++cS)
                        {
                          MD_BOXLOOP(box, i)
                            {
                              fluxSto[MD_IX(i, cS)] = fluxRcl[MD_IX(i, cR)];
                            }
                        }
                    }
                }
              // Shift back
              fluxRcl.shiftHalf(dir, Side::sign(side));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Exchange fluxes once all are known
/** \param[in]  a_fluxSign
 *                      Whether or not to consider the sign change of
 *                      fluxes across MMB boundaries.  Set to true
 *                      (default) for fluxes and false if exchanging
 *                      state information
 *//*-----------------------------------------------------------------*/

void
MMBSingleLevel::exchange(const bool a_fluxSign)
{
  const Interval intv(0, m_lvlSto.nComp() - 1);
  MBAssignFaceOp op(&getCoordSys(),
                    getMBRegions().lvlStoLayout(),
                    getMBRegions().lvlSto2Origin(),
                    1.0,
                    a_fluxSign);
  m_lvlSto.makeItSo(intv, m_lvlSto, m_lvlSto, intv,
                    getMBRegions().lvlRemoteCopier(),
                    op);
}

/*--------------------------------------------------------------------*/
//  Reflux the average.  This updates the divergence in cells.
/** \param[in]  a_RHSlvl
 *                      Usually the right-hand side
 *  \param[out] a_RHSlvl
 *                      Refluxed in cells adjacent to connected block
 *                      interfaces
 *  \param[in]  a_dx    Mesh spacing
 *  \param[in]  a_intv  The interval to average in a_RHSlvl
 *//*-----------------------------------------------------------------*/

void
MMBSingleLevel::refluxAverage(LevelData<FArrayBox>& a_RHSlvl,
                              const RealVect&       a_dx,
                              const Interval&       a_intv)
{
  CH_assert(a_intv.size() <= m_lvlSto.nComp());
  CH_assert(a_RHSlvl.interval().contains(a_intv));

  const DisjointBoxLayout& boxes = a_RHSlvl.disjointBoxLayout();
  CH_assert(boxes.hasLclMMB());  // Block indices available in layout

  for (DataIterator dit = boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = boxes[dit];
      const int idxBlk = boxes.blockIndex(dit);
      const BlockDomain& blockDomain = getCoordSys().problemDomain(idxBlk);
      const Box testBlkBox = grow(blockDomain.domainBox(), -1);
      if (!testBlkBox.contains(disjointBox))
        {
          const FArrayBox& RHSfab = a_RHSlvl[dit];
          const stc::Vector<Vector<Box>, 2*SpaceDim>& lvlLocations =
            getMBRegions().getLvlLocations(dit());
          for (const int dir : EachDir)
            {
              for (const auto side : EachSide)
                {
                  const int idxFace =
                    MultiBlockRegions::FaceTag::indexFace(dir, side);
                  if (lvlLocations[idxFace].size() > 0)
                    {
                      const int sideSign = Side::sign(side);
                      // Get the flux storage location
                      MultiBlockRegions::FaceTag faceTag(dir, side, dit());
                      const DataIndex& didxSto =
                        getMBRegions().getLvlStoDidx(faceTag);
                      const FArrayBox& fluxFab = m_lvlSto[didxSto];
                      // Local flux saved to inner cells.  These overlap the
                      // locations we want to updated.  Remote flux save to
                      // outer cells.
                      for (Box cells : lvlLocations[idxFace])
                        {
                          CH_assert(disjointBox.contains(cells));
                          // Remove half of the local flux from the divergence
                          // and add half of the remote flux
                          // On the high side, a positive flux subtracts from
                          // the state.  On the low side, a positive flux adds
                          // to the state.
                          const Real factor = -sideSign*0.5/a_dx[dir];
                          for (int cS = a_intv.begin(),
                                 cS_end = a_intv.end() + 1, cR = 0;
                               cS != cS_end; ++cS, ++cR)
                            {
                              MD_BOXLOOP(cells, i)
                                {
                                  IntVect ivRmt = MD_GETIV(i);
                                  ivRmt.shift(dir, sideSign);  // To outside
                                  RHSfab[MD_IX(i, cS)] +=
                                    factor*(fluxFab[MD_IV(ivRmt, cR)] -
                                            fluxFab[MD_IX(i, cR)]);
                                }
                            }
                        }
                    }
                }  // Loop over sides
            }  // Loop over directions
        }  // Test for adjacency to block boundary
    }  // Loop over local boxes
}
