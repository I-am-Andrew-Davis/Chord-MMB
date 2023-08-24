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
 * \file MMBSingleLevel.H
 *
 * \brief Implements multiblock-operations relevant to a single level
 *
 * MMB is a chaotic mess.  To achieve sanity, the class in this file is used to
 * interface with the HOMappedGeometry library for operations relevant to a
 * single AMR level.  It also serves as a wrapper in case we want to re-write
 * the MMB implementation in the future.
 *
 *//*+*************************************************************************/

#ifndef _MMBSINGLELEVEL_H_
#define _MMBSINGLELEVEL_H_

//----- Standard Library -----//

//----- Chombo Library -----//

#include "UsingNamespace.H"

//----- Internal -----//

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class LevelGridMetrics;
#include "NamespaceFooter.H"


/*******************************************************************************
 */
///  MMBSingleLevel
/**
 *   \note
 *   <ul>
 *     <li> The CS is obtained through LevelGridMetrics.  This means LGM must be
 *          defined before this class is defined with a grid.  The CS is assumed
 *          to be invariant for the course of a run.
 *   </ul>
 *
 *//*+*************************************************************************/

class MMBSingleLevel
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  MMBSingleLevel(const LevelGridMetrics& a_levelGridMetrics);

  /// Destructor
  ~MMBSingleLevel() = default;

  // Generally prohibit all others
  MMBSingleLevel(const MMBSingleLevel&) = delete;
  MMBSingleLevel(MMBSingleLevel&&) = delete;
  MMBSingleLevel& operator=(const MMBSingleLevel&) = delete;
  MMBSingleLevel& operator=(MMBSingleLevel&&) = delete;

  /// Weak construction
  void define(const int a_maxComp);


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// Retrieve the multiblock coordinate system
  const MultiBlockCoordSys& getCoordSys() const
    { return m_levelGridMetrics.getCoordSys(); }

  /// Retrieve the multiblock regions
  const MultiBlockRegions& getMBRegions() const
    { return m_levelGridMetrics.getMBRegions(); }

  /// Set the flux on a connected block boundary
  void setFlux(const Box&       a_disjointBox,
               const FluxBox&   a_flux,
               const DataIndex& a_didx,
               const Interval&  a_intv);

  /// Exchange fluxes once all are known
  void exchange(const bool a_fluxSign = true);

  /// Reflux the average.  This updates the divergence in cells.
  void refluxAverage(LevelData<FArrayBox>& a_RHSlvl,
                     const RealVect&       a_dx,
                     const Interval&       a_intv);

  /// Apply some operator between the two states
  template <typename F>
  void applyAtBoundary(const DisjointBoxLayout& a_boxes,
                       F&&                      a_op);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  const LevelGridMetrics& m_levelGridMetrics;
                                      ///< Source for multi-block coordinate
                                      ///< system
  BoxLayoutData<FArrayBox> m_lvlSto;  ///< Data on a compact layout supporting
                                      ///< exchanges across mapped multiblock
                                      ///< boundaries
};


/*******************************************************************************
 *
 * Class MMBSingleLevel: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Apply some operator between the two states
/** \tparam F           Operator function object
 *  \param[in]  a_op    Instance of the operator
 *//*-----------------------------------------------------------------*/

template <typename F>
void
MMBSingleLevel::applyAtBoundary(const DisjointBoxLayout& a_boxes,
                                F&&                      a_op)
{
  CH_assert(a_boxes.hasLclMMB());  // Block indices available in layout
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = a_boxes[dit];
      const int idxBlk = a_boxes.blockIndex(dit);
      const BlockDomain& blockDomain = getCoordSys().problemDomain(idxBlk);
      const Box testBlkBox = grow(blockDomain.domainBox(), -1);
      if (!testBlkBox.contains(disjointBox))
        {
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
                      // Get the flux storage location
                      MultiBlockRegions::FaceTag faceTag(dir, side, dit());
                      const DataIndex& didxSto =
                        getMBRegions().getLvlStoDidx(faceTag);
                      const FArrayBox& dataFab = m_lvlSto[didxSto];
                      // Local flux saved to inner cells.  These overlap the
                      // locations we want to updated.  Remote flux save to
                      // outer cells.
                      for (Box cells : lvlLocations[idxFace])
                        {
                          CH_assert(disjointBox.contains(cells));
                          // Do your business
                          a_op(cells, dir, side, dit(), dataFab);
                        }
                    }
                }  // Loop over sides
            }  // Loop over directions
        }  // Test for adjacency to block boundary
    }  // Loop over local boxes
}

#endif  /* ! defined _MMBSINGLELEVEL_H_ */
