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
 * \file TagMethodPhysBox.cpp
 *
 * \brief Member functions for TagMethodPhysBox
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "Box.H"
#include "DataIndex.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"

//----- Internal -----//

#include "TagMethodPhysBox.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "CNSIBC.H"


/*******************************************************************************
 *
 * Class TagMethodPhysBox: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_Xlo   Low physical coordinates of the box
 *  \param[in]  a_Xhi   High physical coordinates of the box
 *  \param[in]  a_fastExlusion
 *                      T - performs a quick test to see if a box
 *                          should be further analysed for having
 *                          cells in the tag box.  For very complex
 *                          geometry, where it is difficult to find
 *                          the minbox, set to false.
 *//*-----------------------------------------------------------------*/

TagMethodPhysBox::TagMethodPhysBox(const RealVect& a_Xlo,
                                   const RealVect& a_Xhi,
                                   const bool      a_fastExclusion)
  :
  m_Xlo(a_Xlo),
  m_Xhi(a_Xhi),
  m_fastExclusion(a_fastExclusion)
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Tag all cells in a box in physical space
/** If m_fastExlcusion is T, a minbox representing the physical space
 *  of a_box is estimated and intersected with the tag box to exclude
 *  boxes from detailed analysis by sampling only corners and face
 *  centers.
 *  \param[out] a_tags  IntVectSet of cells to refine local to
 *                      'a_box'.
 *  \param[in]  a_box   Disjoint box in which cells should be marked
 *                      for refinement
 *  \param[in]  a_didx  Data index of the current box
 *  \param[in]  a_data  Data on the level, both mapped and physical.
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[out] a_gridMetrics
 *                      Only caches should be modified
 *  \param[in]  a_MBCoordSys
 *                      Coordinate system describing the mapping
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
TagMethodPhysBox::operator()(IntVectSet&               a_tags,
                             const Box&                a_box,
                             const DataIndex&          a_didx,
                             const MappedLevelData&    a_data,
                             LevelGridMetrics&         a_gridMetrics,
                             const MultiBlockCoordSys& a_MBCoordSys,
                             const Real                a_time,
                             const int                 a_level)
{
  const int idxBlk = a_gridMetrics.getBoxes().blockIndex(a_didx);
  CH_assert(idxBlk != -1);
  const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(idxBlk));
  // const BlockDomain& blockDomain =
  //   a_gridMetrics.getCoordSys().problemDomain(idxBlk);

  // Do we have to check a_box in detail or can we exclude it?
  if (m_fastExclusion)
    {
      // Try to get a minbox in physical space representing a_box.  Check all
      // corners and face centers.
      const IntVect dims = a_box.size();
      Box corners(IntVect_zero, IntVect_unit);
      RealVect boxXlo(std::numeric_limits<Real>::max());
      RealVect boxXhi(std::numeric_limits<Real>::lowest());
      MD_BOXLOOP(corners, i)
        {
          const RealVect nodeXi(
            a_box.smallEnd() + MD_GETIV(i)*a_box.size());
          const RealVect nodeX = blockCoordSys.realCoord(nodeXi);
          boxXlo = stc::min(boxXlo, nodeX);
          boxXhi = stc::max(boxXhi, nodeX);
        }
      const RealVect halfDims = 0.5*dims;
      const RealVect center = a_box.smallEnd() + halfDims;
      for (const int dir : EachDir)
        {
          for (const auto side : EachSide)
            {
              RealVect nodeXi(center);
              nodeXi[dir] += sign(side)*halfDims[dir];
              const RealVect nodeX = blockCoordSys.realCoord(nodeXi);
              boxXlo = stc::min(boxXlo, nodeX);
              boxXhi = stc::max(boxXhi, nodeX);
            }
        }
      CH_assert(boxXlo <= boxXhi);
      boxXlo.max(m_Xhi);
      boxXhi.min(m_Xlo);
      if (!(m_Xlo <= m_Xhi)) return;
    }
  
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_box, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_box, SpaceDim);   // Physical coordinates
  CRDparam::g_CNSIBC->getCellCoordinates(a_box, XiFab, XFab, blockCoordSys);

  MD_BOXLOOP(a_box, i)
    {
      RealVect X(D_DECL(XFab[MD_IX(i, 0)],
                        XFab[MD_IX(i, 1)],
                        XFab[MD_IX(i, 2)]));
      if (contains(X))
        {
          a_tags |= MD_GETIV(i);
        }
    }
}
