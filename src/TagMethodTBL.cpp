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
 * \file TagMethodTBL.cpp
 *
 * \brief Member functions for TagMethodTBL
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

#include "TagMethodTBL.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "CNSIBC.H"


/*******************************************************************************
 *
 * Class TagMethodTBL: member definitions
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

TagMethodTBL::TagMethodTBL(const int            a_idxBlk,
                           const int            a_dir,
                           const Side::LoHiSide a_side,
                           const int            a_streamwiseDir,
                           const Real           a_delta,
                           const Real           a_ReynoldsNumber)
  :
  m_idxBlk(a_idxBlk),
  m_dir(a_dir),
  m_side(a_side),
  m_streamwiseDir(a_streamwiseDir),
  m_delta(a_delta),
  m_Re(a_ReynoldsNumber)
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
TagMethodTBL::operator()(IntVectSet&               a_tags,
                         const Box&                a_box,
                         const DataIndex&          a_didx,
                         const MappedLevelData&    a_data,
                         LevelGridMetrics&         a_gridMetrics,
                         const MultiBlockCoordSys& a_MBCoordSys,
                         const Real                a_time,
                         const int                 a_level)
{
#if CH_SPACEDIM == 3
  // Get the spanwise direction (assuming there is such a thing)
  int spanwiseDir = 0;
  if ((!m_dir || m_dir == 2) && (!m_streamwiseDir || m_streamwiseDir == 2))
    {
      spanwiseDir = 1;
    }
  else if ((!m_dir || m_dir == 1) && (!m_streamwiseDir || m_streamwiseDir == 1))
    {
      spanwiseDir = 2;
    }
#endif

  // Get the physical coordinates
  FABSTACKTEMP(XiFab, a_box, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_box, SpaceDim);   // Physical coordinates
  const int idxBlk = a_gridMetrics.getBoxes().blockIndex(a_didx);
  CH_assert(idxBlk != -1);
  const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(idxBlk));
  CRDparam::g_CNSIBC->getCellCoordinates(a_box, XiFab, XFab, blockCoordSys);

  // Set up one-dimensional box for computing delta
  Box low1DBox(a_box.smallEnd(), a_box.smallEnd());
  low1DBox.setBig(m_streamwiseDir, a_box.bigEnd()[m_streamwiseDir]);

  // Get the low end of the block
  //**FIXME: This needs to be the actual low side of the block
  const Box lowBox(IntVect_zero, IntVect_unit);
  const RealVect lowXi = RealVect_zero;
  const RealVect lowX = a_gridMetrics.getCoordSys(lowBox)->realCoord(lowXi);

  // Now determine local boundary-layer thickness
  const Real x_star = (m_delta/0.37)*std::pow(m_Re, 0.2);
  const Real c_0 = m_Re/x_star;
  FABSTACKTEMP(delta1D, low1DBox, 1);
  MD_BOXLOOP(low1DBox, i)
    {
      const Real x_equiv = x_star + XFab[MD_IX(i, m_streamwiseDir)]
        - lowX[m_streamwiseDir];
      const Real Re_x = c_0*x_equiv;
      delta1D[MD_IX(i, 0)] = 0.37*x_equiv/std::pow(Re_x, 0.2);
    }

  // Fill the full box with the boundary-layer thickness values
  FABSTACKTEMP(delta, a_box, 1);
  MD_BOXLOOP(low1DBox, i)
    {
      Box normBox(MD_GETIV(i), MD_GETIV(i));
      normBox.setBig(m_dir, a_box.bigEnd()[m_dir]);
#if CH_SPACEDIM == 3
      normBox.setBig(spanwiseDir, a_box.bigEnd()[spanwiseDir]);
#endif
      MD_BOXLOOP(normBox, j)
        {
          delta[MD_IX(j, 0)] = delta1D[MD_IX(i, 0)];
        }
    }

  // If a cell is within the boundary-layer, tag it
  MD_BOXLOOP(a_box, i)
    {
      if ((XFab[MD_IX(i, m_dir)]-lowX[m_dir]) <= delta[MD_IX(i, 0)])
        {
          a_tags |= MD_GETIV(i);
        }
    }
}
