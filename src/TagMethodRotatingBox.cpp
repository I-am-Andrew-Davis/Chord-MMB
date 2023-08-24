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
 * \file TagMethodRotatingBox.cpp
 *
 * \brief Member functions for TagMethodRotatingBox
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

#include "TagMethodRotatingBox.H"
#include "CRDparam.H"


/*******************************************************************************
 *
 * Class TagMethodRotatingBox: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_tagRotatingBox
 *                      Tag any cells in this box.  This box is
 *                      defined on the coarsest (base) level.  The
 *                      tagRotatingBox must be contained in the problem
 *                      domain.
 *//*-----------------------------------------------------------------*/

TagMethodRotatingBox::TagMethodRotatingBox(const Real& a_radius,
                                           const Real& a_cycleTime,
                                           const int   a_boxSize)
  :
  m_radius(a_radius),
  m_cycleTime(a_cycleTime),
  m_boxSize(a_boxSize)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodRotatingBox::~TagMethodRotatingBox()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Tag all cells in baseBox
/** If the baseBox is face-centered, it is first refined to the
 *  current level, grown by 1 in face-centered directions, and then
 *  enclosed cells are found.  This way, one can specify a thin
 *  layer of cells at the resolution of the current level.
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
TagMethodRotatingBox::operator()(IntVectSet&               a_tags,
                                 const Box&                a_box,
                                 const DataIndex&          a_didx,
                                 const MappedLevelData&    a_data,
                                 LevelGridMetrics&         a_gridMetrics,
                                 const MultiBlockCoordSys& a_MBCoordSys,
                                 const Real                a_time,
                                 const int                 a_level)
{
  Real theta = 2*3.14159*(a_time/m_cycleTime);
  Real dx = CRDparam::g_domainLength[0]/CRDparam::g_domainBaseSize[0];
  RealVect length = CRDparam::g_domainLength;
  D_TERM(
    Real xloc = std::cos(theta)*m_radius;,
    Real yloc = std::sin(theta)*m_radius;,);
  D_TERM(
    int xi = int(xloc/dx)+int(length[0]/(2.*dx))-1;,
    int yi = int(yloc/dx)+int(length[1]/(2.*dx))-1;,);
  IntVect center(D_DECL(xi, yi, 0));
  IntVect lo = center - m_boxSize;
  IntVect hi = center + m_boxSize;
  Box tagBox(lo,hi);
  tagBox.refine(CRDparam::g_refFromBase[a_level]);
  if (!tagBox.isEmpty())
    {
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (tagBox.type(dir) == IndexType::NODE)
            {
              tagBox.grow(dir, 1);
            }
        }
      tagBox.enclosedCells();
    }
  tagBox &= a_box;
  a_tags |= tagBox;
}
