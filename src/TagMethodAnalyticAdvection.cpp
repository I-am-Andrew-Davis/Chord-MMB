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
 * \file TagMethodAnalyticAdvection.cpp
 *
 * \brief Member functions for TagMethodAnalyticAdvection
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

#include "TagMethodAnalyticAdvection.H"
#include "CRDparam.H"


/*******************************************************************************
 *
 * Class TagMethodAnalyticAdvection: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

TagMethodAnalyticAdvection::TagMethodAnalyticAdvection(Real a_toggleCycle)
{
  m_toggleCycle = a_toggleCycle;
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodAnalyticAdvection::~TagMethodAnalyticAdvection()
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
TagMethodAnalyticAdvection::operator()(IntVectSet&               a_tags,
                                       const Box&                a_box,
                                       const DataIndex&          a_didx,
                                       const MappedLevelData&    a_data,
                                       LevelGridMetrics&         a_gridMetrics,
                                       const MultiBlockCoordSys& a_MBCoordSys,
                                       const Real                a_time,
                                       const int                 a_level)
{
  // This is specifically for analytically refining the advection of a gaussian
  // profile (see paper ... Guzik2015)
  RealVect p(D_DECL(0.5, 0.5, 0.5));
  RealVect av(D_DECL(1.0, 0.5, 0.0));
  p += a_time*av;
  IntVect offset(D_DECL((int)p[0], (int)p[1], (int)p[2]));
  D_TERM(p[0] -= (Real)offset[0];,
         p[1] -= (Real)offset[1];,
         p[2] -= (Real)offset[2];)
    Real r = 0.0;
  switch (a_level)
    {
    case 0:
      r = 0.35;
      break;
    case 1:
      r = 0.225;
      break;
    }

  int timeCycle = std::floor(a_time/m_toggleCycle);
  // move AMR off by a half period until a specified point
  if (((timeCycle < 1) && (m_toggleCycle > 0)))
    {
      p -= 0.5*RealVect::Unit;
      r /= 4;
    }

  const NewCoordSys *const coordSys = a_gridMetrics.getCoordSys(a_box);
  ShiftIterator shiftIt = a_gridMetrics.getCoordSys().problemDomain(a_box).shiftIterator();
  
  Real dx = 1.0/((Real)a_gridMetrics.getCoordSys().problemDomain(a_box).size(0));

  for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      RealVect Xi = bit();
      Xi += 0.5;
      Xi *= dx;
      const RealVect X = coordSys->realCoord(Xi);
      RealVect len = X - p;
      if (len.vectorLength() <= r)
        {
          a_tags |= bit();
        }
      else
        {
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              len = X;
              len += shiftIt();
              len -= p;
              if (len.vectorLength() <= r)
                {
                  a_tags |= bit();
                  break;
                }
            }
        }
    }
}
