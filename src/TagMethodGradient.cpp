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
 * \file TagMethodGradient.cpp
 *
 * \brief Member functions for TagMethodGradient
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "Box.H"
#include "DataIndex.H"
#include "LoHiCenter.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"

//----- Internal -----//

#include "TagMethodGradient.H"
#include "CRDparam.H"
#include "PolytropicPhysicsF_F.H"
#include "DataTemp.H"
#include "CRDPhysics.H"


/*******************************************************************************
 *
 * Class TagMethodGradient: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
//  FIXME: The current method requires a new TagMethod to be set up for each
//         tag variable However, this is much slower than if one tag method
//         knew all variables, meaning operations like consToPrim would
//         only run once
/** \param[in]  a_comp  Component to use for gradient
 *  \param[in]  a_threshold
 *                      Threshold for refinement.  Only cells with
 *                      gradients above this value will be tagged.
 *  \param[in]  a_usePrimitive
 *                      T - Gradients of primitive variables
 *                      F - Gradients of conservative variables
 *  \param[in]  a_loVal If value is below this number, do not tag
 *  \param[in]  a_restrictBox
 *                      Vector of boxes for which refinement does not occur
 *                      These boxes should only be for the base grid
 *//*-----------------------------------------------------------------*/

TagMethodGradient::TagMethodGradient(const int        a_comp,
                                     const Real       a_threshold,
                                     const bool       a_usePrimitive,
                                     const Real       a_loVal,
                                     std::vector<Box> a_restrictBox)
  :
  m_comp(a_comp),
  m_threshold(a_threshold),
  m_maxGradMag(0.),
  m_usePrimitive(a_usePrimitive),
  m_loVal(a_loVal),
  m_restrictBox(a_restrictBox)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodGradient::~TagMethodGradient()
{
}

/*--------------------------------------------------------------------*/
//  Tag all cells with gradients that exceed a threshold
/** \param[out] a_tags  IntVectSet of cells to refine local to
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
TagMethodGradient::operator()(IntVectSet&               a_tags,
                              const Box&                a_box,
                              const DataIndex&          a_didx,
                              const MappedLevelData&    a_data,
                              LevelGridMetrics&         a_gridMetrics,
                              const MultiBlockCoordSys& a_MBCoordSys,
                              const Real                a_time,
                              const int                 a_level)
{
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box);
  const Box box1 = grow(a_box, 1);
  Box box1Dom(box1);
  box1Dom &= blockDomain;

//--Get the data

  const FArrayBox& Ufab = a_data.getU(1,1)[a_didx];
  const FArrayBox* phiPtr = &Ufab;
  FABSTACKTEMP(Wfab, box1Dom, CRDparam::g_CRDPhysics->numPrimitive());
  if (m_usePrimitive)
    {
      // FIXME: Not all cells in box1Dom in Ufab have data in them
      CRDparam::g_CRDPhysics->consToPrim(Wfab, Ufab, box1Dom, Wfab);
      if (CRDparam::g_CRDPhysics->extraPrimInterval().contains(m_comp))
        {
          CRDparam::g_CRDPhysics->extraPrimitiveState(Wfab, box1Dom);
        }
      phiPtr = &Wfab;
    }

//--Compute the relative gradient

  FABSTACKTEMP(gradFab, a_box, SpaceDim);

  const FArrayBox& phi = *phiPtr;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box box1Dir(a_box);
      box1Dir.grow(dir, 1);
      Box centerCells, entireCells, loCells, hiCells;
      int hasLoCells, hasHiCells;
      loHiCenter(loCells, hasLoCells,
                 hiCells, hasHiCells,
                 centerCells, entireCells,
                 box1Dir, blockDomain, dir);

      FORT_GETRELGRADF(CHF_FRA1(gradFab, dir),
                       CHF_CONST_FRA1(phi, m_comp),
                       CHF_CONST_REAL(m_loVal),
                       CHF_CONST_INT(dir),
                       CHF_BOX(loCells),
                       CHF_CONST_INT(hasLoCells),
                       CHF_BOX(hiCells),
                       CHF_CONST_INT(hasHiCells),
                       CHF_BOX(centerCells));
    }

  FABSTACKTEMP(gradMagFab, a_box, 1);
  FORT_MAGNITUDEF(CHF_FRA1(gradMagFab,0),
                  CHF_CONST_FRA(gradFab),
                  CHF_BOX(a_box));

  m_maxGradMag = 0.;
  MD_ARRAY_RESTRICT(arrGMF, gradMagFab);
  // For when we don't have boxes to restrict refinement
  if (m_restrictBox.size() == 0)
    {
      MD_BOXLOOP(a_box, i)
        {
          const Real diff = arrGMF[MD_IX(i, 0)];
          m_maxGradMag = std::max(m_maxGradMag, diff);
          if (diff >= m_threshold)
            {
              a_tags |= IntVect(D_DECL6(i0 ,i1, i2, i3, i4, i5));
            }
        }
    }
  else
    {
      MD_BOXLOOP(a_box, i)
        {
          const Real diff = arrGMF[MD_IX(i, 0)];
          m_maxGradMag = std::max(m_maxGradMag, diff);
          if (diff >= m_threshold)
            {
              bool nonRestricted = true;
              for (int iv = 0; iv != m_restrictBox.size(); ++iv)
                {
                  Box checkBox = m_restrictBox[iv];
                  checkBox.refine(CRDparam::g_refFromBase[a_level]);
                  if (checkBox.contains(
                       IntVect(D_DECL6(i0, i1, i2, i3, i4, i5))))
                    {
                      nonRestricted = false;
                      continue;
                    }
                }
              // If not inside restricted boxes
              if (nonRestricted)
                {
                  a_tags |= IntVect(D_DECL6(i0 ,i1, i2, i3, i4, i5));
                }
            }
        }
    }
}
