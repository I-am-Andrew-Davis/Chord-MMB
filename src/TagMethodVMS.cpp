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
 * \file TagMethodVMS.cpp
 *
 * \brief Member functions for TagMethodVMS
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

#include "TagMethodVMS.H"
#include "CRDparam.H"
#include "PolytropicPhysicsF_F.H"
#include "DataTemp.H"
#include "CRDPhysics.H"
#include "PatchMappedFunc.H"


/*******************************************************************************
 *
 * Class TagMethodVMS: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_comp  Component to use for value
 *  \param[in]  a_loThreshold
 *                      Lower threshold for refinement.  Only cells with
 *                      values above this value will be tagged.
 *  \param[in]  a_upThreshold
 *                      Upper threshold for refinement.  Only cells with
 *                      values below this value will be tagged.
 *  \param[in]  a_usePrimitive
 *                      T - Values of primitive variables (default)
 *                      F - Values of conservative variables
 *  \param[in]  a_insideThreshold
 *                      T - Tag cells with values inside the thresholds
 *                      F - Tag cells with values outside the thresholds
 *  \param[in]  a_useVorticity
 *                      T - Tag cells based on upThreshold of vorticity
 *                      F - Tag cells based on prim/cons values
 *  \param[in]  a_maxLevel
 *                      Max level to tag with this method (-1 means all
 *                      levels are tagged with this method)
 *//*-----------------------------------------------------------------*/

TagMethodVMS::TagMethodVMS(const int  a_comp,
                           const Real a_loThreshold,
                           const Real a_upThreshold,
                           const bool a_usePrimitive,
                           const bool a_insideThreshold,
                           const bool a_useVorticity,
                           const int  a_maxLevel)
  :
  m_comp(a_comp),
  m_loThreshold(a_loThreshold),
  m_upThreshold(a_upThreshold),
  m_usePrimitive(a_usePrimitive),
  m_insideThreshold(a_insideThreshold),
  m_useVorticity(a_useVorticity),
  m_maxLevel(a_maxLevel)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodVMS::~TagMethodVMS()
{
}

/*--------------------------------------------------------------------*/
//  Tag all cells with values either inside or outside of the thresholds
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
TagMethodVMS::operator()(IntVectSet&               a_tags,
                         const Box&                a_box,
                         const DataIndex&          a_didx,
                         const MappedLevelData&    a_data,
                         LevelGridMetrics&         a_gridMetrics,
                         const MultiBlockCoordSys& a_MBCoordSys,
                         const Real                a_time,
                         const int                 a_level)
{
  // FIXME: (1) Add some support for tagging other turbulence based
  //            parameters such as enstrophy or turbulent dissipation
  //        (2) Include a user-controlled parameter to specify how
  //            many cells wide the projection should be
  //        (3) Since the values are being computed throughout the
  //            box, a cutoff based on the domain-wide sum of the
  //            high-frequency data (and a statistical parameter
  //            based on that data) could be nice to use
  //
  // If the level is above the max level for this tagging, do not tag
  if (a_level >= m_maxLevel && m_maxLevel != -1)
    {
      return;
    }
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box);
  const Box box1 = grow(a_box, 1);
  Box box1Dom(box1);
  box1Dom &= blockDomain;
  int useComp = m_comp;

//--Get the data

  const FArrayBox& Ufab = a_data.getU(2,2)[a_didx];
  const FArrayBox* phiPtr = &Ufab;
  // primative state
  FABSTACKTEMP(Wfab, box1Dom, CRDparam::g_CRDPhysics->numPrimitive());
  // Magnitude of vorticity
  const Box vorticityBox = box1Dom;
  FABSTACKTEMP(magVort, vorticityBox, 1);
  // solve primative vars
  if (m_usePrimitive)
    {
      CRDparam::g_CRDPhysics->consToPrim(Wfab, Ufab, box1Dom, Wfab);
      if (CRDparam::g_CRDPhysics->extraPrimInterval().contains(m_comp))
        {
          CRDparam::g_CRDPhysics->extraPrimitiveState(Wfab, box1Dom);
        }
      phiPtr = &Wfab;
    }
  // solve vorticity magnitude
  else if (m_useVorticity)
    {
      useComp = 0;
      const int endComp = PatchMappedFunc::m_numCurlComps;

      FABSTACKTEMP(vort, vorticityBox, endComp);
      Box boxP = grow(vorticityBox, 1);
      FABSTACKTEMP(vel, boxP, velIntv.size());
      // compute velocity
      PatchMappedFunc::divideVec(boxP,
                                 vel,
                                 Ufab,
                                 velIntv,
                                 rhoIndx,
                                 Interval(0, velIntv.size()-1));
      // compute mapped vorticity
      PatchMappedFunc::curl2OPS(vorticityBox,
                                vort,
                                vel,
                                a_gridMetrics.m_N[a_didx],
                                a_gridMetrics.m_J[a_didx],
                                blockDomain,
                                Interval(0, velIntv.size()-1),
                                a_gridMetrics.dxVect());
          
      // compute magnitude of vorticity
      PatchMappedFunc::magnitude(vorticityBox,
                                 magVort,
                                 0,
                                 vort,
                                 Interval(0, endComp));
      // set the tag pointer
      phiPtr = &magVort;
    }

  // Isolate the highest wavenumbers through projection
  FABSTACKTEMP(sepData, a_box, 1);
  sepData.setVal(0.);
  const FArrayBox& phi = *phiPtr;
  scaleSep(sepData, phi, a_box, a_gridMetrics, a_MBCoordSys);

  // tag values
  MD_ARRAY_RESTRICT(arrPhi, sepData);
  MD_ARRAY_RESTRICT(arrOrig, phi);
  // if both sepData and phi are 0, we need to prevent tagging
  const Real eps = 1e-12;
  // contained inside threshold
  if (m_insideThreshold)
    {
      MD_BOXLOOP(a_box, i)
        {
          const Real phiVal = arrPhi[MD_IX(i, 0)];
          const Real origVal = std::abs(arrOrig[MD_IX(i, useComp)]);
          if (phiVal > m_loThreshold*origVal &&
              phiVal < m_upThreshold*origVal)
            {
              a_tags |= IntVect(D_DECL6(i0, i1, i2, i3, i4, i5));
            }
        }
    }
  // contained outside threshold
  else
    {
      MD_BOXLOOP(a_box, i)
        {
          const Real phiVal = arrPhi[MD_IX(i, 0)];
          const Real origVal = std::abs(arrOrig[MD_IX(i, useComp)]);
          if (phiVal > (m_upThreshold*origVal + eps))
            {
              a_tags |= IntVect(D_DECL6(i0, i1, i2, i3, i4, i5));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Separate scales of variable of interest
/** \param[out] a_sepData  
 *                      Smallest scales of original data
 *  \param[in]  a_origData
 *                      Original, uncoarsened data
 *  \param[in]  a_box   Disjoint box over which data is decomposed
 *  \param[in]  a_didx  Data index of the current box
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[out] a_gridMetrics
 *                      Only caches should be modified
 *  \param[in]  a_MBCoordSys
 *                      Coordinate system describing the mapping
 *//*-----------------------------------------------------------------*/

void
TagMethodVMS::scaleSep(FArrayBox&                a_sepData,
                       const FArrayBox&          a_origData,
                       const Box&                a_box,
                       LevelGridMetrics&         a_gridMetrics,
                       const MultiBlockCoordSys& a_MBCoordSys)
{
  // FIXME: (1) Update this so that it definitely works with mapping
  //        (2) Make the coarsening factor user adjustable
  //        (3) Make the shifting stuff a little more rigorous
  //        (4) Place some reasonable and necessary checks here
  //            (for example, is the blocking factor correct?)
  //
  // Coarsening factor = 2
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box);
  Box coarseBox = a_box;
  coarseBox.coarsen(2);
  const Real volumeScaling = 1./(std::pow(2.,SpaceDim));
  const int c = m_comp;
  // Set aside and initialize containers for data
  FABSTACKTEMP(coarseData, a_box, 1);
  coarseData.setVal(0.);
  Box box1 = grow(a_box,1);
  box1 &= blockDomain;
  FABSTACKTEMP(coarseDataShiftDiag, box1, 1);
  coarseDataShiftDiag.setVal(0.);

  MD_ARRAY_RESTRICT(arrFine, a_origData);
  MD_ARRAY_RESTRICT(arrC, coarseData);
  MD_ARRAY_RESTRICT(arrCD, coarseDataShiftDiag);
  MD_ARRAY_RESTRICT(arrSep, a_sepData);
  MD_BOXLOOP(coarseBox, i)
    {
      D_TERM(
        Real f0 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1,  i2)),c)];
        Real f1 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1,  i2)),c)];,
        Real f2 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1+1,i2)),c)];
        Real f3 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2)),c)];,
        Real f4 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1,  i2+1)),c)];
        Real f5 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1,  i2+1)),c)];
        Real f6 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1+1,i2+1)),c)];
        Real f7 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+1)),c)];);
      Real ave =
        volumeScaling*(D_TERM(f0 + f1,+ f2 + f3,+ f4 + f5 + f6 + f7));
      D_TERM(
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1,  i2)),0)] =
        std::abs(f0 - ave);
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1,  i2)),0)] =
        std::abs(f1 - ave);,
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1+1,i2)),0)] =
        std::abs(f2 - ave);
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2)),0)] =
        std::abs(f3 - ave);,
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1,  i2+1)),0)] =
        std::abs(f4 - ave);
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1,  i2+1)),0)] =
        std::abs(f5 - ave);
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0,  i1+1,i2+1)),0)] =
        std::abs(f6 - ave);
        arrC[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+1)),0)] =
        std::abs(f7 - ave););
    }
  // We shift the boxes around and average each time
  // This way, we can capture the wavenumbers more effectively
  Box coarseBoxShiftDiag = grow(a_box,1);
  BlockDomain coarseBlockDomain = blockDomain;
  coarseBlockDomain.coarsen(2);
  // We must shift the box to avoid a possible Chombo error
  // wherein boxes with odd lower IntVects and even upper IntVects
  // result in an upper IntVect one too far after coarsening
  coarseBoxShiftDiag.shift(IntVect(D_DECL(1,1,1)));
  coarseBoxShiftDiag.coarsen(2);
  coarseBoxShiftDiag &= coarseBlockDomain;
  coarseBoxShiftDiag.shift(IntVect(D_DECL(-1,-1,-1)));
  coarseBoxShiftDiag &= coarseBlockDomain;
  MD_BOXLOOP(coarseBoxShiftDiag, i)
    {
      D_TERM(
        Real f0 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+1)),c)];
        Real f1 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+1,i2+1)),c)];,
        Real f2 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+2,i2+1)),c)];
        Real f3 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+2,i2+1)),c)];,
        Real f4 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+2)),c)];
        Real f5 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+1,i2+2)),c)];
        Real f6 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+2,i2+2)),c)];
        Real f7 =
        arrFine[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+2,i2+2)),c)];);
      Real ave =
        volumeScaling*(D_TERM(f0 + f1,+ f2 + f3,+ f4 + f5 + f6 + f7));
      D_TERM(
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+1)),0)] =
        std::abs(f0 - ave);
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+1,i2+1)),0)] =
        std::abs(f1 - ave);,
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+2,i2+1)),0)] =
        std::abs(f2 - ave);
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+2,i2+1)),0)] =
        std::abs(f3 - ave);,
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+1,i2+2)),0)] =
        std::abs(f4 - ave);
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+1,i2+2)),0)] =
        std::abs(f5 - ave);
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+1,i1+2,i2+2)),0)] =
        std::abs(f6 - ave);
        arrCD[MD_OFFSETIV(i,+,IntVect(D_DECL(i0+2,i1+2,i2+2)),0)] =
        std::abs(f7 - ave););
    }
  MD_BOXLOOP(a_box,i)
    {
      arrSep[MD_IX(i,0)] = arrC[MD_IX(i,0)] + arrCD[MD_IX(i,0)];
    }
}
