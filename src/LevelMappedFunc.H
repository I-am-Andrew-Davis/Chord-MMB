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
 * \file LevelMappedFunc.H
 *
 * \brief Level functions on mapped grids
 *
 *//*+*************************************************************************/

#ifndef _LEVELMAPPEDFUNC_H_
#define _LEVELMAPPEDFUNC_H_

//----- Chombo Library -----//

#include "UsingNamespace.H"

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Interval;
class LevelGridMetrics;
class FArrayBox;
template <typename T> class LevelData;
#include "NamespaceFooter.H"


/*******************************************************************************
 *
 *  Level functions on mapped grids
 *
 ******************************************************************************/

namespace LevelMappedFunc
{
  /// Obtain magnitude of physical space gradient of spatial vector phi
  void magGradient2OPS(LevelData<FArrayBox>&       a_magGradPhi,
                       const int                   a_begCompMGP,
                       const IntVect&              a_numGhost,
                       const LevelData<FArrayBox>& a_phi,
                       const Interval&             a_intvPhi,
                       const LevelGridMetrics&     a_levelGridMetrics);

  /// Obtain physical space divergence of vector phi
  void divergence2OPS(LevelData<FArrayBox>&       a_divPhi,
                      const int                   a_begCompMGP,
                      const IntVect&              a_numGhost,
                      const LevelData<FArrayBox>& a_phi,
                      const Interval&             a_intvPhi,
                      const LevelGridMetrics&     a_levelGridMetrics,
                      const int                   a_normComp = -1);

  /// Obtain physical space curl of vector phi
  void curl2OPS(LevelData<FArrayBox>&       a_curlPhi,
                const int                   a_begCompMGP,
                const IntVect&              a_numGhost,
                const LevelData<FArrayBox>& a_phi,
                const Interval&             a_intvPhi,
                const LevelGridMetrics&     a_levelGridMetrics,
                const int                   a_normComp = -1);
}

#endif  /* ! defined _LEVELMAPPEDFUNC_H_ */
