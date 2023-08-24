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
 * \file LevelMappedFunc.cpp
 *
 * \brief Functions in namespace LevelMappedFunc
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "LevelMappedFunc.H"
#include "PatchMappedFunc.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CNSIBC.H"
#include "DataTemp.H"


/*******************************************************************************
 *
 * Namespace LevelMappedFunc: function definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Obtain magnitude of physical space gradient of spatial vector phi
/** \param[out] a_magGradPhi
 *                      Level data of magnitude of gradient of phi
 *  \param[in]  a_begCompMGP
 *                      Component of 'a_magGradPhi' to store solution in
 *  \param[in]  a_phi   Vector phi
 *  \param[in]  a_intvPhi
 *                      Interval of 'a_phi' for vector.
 *                      Spatial gradients are assumed contiguous
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the current level
 *//*-----------------------------------------------------------------*/

void
LevelMappedFunc::magGradient2OPS(LevelData<FArrayBox>&       a_magGradPhi,
                                 const int                   a_begCompMGP,
                                 const IntVect&              a_numGhost,
                                 const LevelData<FArrayBox>& a_phi,
                                 const Interval&             a_intvPhi,
                                 const LevelGridMetrics&     a_levelGridMetrics)
{
  CH_TIME("LevelMappedfunc::magGradient2OPS");
  CH_assert(a_begCompMGP + a_intvPhi.size() - 1 < a_magGradPhi.nComp());
  CH_assert(a_intvPhi.end() < a_phi.nComp());
  CH_assert(a_phi.ghostVect() >= (a_numGhost + IntVect::Unit));

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_phi.disjointBoxLayout()[dit];
      const BlockDomain& domain =
        a_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box box = grow(disjointBox, a_numGhost);
      FABSTACKTEMP(gradPhiCS, box, a_intvPhi.size()*SpaceDim);
      PatchMappedFunc::gradient2OCS(box,
                                    gradPhiCS,
                                    a_phi[dit],
                                    domain,
                                    a_intvPhi,
                                    a_levelGridMetrics.dxVect());

      FABSTACKTEMP(gradPhiPS, box, a_intvPhi.size()*SpaceDim);
      PatchMappedFunc::gradientCStoPSavg(box,
                                         gradPhiPS,
                                         gradPhiCS,
                                         a_levelGridMetrics.m_N[dit],
                                         a_levelGridMetrics.m_J[dit]);

      // solve the magnitude for each component of phi
      for (int c = 0; c != a_intvPhi.size(); ++c)
        {
          PatchMappedFunc::magnitudeSpace(box,
                                          a_magGradPhi[dit],
                                          c + a_begCompMGP,
                                          gradPhiPS,
                                          c*SpaceDim);
        }
    }
}

/*--------------------------------------------------------------------*/
/// Obtain physical space divergence of vector phi
/** \param[out] a_divPhi
 *                      Level data of magnitude of gradient of phi
 *  \param[in]  a_divMGP
 *                      Component of 'a_magGradPhi' to store solution in
 *  \param[in]  a_phi   Vector phi
 *  \param[in]  a_intvPhi
 *                      Interval of 'a_phi' for vector.
 *                      Spatial gradients are assumed contiguous
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the current level
 *  \param[in]  a_divideComp
 *                      Component to divide by, used for getting
 *                      velocity from momentum. Ignored if negative
 *//*-----------------------------------------------------------------*/
void
LevelMappedFunc::divergence2OPS(LevelData<FArrayBox>&       a_divPhi,
                                const int                   a_begCompDiv,
                                const IntVect&              a_numGhost,
                                const LevelData<FArrayBox>& a_phi,
                                const Interval&             a_intvPhi,
                                const LevelGridMetrics&     a_levelGridMetrics,
                                const int                   a_divideComp)
{
  CH_TIME("LevelMappedfunc::divergence2OPS");
  CH_assert(a_begCompDiv < a_divPhi.nComp());
  CH_assert(a_intvPhi.end() < a_phi.nComp());
  CH_assert(a_phi.ghostVect() >= (a_numGhost + IntVect::Unit));
  CH_assert(a_divideComp < a_phi.nComp());

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_phi.disjointBoxLayout()[dit];
      const BlockDomain& domain =
        a_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box box = grow(disjointBox, a_numGhost);
      const int vecSize = a_intvPhi.size()-1;
      FArrayBox divPhi(Interval(a_begCompDiv, a_begCompDiv), a_divPhi[dit]);
      if (a_divideComp >= 0)
        {
          Box boxP = a_phi[dit].box();
          FABSTACKTEMP(phiMod, boxP, a_intvPhi.size());
          PatchMappedFunc::divideVec(boxP,
                                     phiMod,
                                     a_phi[dit],
                                     a_intvPhi,
                                     a_divideComp,
                                     Interval(0, a_intvPhi.size()-1));
          PatchMappedFunc::divergence2OPS(box,
                                          divPhi,
                                          phiMod,
                                          a_levelGridMetrics.m_N[dit],
                                          a_levelGridMetrics.m_J[dit],
                                          domain,
                                          Interval(0, vecSize),
                                          a_levelGridMetrics.dxVect());
        }
      else
        {
          PatchMappedFunc::divergence2OPS(box,
                                          divPhi,
                                          a_phi[dit],
                                          a_levelGridMetrics.m_N[dit],
                                          a_levelGridMetrics.m_J[dit],
                                          domain,
                                          a_intvPhi,
                                          a_levelGridMetrics.dxVect());
        }
    }
}

/*--------------------------------------------------------------------*/
/// Obtain physical space curl of vector phi
/** \param[out] a_magGradPhi
 *                      Level data of magnitude of gradient of phi
 *  \param[in]  a_begCompCurl
 *                      Component of 'a_curlPhi' to store solution in
 *  \param[in]  a_phi   Vector phi
 *  \param[in]  a_intvPhi
 *                      Interval of 'a_phi' for vector.
 *                      Spatial gradients are assumed contiguous
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the current level
 *  \param[in]  a_divideComp
 *                      Component to divide by, used for getting
 *                      velocity from momentum. Ignored if negative
 *//*-----------------------------------------------------------------*/
void
LevelMappedFunc::curl2OPS(LevelData<FArrayBox>&       a_curlPhi,
                          const int                   a_begCompCurl,
                          const IntVect&              a_numGhost,
                          const LevelData<FArrayBox>& a_phi,
                          const Interval&             a_intvPhi,
                          const LevelGridMetrics&     a_levelGridMetrics,
                          const int                   a_divideComp)
{
  CH_TIME("LevelMappedfunc::curl2OPS");
  CH_assert(a_intvPhi.end() < a_phi.nComp());
  CH_assert(a_phi.ghostVect() >= (a_numGhost + IntVect::Unit));
  CH_assert(a_divideComp < a_phi.nComp());

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = a_phi.disjointBoxLayout()[dit];
      const BlockDomain& domain =
        a_levelGridMetrics.getCoordSys().problemDomain(disjointBox);
      Box box = grow(disjointBox, a_numGhost);
      FArrayBox curlPhi(
        Interval(a_begCompCurl,
                 a_begCompCurl+PatchMappedFunc::m_numCurlComps-1),
        a_curlPhi[dit]);
      
      if (a_divideComp >= 0)
        {
          Box boxP = a_phi[dit].box();
          FABSTACKTEMP(phiMod, boxP, a_intvPhi.size());
          PatchMappedFunc::divideVec(boxP,
                                     phiMod,
                                     a_phi[dit],
                                     a_intvPhi,
                                     a_divideComp,
                                     Interval(0, a_intvPhi.size()-1));
          PatchMappedFunc::curl2OPS(box,
                                    curlPhi,
                                    phiMod,
                                    a_levelGridMetrics.m_N[dit],
                                    a_levelGridMetrics.m_J[dit],
                                    domain,
                                    Interval(0, a_intvPhi.size()-1),
                                    a_levelGridMetrics.dxVect());
        }
      else
        {
          PatchMappedFunc::curl2OPS(box,
                                    curlPhi,
                                    a_phi[dit],
                                    a_levelGridMetrics.m_N[dit],
                                    a_levelGridMetrics.m_J[dit],
                                    domain,
                                    a_intvPhi,
                                    a_levelGridMetrics.dxVect());
        }
    }
}
