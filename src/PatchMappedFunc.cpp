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
 * \file PatchMappedFunc.cpp
 *
 * \brief Functions in namespace PatchMappedFunc
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "LoHiCenter.H"
#include "FluxBox.H"
#include "UnitNormalsF_F.H"

//----- Internal -----//

#include "PatchMappedFunc.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDutil.H"
#include "DataTemp.H"

/*******************************************************************************
 *
 * Namespace PatchMappedFunc: function definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Obtain magnitude of spatial vector phi
//  faster than using magnitude function
/** \param[in]  a_box   Region to compute magnitude
 *  \param[out] a_mag   Magnitude
 *  \param[in]  a_cMag  Component of 'a_mag' to store magnitude
 *  \param[in]  a_phi   Vector phi
 *  \param[in]  a_cPhi  Start component in 'a_phi' for vector.
 *                      Spatial gradients are assumed contiguous in
 *                      components
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::magnitudeSpace(const Box&       a_box,
                                FArrayBox&       a_mag,
                                const int        a_cMag,
                                const FArrayBox& a_phi,
                                const int        a_cPhi)
{
  CH_assert(a_cMag >= 0 && a_cMag < a_mag.nComp());
  CH_assert(a_cPhi >= 0 && (a_cPhi + SpaceDim - 1) < a_phi.nComp());
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_mag.box().contains(a_box));
  
  MD_ARRAY_RESTRICT(magArr, a_mag);
  MD_ARRAY_RESTRICT(phiArr, a_phi);

  MD_BOXLOOP(a_box, i)
    {
      magArr[MD_IX(i, a_cMag)] = sqrt(
        D_TERM6(  std::pow(phiArr[MD_IX(i, a_cPhi + 0)], 2),
                + std::pow(phiArr[MD_IX(i, a_cPhi + 1)], 2),
                + std::pow(phiArr[MD_IX(i, a_cPhi + 2)], 2),
                + std::pow(phiArr[MD_IX(i, a_cPhi + 3)], 2),
                + std::pow(phiArr[MD_IX(i, a_cPhi + 4)], 2),
                + std::pow(phiArr[MD_IX(i, a_cPhi + 5)], 2)));
    }
}

/*--------------------------------------------------------------------*/
//  Obtain magnitude of a vector phi
/** \param[in]  a_box   Region to compute magnitude
 *  \param[out] a_mag   Magnitude
 *  \param[in]  a_cMag  Component of 'a_mag' to store magnitude
 *  \param[in]  a_phi   Vector phi
 *  \param[in]  a_cPhi  Interval of 'a_phi' for vector
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::magnitude(const Box&       a_box,
                           FArrayBox&       a_mag,
                           const int        a_cMag,
                           const FArrayBox& a_phi,
                           const Interval&  a_cPhi)
{
  CH_assert(a_cMag >= 0 && a_cMag < a_mag.nComp());
  CH_assert(a_phi.interval().contains(a_cPhi));
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_mag.box().contains(a_box));

  a_mag.setVal(0.0, a_cMag);

  // sum the square of each element
  for (int c = a_cPhi.begin(), c_end = a_cPhi.end() + 1; c != c_end; ++c)
    {
      MD_BOXLOOP(a_box, i)
        {
          a_mag[MD_IX(i, a_cMag)] += a_phi[MD_IX(i, c)]*a_phi[MD_IX(i, c)];
        }
    }
  // square root of the sum of squares
  MD_BOXLOOP(a_box, i)
    {
      a_mag[MD_IX(i, a_cMag)] = std::sqrt(a_mag[MD_IX(i, a_cMag)]);
    }
}

/*--------------------------------------------------------------------*/
/// divide a vector by an elemnt, happens in place
/// Be careful for division by zero!!
/** \param[in]  a_box   Box to operate over
 *  \param[out] a_phi   vector of phi divided by specified component
 *  \param[in]  a_phi   Field
 *  \param[in]  a_intv  Interval of components in a_phi to devide
 *  \param[in]  a_denominator
 *                      Field that contains element to divide by
 *  \param[in]  a_denomComp
 *                      Component to divide by
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divideVec(const Box&           a_box,
                           FArrayBox&           a_phi,
                           const Interval&      a_intvPhi,
                           const FArrayBox&     a_denominator,
                           const int            a_denomComp)
{
  // FIXME: Inconsistent treatment of Interval between
  //        this function and PatchMappedFunc::magnitude (?)
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_phi.interval().contains(a_intvPhi));
  
  MD_ARRAY_RESTRICT(phi, a_phi);
  MD_ARRAY_RESTRICT(denom, a_denominator);
  for (int comp = a_intvPhi.begin(), comp_end = a_intvPhi.end() + 1;
       comp != comp_end; ++comp)
    {
      MD_BOXLOOP(a_box, i)
        {
          phi[MD_IX(i, comp)] /= denom[MD_IX(i, a_denomComp)];
        }
    }
}

/*--------------------------------------------------------------------*/
/// divide a vector by an elemnt
/// Be careful for division by zero!!
/** \param[in]  a_box   Box to operate over
 *  \param[out] a_phiDiv
 *                      vector of phi divided by specified component
 *  \param[in]  a_phi   Field
 *  \param[in]  a_intv  Interval of components in a_phi to divide
 *  \param[in]  a_denominator
 *                      Field that contains element to divide by
 *  \param[in]  a_denomComp
 *                      Component to divide by
 *  \param[in]  a_intvDest
 *                      Interval to store in 'a_phiDiv'
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divideVec(const Box&           a_box,
                           FArrayBox&           a_phiDiv,
                           const FArrayBox&     a_phi,
                           const Interval&      a_intvPhi,
                           const FArrayBox&     a_denominator,
                           const int            a_denomComp,
                           const Interval&      a_intvDest)
{
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_phiDiv.box().contains(a_box));
  CH_assert(a_phi.interval().contains(a_intvPhi));
  CH_assert(a_phiDiv.interval().contains(a_intvDest));
  CH_assert(a_intvPhi.size() == a_intvDest.size());
  
  for (int c = a_intvPhi.begin(), c_end = a_intvPhi.end() + 1,
         n = a_intvDest.begin(); c != c_end; ++c, ++n)
    {
      MD_BOXLOOP(a_box, i)
        {
          a_phiDiv[MD_IX(i, n)] = a_phi[MD_IX(i, c)]/
            a_denominator[MD_IX(i, a_denomComp)];
        }
    }
}

/*--------------------------------------------------------------------*/
/// divide a vector by an elemnt
/// Be careful for division by zero!!
/** \param[in]  a_box   Box to operate over
 *  \param[out] a_phiDiv
 *                      vector of phi divided by specified component
 *  \param[in]  a_phi   Field
 *  \param[in]  a_intv  Interval of components in a_phi to divide
 *  \param[in]  a_denomComp
 *                      Component to divide by
 *  \param[in]  a_intvDest
 *                      Interval to store in 'a_phiDiv'
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divideVec(const Box&           a_box,
                           FArrayBox&           a_phiDiv,
                           const FArrayBox&     a_phi,
                           const Interval&      a_intvPhi,
                           const int            a_denomComp,
                           const Interval&      a_intvDest)
{
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_phiDiv.box().contains(a_box));
  CH_assert(a_phi.interval().contains(a_intvPhi));
  CH_assert(a_phiDiv.interval().contains(a_intvDest));
  CH_assert(a_intvPhi.size() == a_intvDest.size());
  
  MD_ARRAY_RESTRICT(phiDiv, a_phiDiv);
  MD_ARRAY_RESTRICT(phi, a_phi);
  for (int c = a_intvPhi.begin(), n = a_intvDest.begin();
       c != a_intvPhi.end()+1;
       ++c, ++n)
    {
      MD_BOXLOOP(a_box, i)
        {
          phiDiv[MD_IX(i, n)] = phi[MD_IX(i, c)]/phi[MD_IX(i, a_denomComp)];
        }
    }
}

/*--------------------------------------------------------------------*/
/// divide by an elemnt, for getting velocity from momentum
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_phi   output normalized interval
 *  \param[in]  a_phi   field containing vector to scale
 *  \param[in]  a_intv  Interval of components in a_phi to normalize
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::normalize(const Box&           a_box,
                           FArrayBox&           a_phi,
                           const Interval&      a_intv)
{
  CH_assert(a_phi.box().contains(a_box));
  CH_assert(a_phi.interval().contains(a_intv));
  
  // get the magnitude of the vector
  FABSTACKTEMP(magFab, a_box, 1);
  const int magComp = 0;
  if (a_intv.size() == SpaceDim)
    {
      // this version probably is a little faster if we can use it
      PatchMappedFunc::magnitudeSpace(a_box, magFab, magComp, a_phi, a_intv.begin());
    }
  else
    {
      PatchMappedFunc::magnitude(a_box, magFab, magComp, a_phi, a_intv); 
    }

  // divide each vector by the length
  divideVec(a_box, a_phi, a_intv, magFab, magComp);
}

/*--------------------------------------------------------------------*/
//  Obtain N^T in a cell with components that are contiguous in
//  memory.
/** Both N and N^T are assumed to have column-major order
 *  \param[in]  a_box   Cell box on which to find N^T
 *  \param[out] a_NtCtg N^T in the cells
 *  \param[in]  a_N     Metrics on the faces
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::cellNtCtg(const Box&                             a_box,
                           CHArray<Real, SpaceDim+1, ArRangeCol>& a_NtCtg,
                           const FluxBox&                         a_N)
{
  CH_assert(a_N.box().contains(a_box));
  
  MD_ARRAY_RESTRICT_CHARRAY(NtCtgArr, a_NtCtg);
  for (int xi = 0; xi != SpaceDim; ++xi)
    {
      const FArrayBox& NfaceDir = a_N[xi];
      MD_ARRAY_RESTRICT(NfaceDirArr, NfaceDir);
      const int MD_ID(ii, xi);
      MD_BOXLOOP(a_box, i)
        {
          for (int x = 0; x != SpaceDim; ++x)
            {
              NtCtgArr[MD_CIX(getTensorCompTranspose(x, xi), i)] =
                0.5*(NfaceDirArr[MD_IX(i, getTensorComp(x, xi))] +
                     NfaceDirArr[MD_OFFSETIX(i,+,ii, getTensorComp(x, xi))]);
              
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Obtain N^T in a cell with components that are contiguous in
//  memory.
/** Both N and N^T are assumed to have column-major order
 *  \param[in]  a_box   Cell box on which to find N^T
 *  \param[out] a_Ncell Cell averaged N
 *  \param[in]  a_N     Metrics on the faces
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::cellN(const Box&     a_box,
                       FArrayBox&     a_Ncell,
                       const FluxBox& a_N)
{
  CH_assert(a_N.box().contains(a_box));
  CH_assert(a_Ncell.box().contains(a_box));
  CH_assert(a_Ncell.nComp() >= m_dyadSize);
  
  MD_ARRAY_RESTRICT(NcellArr, a_Ncell);
  // loop over face directions
  for (int xi = 0; xi != SpaceDim; ++xi)
    {
      // get the face and direction
      const FArrayBox& NfaceDir = a_N[xi];
      MD_ARRAY_RESTRICT(NfaceDirArr, NfaceDir);
      const int MD_ID(ii, xi);
      MD_BOXLOOP(a_box, i)
        {
          // loop over columns of N
          for (int x = 0; x != SpaceDim; ++x)
            {
              // average face components of cell averaged value
              int compN = getTensorComp(x, xi);
              NcellArr[MD_IX(i, compN)] =
                0.5*(NfaceDirArr[MD_IX(i, compN)] +
                     NfaceDirArr[MD_OFFSETIX(i,+,ii, compN)]);
              
            }
        }
    }
}

/// get unit normals
void
PatchMappedFunc::getUnitNormals(FArrayBox& a_unitNormals,
                                const FArrayBox& a_faceMetrics,
                                const IntVect& a_normalComps,
                                const int& a_dir,
                                const Box& a_box)
{
  CH_assert(a_unitNormals.contains(a_box));
  CH_assert(a_unitNormals.nComp() == SpaceDim*SpaceDim);
  CH_assert(a_faceMetrics.contains(a_box));
  FORT_GETUNITNORMALS(CHF_FRA(a_unitNormals),
                      CHF_CONST_FRA(a_faceMetrics),
                      CHF_CONST_INTVECT(a_normalComps),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
// Transform a vector in computational space to one in physical space
/** \param[out] a_U     The vector data in physical space
 *  \param[in]  a_U     The vector data in computational space
 *  \param[in]  a_unitNormals
 *                      unit normal data
 *  \param[in]  a_box   Box to operate on
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::forwardTransform(FArrayBox& a_U,
                                  const FArrayBox& a_unitNormals,
                                  const Box& a_box)
{
  CH_assert(a_U.contains(a_box));
  CH_assert(a_U.nComp() >= SpaceDim);
  CH_assert(a_unitNormals.contains(a_box));
  CH_assert(a_unitNormals.nComp() == SpaceDim*SpaceDim);
  FORT_FORWARDTRANSFORMF(CHF_FRA(a_U),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
// Transform a vector in computational space to one in physical space
/** \param[out] a_Uout    The vector data in physical space
 *  \param[in]  a_outComp The start component for the vector data
 *  \param[in]  a_Uin     The vector data in computational space
 *  \param[in]  a_inComp  The start component for the vector data
 *  \param[in]  a_unitNormals
 *                      unit normal data
 *  \param[in]  a_box   Box to operate on
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::forwardTransform(FArrayBox& a_Uout,
                                  const int& a_outComp,
                                  const FArrayBox& a_Uin,
                                  const int& a_inComp,
                                  const FArrayBox& a_unitNormals,
                                  const Box& a_box)
{
  CH_assert(a_Uout.contains(a_box));
  CH_assert(a_Uout.nComp() >= SpaceDim);
  CH_assert(a_Uin.contains(a_box));
  CH_assert(a_Uin.nComp() >= SpaceDim);
  CH_assert(a_unitNormals.contains(a_box));
  CH_assert(a_unitNormals.nComp() == SpaceDim*SpaceDim);
  FORT_FORWARDTRANSFORMGENF(CHF_FRA(a_Uout),
                            CHF_CONST_INT(a_outComp),
                            CHF_CONST_FRA(a_Uin),
                            CHF_CONST_INT(a_inComp),
                            CHF_CONST_FRA(a_unitNormals),
                            CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Transform a vector in physical space to computational physical space
/** \param[out] a_U     The vector data in computational space
 *  \param[in]  a_U     The vector data in physical space
 *  \param[in]  a_unitNormals
 *                      unit normal data
 *  \param[in]  a_box   Box to operate on
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::reverseTransform(FArrayBox& a_U,
                                  const FArrayBox& a_unitNormals,
                                  const Box& a_box)
{
  CH_assert(a_U.contains(a_box));
  CH_assert(a_U.nComp() >= SpaceDim);
  CH_assert(a_unitNormals.contains(a_box));
  CH_assert(a_unitNormals.nComp() == SpaceDim*SpaceDim);
  FORT_REVERSETRANSFORMF(CHF_FRA(a_U),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Transform a vector in physical space to computational physical space
/** \param[out] a_Uout    The vector data in computational space
 *  \param[in]  a_outComp The start component for the vector data
 *  \param[in]  a_Uin     The vector data in physical space
 *  \param[in]  a_inComp  The start component for the vector data
 *  \param[in]  a_unitNormals
 *                      unit normal data
 *  \param[in]  a_box   Box to operate on
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::reverseTransform(FArrayBox& a_Uout,
                                  const int& a_outComp,
                                  const FArrayBox& a_Uin,
                                  const int& a_inComp,
                                  const FArrayBox& a_unitNormals,
                                  const Box& a_box)
{
  CH_assert(a_Uout.contains(a_box));
  CH_assert(a_Uout.nComp() >= SpaceDim);
  CH_assert(a_Uin.contains(a_box));
  CH_assert(a_Uin.nComp() >= SpaceDim);
  CH_assert(a_unitNormals.contains(a_box));
  CH_assert(a_unitNormals.nComp() == SpaceDim*SpaceDim);
  FORT_REVERSETRANSFORMGENF(CHF_FRA(a_Uout),
                            CHF_CONST_INT(a_outComp),
                            CHF_CONST_FRA(a_Uin),
                            CHF_CONST_INT(a_inComp),
                            CHF_CONST_FRA(a_unitNormals),
                            CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Compute gradients of phi in computational space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhi
 *                      Gradients of phi stored in components
 *                      [0, a_intv.size()*SpaceDim - 1] with order
 *                      (gradient direction, component)
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *  \param[in]  a_dxi   Computational mesh spacing
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradient2OCS(const Box&           a_box,
                              FArrayBox&           a_gradPhi,
                              const FArrayBox&     a_phi,
                              const ProblemDomain& a_problemDomain,
                              const Interval&      a_intv,
                              const RealVect&      a_dxi)
{
  const int phiNComp = a_phi.nComp();
  CH_assert(a_intv.begin() >= 0 && a_intv.begin() < phiNComp);
  CH_assert(a_intv.end() < phiNComp);
  CH_assert(a_gradPhi.nComp() >= a_intv.size()*SpaceDim);

  // MD Arrays for calculations in the boxes
  MD_ARRAY_RESTRICT(gradPhiArr, a_gradPhi);
  MD_ARRAY_RESTRICT(phiArr, a_phi);

  int cGradPhi = 0;
  // For each component in phi
  for (int cPhi = a_intv.begin(); cPhi <= a_intv.end(); ++cPhi)
    {

//--Get the gradient of phi in computational space

      for (int dir = 0; dir != SpaceDim; ++dir, ++cGradPhi)
        {                    
          Box box1Dir(a_box);
          box1Dir.grow(dir, 1);
          Box centerCells, entireCells, loCells, hiCells;
          int hasLoCells, hasHiCells;
          loHiCenter(loCells, hasLoCells,
                     hiCells, hasHiCells,
                     centerCells, entireCells,
                     box1Dir, a_problemDomain, dir);

          const int MD_ID(ii, dir);
          const Real factor = 1./a_dxi[dir];

          MD_BOXLOOP(centerCells, i)
            {
              gradPhiArr[MD_IX(i, cGradPhi)] =
                0.5*factor*(phiArr[MD_OFFSETIX(i,+,ii, cPhi)] -
                            phiArr[MD_OFFSETIX(i,-,ii, cPhi)]);
            }

          if (hasLoCells)
            {
              MD_BOXLOOP(loCells, i)
                {
                  gradPhiArr[MD_IX(i, cGradPhi)] =
                    -0.5*factor*(phiArr[MD_OFFSETIX(i,+,2*ii, cPhi)] 
                                 - 4.0*phiArr[MD_OFFSETIX(i,+,ii, cPhi)]
                                 + 3.0*phiArr[MD_IX(i, cPhi)]);
                }
            }

          if (hasHiCells)
            {
              MD_BOXLOOP(hiCells, i)
                {
                  gradPhiArr[MD_IX(i, cGradPhi)] =
                    0.5*factor*(phiArr[MD_OFFSETIX(i,-,2*ii, cPhi)] 
                                - 4.0*phiArr[MD_OFFSETIX(i,-,ii, cPhi)]
                                + 3.0*phiArr[MD_IX(i, cPhi)]);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute undivided gradients of phi in computational space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhi
 *                      Gradients of phi stored in components
 *                      [0, a_intv.size()*SpaceDim - 1] with order
 *                      (gradient direction, component)
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::undividedGrad2OCS(const Box&         a_box,
                                   FArrayBox&         a_gradPhi,
                                   const FArrayBox&   a_phi,
                                   const BlockDomain& a_problemDomain,
                                   const Interval&    a_intv)
{
  const int phiNComp = a_phi.nComp();
  CH_assert(a_intv.begin() >= 0 && a_intv.begin() < phiNComp);
  CH_assert(a_intv.end() < phiNComp);
  CH_assert(a_gradPhi.nComp() >= a_intv.size()*SpaceDim);

  // MD Arrays for calculations in the boxes
  MD_ARRAY_RESTRICT(gradPhiArr, a_gradPhi);
  MD_ARRAY_RESTRICT(phiArr, a_phi);

  int cGradPhi = 0;
  // For each component in phi
  for (int cPhi = a_intv.begin(); cPhi <= a_intv.end(); ++cPhi)
    {

//--Get the gradient of phi in computational space

      for (int dir = 0; dir != SpaceDim; ++dir, ++cGradPhi)
        {                    
          Box box1Dir(a_box);
          box1Dir.grow(dir, 1);
          Box centerCells, entireCells, loCells, hiCells;
          int hasLoCells, hasHiCells;
          loHiCenter(loCells, hasLoCells,
                     hiCells, hasHiCells,
                     centerCells, entireCells,
                     box1Dir, a_problemDomain, dir);

          const int MD_ID(ii, dir);

          MD_BOXLOOP(centerCells, i)
            {
              gradPhiArr[MD_IX(i, cGradPhi)] =
                0.5*(phiArr[MD_OFFSETIX(i,+,ii, cPhi)] -
                     phiArr[MD_OFFSETIX(i,-,ii, cPhi)]);
            }

          if (hasLoCells)
            {
              MD_BOXLOOP(loCells, i)
                {
                  gradPhiArr[MD_IX(i, cGradPhi)] =
                    -0.5*(phiArr[MD_OFFSETIX(i,+,2*ii, cPhi)] 
                          - 4.0*phiArr[MD_OFFSETIX(i,+,ii, cPhi)]
                          + 3.0*phiArr[MD_IX(i, cPhi)]);
                }
            }

          if (hasHiCells)
            {
              MD_BOXLOOP(hiCells, i)
                {
                  gradPhiArr[MD_IX(i, cGradPhi)] =
                    0.5*(phiArr[MD_OFFSETIX(i,-,2*ii, cPhi)] 
                         - 4.0*phiArr[MD_OFFSETIX(i,-,ii, cPhi)]
                         + 3.0*phiArr[MD_IX(i, cPhi)]);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Comp space face-avgd face-normal 1st-deriv from cell-avgd vals, codim-1
/** \param[out] a_faceAvgNormalDeriv
 *                      Fab of face-averaged face-normal first derivatives
 *                      computed using cell-averaged values
 *  \param[in]  a_cellAvgVal
 *                      Fab of cell-averaged values
 *  \param[in]  a_valComp
 *                      Component in a_cellAvgVal that
 *                      will be used to compute normal derivative
 *  \param[in]  a_derivComp
 *                      Component in a_faceAvgNormalDeriv
 *                      where normal derivative will be placed
 *  \param[in]  a_problemDomain
 *                      Computational space domain for this block
 *  \param[in]  a_box   Face box over which normal derivative
 *                      will be computed
 *  \param[in]  a_dxi   Vector of computational mesh spacing
 *  \param[in]  a_dir   Face normal direction
 *  \param[in]  a_fourthOrder
 *                      If true, derivative is computed to 4th-order
 *                      accuracy, if false, to 2nd-order accuracy
 *  \param[in]  a_interior
 *                      If true, it is assumed centered operations
 *                      can be used everywhere and no bounds will be
 *                      checked, if false, one-sided stencils will
 *                      be used where necessary
 *  \param[in]  a_unitStrideDim
 *                      This may be important later
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(
  FArrayBox&         a_faceAvgNormalDeriv,
  const FArrayBox&   a_cellAvgVal,
  const int          a_valComp,
  const int          a_derivComp,
  const BlockDomain& a_problemDomain,
  const Box&         a_box,
  const RealVect&    a_dxi,
  const int          a_dir,
  const bool         a_fourthOrder,
  const bool         a_interior,
  const bool         a_unitStrideDim)
{
  CH_TIME("PatchMappedFunc::faceAvgNormDerivFromCellAvgCS");
  CH_assert(a_cellAvgVal.interval().contains(a_valComp));
  CH_assert(a_faceAvgNormalDeriv.interval().contains(a_derivComp));
  CH_assert(a_dir >= 0 && a_dir <= SpaceDim);
  if (!a_interior) // if using boundary-aware interp, better be in domain!!!
    {
      CH_assert(a_problemDomain.contains(a_box));
    }
  CH_assert(a_faceAvgNormalDeriv.box().contains(a_box));

  Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox, entireBox;
  int hasLo, hasHi, hasLowest, hasHighest;
  if (a_fourthOrder)
    {
      hasLowest = 0;
      hasHighest = 0;
      // create the boxes necessary for the fourth-order scheme and test
      // them thoroughly (try to catch as many issues as possible)
      Box inBox = a_box;
      inBox.growHi(a_dir, 1); // prepare to turn into a cell-box
      inBox.enclosedCells(a_dir); // turn into a cell-box
      inBox.growLo(a_dir, 1); // now this is equivalent to a second-order box
      inBox.grow(a_dir, 1); // now this is equivalent to a fourth-order box
      loHiCenterFace4(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                      centerBox, innerCenterBox, entireBox, inBox,
                      a_problemDomain, a_dir); //**FIXME: check if correct
      innerCenterBox &= a_box; // make sure this doesn't extend past a_box
      Box innerCenterBoxTest = innerCenterBox;
      innerCenterBoxTest.enclosedCells(a_dir); // turn this into a cell-box
      innerCenterBoxTest.grow(a_dir, 2); // check if cntrd stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(innerCenterBoxTest));
      if (hasLo)
        {
          // check to see if loBox is needed; loHiCenterFace4 does not check
          // this for us and will always give both loBox and nextLoBox
          if (a_box.contains(loBox))
            {
              hasLowest = 1; // this a_box needs loBox
              Box loBoxTest = loBox;
              loBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
              loBoxTest.enclosedCells(a_dir); // turn this into a cell-box
              loBoxTest.growHi(a_dir, 3); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(loBoxTest));
            }
          Box nextLoBoxTest = nextLoBox;
          nextLoBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
          nextLoBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          nextLoBoxTest.growLo(a_dir, 1); // extend to lo-boundary
          nextLoBoxTest.growHi(a_dir, 2); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextLoBoxTest));
        }
      if (hasHi)
        {
          // check to see if hiBox is needed; loHiCenterFace4 does not check
          // this for us and will always give both hiBox and nextHiBox
          if (a_box.contains(hiBox))
            {
              hasHighest = 1; // this a_box needs hiBox
              Box hiBoxTest = hiBox;
              hiBoxTest.growLo(a_dir, 1); // prepare to change it to a cell-box
              hiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
              hiBoxTest.growLo(a_dir, 3); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
            }
          Box nextHiBoxTest = nextHiBox;
          nextHiBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
          nextHiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          nextHiBoxTest.growLo(a_dir, 3); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextHiBoxTest));
        }
    }
  else
    {
      Box testBox = a_box;
      bool testBoxIsCellBox = testBox.ixType().test(a_dir);
      if (testBoxIsCellBox) // test if it's a face-box in direction a_dir
        {
          testBox.growHi(a_dir, 1);
          testBox.enclosedCells(a_dir);
        }
      // create the boxes necessary for the second-order scheme
      Box inBox = testBox;
      inBox.grow(a_dir, 1);
      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox,
                     entireBox, inBox, a_problemDomain, a_dir); //**FIXME
      centerBox &= a_box; //**FIXME: we really shouldn't need this
      if (hasLo && !(a_box.contains(loBox)))
        {
          hasLo = 0; //**FIXME: we really shouldn't need this
        }
      if (hasHi && !(a_box.contains(hiBox)))
        {
          hasHi = 0; //**FIXME: we really shouldn't need this
        }
      if (hasLo)
        {
          Box loBoxTest = loBox;
          loBoxTest.growHi(a_dir, 1); // prepare to turn into a cell-box
          loBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          loBoxTest.growHi(a_dir, 1); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(loBoxTest));
        }
      if (hasHi)
        {
          Box hiBoxTest = hiBox;
          hiBoxTest.growLo(a_dir, 1); // prepare to turn this into a cell-box
          hiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          hiBoxTest.growLo(a_dir, 1); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
        }
      Box centerBoxTest = centerBox;
      centerBoxTest.enclosedCells(a_dir); // turn this into a cell-box
      centerBoxTest.grow(a_dir, 1); // check if centered stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(centerBoxTest));
    }

  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  if (a_interior && a_fourthOrder) // fourth-order interior scheme
    {
      const Real c_0 = 1./(12.*a_dxi[a_dir]);
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
            c_0*(u_jmm - 15.*u_jm + 15.*u_j - u_jp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
    }
  else if (a_interior) // second-order interior scheme
    {
      const Real c_0 = 1./(a_dxi[a_dir]);
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_j  = a_cellAvgVal[MD_IX(i, a_valComp)];
          a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] = c_0*(u_j - u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
    }
  else if (a_fourthOrder) // fourth-order boundary scheme
    {
      const Real c_0 = 1./(12.*a_dxi[a_dir]);
      MD_BOXLOOP(innerCenterBox, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
            c_0*(u_jmm - 15.*u_jm + 15.*u_j - u_jp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
      if (hasLowest)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_valComp)];
              const Real u_j3p = a_cellAvgVal[MD_OFFSETIX(i,+,3*o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
                c_0*(-35.*u_j + 69.*u_jp - 45.*u_jpp + 11.*u_j3p);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(u_j3p < hiTol && u_j3p > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasLo)
        {
          MD_BOXLOOP(nextLoBox, i)
            {
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
                c_0*(-11.*u_jm + 9.*u_j + 3.*u_jp - u_jpp);

              // Check to make sure values are valid
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHighest)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_valComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_valComp)];
              const Real u_j4m = a_cellAvgVal[MD_OFFSETIX(i,-,4*o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
                c_0*(35.*u_jm - 69.*u_jmm + 45.*u_j3m - 11.*u_j4m);

              // Check to make sure values are valid
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(u_j4m < hiTol && u_j4m > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(nextHiBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_valComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] =
                c_0*(11.*u_j - 9.*u_jm - 3.*u_jmm + u_j3m);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
    }
  else // second-order boundary scheme
    {
      const Real c_0 = 1./(a_dxi[a_dir]);
      MD_BOXLOOP(centerBox, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_j  = a_cellAvgVal[MD_IX(i, a_valComp)];
          a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] = c_0*(u_j - u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j  = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] = c_0*(-u_j + u_jp);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_jmm =
                a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
              a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] = c_0*(-u_jmm + u_jm);

              // Check to make sure values are valid
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_faceAvgNormalDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Comp space cell-avgd 1st-deriv from cell-avgd vals
/** \param[out] a_cellAvgDeriv
 *                      Fab of cell-averaged first derivatives
 *                      computed using cell-averaged values
 *  \param[in]  a_cellAvgVal
 *                      Fab of cell-averaged values
 *  \param[in]  a_valComp
 *                      Component in a_cellAvgVal that will be used
 *                      to compute normal derivative
 *  \param[in]  a_derivComp
 *                      Component in a_cellAvgDeriv where derivative
 *                      will be placed
 *  \param[in]  a_problemDomain
 *                      Computational space domain for this block
 *  \param[in]  a_box   Cell box over which derivative will be computed
 *  \param[in]  a_dxi   Vector of computational mesh spacing
 *  \param[in]  a_dir   Derivative direction
 *  \param[in]  a_fourthOrder
 *                      If true, derivative is computed to 4th-order
 *                      accuracy, if false, to 2nd-order accuracy
 *  \param[in]  a_interior
 *                      If true, it is assumed centered operations
 *                      can be used everywhere and no bounds will be
 *                      checked, if false, one-sided stencils will
 *                      be used where necessary
 *  \param[in]  a_unitStrideDim
 *                      This may be important later
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::cellAvgDerivFromCellAvgCS(
  FArrayBox&         a_cellAvgDeriv,
  const FArrayBox&   a_cellAvgVal,
  const int          a_valComp,
  const int          a_derivComp,
  const BlockDomain& a_problemDomain,
  const Box&         a_box,
  const RealVect&    a_dxi,
  const int          a_dir,
  const bool         a_fourthOrder,
  const bool         a_interior,
  const bool         a_unitStrideDim)
{
  CH_TIME("PatchMappedFunc::cellAvgDerivFromCellAvgCS");
  CH_assert(a_cellAvgVal.interval().contains(a_valComp));
  CH_assert(a_cellAvgDeriv.interval().contains(a_derivComp));
  CH_assert(a_dir >= 0 && a_dir <= SpaceDim);
  if (!a_interior) // if using boundary-aware interp, better be in domain!!!
    {
      CH_assert(a_problemDomain.contains(a_box));
    }
  CH_assert(a_cellAvgDeriv.box().contains(a_box));

  Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox, entireBox;
  int hasLo, hasHi, hasLowest, hasHighest;
  if (a_fourthOrder)
    {
      hasLowest = 0;
      hasHighest = 0;
      // create the boxes necessary for the fourth-order scheme and test
      // them thoroughly (try to catch as many issues as possible)
      Box inBox = a_box;
      inBox.grow(a_dir, 1);
      loHiCenter5(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                  centerBox, innerCenterBox, entireBox, inBox,
                  a_problemDomain, a_dir);
      Box innerCenterBoxTest = innerCenterBox;
      innerCenterBoxTest.grow(a_dir, 2); // check if cntrd stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(innerCenterBoxTest));
      if (hasLo)
        {
          // check to see if loBox is needed; loHiCenter5 does not check
          // this for us and will always give both loBox and nextLoBox
          if (a_box.contains(loBox))
            {
              hasLowest = 1; // this a_box needs loBox
              Box loBoxTest = loBox;
              loBoxTest.growHi(a_dir, 4); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(loBoxTest));
            }
          Box nextLoBoxTest = nextLoBox;
          nextLoBoxTest.growLo(a_dir, 1); // extend to lo-boundary
          nextLoBoxTest.growHi(a_dir, 3); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextLoBoxTest));
        }
      if (hasHi)
        {
          // check to see if hiBox is needed; loHiCenterFace4 does not check
          // this for us and will always give both hiBox and nextHiBox
          if (a_box.contains(hiBox))
            {
              hasHighest = 1; // this a_box needs hiBox
              Box hiBoxTest = hiBox;
              hiBoxTest.growLo(a_dir, 4); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
            }
          Box nextHiBoxTest = nextHiBox;
          nextHiBoxTest.growHi(a_dir, 1); // extend to hi-boundary
          nextHiBoxTest.growLo(a_dir, 3); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextHiBoxTest));
        }
    }
  else
    {
      // create the boxes necessary for the second-order scheme
      Box inBox = a_box;
      inBox.grow(a_dir, 1);
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox,
                 entireBox, inBox, a_problemDomain, a_dir);
      if (hasLo)
        {
          Box loBoxTest = loBox;
          loBoxTest.growHi(a_dir, 2); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(loBoxTest));
        }
      if (hasHi)
        {
          Box hiBoxTest = hiBox;
          hiBoxTest.growLo(a_dir, 2); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
        }
      Box centerBoxTest = centerBox;
      centerBoxTest.grow(a_dir, 1); // check if centered stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(centerBoxTest));
    }

  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  if (a_interior && a_fourthOrder) // fourth-order interior scheme
    {
      const Real c_0 = 1./(12.*a_dxi[a_dir]);
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i, +, 2*o, a_valComp)];
          a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
            c_0*(u_jmm - 8.*u_jm + 8.*u_jp - u_jpp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(u_jpp < hiTol && u_jpp > loTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
    }
  else if (a_interior) // second-order interior scheme
    {
      const Real c_0 = 1./(2.*a_dxi[a_dir]);
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_jp = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          a_cellAvgDeriv[MD_IX(i, a_derivComp)] = c_0*(u_jp - u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
    }
  else if (a_fourthOrder) // fourth-order boundary scheme
    {
      const Real c_0 = 1./(12.*a_dxi[a_dir]);
      MD_BOXLOOP(innerCenterBox, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i, +, 2*o, a_valComp)];
          a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
            c_0*(u_jmm - 8.*u_jm + 8.*u_jp - u_jpp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(u_jpp < hiTol && u_jpp > loTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
      if (hasLowest)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_valComp)];
              const Real u_j3p = a_cellAvgVal[MD_OFFSETIX(i,+,3*o, a_valComp)];
              const Real u_j4p = a_cellAvgVal[MD_OFFSETIX(i,+,4*o, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(-25.*u_j + 48.*u_jp - 36.*u_jpp + 16.*u_j3p - 3.*u_j4p);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(u_j3p < hiTol && u_j3p > loTol);
              CH_assert(u_j4p < hiTol && u_j4p > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasLo)
        {
          MD_BOXLOOP(nextLoBox, i)
            {
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_valComp)];
              const Real u_j3p = a_cellAvgVal[MD_OFFSETIX(i,+,3*o, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(-3.*u_jm - 10.*u_j + 18.*u_jp - 6.*u_jpp + u_j3p);

              // Check to make sure values are valid
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(u_j3p < hiTol && u_j3p > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHighest)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_valComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_valComp)];
              const Real u_j4m = a_cellAvgVal[MD_OFFSETIX(i,-,4*o, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(25.*u_j - 48.*u_jm + 36.*u_jmm - 16.*u_j3m + 3.*u_j4m);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(u_j4m < hiTol && u_j4m > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(nextHiBox, i)
            {
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_valComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_valComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(3.*u_jp + 10.*u_j - 18.*u_jm + 6.*u_jmm - u_j3m);

              // Check to make sure values are valid
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
    }
  else // second-order boundary scheme
    {
      const Real c_0 = 1./(2.*a_dxi[a_dir]);
      MD_BOXLOOP(centerBox, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
          const Real u_jp = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_valComp)];
          a_cellAvgDeriv[MD_IX(i, a_derivComp)] = c_0*(u_jp - u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
          CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
        }
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_valComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(-3.*u_j + 4.*u_jp - u_jpp);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_jmm =
                a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_valComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_valComp)];
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_valComp)];
              a_cellAvgDeriv[MD_IX(i, a_derivComp)] =
                c_0*(3.*u_j - 4.*u_jm + u_jmm);

              // Check to make sure values are valid
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] < hiTol);
              CH_assert(a_cellAvgDeriv[MD_IX(i, a_derivComp)] > loTol);
            }
        }
    }
}


/*--------------------------------------------------------------------*/
//  Comp space face-avgd value from cell-avgd vals, codim-1 version
/** \param[out] a_faceAvgVal
 *                      Fab of face-averaged values
 *                      computed using cell-averaged values
 *  \param[in]  a_cellAvgVal
 *                      Fab of cell-averaged values
 *  \param[in]  a_cellComp
 *                      Component in a_cellAvgVal that
 *                      will be used to compute a_faceAvgVal
 *  \param[in]  a_faceComp
 *                      Component in a_faceAvgVal
 *                      where face-avgd value will be placed
 *  \param[in]  a_problemDomain
 *                      Computational space domain for this block
 *  \param[in]  a_box   Face box over which face values
 *                      will be computed
 *  \param[in]  a_dir   Face normal direction
 *  \param[in]  a_order Order of accuracy of face interpolation
 *  \param[in]  a_interior
 *                      If true, it is assumed centered operations
 *                      can be used everywhere and no bounds will be
 *                      checked, if false, one-sided stencils will
 *                      be used where necessary
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::faceAvgValFromCellAvgCS(FArrayBox&         a_faceAvgVal,
                                         const FArrayBox&   a_cellAvgVal,
                                         const int          a_cellComp,
                                         const int          a_faceComp,
                                         const BlockDomain& a_problemDomain,
                                         const Box&         a_box,
                                         const int          a_dir,
                                         const int          a_order,
                                         const bool         a_interior)
{
  CH_TIME("PatchMappedFunc::faceAvgValFromCellAvgCS");
  CH_assert(a_cellAvgVal.interval().contains(a_cellComp));
  CH_assert(a_faceAvgVal.interval().contains(a_faceComp));
  CH_assert(a_dir >= 0 && a_dir <= SpaceDim);
  if (!a_interior) // if using boundary-aware interp, better be in domain!!!
    {
      CH_assert(a_problemDomain.contains(a_box));
    }
  CH_assert(a_faceAvgVal.box().contains(a_box));

  //**FIXME: add 1st, 3rd, and 5th-order interpolations

  Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox, entireBox;
  int hasLo, hasHi, hasLowest, hasHighest;
  if (a_order == 4) // 4th-order
    {
      hasLowest = 0;
      hasHighest = 0;
      // create the boxes necessary for the fourth-order scheme and test
      // them thoroughly (try to catch as many issues as possible)
      Box inBox = a_box;
      inBox.growHi(a_dir, 1); // prepare to turn into a cell-box
      inBox.enclosedCells(a_dir); // turn into a cell-box
      inBox.growLo(a_dir, 1); // now this is equivalent to a second-order box
      inBox.grow(a_dir, 1); // now this is equivalent to a fourth-order box
      loHiCenterFace4(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                      centerBox, innerCenterBox, entireBox, inBox,
                      a_problemDomain, a_dir); //**FIXME: check if correct
      innerCenterBox &= a_box; // make sure this doesn't extend past a_box
      Box innerCenterBoxTest = innerCenterBox;
      innerCenterBoxTest.enclosedCells(a_dir); // turn this into a cell-box
      innerCenterBoxTest.grow(a_dir, 2); // check if cntrd stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(innerCenterBoxTest));
      if (hasLo)
        {
          // check to see if loBox is needed; loHiCenterFace4 does not check
          // this for us and will always give both loBox and nextLoBox
          if (a_box.contains(loBox))
            {
              hasLowest = 1; // this a_box needs loBox
              Box loBoxTest = loBox;
              loBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
              loBoxTest.enclosedCells(a_dir); // turn this into a cell-box
              loBoxTest.growHi(a_dir, 3); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(loBoxTest));
            }
          Box nextLoBoxTest = nextLoBox;
          nextLoBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
          nextLoBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          nextLoBoxTest.growLo(a_dir, 1); // extend to lo-boundary
          nextLoBoxTest.growHi(a_dir, 2); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextLoBoxTest));
        }
      if (hasHi)
        {
          // check to see if hiBox is needed; loHiCenterFace4 does not check
          // this for us and will always give both hiBox and nextHiBox
          if (a_box.contains(hiBox))
            {
              hasHighest = 1; // this a_box needs hiBox
              Box hiBoxTest = hiBox;
              hiBoxTest.growLo(a_dir, 1); // prepare to change it to a cell-box
              hiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
              hiBoxTest.growLo(a_dir, 3); // see if biased stencil fits in data
              CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
            }
          Box nextHiBoxTest = nextHiBox;
          nextHiBoxTest.growHi(a_dir, 1); // prepare to change it to a cell-box
          nextHiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          nextHiBoxTest.growLo(a_dir, 3); // see if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(nextHiBoxTest));
        }
    }
  else if (a_order == 2) // 2nd-order
    {
      Box testBox = a_box;
      bool testBoxIsCellBox = testBox.ixType().test(a_dir);
      if (testBoxIsCellBox) // test if it's a face-box in direction a_dir
        {
          testBox.growHi(a_dir, 1);
          testBox.enclosedCells(a_dir);
        }
      // create the boxes necessary for the second-order scheme
      Box inBox = testBox;
      inBox.grow(a_dir, 1);
      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox,
                     entireBox, inBox, a_problemDomain, a_dir); //**FIXME
      centerBox &= a_box; //**FIXME: we really shouldn't need this
      if (hasLo && !(a_box.contains(loBox)))
        {
          hasLo = 0; //**FIXME: we really shouldn't need this
        }
      if (hasHi && !(a_box.contains(hiBox)))
        {
          hasHi = 0; //**FIXME: we really shouldn't need this
        }
      if (hasLo)
        {
          Box loBoxTest = loBox;
          loBoxTest.growHi(a_dir, 1); // prepare to turn into a cell-box
          loBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          loBoxTest.growHi(a_dir, 1); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(loBoxTest));
        }
      if (hasHi)
        {
          Box hiBoxTest = hiBox;
          hiBoxTest.growLo(a_dir, 1); // prepare to turn this into a cell-box
          hiBoxTest.enclosedCells(a_dir); // turn this into a cell-box
          hiBoxTest.growLo(a_dir, 1); // check if biased stencil fits in data
          CH_assert(a_cellAvgVal.box().contains(hiBoxTest));
        }
      Box centerBoxTest = centerBox;
      centerBoxTest.enclosedCells(a_dir); // turn this into a cell-box
      centerBoxTest.grow(a_dir, 1); // check if centered stencil fits in data
      CH_assert(a_cellAvgVal.box().contains(centerBoxTest));
    }

  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  if (a_interior && (a_order == 4)) // fourth-order interior scheme
    {
      const Real c_0 = 1./12.;
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_cellComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_cellComp)];
          const Real u_j   = a_cellAvgVal[MD_IX(i, a_cellComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_cellComp)];
          a_faceAvgVal[MD_IX(i, a_faceComp)] =
            c_0*(-u_jmm + 7.*u_jm + 7.*u_j - u_jp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
        }
    }
  else if (a_interior && (a_order == 2)) // second-order interior scheme
    {
      const Real c_0 = 0.5;
      MD_BOXLOOP(a_box, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_cellComp)];
          const Real u_j  = a_cellAvgVal[MD_IX(i, a_cellComp)];
          a_faceAvgVal[MD_IX(i, a_faceComp)] = c_0*(u_j + u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
        }
    }
  else if (!a_interior && (a_order == 4)) // fourth-order boundary scheme
    {
      const Real c_0 = 1./12.;
      MD_BOXLOOP(innerCenterBox, i)
        {
          const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i, -, 2*o, a_cellComp)];
          const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_cellComp)];
          const Real u_j   = a_cellAvgVal[MD_IX(i, a_cellComp)];
          const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_cellComp)];
          a_faceAvgVal[MD_IX(i, a_faceComp)] =
            c_0*(-u_jmm + 7.*u_jm + 7.*u_j - u_jp);

          // Check to make sure values are valid
          CH_assert(u_jmm < hiTol && u_jmm > loTol);
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(u_jp < hiTol && u_jp > loTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
        }
      if (hasLowest)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_cellComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_cellComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_cellComp)];
              const Real u_j3p = a_cellAvgVal[MD_OFFSETIX(i,+,3*o, a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] =
                c_0*(25.*u_j - 23.*u_jp + 13.*u_jpp - 3.*u_j3p);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(u_j3p < hiTol && u_j3p > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
      if (hasLo)
        {
          MD_BOXLOOP(nextLoBox, i)
            {
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_cellComp)];
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_cellComp)];
              const Real u_jp  = a_cellAvgVal[MD_OFFSETIX(i,+,o, a_cellComp)];
              const Real u_jpp = a_cellAvgVal[MD_OFFSETIX(i,+,2*o, a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] =
                c_0*(3.*u_jm + 13.*u_j - 5.*u_jp + u_jpp);

              // Check to make sure values are valid
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(u_jpp < hiTol && u_jpp > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
      if (hasHighest)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_cellComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_cellComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_cellComp)];
              const Real u_j4m = a_cellAvgVal[MD_OFFSETIX(i,-,4*o, a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] =
                c_0*(25.*u_jm - 23.*u_jmm + 13.*u_j3m - 3.*u_j4m);

              // Check to make sure values are valid
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(u_j4m < hiTol && u_j4m > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(nextHiBox, i)
            {
              const Real u_j   = a_cellAvgVal[MD_IX(i, a_cellComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o, a_cellComp)];
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o, a_cellComp)];
              const Real u_j3m = a_cellAvgVal[MD_OFFSETIX(i,-,3*o, a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] =
                c_0*(3.*u_j + 13.*u_jm - 5.*u_jmm + u_j3m);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_j3m < hiTol && u_j3m > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
    }
  else if (!a_interior && (a_order == 2)) // second-order boundary scheme
    {
      const Real c_0 = 0.5;
      MD_BOXLOOP(centerBox, i)
        {
          const Real u_jm = a_cellAvgVal[MD_OFFSETIX(i, -, o, a_cellComp)];
          const Real u_j  = a_cellAvgVal[MD_IX(i, a_cellComp)];
          a_faceAvgVal[MD_IX(i, a_faceComp)] = c_0*(u_j + u_jm);

          // Check to make sure values are valid
          CH_assert(u_jm < hiTol && u_jm > loTol);
          CH_assert(u_j < hiTol && u_j > loTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
          CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
        }
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real u_j  = a_cellAvgVal[MD_IX(i, a_cellComp)];
              const Real u_jp = a_cellAvgVal[MD_OFFSETIX(i, +, o, a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] = c_0*(3.*u_j - u_jp);

              // Check to make sure values are valid
              CH_assert(u_j < hiTol && u_j > loTol);
              CH_assert(u_jp < hiTol && u_jp > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real u_jmm = a_cellAvgVal[MD_OFFSETIX(i,-,2*o,a_cellComp)];
              const Real u_jm  = a_cellAvgVal[MD_OFFSETIX(i,-,o,a_cellComp)];
              a_faceAvgVal[MD_IX(i, a_faceComp)] = c_0*(3.*u_jm - u_jmm);

              // Check to make sure values are valid
              CH_assert(u_jmm < hiTol && u_jmm > loTol);
              CH_assert(u_jm < hiTol && u_jm > loTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] < hiTol);
              CH_assert(a_faceAvgVal[MD_IX(i, a_faceComp)] > loTol);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Comp space cell-avgd grad from face-avgd vals, exact calculation
/** \param[out] a_cellAvgGrad
 *                      Fab of cell-averaged derivatives in a_dir
 *                      computed using a_faceAvgVal
 *  \param[in]  a_faceAvgVal
 *                      FluxBox of face-averaged values
 *  \param[in]  a_inputIntv
 *                      Interval of components in a_faceAvgVal for
 *                      which to compute gradients
 *  \param[in]  a_outputStartComp
 *                      Component in a_cellAvgGrad at which to start
 *                      writing the gradients
 *  \param[in]  a_box   Cell box over which cell-avgd derivatives
 *                      will be computed
 *  \param[in]  a_dxi   Vector of computational mesh spacing
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::cellAvgGradFromFaceAvgCS(FArrayBox&      a_cellAvgGrad,
                                          const FluxBox&  a_faceAvgVal,
                                          const Interval& a_inputIntv,
                                          const int       a_outputStartComp,
                                          const Box&      a_box,
                                          const RealVect& a_dxi)
{
  CH_TIME("PatchMappedFunc::cellAvgGradFromFaceAvgCS");
  const int inBeginComp = a_inputIntv.begin();
  const int numComps = a_inputIntv.size();
  CH_assert(a_cellAvgGrad.interval().contains(a_outputStartComp));
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      CH_assert(a_faceAvgVal[dir].interval().contains(a_inputIntv));
      for (int comp = 0; comp != numComps; ++comp)
        {
          int inputComp = inBeginComp + comp;
          int derivComp = a_outputStartComp + dir + comp*SpaceDim;
          PatchMappedFunc::cellAvgDerivFromFaceAvgCS(a_cellAvgGrad,
                                                     a_faceAvgVal[dir],
                                                     inputComp,
                                                     derivComp,
                                                     a_box,
                                                     a_dxi,
                                                     dir);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Comp space cell-avgd derivative from face-avgd vals, exact calculation
/** \param[out] a_cellAvgGrad
 *                      Fab of cell-averaged derivatives in a_dir
 *                      computed using a_faceAvgVal
 *  \param[in]  a_faceAvgVal
 *                      Fab of face-averaged values
 *  \param[in]  a_inputComp
 *                      Component in a_faceAvgVal that
 *                      will be used to compute a_cellAvgGrad
 *  \param[in]  a_outputComp
 *                      Component in a_cellAvgGrad
 *                      where cell-avgd derivative will be placed
 *  \param[in]  a_box   Cell box over which cell-avgd derivatives
 *                      will be computed
 *  \param[in]  a_dxi   Vector of computational mesh spacing
 *  \param[in]  a_dir   Face normal direction
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::cellAvgDerivFromFaceAvgCS(FArrayBox&       a_cellAvgGrad,
                                           const FArrayBox& a_faceAvgVal,
                                           const int&       a_inputComp,
                                           const int&       a_outputComp,
                                           const Box&       a_box,
                                           const RealVect&  a_dxi,
                                           const int        a_dir)
{
  CH_TIME("PatchMappedFunc::cellAvgDerivFromFaceAvgCS");
  CH_assert(a_faceAvgVal.interval().contains(a_inputComp));
  CH_assert(a_cellAvgGrad.interval().contains(a_outputComp));
  CH_assert(a_dir >= 0 && a_dir <= SpaceDim);
  CH_assert(a_cellAvgGrad.box().contains(a_box));
  Box faceAvgValBox = a_faceAvgVal.box();
  bool faceAvgValBoxIsFaceBox = faceAvgValBox.ixType().test(a_dir);
  if (!faceAvgValBoxIsFaceBox) // test if it's a cell-box in direction a_dir
    {
      faceAvgValBox.growHi(a_dir, -1); // if it is, just shrink hi-side
    }
  else // if it's a face-box, just convert it to the correct cell-box
    {
      faceAvgValBox.enclosedCells(a_dir);
    }
  CH_assert(faceAvgValBox.contains(a_box));
  const Real c_0 = 1./a_dxi[a_dir];
  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  MD_BOXLOOP(a_box, i)
    {
      const Real u_j  = a_faceAvgVal[MD_IX(i, a_inputComp)];
      const Real u_jp = a_faceAvgVal[MD_OFFSETIX(i, +, o, a_inputComp)];
      a_cellAvgGrad[MD_IX(i, a_outputComp)] = c_0*(u_jp - u_j);

      // check that valid values are used when computing a_cellAvgGrad
      CH_assert(u_j < hiTol && u_j > loTol);
      CH_assert(u_jp < hiTol && u_jp > loTol);
      CH_assert(a_cellAvgGrad[MD_IX(i, a_outputComp)] < hiTol);
      CH_assert(a_cellAvgGrad[MD_IX(i, a_outputComp)] > loTol);
    }
}

/*--------------------------------------------------------------------*/
/// Convert gradients of phi in computational space to physical space
///   a_gradPhiCS, and a_NtJ must be of same type (ie. point, face avg, ...)
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhiPS
 *                      Gradients of phi in physical space with
 *                      derivatives of each component store contiguous
 *                      ie. [d_u/d_\xi, d_u/d_\eta, d_v/d_\xi, ...]
 *  \param[in]  a_gradPhiCS
 *                      Gradients of phi in computational space, with
 *                      same storage as a_gradPhiPS.
 *  \param[in]  a_NtJ   Grid metrics
 *  \param[in]  a_intvPS
 *                      Interval of computational gradient to transform
 *                      By default this is the entire FAB
 *  \param[in]  a_intvCS
 *                      Interval to write physical space gradient to
 *                      By default this is the entire FAB
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradientCStoPS(const Box&       a_box,
                                FArrayBox&       a_gradPhiPS,
                                const FArrayBox& a_gradPhiCS,
                                const FArrayBox& a_NtJ,
                                const Interval&  a_intvPS,
                                const Interval&  a_intvCS)
{
  const int numCompPhi = a_intvCS.size()/SpaceDim;
  CH_assert(a_intvCS.size() == a_intvPS.size());
  CH_assert(a_gradPhiPS.interval().contains(a_intvPS));
  CH_assert(a_gradPhiCS.interval().contains(a_intvCS));
  CH_assert(a_gradPhiPS.box().contains(a_box));
  CH_assert(a_gradPhiCS.box().contains(a_box));
  CH_assert(a_NtJ.box().contains(a_box));

  // Arrays are stored in row-major order
  // gradient arrays are arranged in the order of phi components
  //   ex) (d phi_0/d xi_0), (d phi_0/d xi_1), (d phi_1/d xi_0), (d phi_1/d xi_1)
  // transformation arrays are arranged in the order of x components
  //   ex) (d xi_0/d x_0), (d xi_0/d x_1), (d xi_1/d x_0), (d xi_1/d x_1)
  
  // MD Arrays for calculations in the boxes
  MD_ARRAY_RESTRICT(gradPhiPSArr, a_gradPhiPS);
  MD_ARRAY_RESTRICT(gradPhiCSArr, a_gradPhiCS);
  MD_ARRAY_RESTRICT(NtJArr, a_NtJ);
  // loop over the components of phi
  for (int cPhi = 0; cPhi != numCompPhi; ++cPhi)
    {
      // loop over the x-y-z direction for gradients in physical space 
      for (int xDir = 0; xDir != SpaceDim; ++xDir)
        {
          const int cGradPhi = a_intvPS.begin() + getGradientComp(cPhi, xDir);
          MD_BOXLOOP(a_box, i)
            {
              // unrolled row of gradPhiCS times column of Nt/J
              gradPhiPSArr[MD_IX(i, cGradPhi)] = (
                D_TERM6(  NtJArr[MD_IX(i, getTensorComp(xDir,0))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,0))],
                        + NtJArr[MD_IX(i, getTensorComp(xDir,1))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,1))],
                        + NtJArr[MD_IX(i, getTensorComp(xDir,2))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,2))],
                        + NtJArr[MD_IX(i, getTensorComp(xDir,3))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,3))],
                        + NtJArr[MD_IX(i, getTensorComp(xDir,4))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,4))],
                        + NtJArr[MD_IX(i, getTensorComp(xDir,5))]
                        * gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,5))]));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
/// Convert gradients of phi in computational space to physical space
///   a_gradPhiCS, and a_NtJ must be of same type (ie. point, face avg, ...)
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhiPS
 *                      Gradients of phi in physical space with
 *                      derivatives of each component store contiguous
 *                      ie. [d_u/d_\xi, d_u/d_\eta, d_v/d_\xi, ...]
 *  \param[in]  a_gradPhiCS
 *                      Gradients of phi in computational space, with
 *                      same storage as a_gradPhiPS.
 *  \param[in]  a_NtJ   Grid metrics
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradientCStoPS(const Box&       a_box,
                                FArrayBox&       a_gradPhiPS,
                                const FArrayBox& a_gradPhiCS,
                                const FArrayBox& a_NtJ)
{
  PatchMappedFunc::gradientCStoPS(a_box,
                                  a_gradPhiPS,
                                  a_gradPhiCS,
                                  a_NtJ,
                                  a_gradPhiPS.interval(),
                                  a_gradPhiCS.interval());
}

/*--------------------------------------------------------------------*/
/// Convert gradients of phi in computational space to physical space
///   a_gradPhiCS, and a_N, and a_J must be of same type (ie. point, face avg, ...)
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhiPS
 *                      Gradients of phi in physical space with
 *                      derivatives of each component store contiguous
 *                      ie. [d_u/d_\xi, d_u/d_\eta, d_v/d_\xi, ...]
 *  \param[in]  a_gradPhiCS
 *                      Gradients of phi in computational space, with
 *                      same storage as a_gradPhiPS.
 *  \param[in]  a_Nt    Grid transformation matrix
 *  \param[in]  a_J     Grid metric jabobian
 *  \param[in]  a_intvPS
 *                      Interval of computational gradient to transform
 *                      By default this is the entire FAB
 *  \param[in]  a_intvCS
 *                      Interval to write physical space gradient to
 *                      By default this is the entire FAB
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradientCStoPS(const Box&       a_box,
                                FArrayBox&       a_gradPhiPS,
                                const FArrayBox& a_gradPhiCS,
                                const FArrayBox& a_Nt,
                                const FArrayBox& a_J,
                                const Interval&  a_intvPS,
                                const Interval&  a_intvCS)
{
  // since J is a scalar, pull it out of transformation
  const Interval gradPhiPSintv(0, m_dyadSize-1);
  PatchMappedFunc::gradientCStoPS(a_box,
                                  a_gradPhiPS,
                                  a_gradPhiCS,
                                  a_Nt,
                                  a_intvPS,
                                  a_intvCS);
  const int Jcomp = 0;
  PatchMappedFunc::divideVec(a_box,
                             a_gradPhiPS,
                             a_intvPS,
                             a_J,
                             Jcomp);
}

/*--------------------------------------------------------------------*/
/// Convert gradients of phi in computational space to physical space
///   a_gradPhiCS, and a_NtJ must be of same type (ie. point, face avg, ...)
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhiPS
 *                      Gradients of phi in physical space with
 *                      derivatives of each component store contiguous
 *                      ie. [d_u/d_\xi, d_u/d_\eta, d_v/d_\xi, ...]
 *  \param[in]  a_gradPhiCS
 *                      Gradients of phi in computational space, with
 *                      same storage as a_gradPhiPS.
 *  \param[in]  a_NtJ   Grid metrics
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradientCStoPS(const Box&       a_box,
                                FArrayBox&       a_gradPhiPS,
                                const FArrayBox& a_gradPhiCS,
                                const FArrayBox& a_N,
                                const FArrayBox& a_J)
{
  PatchMappedFunc::gradientCStoPS(a_box,
                                  a_gradPhiPS,
                                  a_gradPhiCS,
                                  a_N,
                                  a_J,
                                  a_gradPhiPS.interval(),
                                  a_gradPhiCS.interval());
}

/*--------------------------------------------------------------------*/
//  Convert gradients of phi in computational space to physical space
//    This for cell averaged values only!
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhiPS
 *                      Gradients of phi in physical space, stored in
 *                      components [0, a_intv.size()*SpaceDim - 1]
 *                      with order (gradient direction, component)
 *  \param[in] a_gradPhiCS
 *                      Gradients of phi in computational space, with
 *                      same storage as a_gradPhiPS.
 *  \param[in]  a_N     Metrics on the faces
 *  \param[in]  a_J     Metrics Jacobian in the cells
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradientCStoPSavg(const Box&           a_box,
                                   FArrayBox&           a_gradPhiPS,
                                   const FArrayBox&     a_gradPhiCS,
                                   const FluxBox&       a_N,
                                   const FArrayBox&     a_J)
{
  const int numCompPhi = a_gradPhiCS.nComp()/SpaceDim;
  CH_assert(numCompPhi*SpaceDim == a_gradPhiCS.nComp());
  CH_assert(a_gradPhiPS.nComp() >= a_gradPhiCS.nComp());
  CH_assert(a_gradPhiPS.box().contains(a_box));
  CH_assert(a_gradPhiCS.box().contains(a_box));
  CH_assert(a_N.box().contains(a_box));
  CH_assert(a_J.box().contains(a_box));

//--Store a contiguous N^T for each direction in the cells.  Both N and N^T
//--have column-major order.  We store N^T because we will be accessing the
//--rows of N.

  const int numNComp = m_dyadSize;
  CHARRAYNBSTACKTEMP(NtCtg, Real, SpaceDim+1, ArRangeCol, numNComp, a_box);
  cellNtCtg(a_box, NtCtg, a_N);

  // MD Arrays for calculations in the boxes
  MD_ARRAY_RESTRICT(gradPhiPSArr, a_gradPhiPS);
  MD_ARRAY_RESTRICT(gradPhiCSArr, a_gradPhiCS);
  MD_ARRAY_RESTRICT_CHARRAY(NtCtgArr, NtCtg);
  MD_ARRAY_RESTRICT(JArr, a_J);
  
  // loop over components of phi
  for (int cPhi = 0; cPhi != numCompPhi; ++cPhi)
    {
      // loop over the x-y-z directions of the physical space gradient
      for (int xDir = 0; xDir != SpaceDim; ++xDir)
        {
          MD_BOXLOOP(a_box, i)
            {
              gradPhiPSArr[MD_IX(i, getGradientComp(cPhi, xDir))] = (
                // unrolled row of gradPhiCS times column of Nt/J
                D_TERM6(  NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,0), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,0))],
                        + NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,1), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,1))],
                        + NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,2), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,2))],
                        + NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,3), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,3))],
                        + NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,4), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,4))],
                        +  NtCtgArr[MD_CIX(getTensorCompTranspose(xDir,5), i)]*
                          gradPhiCSArr[MD_IX(i, getGradientComp(cPhi,5))])
                )/JArr[MD_IX(i, 0)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
/// Compute gradient of phi in physical space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhi
 *                      Gradients of phi stored in components
 *                      [0, a_intv.size()*SpaceDim - 1] with order
 *                      (gradient direction, component)
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *  \param[in]  a_dxi   Computational mesh spacing
 *  \param[in]  a_N     Metrics on the faces
 *  \param[in]  a_J     Metrics Jacobian for cell averages
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradient2OPS(const Box&           a_box,
                              FArrayBox&           a_gradPhi,
                              const FArrayBox&     a_phi,
                              const FluxBox&       a_N,
                              const FArrayBox&     a_J,
                              const ProblemDomain& a_problemDomain,
                              const Interval&      a_intv,
                              const RealVect&      a_dxi)
{
  // solve the computational gradient to second order
  FABSTACKTEMP(gradPhiCS, a_box, a_intv.size()*SpaceDim);
  PatchMappedFunc::gradient2OCS(a_box,
                                gradPhiCS,
                                a_phi,
                                a_problemDomain,
                                a_intv,
                                a_dxi);
  // transform to physical space
  FABSTACKTEMP(Ncell, a_box, m_dyadSize);
  PatchMappedFunc::cellN(a_box,
                         Ncell,
                         a_N);
  PatchMappedFunc::gradientCStoPS(a_box,
                                  a_gradPhi,
                                  gradPhiCS,
                                  Ncell,
                                  a_J);
  // PatchMappedFunc::gradientCStoPSavg(a_box,
  //                                    a_gradPhi,
  //                                    gradPhiCS,
  //                                    a_N,
  //                                    a_J);
}

/*--------------------------------------------------------------------*/
/// Compute gradient of phi in physical space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhi
 *                      Gradients of phi stored in components
 *                      [0, a_intv.size()*SpaceDim - 1] with order
 *                      (gradient direction, component)
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *  \param[in]  a_dxi   Computational mesh spacing
 *  \param[in]  a_NtJ   Grid metrics
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::gradient2OPS(const Box&           a_box,
                              FArrayBox&           a_gradPhi,
                              const FArrayBox&     a_phi,
                              const FArrayBox&     a_NtJ,
                              const ProblemDomain& a_problemDomain,
                              const Interval&      a_intv,
                              const RealVect&      a_dxi)
{
  // solve the computational gradient to second order
  FABSTACKTEMP(gradPhiCS, a_box, a_intv.size()*SpaceDim);
  PatchMappedFunc::gradient2OCS(a_box,
                                gradPhiCS,
                                a_phi,
                                a_problemDomain,
                                a_intv,
                                a_dxi);
  // solve the physical gradient to second order
  PatchMappedFunc::gradientCStoPS(a_box,
                                  a_gradPhi,
                                  gradPhiCS,
                                  a_NtJ,
                                  a_gradPhi.interval(),
                                  gradPhiCS.interval());
}

/*--------------------------------------------------------------------*/
/// Compute divergence of phi in physical space given the gradient in
///  computational space and the grid metrics
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_divPhi
 *                      Divergence of phi stored in component 0
 *  \param[in]  a_gradPhiCS
 *                      gradient in computational space
 *  \param[in]  a_NtJ   grid metrics
 *  \param[in]  a_intvGrad
 *                      interval for gradient components
 *  \param[in]  a_compDiv
 *                      Component to save the divergence to
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divergencePS(const Box&           a_box,
                              FArrayBox&           a_divPhi,
                              const FArrayBox&     a_gradPhiCS,
                              const FArrayBox&     a_NtJ,
                              const Interval&      a_intvGrad,
                              const int            a_compDiv)
{
  // !!CHECK ME, function is untested!!
  CH_assert(a_intvGrad.size() == m_dyadSize);
  CH_assert(a_divPhi.interval().contains(a_compDiv));
  CH_assert(a_gradPhiCS.interval().contains(a_intvGrad));
  CH_assert(a_gradPhiCS.box().contains(a_box));
  CH_assert(a_divPhi.box().contains(a_box));
  CH_assert(a_NtJ.box().contains(a_box));

  // Solve the curl of phi in physical space
  a_divPhi.setVal(0., a_box, a_compDiv);
  // Solve the inner product (N/J):(\nabla_\xi phi)
  // MD Arrays for calculations in the boxes
  MD_ARRAY_RESTRICT(divPhiPSarr, a_divPhi);
  MD_ARRAY_RESTRICT(gradPhiCSarr, a_gradPhiCS);
  MD_ARRAY_RESTRICT(NtJarr, a_NtJ);
  // Get the gradient of phi in physical space
  // loop over the x-y-z direction in physical space
  for (int xDir = 0; xDir != SpaceDim; ++xDir)
    {
      MD_BOXLOOP(a_box, i)
        {
          // unroll the product of a column Nt/J by column of \nabla_\xi u 
          divPhiPSarr[MD_IX(i, a_compDiv)] +=
            D_TERM6(  NtJarr[MD_IX(i, getTensorComp(xDir,0))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,0))],
                      + NtJarr[MD_IX(i, getTensorComp(xDir,1))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,1))],
                      + NtJarr[MD_IX(i, getTensorComp(xDir,2))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,2))],
                      + NtJarr[MD_IX(i, getTensorComp(xDir,3))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,3))],
                      + NtJarr[MD_IX(i, getTensorComp(xDir,4))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,4))],
                      + NtJarr[MD_IX(i, getTensorComp(xDir,5))]
                      * gradPhiCSarr[MD_IX(i, getGradientComp(xDir,5))]);
        }
    }
}

/*--------------------------------------------------------------------*/
/// Compute divergence of phi in physical space given the gradient in
///  computational space and the grid metrics
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_divPhi
 *                      Divergence of phi stored in component 0
 *  \param[in]  a_gradPhiCS
 *                      gradient in computational space
 *  \param[in]  a_N     grid transformation matrix
 *  \param[in]  a_J     grid metric Jacobian
 *  \param[in]  a_intvGrad
 *                      interval for gradient components
 *  \param[in]  a_compDiv
 *                      Component to save the divergence to
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divergencePS(const Box&           a_box,
                              FArrayBox&           a_divPhi,
                              const FArrayBox&     a_gradPhiCS,
                              const FArrayBox&     a_N,
                              const FArrayBox&     a_J,
                              const Interval&      a_intvGrad,
                              const int            a_compDiv)
{
  // The division by J can be pulled out of the transformation
  PatchMappedFunc::divergencePS(a_box,
                                a_divPhi,
                                a_gradPhiCS,
                                a_N,
                                a_intvGrad,
                                a_compDiv);
  const int Jcomp = 0;
  PatchMappedFunc::divideVec(a_box,
                             a_divPhi,
                             Interval(a_compDiv, a_compDiv),
                             a_J,
                             Jcomp);
}

/*--------------------------------------------------------------------*/
/// Compute divergence of phi in physical space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_divPhi
 *                      Divergence of phi stored in component 0
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *                      Divergence is vector op, so requires a
 *                      SpaceDim sized vector to work on
 *  \param[in]  a_dxi   Computational mesh spacing
 *  \param[in]  a_N     Metrics on the faces
 *  \param[in]  a_J     Metrics Jacobian in the cells
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::divergence2OPS(const Box&           a_box,
                                FArrayBox&           a_divPhi,
                                const FArrayBox&     a_phi,
                                const FluxBox&       a_N,
                                const FArrayBox&     a_J,
                                const ProblemDomain& a_problemDomain,
                                const Interval&      a_intv,
                                const RealVect&      a_dxi)
{
  CH_assert(a_phi.nComp() >= SpaceDim);
  CH_assert(a_intv.size() == SpaceDim);
  CH_assert(a_divPhi.nComp() >= 1);
  CH_assert(a_divPhi.box().contains(a_box));
  CH_assert(a_N.box().contains(a_box));
  CH_assert(a_J.box().contains(a_box));
  
  // Solve \nabla_{x} Phi = (N/J) : (\nabla_{\xi} Phi)

  // solve the computational gradient to second order
  FABSTACKTEMP(gradPhiCS, a_box, a_intv.size()*SpaceDim);
  PatchMappedFunc::gradient2OCS(a_box,
                                gradPhiCS,
                                a_phi,
                                a_problemDomain,
                                a_intv,
                                a_dxi);
  // transform to physical space, and solve divergence
  FABSTACKTEMP(Ncell, a_box, m_dyadSize);
  PatchMappedFunc::cellN(a_box,
                         Ncell,
                         a_N);
  PatchMappedFunc::divergencePS(a_box,
                                a_divPhi,
                                gradPhiCS,
                                Ncell,
                                a_J);
}

// /*--------------------------------------------------------------------*/
// /// Compute curl of phi in physical space
// /** \param[in]  a_box   Cell box on which to compute the curl
//  *  \param[out] a_curlPhi
//  *                      Curl of phi stored in components
//  *                      [0, a_intv.size()] with order x, y, z for 3D
//  *                      or component 0 for 2D
//  *  \param[in]  a_gradPhiPS
//  *                      Gradient in computational space
//  *  \param[in]  a_intvGrad
//  *                      Interval for the gradient components
//  *  \param[in]  a_intvCurl
//  *                      Interval for curl components to be written
//  *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::curlPS(const Box&           a_box,
                        FArrayBox&           a_curlPhi,
                        const FArrayBox&     a_gradPhiPS,
                        const Interval&      a_intvGrad,
                        const Interval&      a_intvCurl)
{
  // curl only exists in 3D and 2D
  if (!((SpaceDim==2) || (SpaceDim==3)))
    {
      return;
    }
  
  CH_assert(a_gradPhiPS.interval().contains(a_intvGrad));
  CH_assert(a_curlPhi.interval().contains(a_intvCurl));
  CH_assert(a_intvCurl.size() == m_numCurlComps);
  CH_assert(a_intvGrad.size() == m_dyadSize);
  CH_assert(a_curlPhi.box().contains(a_box));
  CH_assert(a_gradPhiPS.box().contains(a_box));

  // Solve the curl of phi in physical space
  a_curlPhi.setVal(0., a_box, 0, m_numCurlComps);
  // loop components of the vector phi
  for (int phiComp = 0; phiComp != SpaceDim; ++phiComp)
    {
      // loop over x direction derivatives
      for (int xDir = 0; xDir != SpaceDim; ++xDir)
        {
          // skip aligned terms d_u[i]/d_x[i]
          if (phiComp == xDir) continue;
          // get the curl component to add element to
          int curlComp = (SpaceDim==3) ? SpaceDim - (phiComp + xDir) : 0;
          // sign of the element to add
          Real sign = (phiComp < xDir) ? -1 : 1;
          if (curlComp == 1) sign *= -1;
          // add the gradient element to the curl vector
          a_curlPhi.plus(a_gradPhiPS, sign, phiComp*SpaceDim + xDir, curlComp);
        }
    }
}

/*--------------------------------------------------------------------*/
/// Compute curl of phi in physical space to second order
/** \param[in]  a_box   Cell box on which to compute gradient
 *  \param[out] a_gradPhi
 *                      Gradients of phi stored in components
 *                      [0, a_intv.size()*SpaceDim - 1] with order
 *                      (gradient direction, component)
 *  \param[in]  a_phi   Field
 *  \param[in]  a_problemDomain
 *                      Problem domain
 *  \param[in]  a_intv  Interval of components in a_phi to consider
 *                      Curl is a vector op, so requires a
 *                      SpaceDim sized vector to work on
 *  \param[in]  a_dxi   Computational mesh spacing
 *  \param[in]  a_N     Metrics on the faces
 *  \param[in]  a_J     Metrics Jacobian in the cells
 *//*-----------------------------------------------------------------*/

void
PatchMappedFunc::curl2OPS(const Box&           a_box,
                          FArrayBox&           a_curlPhi,
                          const FArrayBox&     a_phi,
                          const FluxBox&       a_N,
                          const FArrayBox&     a_J,
                          const ProblemDomain& a_problemDomain,
                          const Interval&      a_intv,
                          const RealVect&      a_dxi)
{
  // curl only exists in 3D and 2D
  if (!((SpaceDim==2) || (SpaceDim==3)))
    {
      MayDay::Error("Curl requested in 1D");
      return;
    }
  
  // get \nabla_{\x} Phi, all components of this are used for computing curl
  FABSTACKTEMP(gradPhiPS, a_box, a_intv.size()*SpaceDim);
  PatchMappedFunc::gradient2OPS(a_box,
                                gradPhiPS,
                                a_phi,
                                a_N,
                                a_J,
                                a_problemDomain,
                                a_intv,
                                a_dxi);
  // Solve \nabla_{x} \cross Phi
  PatchMappedFunc::curlPS(a_box,
                          a_curlPhi,
                          gradPhiPS);
}

/*--------------------------------------------------------------------*/
/// Time-averaged velocity and Reynolds stress-tensor
/** \param[out] a_timeAvgData
 *                      FAB of time-averaged quantities
 *  \param[in]  a_WfaceAvgFxb
 *                      FluxBox of face-averaged primitive state
 *  \param[in]  a_WfacePntFxb
 *                      FluxBox of face-centered primitive state
 *  \param[in]  a_facePntVSTFxb
 *                      FluxBox of face-centered viscous stress tensor
 *  \param[in]  a_domain
 *                      Block domain on the current level
 *  \param[in]  a_box   Box over which to time-average
 *  \param[in]  a_dt    Time-step size on current level
 *  \param[in]  a_timeOld
 *                      Time at beginning of current time-step on
 *                      current level.
 *  \param[in]  a_timeAfterInit
 *                      Total time after starting to filter
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::timeAvgData(FArrayBox&           a_timeAvgData,
                             const FArrayBox&     a_WcellAvgFab,
                             const FluxBox&       a_WfaceAvgFxb,
                             const FluxBox&       a_WfacePntFxb,
                             const FluxBox&       a_facePntVSTFxb,
                             const FluxBox&       a_sgsMomentumFxb,
                             const ProblemDomain& a_domain,
                             const Box&           a_box,
                             const Real           a_dt,
                             const Real           a_timeOld,
                             const Real           a_timeAfterInit)
{
  CH_TIME("LevelMappedfunc::timeAvgData");

  if (!CRDparam::g_plotTimeAvgTurb) { return; }

  //**FIXME: All of the following plots low-faces as cell-centered data
  //         so that it can be printed in the standard plot file. However,
  //         this means that high-faces will be excluded from the plot file.
  //         For high wall boundaries, this is an issue.

  const int cRho = 0;
  const int cVel = 1; // Proper format would be to call velocityInterval() but
                      // we can't access that here
  const int cPress = 1 + SpaceDim;
  // Setup the necessary face-averaged values right here
  // 1) Compute uu_total = uu + tau_sgs
  Box box1Dom = grow(a_box, 1);
  box1Dom &= a_domain;
  FLUXBOXSTACKTEMP(facePntVelStress, box1Dom, SpaceDim*SpaceDim);
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      FArrayBox& facePntVelStressFab = facePntVelStress[faceDir];
      const FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
      const FArrayBox& facePntSGSFab = a_sgsMomentumFxb[faceDir];
      Box faceBox = grow(a_box, 1);
      faceBox.grow(faceDir, -1);
      faceBox.surroundingNodes(faceDir);
      faceBox &= a_domain;
      for (int velComp = 0; velComp != SpaceDim; ++velComp)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              int stressIndx = velComp + dir*SpaceDim;
              MD_BOXLOOP(faceBox, i)
                {
                  const Real rho = WfacePntFab[MD_IX(i, cRho)];
                  const Real vel1 = WfacePntFab[MD_IX(i, cVel+velComp)];
                  const Real vel2 = WfacePntFab[MD_IX(i, cVel+dir)];
                  const Real tauSGS = facePntSGSFab[MD_IX(i, stressIndx)]/rho;
                  facePntVelStressFab[MD_IX(i, stressIndx)] =
                    vel1*vel2 + tauSGS;
                }
            }
        }
    }
  // 2) Compute <uu_total> from the convolution
  FLUXBOXSTACKTEMP(faceAvgVelStress, a_box, SpaceDim*SpaceDim);
  // Perform the convolution right here
  CRDutil::convolveFace(faceAvgVelStress, facePntVelStress,
                        a_box, a_domain,
                        Interval(0,SpaceDim*SpaceDim-1), 4, false, false);

  // 3) Compute <SGS_stress> from the convolution
  FLUXBOXSTACKTEMP(facePntSGSStress, box1Dom, SpaceDim*SpaceDim);
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      FArrayBox& facePntSGSStressFab = facePntSGSStress[faceDir];
      const FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
      const FArrayBox& facePntSGSFab = a_sgsMomentumFxb[faceDir];
      Box faceBox = grow(a_box, 1);
      faceBox.grow(faceDir, -1);
      faceBox.surroundingNodes(faceDir);
      faceBox &= a_domain;
      for (int velComp = 0; velComp != SpaceDim; ++velComp)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              int stressIndx = velComp + dir*SpaceDim;
              MD_BOXLOOP(faceBox, i)
                {
                  const Real rho = WfacePntFab[MD_IX(i, cRho)];
                  const Real tauSGS = facePntSGSFab[MD_IX(i, stressIndx)]/rho;
                  facePntSGSStressFab[MD_IX(i, stressIndx)] = tauSGS;
                }
            }
        }
    }
  FLUXBOXSTACKTEMP(faceAvgSGSStress, a_box, SpaceDim*SpaceDim);
  // Perform the convolution right here
  CRDutil::convolveFace(faceAvgSGSStress, facePntSGSStress,
                        a_box, a_domain,
                        Interval(0,SpaceDim*SpaceDim-1), 4, false, false);
  // 4) Compute <rho*u> and <rho*u*u>
  FLUXBOXSTACKTEMP(facePntMomentumTerms, box1Dom, 2);
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      FArrayBox& facePntMomentumTermsFab = facePntMomentumTerms[faceDir];
      const FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
      Box faceBox = grow(a_box, 1);
      faceBox.grow(faceDir, -1);
      faceBox.surroundingNodes(faceDir);
      faceBox &= a_domain;
      MD_BOXLOOP(faceBox, i)
        {
          const Real rho = WfacePntFab[MD_IX(i, cRho)];
          const Real uVel = WfacePntFab[MD_IX(i, cVel)];
          facePntMomentumTermsFab[MD_IX(i, 0)] = rho*uVel;
          facePntMomentumTermsFab[MD_IX(i, 1)] = rho*uVel*uVel;
        }
    }
  FLUXBOXSTACKTEMP(faceAvgMomentumTerms, a_box, 2);
  // Perform the convolution right here
  CRDutil::convolveFace(faceAvgMomentumTerms, facePntMomentumTerms,
                        a_box, a_domain, Interval(0,1), 4, false, false);

  // If this is the first time-step, we need to initialize the time-avg
  // to the beginning instantaneous time value
  if ((a_timeOld == 0.) || (a_timeAfterInit == 0))
    {
      int comp = 0;
      // 1) Wall shear-stress from extrapolation
      //**FIXME: This definitely needs some work to be generalized
      const int normalDir = 1;
      const FArrayBox& faceVSTFab = a_facePntVSTFxb[normalDir];
      const FArrayBox& WfaceAvgYDirFab = a_WfaceAvgFxb[1];
      MD_BOXLOOP(a_box, i)
        {
          Real stressSqrd = 0.;
          for (int j = 0; j != SpaceDim; ++j)
            {
              if (j != normalDir)
                {
                  int stressIndx = j + normalDir*SpaceDim;
                  stressSqrd = faceVSTFab[MD_IX(i, stressIndx)]*
                    faceVSTFab[MD_IX(i, stressIndx)];
                }
            }
          const Real rho = WfaceAvgYDirFab[MD_IX(i, cRho)];
          a_timeAvgData[MD_IX(i, comp)] = std::sqrt(stressSqrd)/rho;
        }
      ++comp;
      // 2) Wall shear-stress from LES model
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          //**FIXME: Does not work when using species/thermally perfect
          const int etaIndx = 1  + SpaceDim + 1 + 1   + 1;
          //                  rho  velocity   p   sgsKE x-eta
          const Real eps2 = 1.e-40;
          MD_BOXLOOP(a_box, i)
            {
              Real eta_0 = WfaceAvgYDirFab[MD_IX(i, etaIndx)];
              RealVect surfVelVect(D_DECL(WfaceAvgYDirFab[MD_IX(i, cVel)],
                                          WfaceAvgYDirFab[MD_IX(i, cVel+1)],
                                          WfaceAvgYDirFab[MD_IX(i, cVel+2)]));
              Real signXVel = 1.;
              if (surfVelVect[0] < 0.) { signXVel = -1.; }
              // Project the velocity onto a 2D plane
              Real streamwiseVel =
                signXVel*std::sqrt(surfVelVect[0]*surfVelVect[0]
                                   + surfVelVect[1]*surfVelVect[1]);
              Real velMag = surfVelVect.vectorLength();
              Real unitStreamwiseVelComp = streamwiseVel/(velMag+eps2);
              a_timeAvgData[MD_IX(i, comp)] = eta_0*unitStreamwiseVelComp;
            }
          ++comp;
        }
      // 3) Sharp space-time-avg of face-averaged velocity on all faces
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const FArrayBox& WfaceAvgFab = a_WfaceAvgFxb[dir];
          for (int j = 0; j != SpaceDim; ++j)
            {
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] = WfaceAvgFab[MD_IX(i, cVel+j)];
                }
              ++comp;
            }
        }
      // 4) Sharp space-time-avg of face-centered velocity on all faces
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const FArrayBox& WfacePntFab = a_WfacePntFxb[dir];
          for (int j = 0; j != SpaceDim; ++j)
            {
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] = WfacePntFab[MD_IX(i, cVel+j)];
                }
              ++comp;
            }
        }
      // 5) Sharp space-time-avg of face-avg stress tensor on all faces
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          const FArrayBox& faceAvgVelStressFab = faceAvgVelStress[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int j = dir; j != SpaceDim; ++j)
                {
                  const int cStress = j + dir*SpaceDim;
                  MD_BOXLOOP(a_box, i)
                    {
                      a_timeAvgData[MD_IX(i, comp)] =
                        faceAvgVelStressFab[MD_IX(i, cStress)];
                    }
                  ++comp;
                }
            }
        }
      // 6) Sharp space-time-avg of face-centered stress tensor on all faces
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          const FArrayBox& facePntVelStressFab = facePntVelStress[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int j = dir; j != SpaceDim; ++j)
                {
                  const int cStress = j + dir*SpaceDim;
                  MD_BOXLOOP(a_box, i)
                    {
                      a_timeAvgData[MD_IX(i, comp)] =
                        facePntVelStressFab[MD_IX(i, cStress)];
                    }
                  ++comp;
                }
            }
        }
      // 7) Sharp space-time-avg of density and pressure in cells
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] = a_WcellAvgFab[MD_IX(i, cRho)];
        }
      ++comp;
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] = a_WcellAvgFab[MD_IX(i, cPress)];
        }
      ++comp;
      // 8) Sharp space-time-avg of rho*u and rho*u*u on x-faces
      const FArrayBox& faceAvgMomentumTermsXFab = faceAvgMomentumTerms[0];
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] = faceAvgMomentumTermsXFab[MD_IX(i, 0)];
        }
      ++comp;
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] = faceAvgMomentumTermsXFab[MD_IX(i, 1)];
        }
      ++comp;
      // 9) Sharp space-time-avg of eta_0^2
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          //**FIXME: Does not work when using species/thermally perfect
          // DEBUG ONLY -- removing sgsKE comp and x-eta comp
          const int etaIndx = 1  + SpaceDim + 1   + 1    + 1;
          //                  rho  velocity   p   sgsKE x-eta
          MD_BOXLOOP(a_box, i)
            {
              Real eta_0 = WfaceAvgYDirFab[MD_IX(i, etaIndx)];
              a_timeAvgData[MD_IX(i, comp)] = eta_0*eta_0;
            }
          ++comp;
        }
      // 10) Sharp space-time-avg of SGS stress tensor on x-faces
      const FArrayBox& faceAvgSGSStressFab = faceAvgSGSStress[0];
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          for (int j = dir; j != SpaceDim; ++j)
            {
              const int cStress = j + dir*SpaceDim;
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] =
                    faceAvgSGSStressFab[MD_IX(i, cStress)];
                }
              ++comp;
            }
        }
    }
  else // If this is after the first time-step, start averaging
    {
      // Parameter for sharp time-filter
      const Real a_0 = 1./(a_timeAfterInit + a_dt);

      int comp = 0;
      // 1) Wall shear-stress from extrapolation
      //**FIXME: This definitely needs some work to be generalized
      const int normalDir = 1;
      const FArrayBox& faceVSTFab = a_facePntVSTFxb[normalDir];
      const FArrayBox& WfaceAvgYDirFab = a_WfaceAvgFxb[1];
      MD_BOXLOOP(a_box, i)
        {
          Real stressSqrd = 0.;
          for (int j = 0; j != SpaceDim; ++j)
            {
              if (j != normalDir)
                {
                  int stressIndx = j + normalDir*SpaceDim;
                  stressSqrd = faceVSTFab[MD_IX(i, stressIndx)]*
                    faceVSTFab[MD_IX(i, stressIndx)];
                }
            }
          const Real rho = WfaceAvgYDirFab[MD_IX(i, cRho)];
          a_timeAvgData[MD_IX(i, comp)] =
            a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                 + a_dt*std::sqrt(stressSqrd)/rho);
        }
      ++comp;
      // 2) Wall shear-stress from LES model
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          //**FIXME: Does not work when using species/thermally perfect
          const int etaIndx = 1  + SpaceDim + 1 + 1   + 1;
          //                  rho  velocity   p   sgsKE x-eta
          MD_BOXLOOP(a_box, i)
            {
              RealVect surfVelVect(D_DECL(WfaceAvgYDirFab[MD_IX(i, cVel)],
                                          WfaceAvgYDirFab[MD_IX(i, cVel+1)],
                                          WfaceAvgYDirFab[MD_IX(i, cVel+2)]));
              Real signXVel = 1.;
              if (surfVelVect[0] < 0.) { signXVel = -1.; }
              // Project the velocity onto a 2D plane
              Real streamwiseVel =
                signXVel*std::sqrt(surfVelVect[0]*surfVelVect[0]
                                   + surfVelVect[1]*surfVelVect[1]);
              Real velMag = surfVelVect.vectorLength();
              Real unitStreamwiseVelComp = streamwiseVel/velMag;
              Real eta_0 =
                unitStreamwiseVelComp*WfaceAvgYDirFab[MD_IX(i, etaIndx)];
              a_timeAvgData[MD_IX(i, comp)] =
                a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                     + a_dt*eta_0);
            }
          ++comp;
        }
      // 3) Sharp space-time-avg of face-averaged velocity on all faces
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const FArrayBox& WfaceAvgFab = a_WfaceAvgFxb[dir];
          for (int j = 0; j != SpaceDim; ++j)
            {
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] =
                    a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                         + a_dt*WfaceAvgFab[MD_IX(i, cVel+j)]);
                }
              ++comp;
            }
        }
      // 4) Sharp space-time-avg of face-centered velocity on all faces
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          const FArrayBox& WfacePntFab = a_WfacePntFxb[dir];
          for (int j = 0; j != SpaceDim; ++j)
            {
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] =
                    a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                         + a_dt*WfacePntFab[MD_IX(i, cVel+j)]);
                }
              ++comp;
            }
        }
      // 5) Sharp space-time-avg of face-avg stress tensor on all faces
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          const FArrayBox& faceAvgVelStressFab = faceAvgVelStress[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int j = dir; j != SpaceDim; ++j)
                {
                  const int cStress = j + dir*SpaceDim;
                  MD_BOXLOOP(a_box, i)
                    {
                      a_timeAvgData[MD_IX(i, comp)] =
                        a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                             + a_dt*faceAvgVelStressFab[MD_IX(i, cStress)]);
                    }
                  ++comp;
                }
            }
        }
      // 6) Sharp space-time-avg of face-centered stress tensor on all faces
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          const FArrayBox& facePntVelStressFab = facePntVelStress[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int j = dir; j != SpaceDim; ++j)
                {
                  const int cStress = j + dir*SpaceDim;
                  MD_BOXLOOP(a_box, i)
                    {
                      a_timeAvgData[MD_IX(i, comp)] =
                        a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                             + a_dt*facePntVelStressFab[MD_IX(i, cStress)]);
                    }
                  ++comp;
                }
            }
        }
      // 7) Sharp space-time-avg of density and pressure in cells
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] =
            a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                 + a_dt*a_WcellAvgFab[MD_IX(i, cRho)]);
        }
      ++comp;
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] =
            a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                 + a_dt*a_WcellAvgFab[MD_IX(i, cPress)]);
        }
      ++comp;
      // 8) Sharp space-time-avg of rho*u and rho*u*u on x-faces
      const FArrayBox& faceAvgMomentumTermsXFab = faceAvgMomentumTerms[0];
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] =
            a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                 + a_dt*faceAvgMomentumTermsXFab[MD_IX(i, 0)]);
        }
      ++comp;
      MD_BOXLOOP(a_box, i)
        {
          a_timeAvgData[MD_IX(i, comp)] =
            a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                 + a_dt*faceAvgMomentumTermsXFab[MD_IX(i, 1)]);
        }
      ++comp;
      // 9) Sharp space-time-avg of eta_0^2
      if (CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)
        {
          //**FIXME: Does not work when using species/thermally perfect
          // DEBUG ONLY -- removing sgsKE comp and x-eta comp
          const int etaIndx = 1  + SpaceDim + 1 + 1   + 1;
          //                  rho  velocity   p   sgsKE x-eta
          MD_BOXLOOP(a_box, i)
            {
              Real eta_0 = WfaceAvgYDirFab[MD_IX(i, etaIndx)];
              a_timeAvgData[MD_IX(i, comp)] =
                a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                     + a_dt*eta_0*eta_0);
            }
          ++comp;
        }
      // 10) Sharp space-time-avg of SGS stress tensor on x-faces
      const FArrayBox& faceAvgSGSStressFab = faceAvgSGSStress[0];
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          for (int j = dir; j != SpaceDim; ++j)
            {
              const int cStress = j + dir*SpaceDim;
              MD_BOXLOOP(a_box, i)
                {
                  a_timeAvgData[MD_IX(i, comp)] =
                    a_0*(a_timeAfterInit*a_timeAvgData[MD_IX(i, comp)]
                         + a_dt*faceAvgSGSStressFab[MD_IX(i, cStress)]);
                }
              ++comp;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
/// Time-averaged velocity and Reynolds stress-tensor
/** \param[out] a_timeAvgData
 *                      FAB of time-averaged quantities
 *  \param[in]  a_WfaceAvgFxb
 *                      FluxBox of face-averaged primitive state
 *  \param[in]  a_WfacePntFxb
 *                      FluxBox of face-centered primitive state
 *  \param[in]  a_facePntVSTFxb
 *                      FluxBox of face-centered viscous stress tensor
 *  \param[in]  a_domain
 *                      Block domain on the current level
 *  \param[in]  a_box   Box over which to time-average
 *  \param[in]  a_dt    Time-step size on current level
 *  \param[in]  a_timeOld
 *                      Time at beginning of current time-step on
 *                      current level.
 *  \param[in]  a_timeAfterInit
 *                      Total time after starting to filter
 *//*-----------------------------------------------------------------*/
void
PatchMappedFunc::timeAvgData(FluxBox&             a_timeAvgDataFxb,
                             const FluxBox&       a_instantaneousDataFxb,
                             const ProblemDomain& a_domain,
                             const Box&           a_box,
                             const Real           a_dt,
                             const Real           a_timeOld,
                             const Real           a_timeAfterInit)
{
  CH_TIME("LevelMappedfunc::timeAvgData");

  if (!CRDparam::g_plotLoFaceAvgTimeAvgComps) { return; }

  // If this is the first time-step, we need to initialize the time-avg
  // to the beginning instantaneous time value
  if ((a_timeOld == 0.) || (a_timeAfterInit == 0))
    {
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          FArrayBox& timeAvgDataFab = a_timeAvgDataFxb[faceDir];
          const FArrayBox& instantaneousDataFab =
            a_instantaneousDataFxb[faceDir];
          Box faceBox = a_box;
          faceBox.surroundingNodes(faceDir);
          int numComps = timeAvgDataFab.nComp();
          for (int comp = 0; comp != numComps; ++comp)
            {
              MD_BOXLOOP(faceBox, i)
                {
                  timeAvgDataFab[MD_IX(i, comp)] =
                    instantaneousDataFab[MD_IX(i, comp)];
                }
            }
        }
    }
  else // If this is after the first time-step, start averaging
    {
      // Parameter for sharp time-filter
      const Real a_0 = 1./(a_timeAfterInit + a_dt);

      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          FArrayBox& timeAvgDataFab = a_timeAvgDataFxb[faceDir];
          const FArrayBox& instantaneousDataFab =
            a_instantaneousDataFxb[faceDir];
          Box faceBox = a_box;
          faceBox.surroundingNodes(faceDir);
          int numComps = timeAvgDataFab.nComp();
          for (int comp = 0; comp != numComps; ++comp)
            {
              MD_BOXLOOP(faceBox, i)
                {
                  timeAvgDataFab[MD_IX(i, comp)] =
                    a_0*(a_timeAfterInit*timeAvgDataFab[MD_IX(i, comp)]
                         + a_dt*(instantaneousDataFab[MD_IX(i, comp)]));
                }
            }
        }
    }
}
