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
 * \file PatchMappedFunc.H
 *
 * \brief Patch functions on mapped grids
 *
 *//*+*************************************************************************/

#ifndef _PATCHMAPPEDFUNC_H_
#define _PATCHMAPPEDFUNC_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "CHArray_fwd.H"
#include "UsingNamespace.H"

//----- Internal -----//

//--Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class ProblemDomain;
class FArrayBox;
class FluxBox;
class Interval;
#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"
namespace stc
{

template <typename T,
          array_size_type N,
          template<typename, array_size_type> class Impl>
class VectorImpl;

}
#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"
using RealVect = stc::VectorImpl<Real, SpaceDim, stc::DefaultImpl>;
#include "NamespaceFooter.H"


/*******************************************************************************
 *
 * Patch functions on mapped grids
 *
 ******************************************************************************/

namespace PatchMappedFunc
{
  /*****************************************************************************
   * namespace constants
   ****************************************************************************/

  // size of a dyad
  constexpr int m_dyadSize = SpaceDim*SpaceDim;

  // number of compents in curl
  constexpr int m_numCurlComps = (SpaceDim==3) ? 3 : 1;
  
  /*****************************************************************************
   * Generic vector functions 
   ****************************************************************************/
  
  /// element of a tensor, in row major order
  /// FIXME: There should really be a global version of this somewhere and
  ///        maybe there, but it is poorly documented and not frequently used
  ///        This matches LGM getNcomponent
  inline int getTensorComp(const int a_row, const int a_col)
  {
    return (SpaceDim*a_col) + a_row;
  }

  /// transpose of a tensor, or column major order
  ///   This is how gradients are stored for a number of
  ///   practical reasons both in the math and code.
  inline int getTensorCompTranspose(const int a_row, const int a_col)
  {
    return PatchMappedFunc::getTensorComp(a_col, a_row);
  }

  /// Same as above, but gives the gradient its own function to reduce confusion
  ///  return component for d_u[a_comp]/d_x[a_deriv]
  inline int getGradientComp(const int a_comp, const int a_deriv)
  {
    return PatchMappedFunc::getTensorComp(a_deriv, a_comp);
  }
  
  /// Obtain magnitude of spatial vector phi
  void magnitudeSpace(const Box&       a_box,
                      FArrayBox&       a_mag,
                      const int        a_cMag,
                      const FArrayBox& a_phi,
                      const int        a_cPhi);

  /// Obtain magnitude of arbitrary sized vector phi
  void magnitude(const Box&       a_box,
                 FArrayBox&       a_mag,
                 const int        a_cMag,
                 const FArrayBox& a_phi,
                 const Interval&  a_cPhi);

  /// divide a vector by scalar, element by element. Happens in place
  void divideVec(const Box&           a_box,
                 FArrayBox&           a_phi,
                 const Interval&      a_intvPhi,
                 const FArrayBox&     a_denominator,
                 const int            a_denomComp);
  
  /// divide a vector by scalar, element by element
  void divideVec(const Box&           a_box,
                 FArrayBox&           a_phiDiv,
                 const FArrayBox&     a_phi,
                 const Interval&      a_intvPhi,
                 const FArrayBox&     a_denominator,
                 const int            a_denomComp,
                 const Interval&      a_intvDest);

  /// divide a vector by scalar, element by element
  void divideVec(const Box&           a_box,
                 FArrayBox&           a_phiDiv,
                 const FArrayBox&     a_phi,
                 const Interval&      a_intvPhi,
                 const int            a_denomComp,
                 const Interval&      a_intvDest);

  /// scale a vector to be of unit length
  void normalize(const Box&           a_box,
                 FArrayBox&           a_Phi,
                 const Interval&      a_intv);

  /*****************************************************************************
   * Get forms of the mapping variables
   ****************************************************************************/

  // FIXME:: These functions probably should not exist in their current from but
  //         instead should have some way to get cell averaged mapping from LGM
  
  /// Obtain N^T in a cell with components that are contiguous in memory.
  /// This is a second order average
  /// Although probably more effecient to use this version, it makes
  // generalizing code challenging so I try to avoid it
  void cellNtCtg(const Box&                             a_box,
                 CHArray<Real, SpaceDim+1, ArRangeCol>& a_NtCtg,
                 const FluxBox&                         a_N);

  // get N in a cell using a second order average from face values
  void cellN(const Box&     a_box,
             FArrayBox&     a_Ncell,
             const FluxBox& a_N);

  /*****************************************************************************
   * Do mapping transformations
   ****************************************************************************/
  /// get unit normals
  void getUnitNormals(FArrayBox& a_normals,
                      const FArrayBox& a_faceMetrics,
                      const IntVect& a_normalComps,
                      const int& a_dir,
                      const Box& a_box);

  /// Transform a vector in computational space to one in physical space
  void forwardTransform(FArrayBox& a_U,
                        const FArrayBox& a_unitNormals,
                        const Box& a_box);

  /// Transform a vector in computational space to one in physical space
  void forwardTransform(FArrayBox& a_Uout,
                        const int& a_outComp,
                        const FArrayBox& a_Uin,
                        const int& a_inComp,
                        const FArrayBox& a_unitNormals,
                        const Box& a_box);

  /// Transform a vector in physical space to computational physical space
  void reverseTransform(FArrayBox& a_U,
                        const FArrayBox& a_unitNormals,
                        const Box& a_box);

  /// Transform a vector in physical space to computational physical space
  void reverseTransform(FArrayBox& a_Uout,
                        const int& a_outComp,
                        const FArrayBox& a_Uin,
                        const int& a_inComp,
                        const FArrayBox& a_unitNormals,
                        const Box& a_box);
  
  /*****************************************************************************
   * Calculate computational space gradients
   ****************************************************************************/
  
  /// Compute gradients of phi in computational space to second order
  ///  one sided on domain boundaries, centered elsewhere
  void gradient2OCS(const Box&           a_box,
                    FArrayBox&           a_gradPhi,
                    const FArrayBox&     a_phi,
                    const ProblemDomain& a_problemDomain,
                    const Interval&      a_intv,
                    const RealVect&      a_dxi);

  /// Compute undivided gradients of phi in computational space to second order
  ///  one sided on domain boundaries, centered elsewhere
  void undividedGrad2OCS(const Box&         a_box,
                         FArrayBox&         a_gradPhi,
                         const FArrayBox&   a_phi,
                         const BlockDomain& a_problemDomain,
                         const Interval&    a_intv);

  /// Comp space face-avgd face-normal 1st-deriv from cell-avgd vals, codim-1
  void faceAvgNormDerivFromCellAvgCS(FArrayBox&         a_faceAvgNormalDeriv,
                                     const FArrayBox&   a_cellAvgVal,
                                     const int          a_valComp,
                                     const int          a_derivComp,
                                     const BlockDomain& a_problemDomain,
                                     const Box&         a_box,
                                     const RealVect&    a_dxi,
                                     const int          a_dir,
                                     const bool         a_fourthOrder   = true,
                                     const bool         a_interior      = true,
                                     const bool         a_unitStrideDim = true);

  /// Comp space cell-avgd 1st-deriv from cell-avgd vals
  void cellAvgDerivFromCellAvgCS(FArrayBox&         a_cellAvgDeriv,
                                 const FArrayBox&   a_cellAvgVal,
                                 const int          a_valComp,
                                 const int          a_derivComp,
                                 const BlockDomain& a_problemDomain,
                                 const Box&         a_box,
                                 const RealVect&    a_dxi,
                                 const int          a_dir,
                                 const bool         a_fourthOrder   = true,
                                 const bool         a_interior      = true,
                                 const bool         a_unitStrideDim = true);

  /// Comp space face-avgd value from cell-avgd vals, codim-1 version
  void faceAvgValFromCellAvgCS(FArrayBox&         a_faceAvgVal,
                               const FArrayBox&   a_cellAvgVal,
                               const int          a_cellComp,
                               const int          a_faceComp,
                               const BlockDomain& a_problemDomain,
                               const Box&         a_box,
                               const int          a_dir,
                               const int          a_order    = 4,
                               const bool         a_interior = true);

  /// Comp space cell-avgd grad from face-avgd vals, exact calculation
  void cellAvgGradFromFaceAvgCS(FArrayBox&      a_cellAvgGrad,
                                const FluxBox&  a_faceAvgVal,
                                const Interval& a_inputIntv,
                                const int       a_outputStartComp,
                                const Box&      a_box,
                                const RealVect& a_dxi);

  /// Comp space cell-avgd derivative from face-avgd vals, exact calculation
  void cellAvgDerivFromFaceAvgCS(FArrayBox&       a_cellAvgGrad,
                                 const FArrayBox& a_faceAvgVal,
                                 const int&       a_inputComp,
                                 const int&       a_outputComp,
                                 const Box&       a_box,
                                 const RealVect&  a_dxi,
                                 const int        a_dir);

  /*****************************************************************************
   * Calculate physical space gradients
   *   These transformation look like
   *   $\nabla_{x} \vec{u}  =  (\nabla_{\xi} \vec{u})(N^T/J)$
   ****************************************************************************/

  /// Convert gradients of phi in computational space to physical space
  ///   a_gradPhiCS, and a_NtJ must be of same type (ie. point, face avg, ...)
  ///   this is the primary gradient transformation
  void gradientCStoPS(const Box&       a_box,
                      FArrayBox&       a_gradPhiPS,
                      const FArrayBox& a_gradPhiCS,
                      const FArrayBox& a_NtJ,
                      const Interval&  a_intvPS,
                      const Interval&  a_intvCS);

  /// Convert gradients of phi in computational space to physical space
  /// for all components
  void gradientCStoPS(const Box&       a_box,
                      FArrayBox&       a_gradPhiPS,
                      const FArrayBox& a_gradPhiCS,
                      const FArrayBox& a_NtJ);
  
  /// Convert gradients of phi in computational space to physical space
  ///   a_gradPhiCS, and a_N, and a_J must be of same type (ie. point, face avg, ...)
  void gradientCStoPS(const Box&       a_box,
                      FArrayBox&       a_gradPhiPS,
                      const FArrayBox& a_gradPhiCS,
                      const FArrayBox& a_N,
                      const FArrayBox& a_J,
                      const Interval&  a_intvPS,
                      const Interval&  a_intvCS);

  /// Convert gradients of phi in computational space to physical space
  /// for all components
  void gradientCStoPS(const Box&       a_box,
                      FArrayBox&       a_gradPhiPS,
                      const FArrayBox& a_gradPhiCS,
                      const FArrayBox& a_N,
                      const FArrayBox& a_J);
  
  /// Convert gradients of phi in computational space to physical space
  ///  for cell averaged quantities
  void gradientCStoPSavg(const Box&       a_box,
                         FArrayBox&       a_gradPhiPS,
                         const FArrayBox& a_gradPhiCS,
                         const FluxBox&   a_N,
                         const FArrayBox& a_J);
  
  /// Compute gradient of phi in physical space to second order
  void gradient2OPS(const Box&           a_box,
                    FArrayBox&           a_gradPhi,
                    const FArrayBox&     a_phi,
                    const FluxBox&       a_N,
                    const FArrayBox&     a_J,
                    const ProblemDomain& a_problemDomain,
                    const Interval&      a_intv,
                    const RealVect&      a_dxi);

  /// Compute gradient of phi in physical space to second order
  void gradient2OPS(const Box&           a_box,
                    FArrayBox&           a_gradPhi,
                    const FArrayBox&     a_phi,
                    const FArrayBox&     a_NtJ,
                    const ProblemDomain& a_problemDomain,
                    const Interval&      a_intv,
                    const RealVect&      a_dxi);
  
  /*****************************************************************************
   * Calculate physical space divergence
   *   These transformation look like
   *   $\nabla_{x} \cdot \vec{u}  =  (N/J) : (\nabla_{\xi} \vec{u})$
   ****************************************************************************/

  /// Divergence of phi in physical space, given the computational space gradient
  ///  this is the primary divergence operation
  void divergencePS(const Box&       a_box,
                    FArrayBox&       a_divPhi,
                    const FArrayBox& a_gradPhiCS,
                    const FArrayBox& a_NtJ,
                    const Interval&  a_intvGrad = Interval(0, m_dyadSize-1),
                    const int        a_compDiv = 0);

  /// Divergence of phi in physical space, given the computational space gradient
  void divergencePS(const Box&       a_box,
                    FArrayBox&       a_divPhi,
                    const FArrayBox& a_gradPhiCS,
                    const FArrayBox& a_N,
                    const FArrayBox& a_J,
                    const Interval&  a_intvGrad = Interval(0, m_dyadSize-1),
                    const int        a_compDiv = 0);
  
  /// Compute divergence of phi in physical space to second order
  void divergence2OPS(const Box&           a_box,
                      FArrayBox&           a_divPhi,
                      const FArrayBox&     a_phi,
                      const FluxBox&       a_N,
                      const FArrayBox&     a_J,
                      const ProblemDomain& a_problemDomain,
                      const Interval&      a_intv,
                      const RealVect&      a_dxi);

  /*****************************************************************************
   * Calculate physical space curl
   *   These do not have nice transformation from computational space
   ****************************************************************************/

  /// Curl of phi in physical space, given the physical space gradient
  ///   this is the primary curl operation
  void curlPS(const Box&       a_box,
              FArrayBox&       a_curlPhi,
              const FArrayBox& a_gradPhiPS,
              const Interval&  a_intvGrad = Interval(0, m_dyadSize-1),
              const Interval&  a_intvCurl = Interval(0, m_numCurlComps-1));

  /// Curl of phi in physical space, given the computational space gradient
  void curlPS(const Box&       a_box,
              FArrayBox&       a_curlPhi,
              const FArrayBox& a_gradPhiCS,
              const FArrayBox& a_NtJ,
              const Interval&  a_intvGrad = Interval(0, m_dyadSize-1),
              const Interval&  a_intvCurl = Interval(0, m_numCurlComps-1));
  
  /// Compute curl of phi in physical space to second order
  void curl2OPS(const Box&           a_box,
                FArrayBox&           a_curlPhi,
                const FArrayBox&     a_phi,
                const FluxBox&       a_N,
                const FArrayBox&     a_J,
                const ProblemDomain& a_problemDomain,
                const Interval&      a_intv,
                const RealVect&      a_dxi);

  /// Time-averaged velocity and Reynolds stress-tensor
  void timeAvgData(FArrayBox&           a_timeAvgData,
                   const FArrayBox&     a_WcellAvgFab,
                   const FluxBox&       a_WfaceAvgFxb,
                   const FluxBox&       a_WfacePntFxb,
                   const FluxBox&       a_facePntVSTFxb,
                   const FluxBox&       a_sgsMomentumFxb,
                   const ProblemDomain& a_domain,
                   const Box&           a_box,
                   const Real           a_dt,
                   const Real           a_timeOld,
                   const Real           a_timeAfterInit);

  /// Time-averaged flux-box
  void timeAvgData(FluxBox&             a_timeAvgDataFxb,
                   const FluxBox&       a_instantaneousDataFxb,
                   const ProblemDomain& a_domain,
                   const Box&           a_box,
                   const Real           a_dt,
                   const Real           a_timeOld,
                   const Real           a_timeAfterInit);
}

#endif  /* ! defined _PATCHMAPPEDFUNC_H_ */
