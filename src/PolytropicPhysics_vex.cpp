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
 * \file PolytropicPhysics_vex.cpp
 *
 * \brief Vector implementations of kernels in PolytropicPhysics
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "BaseFabMacros.H"
#include "CH_Timer.H"

//----- Internal -----//

#include "VEx.H"
#include "CRDparam.H"
#include "LGintegrator.H"

//----- VEx -----//

#ifdef CRD_USE_VEX
#include "SSETypes.H"
#include "SSESupport.H"
#endif

/*
 * In tests, O means ordered (false if either argument is NAN) and Q means
 * non-signalling.  All tests ar OQ in this routine.
 */
#ifdef CRD_USE_VEX
void riemann_vex(FArrayBox&       a_Wstar,
                 const FArrayBox& a_Wleft,
                 const FArrayBox& a_Wright,
                 const int        a_dir,
                 const Box&       a_box)
{
  CH_TIMELEAF("VEX::riemann_vex");
  CH_assert(WRHO == 0);
  D_TERM(const int inorm = WVELX + a_dir;,
         const int itan1 = WVELX + ((a_dir + 1) % SpaceDim);,
         const int itan2 = WVELX + ((a_dir + 2) % SpaceDim);)

  const int pencilSize   = a_box.size(0);
  const int vecPacked    = pencilSize/VecSz_r;
  const int vecRmdr      = pencilSize - vecPacked*VecSz_r;
  const int i0EndRPacked = a_box.smallEnd(0) + vecPacked*VecSz_r;

  const __mvr gamma  = _mm_vr(set1)(CRDparam::g_gamma);
  const __mvr small  = _mm_vr(set1)(CRDparam::g_small);
  const __mvr smallp = _mm_vr(set1)(CRDparam::g_smallp);
  const __mvr smallr = _mm_vr(set1)(CRDparam::g_smallr);

  // Component strides
  const ptrdiff_t WleftCmpStride =
    &a_Wleft(a_Wleft.smallEnd(), WRHO+1) - &a_Wleft(a_Wleft.smallEnd(), WRHO);
  const ptrdiff_t WrightCmpStride =
    &a_Wright(a_Wright.smallEnd(), WRHO+1) - &a_Wright(a_Wright.smallEnd(), WRHO);
  const ptrdiff_t WstarCmpStride =
    &a_Wstar(a_Wstar.smallEnd(), WRHO+1) - &a_Wstar(a_Wstar.smallEnd(), WRHO);

  // Set the load/store mask for the last vector
  CHvi_t lsmask;
  lsmask.m = _mm_si(setzero)();
  for (int i = 0; i != vecRmdr; ++i)
    {
      lsmask.u[i] = itrue;
    }

  ForAllBPencil(a_box)
    IntVect iv(D_DECL(i0BegR, i1R, i2R));
    const Real *Wleft_ptr = &a_Wleft(iv, WRHO);
    const Real *Wright_ptr = &a_Wright(iv, WRHO);
    Real *Wstar_ptr = &a_Wstar(iv, WRHO);
    for (int i0R = i0BegR; i0R < i0EndR; i0R += VecSz_r)
      {
        iv[0] = i0R;
        __mvr pl, rhol, D_DECL(unl, utan1l, utan2l);
        __mvr pr, rhor, D_DECL(unr, utan1r, utan2r);

//--Load

        if (i0R < i0EndRPacked)  // Can fully load the vector
          {
                   rhol   = _mm_vr(loadu)(Wleft_ptr);
            D_TERM(unl    = _mm_vr(loadu)(Wleft_ptr + WleftCmpStride*inorm);,
                   utan1l = _mm_vr(loadu)(Wleft_ptr + WleftCmpStride*itan1);,
                   utan2l = _mm_vr(loadu)(Wleft_ptr + WleftCmpStride*itan2);)
                   pl     = _mm_vr(loadu)(Wleft_ptr + WleftCmpStride*WPRES);

                   rhor   = _mm_vr(loadu)(Wright_ptr);
            D_TERM(unr    = _mm_vr(loadu)(Wright_ptr + WrightCmpStride*inorm);,
                   utan1r = _mm_vr(loadu)(Wright_ptr + WrightCmpStride*itan1);,
                   utan2r = _mm_vr(loadu)(Wright_ptr + WrightCmpStride*itan2);)
                   pr     = _mm_vr(loadu)(Wright_ptr + WrightCmpStride*WPRES);
          }
        else                   // Partial loading of vector
          {
                   rhol   = _mm_vr(maskload)(Wleft_ptr                       , lsmask.m);
            D_TERM(unl    = _mm_vr(maskload)(Wleft_ptr + WleftCmpStride*inorm, lsmask.m);,
                   utan1l = _mm_vr(maskload)(Wleft_ptr + WleftCmpStride*itan1, lsmask.m);,
                   utan2l = _mm_vr(maskload)(Wleft_ptr + WleftCmpStride*itan2, lsmask.m);)
                   pl     = _mm_vr(maskload)(Wleft_ptr + WleftCmpStride*WPRES, lsmask.m);

                   rhor   = _mm_vr(maskload)(Wright_ptr                        , lsmask.m);
            D_TERM(unr    = _mm_vr(maskload)(Wright_ptr + WrightCmpStride*inorm, lsmask.m);,
                   utan1r = _mm_vr(maskload)(Wright_ptr + WrightCmpStride*itan1, lsmask.m);,
                   utan2r = _mm_vr(maskload)(Wright_ptr + WrightCmpStride*itan2, lsmask.m);)
                   pr     = _mm_vr(maskload)(Wright_ptr + WrightCmpStride*WPRES, lsmask.m);
          }
        // Debug example - write out pl in cell (9,15)
        // if ((a_dir == 0) && ((i0R-9)*(i0R+VecSz_r-1-9) <= 0) && (i1R == 15))
        //   {
        //     CHvr_t cvr;
        //     cvr.m = pl;
        //     std::cout << cvr.f[9 - i0R] << std::endl;
        //   }

//--Compute

        pl = _mm_vr(max)(smallp, pl);
        pr = _mm_vr(max)(smallp, pr);
        rhol = _mm_vr(max)(smallr, rhol);
        rhor = _mm_vr(max)(smallr, rhor);

        const __mvr cl = _mm_vr(sqrt)(gamma*pl/rhol);
        const __mvr cr = _mm_vr(sqrt)(gamma*pr/rhor);

        const __mvr wl = rhol*cl;
        const __mvr wr = rhor*cr;

        const __mvr pstar = (wr*pl + wl*pr + wl*wr*(unl - unr))/(wl + wr);
        const __mvr ustar = (wl*unl + wr*unr + pl - pr)/(wl + wr);

        __mvr test_vr;  // conditional flag in full word

        // test: ustar > 0
        test_vr = _mm_vr(cmp)(ustar, _mm_vr(setzero)(), _CMP_GT_OQ);
               const __mvr ro    = _mm_vr(blendv)(rhor  , rhol  , test_vr);
               const __mvr po    = _mm_vr(blendv)(pr    , pl    , test_vr);
        D_TERM(const __mvr uno   = _mm_vr(blendv)(unr   , unl   , test_vr);,
               const __mvr utan1 = _mm_vr(blendv)(utan1r, utan1l, test_vr);,
               const __mvr utan2 = _mm_vr(blendv)(utan2r, utan2l, test_vr);)
               const __mvr co    = _mm_vr(blendv)(cr    , cl    , test_vr);
        const __mvr sgnm = _mm_vr(blendv)(minusone_vr, one_vr, test_vr);

        __mvr rstar = ro + (pstar - po)/pow_vr_si<2>(co);
        rstar = _mm_vr(max)(smallr, rstar);

        const __mvr cstar = _mm_vr(sqrt)(abs_vr(gamma*pstar/rstar));
        const __mvr wstar = _mm_vr(set1)(0.5)*(cstar*rstar + co*ro);
        
        __mvr spout = co - sgnm*uno;
        __mvr spin  = cstar - sgnm*ustar;

        const __mvr ushock = wstar/rstar - sgnm*ustar;

        // test: pstar > po
        test_vr = _mm_vr(cmp)(pstar, po, _CMP_GT_OQ);
        // if (_mm_vr(movemask)(test_vr))  // No apparent benefit from this
        //   {
        spout = _mm_vr(blendv)(spout, ushock, test_vr);
        spin  = _mm_vr(blendv)(spin , ushock, test_vr);
        //   }

        __mvr frac = _mm_vr(set1)(0.5)*
          (one_vr + (spout + spin)/_mm_vr(max)(spout-spin, small));
        frac = _mm_vr(max)(_mm_vr(setzero)(),_mm_vr(min)(one_vr, frac));

        __mvr rsto   = ro  + frac*(rstar - ro);
        __mvr unosto = uno + frac*(ustar - uno);
        __mvr psto   = po  + frac*(pstar - po);

        // test: spout <= 0
        test_vr = _mm_vr(cmp)(spout, _mm_vr(setzero)(), _CMP_LE_OQ);
        rsto   = _mm_vr(blendv)(rsto  , ro , test_vr);
        unosto = _mm_vr(blendv)(unosto, uno, test_vr);
        psto   = _mm_vr(blendv)(psto  , po , test_vr);

        // test: spin > 0
        test_vr = _mm_vr(cmp)(spin, _mm_vr(setzero)(), _CMP_GT_OQ);
        rsto   = _mm_vr(blendv)(rsto  , rstar, test_vr);
        unosto = _mm_vr(blendv)(unosto, ustar, test_vr);
        psto   = _mm_vr(blendv)(psto  , pstar, test_vr);

//--Store

        if (i0R < i0EndRPacked)  // Can fully store the vector
          {
                   _mm_vr(storeu)(Wstar_ptr                       , rsto);
            D_TERM(_mm_vr(storeu)(Wstar_ptr + WstarCmpStride*inorm, unosto);,
                   _mm_vr(storeu)(Wstar_ptr + WstarCmpStride*itan1, utan1);,
                   _mm_vr(storeu)(Wstar_ptr + WstarCmpStride*itan2, utan2);)
                   _mm_vr(storeu)(Wstar_ptr + WstarCmpStride*WPRES, psto);
          }
        else                   // Partial storing of vector
          {
                   _mm_vr(maskstore)(Wstar_ptr                       , lsmask.m, rsto);
            D_TERM(_mm_vr(maskstore)(Wstar_ptr + WstarCmpStride*inorm, lsmask.m, unosto);,
                   _mm_vr(maskstore)(Wstar_ptr + WstarCmpStride*itan1, lsmask.m, utan1);,
                   _mm_vr(maskstore)(Wstar_ptr + WstarCmpStride*itan2, lsmask.m, utan2);)
                   _mm_vr(maskstore)(Wstar_ptr + WstarCmpStride*WPRES, lsmask.m, psto);
          }
        Wleft_ptr  += VecSz_r;
        Wright_ptr += VecSz_r;
        Wstar_ptr  += VecSz_r;
      }
  EndForPencil
}
#endif  /* defined CRD_USE_VEX */
