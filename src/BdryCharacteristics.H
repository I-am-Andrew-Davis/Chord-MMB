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
 * \file BdryCharacteristics.H
 *
 * \brief Contains the functions for characteristic boundaries
 *//*+*************************************************************************/

#ifndef _BDRYCHARACTERISTICS_H_
#define _BDRYCHARACTERISTICS_H_

//----- Standard Library -----//

#include <string>
#include <vector>

//----- Chombo Library -----//

#include "REAL.H"

//----- Internal -----//

#include "CRDparam.H"

/*==============================================================================
 * Characteristic functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Functions that are used to solve the characteristics
//  FIXME: Functions do not exist for resolving outflow-outflow-outflow
//         or any combination of outflow or inflows in 3D
/** General nomenclature for these functions.
 *  \param[out] a_dqdx  FAB of characteristic derivatives
 *  \param[in]  a_gammaFab
 *                      FAB of gamma values
 *  \param[in]  a_Wface FAB of prescribed inlet values
 *  \param[in]  a_Wcell FAB of cell-averaged values
 *  \param[in]  a_d1FabN
 *                      FAB of normal gradients
 *  \param[in]  a_d1FabT
 *                      FAB of transverse gradients
 *  \param[in]  a_box   Box to solve over, should be first interior box
 *  \param[in]  a_dir   Direction normal to boundary
 *  \param[in]  a_etaM  Smoothing coefficient for u and T inlets
 *  \param[in]  a_etaCN Smoothing coefficient for all mass fraction inlets
 *//*-----------------------------------------------------------------*/

void
CNSCBCInflow(FArrayBox&            a_dqdx,
             const FArrayBox&      a_Wcell,
             const FArrayBox&      a_Wface,
             const FArrayBox&      a_d1FabN,
             const FArrayBox&      a_d1FabT,
             const FArrayBox&      a_gammaFab,
             const Box&            a_box,
             const Real&           a_etaM,
             const Real&           a_etaCN,
             const Real&           a_beta,
             const int             a_dir,
             const Side::LoHiSide& a_side)
{
  const int lohiSign = sign(a_side);
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive(); (void)numWcomp;
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int speciesComp = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int turbWVar = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const Real length = CRDparam::g_domainLength[a_dir];
  CH_assert(a_Wcell.box().contains(a_box));
  std::vector<Real> spec(numSpecies);
  std::vector<Real> inletSpec(numSpecies);
  std::vector<Real> tauv(SpaceDim-1, 0.);
  std::vector<Real> Rn(numSpecies);
  Real origR = CRDparam::g_R;
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      origR = 0.;
      for(int comp = 0; comp != numSpecies; ++comp)
        {
          Rn[comp] = CRDparam::g_CRDPhysics->speciesGasConstant(comp);
        }
    }
  D_TERM(const int idnorm = a_dir;,
         const int idtan1 = (a_dir + 1) % SpaceDim;,
         const int idtan2 = (a_dir + 2) % SpaceDim;);
  D_TERM(const int inorm = velIndx + idnorm;,
         const int itan1 = velIndx + idtan1;,
         const int itan2 = velIndx + idtan2;);
  MD_ARRAY_RESTRICT(arrdqdx, a_dqdx);
  MD_ARRAY_RESTRICT(arrWface, a_Wface);
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_ARRAY_RESTRICT(arrd1FabT, a_d1FabT);
  MD_ARRAY_RESTRICT(arrd1FabN, a_d1FabN);
  MD_ARRAY_RESTRICT(arrGamma, a_gammaFab);
  MD_BOXLOOP(a_box, i)
    {
      D_TERM(Real uu = arrWcell[MD_IX(i, inorm)];,
             Real vv = arrWcell[MD_IX(i, itan1)];,
             Real ww = arrWcell[MD_IX(i, itan2)];);
      D_TERM(Real inletunorm = arrWface[MD_IX(i, inorm)];,
             Real inletutan1 = arrWface[MD_IX(i, itan1)];,
             Real inletutan2 = arrWface[MD_IX(i, itan2)];);
      D_TERM(Real duin = uu - inletunorm;,
             Real dvin = vv - inletutan1;,
             Real dwin = ww - inletutan2;);
      Real rho = arrWcell[MD_IX(i, rhoIndx)];
      Real p = arrWcell[MD_IX(i, presIndx)];
      Real T = arrWcell[MD_IX(i, tempIndx)];
      Real gamma = arrGamma[MD_IX(i, 0)];
      Real Tinlet = arrWface[MD_IX(i, tempIndx)];
      Real dpdx = arrd1FabN[MD_IX(i, presIndx)];
      Real dudx = arrd1FabN[MD_IX(i, inorm)];
      Real dTin = T - Tinlet;
      Real Rval = origR;
      Real Rinlet = origR;
      for(int comp = 0; comp != numSpecies; ++comp)
        {
          int cncomp = comp + speciesComp;
          spec[comp] = arrWcell[MD_IX(i, cncomp)];
          inletSpec[comp] = arrWface[MD_IX(i, cncomp)];
          Rinlet += inletSpec[comp]*Rn[comp];
          Rval += spec[comp]*Rn[comp];
        }
      Real cs2 = p*gamma/rho;
      Real cs = sqrt(cs2);
      Real csin2 = Tinlet*Rinlet*gamma;
      Real csin = sqrt(csin2);
      Real M = (D_TERM(std::abs(inletunorm),+std::abs(inletutan1),
                       +std::abs(inletutan2)))/csin;
      Real lambda2 = inletunorm;
      Real lambdain = inletunorm - lohiSign*csin;
      Real betain = a_beta*M;
      Real etaCN = a_etaCN/length;
      Real etaM = a_etaM/length;
      if (M >= 1.)
        {
          M = 0.;
        }
      for(int comp = 0; comp != numSpecies; ++comp)
        {
          int cncomp = comp + speciesComp;
          Real taucn = D_TERM(0.,+vv*arrd1FabT[MD_IX(i, cncomp)],
                              +ww*arrd1FabT[MD_IX(i, numWcomp + cncomp)]);
          Real cndiff = spec[comp] - inletSpec[comp];
          arrdqdx[MD_IX(i, cncomp)] = etaCN*csin*cndiff/lambda2 - betain*taucn;
        }
      // Add any turbulent variables
      for(int tComp = 0; tComp != numTurbVar; ++tComp)
        {
          int turbComp = tComp + turbWVar;
          Real nutCell = arrWcell[MD_IX(i, turbComp)];
          Real nutFace = arrWface[MD_IX(i, turbComp)];
          arrdqdx[MD_IX(i, turbComp)] = etaCN*(nutCell - nutFace)/lambda2;
        }
      Real tauin = 0.;
      Real tau2 = 0.;
      // FIXME: Ignore transverse aspect of inflow
      // for(int tdir = 0; tdir != SpaceDim - 1; ++tdir)
      //   {
      //     int scomp = tdir*numWcomp;
      //     int vtcomp = velIndx + ((a_dir+1+tdir)%SpaceDim);
      //     Real vc = arrWcell[MD_IX(i, vtcomp)];
      //     Real drdy = arrd1FabT[MD_IX(i, scomp+rhoIndx)];
      //     Real dpdy = arrd1FabT[MD_IX(i, scomp+presIndx)];
      //     Real dudy = arrd1FabT[MD_IX(i, scomp+inorm)];
      //     Real dvtdy = arrd1FabT[MD_IX(i, scomp+vtcomp)];
      //     tauin += vc*(dpdy-lohiSign*rho*cs*dudy)+gamma*p*dvtdy;
      //     tau2 += csin2*(vc*drdy+rho*dvtdy)-vc*dpdy-gamma*p*dvtdy;
      //     tauv[tdir] = vc*dvtdy + dpdy/rho;
      //   }
      Real Lin = -lohiSign*etaM*rho*csin2*(1. - M*M)*duin - betain*tauin;
      Real Lout = (dpdx + lohiSign*rho*cs*dudx);
      D_TERM(
        Real L2 = -etaM*rho*csin*Rval*dTin - betain*tau2;,
        Real L3 = etaM*dvin - betain*tauv[0];,
        Real L4 = etaM*dwin - betain*tauv[1];);
      arrdqdx[MD_IX(i, rhoIndx)] = (L2/lambda2 + (Lout + Lin/lambdain)/2.)/cs2;
      D_TERM(arrdqdx[MD_IX(i, inorm)] = lohiSign/(2.*csin*rho)*
             (Lout - Lin/lambdain);,
             arrdqdx[MD_IX(i, itan1)] = L3/lambda2;,
             arrdqdx[MD_IX(i, itan2)] = L4/lambda2;);
      arrdqdx[MD_IX(i, presIndx)] = (Lin/lambdain + Lout)/2.;
    }
}

void
CNSCBCOutflow(FArrayBox&            a_dqdx,
              const FArrayBox&      a_Wcell,
              const FArrayBox&      a_Wface,
              const FArrayBox&      a_d1FabN,
              const FArrayBox&      a_d1FabT,
              const FArrayBox&      a_gammaFab,
              const Box&            a_box,
              const Real&           a_sigma,
              const Real&           a_beta,
              const int             a_dir,
              const Side::LoHiSide& a_side)
{
  const int lohiSign = sign(a_side);
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int speciesComp = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int turbWVar = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const Real length = CRDparam::g_domainLength[a_dir];
  CH_assert(a_Wcell.box().contains(a_box));
  std::vector<Real> tauv(SpaceDim-1);
  D_TERM(const int idnorm = a_dir;,
         const int idtan1 = (a_dir + 1) % SpaceDim;,
         const int idtan2 = (a_dir + 2) % SpaceDim;);
  D_TERM(const int inorm = velIndx + idnorm;,
         const int itan1 = velIndx + idtan1;,
         const int itan2 = velIndx + idtan2;);
  MD_ARRAY_RESTRICT(arrdqdx, a_dqdx);
  MD_ARRAY_RESTRICT(arrWface, a_Wface);
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_ARRAY_RESTRICT(arrd1FabT, a_d1FabT);
  MD_ARRAY_RESTRICT(arrd1FabN, a_d1FabN);
  MD_ARRAY_RESTRICT(arrGamma, a_gammaFab);
  MD_BOXLOOP(a_box, i)
    {
      D_TERM(Real uu = arrWcell[MD_IX(i, inorm)];,
             Real vv = arrWcell[MD_IX(i, itan1)];,
             Real ww = arrWcell[MD_IX(i, itan2)];);
      Real rho = arrWcell[MD_IX(i, rhoIndx)];
      Real p = arrWcell[MD_IX(i, presIndx)];
      Real gamma = arrGamma[MD_IX(i, 0)];
      Real pOut = arrWface[MD_IX(i, presIndx)];
      Real diffp = p - pOut;
      Real dpdx = arrd1FabN[MD_IX(i, presIndx)];
      Real drdx = arrd1FabN[MD_IX(i, rhoIndx)];
      Real dudx = arrd1FabN[MD_IX(i, inorm)];
      for(int comp = 0; comp != numSpecies; ++comp)
        {
          int cncomp = comp + speciesComp;
          arrdqdx[MD_IX(i, cncomp)] = arrd1FabN[MD_IX(i, cncomp)];
        }
      // Add any turbulent variables
      for(int tComp = 0; tComp != numTurbVar; ++tComp)
        {
          int turbComp = tComp + turbWVar;
          arrdqdx[MD_IX(i, turbComp)] = arrd1FabN[MD_IX(i, turbComp)];
        }
      Real cs2 = p*gamma/rho;
      Real cs = sqrt(cs2);
      Real lambdain = uu - lohiSign*cs;
      // Fixed the mach number calculation here.
      Real M = std::sqrt(D_TERM(uu*uu,+(vv*vv),+(ww*ww)))/cs;
      Real betaout = a_beta*(1. - M);
      if (M >= 1.)
        {
          betaout = 0.;
        }
      Real taupr = 0.;
      for(int tdir = 0; tdir != SpaceDim - 1; ++tdir)
        {
          int scomp = tdir*numWcomp;
          int vtcomp = velIndx + ((a_dir+1+tdir)%SpaceDim);
          Real vc = arrWcell[MD_IX(i, vtcomp)];
          Real dpdy = arrd1FabT[MD_IX(i, scomp+presIndx)];
          Real dudy = arrd1FabT[MD_IX(i, scomp+inorm)];
          Real dvtdy = arrd1FabT[MD_IX(i, scomp+vtcomp)];
          taupr += vc*(dpdy - lohiSign*rho*cs*dudy)+gamma*p*dvtdy;
        }
      Real Kval = a_sigma*cs*(1. - M*M)/length;
      Real Lin_lin = (Kval*diffp - betaout*taupr)/lambdain;
      Real Lout = (dpdx + lohiSign*rho*cs*dudx);
      D_TERM(Real L2 = (cs2*drdx - dpdx);,
             Real L3 = arrd1FabN[MD_IX(i, itan1)];,
             Real L4 = arrd1FabN[MD_IX(i, itan2)];);
      arrdqdx[MD_IX(i, rhoIndx)] = (L2 + (Lout + Lin_lin)/2.)/cs2;
      D_TERM(arrdqdx[MD_IX(i, inorm)] = lohiSign/(2.*cs*rho)*(Lout-Lin_lin);,
             arrdqdx[MD_IX(i, itan1)] = L3;,
             arrdqdx[MD_IX(i, itan2)] = L4;);
      arrdqdx[MD_IX(i, presIndx)] = (Lin_lin + Lout)/2.;
    }
}

void
CNSCBCOutflowOutflow(FArrayBox&            a_Wcell,
                     FArrayBox&            a_Wface1,
                     FArrayBox&            a_Wface2,
                     const FArrayBox&      a_d1FabN1,
                     const FArrayBox&      a_d1FabN2,
                     const FArrayBox&      a_gammaFab,
                     const Box&            a_box,
                     const Real&           a_sigma,
                     const Real&           a_beta,
                     const int             a_dir1,
                     const int             a_dir2,
                     const Side::LoHiSide& a_side1,
                     const Side::LoHiSide& a_side2,
                     const RealVect&       a_dxVect)
{
  Real sigma = a_sigma;
  const int numSpecies = CRDparam::g_numSpecies;
  Interval velInt = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx = velInt.begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  Interval specInt = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int speciesIndx = specInt.begin();
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int turbWVar = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int numWpcomp = CRDparam::g_CRDPhysics->numNativePrimitive();
#if CH_SPACEDIM==3
  const int dir3 = SpaceDim - (a_dir1 + a_dir2);
  const int uTanDir = velIndx + dir3;
#endif
  const Real dxN1 = a_dxVect[a_dir1];
  const Real dxN2 = a_dxVect[a_dir2];
  const Real Lx1 = CRDparam::g_domainLength[a_dir1];
  const Real Lx2 = CRDparam::g_domainLength[a_dir2];
  std::vector<Real> cn(numSpecies);
  std::vector<Real> dqdx1(numWpcomp);
  std::vector<Real> dqdx2(numWpcomp);
  const int uNormIndx1 = velIndx + a_dir1;
  const int uNormIndx2 = velIndx + a_dir2;
  const int lohiSign1 = sign(a_side1);
  const int lohiSign2 = sign(a_side2);
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_ARRAY_RESTRICT(arrGamma, a_gammaFab);
  MD_ARRAY_RESTRICT(arrWface1, a_Wface1);
  MD_ARRAY_RESTRICT(arrWface2, a_Wface2);
  MD_ARRAY_RESTRICT(arrd1FabN1, a_d1FabN1);
  MD_ARRAY_RESTRICT(arrd1FabN2, a_d1FabN2);
  const int MD_ID(o, a_dir1);
  const int MD_ID(q, a_dir2);
  MD_BOXLOOP(a_box, i)
    {
      Real p = arrWcell[MD_IX(i, presIndx)];
      // Target pressure, should be the same on both faces
      Real pOut = arrWface1[MD_IX(i, presIndx)];
      Real diffp = (p - pOut);
      Real rho = arrWcell[MD_IX(i, rhoIndx)];
      D_TERM(Real uu = arrWcell[MD_IX(i, uNormIndx1)];,
             Real vv = arrWcell[MD_IX(i, uNormIndx2)];,
             Real ww = arrWcell[MD_IX(i, dir3)];);
      Real dpdx1 = arrd1FabN1[MD_IX(i, presIndx)];
      Real drdx1 = arrd1FabN1[MD_IX(i, rhoIndx)];
      Real dudx1 = arrd1FabN1[MD_IX(i, uNormIndx1)];
      Real dvdx1 = arrd1FabN1[MD_IX(i, uNormIndx2)];
      Real dpdx2 = arrd1FabN2[MD_IX(i, presIndx)];
      Real drdx2 = arrd1FabN2[MD_IX(i, rhoIndx)];
      Real dudx2 = arrd1FabN2[MD_IX(i, uNormIndx2)];
      Real dvdx2 = arrd1FabN2[MD_IX(i, uNormIndx1)];
      Real gamma = arrGamma[MD_IX(i, 0)];
      Real cs = sqrt(p*gamma/rho);
      Real cs2 = p*gamma/rho;
      Real Ms = (D_TERM(uu, + vv, + ww))/cs;
      Real beta1 = a_beta*(1. - Ms);
      Real factor1 = 1./(beta1*beta1 - 4.);
      Real lam1PR = uu + lohiSign1*cs;
      Real lam2PR = vv + lohiSign2*cs;
      Real lam1Phi = uu - lohiSign1*cs;
      Real lam2Phi = vv - lohiSign2*cs;
      Real L2 = (cs2*drdx1 - dpdx1);
      Real M3 = (cs2*drdx2 - dpdx2);
      Real LPStar = lam1PR*(dpdx1 + lohiSign1*rho*cs*dudx1);
      Real MPStar = lam2PR*(dpdx2 + lohiSign2*rho*cs*dudx2);
      Real DL = beta1*(LPStar/2. - lohiSign2*rho*cs*uu*dvdx1);
      Real DM = beta1*(MPStar/2. - lohiSign1*rho*cs*vv*dvdx2);
      Real dTerm = sigma*cs*(1. - Ms*Ms)*diffp;
      Real LPR = ((2.*Lx1*beta1 - 4.*Lx2)*dTerm/(Lx1*Lx2)
                  - 2.*DL*beta1 + 4.*DM)*factor1;
      Real MPR = ((2.*Lx2*beta1 - 4.*Lx1)*dTerm/(Lx1*Lx2)
                  - 2.*DM*beta1 + 4.*DL)*factor1;
      dqdx1[rhoIndx] = (L2+(LPStar/lam1PR + LPR/lam1Phi)/2.)/cs2;
      dqdx2[rhoIndx] = (M3+(MPStar/lam2PR + MPR/lam2Phi)/2.)/cs2;
      dqdx1[presIndx] = (LPStar/lam1PR + LPR/lam1Phi)/2.;
      dqdx2[presIndx] = (MPStar/lam2PR + MPR/lam2Phi)/2.;
      D_TERM(dqdx1[uNormIndx1] = lohiSign1/(2.*rho*cs)*
             (LPStar/lam1PR - LPR/lam1Phi);
             dqdx2[uNormIndx2] = lohiSign2/(2.*rho*cs)*
             (MPStar/lam2PR - MPR/lam2Phi);,
             dqdx1[uNormIndx2] = dvdx1;
             dqdx2[uNormIndx1] = dvdx2;,
             // FIXME: Add third dimension
             dqdx1[uTanDir] = 0.;
             dqdx2[uTanDir] = 0.;);
      for(int spec = 0; spec != numSpecies; ++spec)
        {
          const int wIndx = speciesIndx + spec;
          dqdx1[wIndx] = arrd1FabN1[MD_IX(i, wIndx)];
          dqdx2[wIndx] = arrd1FabN2[MD_IX(i, wIndx)];
        }
      // Add any turbulent variables
      for(int tComp = 0; tComp != numTurbVar; ++tComp)
        {
          int turbComp = tComp + turbWVar;
          dqdx1[turbComp] = arrd1FabN1[MD_IX(i, turbComp)];
          dqdx2[turbComp] = arrd1FabN2[MD_IX(i, turbComp)];
        }
      // Calculate boundary face values and ghost cell values
      for(int comp = 0; comp != numWpcomp; ++comp)
        {
          Real Wj = arrWcell[MD_IX(i,comp)];
          Real Wjp = arrWcell[MD_OFFSETIX(i,-,lohiSign1*o,comp)];
          Real Wjpp = arrWcell[MD_OFFSETIX(i,-,lohiSign1*2*o,comp)];
          Real Wjm = -(Wjpp - 6.*Wjp+3.*Wj - lohiSign1*6.*dxN1*dqdx1[comp])/2.;
          Real Wjmm = -2.*Wjp - 3.*Wj + 6.*Wjm
            - lohiSign1*6.*dxN1*dqdx1[comp];
          arrWcell[MD_OFFSETIX(i,+,lohiSign1*o,comp)] = Wjm;
          arrWcell[MD_OFFSETIX(i,+,lohiSign1*2*o,comp)] = Wjmm;
          arrWface1[MD_IX(i,comp)] = (7.*(Wjm + Wj) - (Wjmm + Wjp))/12.;
        }
      for(int comp = 0; comp != numWpcomp; ++comp)
        {
          Real Wj = arrWcell[MD_IX(i,comp)];
          Real Wjp = arrWcell[MD_OFFSETIX(i,-,lohiSign2*q,comp)];
          Real Wjpp =
            arrWcell[MD_OFFSETIX(i,-,lohiSign2*2*q,comp)];
          Real Wjm = -(Wjpp - 6.*Wjp + 3.*Wj
                       - lohiSign2*6.*dxN2*dqdx2[comp])/2.;
          Real Wjmm = -2.*Wjp - 3.*Wj + 6.*Wjm
            - lohiSign2*6.*dxN2*dqdx2[comp];
          arrWcell[MD_OFFSETIX(i,+,lohiSign2*q,comp)] = Wjm;
          arrWcell[MD_OFFSETIX(i,+,lohiSign2*2*q,comp)] = Wjmm;
          arrWface2[MD_IX(i,comp)] = (7.*(Wjm + Wj) - (Wjmm + Wjp))/12.;
        }
    }
}

void
CNSCBCInflowOutflow(FArrayBox&                    a_Wcell,
                    FArrayBox&                    a_Wface1,
                    FArrayBox&                    a_Wface2,
                    const FArrayBox&              a_d1FabN1,
                    const FArrayBox&              a_d1FabN2,
                    const FArrayBox&              a_gammaFab,
                    const Box&                    a_box,
                    const Real&                   a_sigma,
                    const Real&                   a_etaMax,
                    const Real&                   a_etaCN,
                    const int                     a_dir1,
                    const int                     a_dir2,
                    const Side::LoHiSide&         a_side1,
                    const Side::LoHiSide&         a_side2,
                    const CRDparam::DomainBCType& a_domT1,
                    const CRDparam::DomainBCType& a_domT2,
                    const RealVect&               a_dxVect)
{
  const int numSpecies = CRDparam::g_numSpecies;
  Interval velInt = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx = velInt.begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  Interval specInt = CRDparam::g_CRDPhysics->speciesPrimInterval();
  const int speciesIndx = specInt.begin();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int numWpcomp = CRDparam::g_CRDPhysics->numNativePrimitive();
  const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
  const int turbWVar = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
#if CH_SPACEDIM==3
  const int dir3 = SpaceDim - (a_dir1 + a_dir2);
  const int uTanDir = velIndx + dir3;
#endif
  const Real dxN1 = a_dxVect[a_dir1];
  const Real dxN2 = a_dxVect[a_dir2];
  const Real Lx1 = CRDparam::g_domainLength[a_dir1];
  const Real Lx2 = CRDparam::g_domainLength[a_dir2];
  std::vector<Real> cn(numSpecies);
  std::vector<Real> dqdx1(numWpcomp);
  std::vector<Real> dqdx2(numWpcomp);
  const int lohiSign1 = sign(a_side1);
  const int lohiSign2 = sign(a_side2);
  int inflowDir = a_dir1, outflowDir = a_dir2;
  Real LxIn = Lx1;
  Real dxIn = dxN1, dxOut = dxN2;
  int inSign = lohiSign1, outSign = lohiSign2;
  FArrayBox* inflowFab; FArrayBox* outflowFab;
  const FArrayBox* d1Inflow; const FArrayBox* d1Outflow;
  // Assign inflow or outflow conditions according to directions
  if (a_domT1 & CRDparam::DomainBCTypeCNSCBCInflow)
    {
      inflowFab = &a_Wface1; outflowFab = &a_Wface2;
      d1Inflow = &a_d1FabN1; d1Outflow = &a_d1FabN2;
    }
  else
    {
      inflowFab = &a_Wface2; outflowFab = &a_Wface1;
      d1Inflow = &a_d1FabN2; d1Outflow = &a_d1FabN1;
      inflowDir = a_dir2; outflowDir = a_dir1;
      LxIn = Lx2;
      dxIn = dxN2; dxOut = dxN1;
      inSign = lohiSign2; outSign = lohiSign1;
    }
  MD_ARRAY_RESTRICT(arrWcell, a_Wcell);
  MD_ARRAY_RESTRICT(arrGamma, a_gammaFab);
  MD_ARRAY_RESTRICT(arrInflow, *inflowFab);
  MD_ARRAY_RESTRICT(arrOutflow, *outflowFab);
  MD_ARRAY_RESTRICT(arrd1FabIn, *d1Inflow);
  MD_ARRAY_RESTRICT(arrd1FabOut, *d1Outflow);
  const int MD_ID(o, inflowDir);
  const int MD_ID(q, outflowDir);
  const int uNormIn = velIndx + inflowDir;
  const int uNormOut = velIndx + outflowDir;
  std::vector<Real> Rn(numSpecies);
  for(int spec = 0; spec != numSpecies; ++spec)
    {
      Rn[spec] = CRDparam::g_CRDPhysics->speciesGasConstant(spec);
    }
  std::vector<Real> inletCn(numSpecies);
  MD_BOXLOOP(a_box, i)
    {
      Real p = arrWcell[MD_IX(i, presIndx)];
      Real T = arrWcell[MD_IX(i, tempIndx)];
      Real inletT = arrInflow[MD_IX(i, tempIndx)];
      Real rho = arrWcell[MD_IX(i, rhoIndx)];
      D_TERM(Real uu = arrWcell[MD_IX(i, uNormIn)];,
             Real vv = arrWcell[MD_IX(i, uNormOut)];,
             Real ww = arrWcell[MD_IX(i, dir3)]);
      Real inletU = arrInflow[MD_IX(i, uNormIn)];
      Real inletV = arrInflow[MD_IX(i, uNormOut)];
      Real dpdx1 = arrd1FabIn[MD_IX(i, presIndx)];
      Real dudx1 = arrd1FabIn[MD_IX(i, uNormIn)];
      Real dpdx2 = arrd1FabOut[MD_IX(i, presIndx)];
      Real drdx2 = arrd1FabOut[MD_IX(i, rhoIndx)];
      Real dudx2 = arrd1FabOut[MD_IX(i, uNormOut)];
      Real dvdx2 = arrd1FabOut[MD_IX(i, uNormIn)];
      Real gamma = arrGamma[MD_IX(i, 0)];
      Real Rval = 0.;
      Real inletR = 0.;
      for(int spec = 0; spec != numSpecies; ++spec)
        {
          const int wIndx = speciesIndx + spec;
          inletCn[spec] = arrInflow[MD_IX(i, wIndx)];
          cn[spec] = arrWcell[MD_IX(i, wIndx)];
          Rval += Rn[spec]*cn[spec];
          inletR += Rn[spec]*inletCn[spec];
        }
      Real cs2 = p*gamma/rho;
      Real cs = sqrt(cs2);
      Real Ms = (D_TERM(inletU, + vv, + ww))/cs;
      Real etaCN = a_etaCN;
      Real etaM = a_etaMax;
      Real lamInPR = inletU + inSign*cs;
      Real lamOutPR = vv + outSign*cs;
      Real lamInPhi = inletU - inSign*cs;
      Real lamOutPhi = vv - outSign*cs;
      // FIXME: Make these transverse terms non-zero for 3D
      Real Tau1 = 0., Tau2 = 0., Tau3 = 0., Tau5 = 0.,Tau6 = 0.;
      Real M2 = dvdx2;
      Real M3 = (cs2*drdx2 - dpdx2);
      Real L2 = -etaM*rho*cs*Rval/LxIn*(T - inletT) - vv*M3 + (cs2*Tau1 - Tau5);
      Real MPStar = lamOutPR*(dpdx2 + outSign*rho*cs*dudx2);
      Real LPStar = lamInPR*(dpdx1 + inSign*rho*cs*dudx1);
      Real LPR = -inSign*etaM*rho*cs2*(1. - Ms*Ms)/LxIn*
        (uu - inletU) - MPStar/2. + Tau5 + inSign*rho*cs*(Tau2-vv*M2);
      // Compatibility condition
      Real MPR = 0.;
      Real L3 = etaM*cs/LxIn*(vv - inletV) + outSign/(2.*rho*cs)*MPStar + Tau3;
      dqdx1[rhoIndx] = (L2/inletU + (LPStar/lamInPR+LPR/lamInPhi)/2.)/cs2;
      dqdx2[rhoIndx] = (M3+(MPStar/lamOutPR + MPR/lamOutPhi)/2.)/cs2;
      dqdx1[presIndx] = (LPStar/lamInPR + LPR/lamInPhi)/2.;
      dqdx2[presIndx] = (MPStar/lamOutPR + MPR/lamOutPhi)/2.;
      D_TERM(dqdx1[uNormIn] = inSign/(2.*rho*cs)*
             (LPStar/lamInPR - LPR/lamInPhi);
             dqdx2[uNormOut] = outSign/(2.*rho*cs)*
             (MPStar/lamOutPR - MPR/lamOutPhi);,
             dqdx1[uNormOut] = L3/inletU;
             dqdx2[uNormIn] = dvdx2;,
             // FIXME: update this, add 3rd dimension
             dqdx1[uTanDir] = 0.;
             dqdx2[uTanDir] = 0.;);
      for(int spec = 0; spec != numSpecies; ++spec)
        {
          const int wIndx = speciesIndx + spec;
          Real dcndx2 = arrd1FabOut[MD_IX(i, wIndx)];
          dqdx1[wIndx] = (etaCN*cs/LxIn*(cn[spec] - inletCn[spec])
                          - vv*dcndx2 + Tau6)/inletU;
        }
      // Add any turbulent variables
      for(int tComp = 0; tComp != numTurbVar; ++tComp)
        {
          int turbComp = tComp + turbWVar;
          Real dnutdx2 = arrd1FabOut[MD_IX(i, turbComp)];
          Real cellNut = arrWcell[MD_IX(i, turbComp)];
          Real inletNut = arrInflow[MD_IX(i, turbComp)];
          dqdx1[turbComp] = (etaCN*cs/LxIn*(cellNut - inletNut)
                             - vv*dnutdx2 + Tau6)/inletU;
        }
      // Calculate boundary face values and ghost cell values
      for(int comp = 0; comp != numWpcomp; ++comp)
        {
          Real Wj = arrWcell[MD_IX(i,comp)];
          Real Wjp = arrWcell[MD_OFFSETIX(i,-,inSign*o,comp)];
          Real Wjpp = arrWcell[MD_OFFSETIX(i,-,inSign*2*o,comp)];
          Real Wjm = -(Wjpp - 6.*Wjp + 3.*Wj
                       - inSign*6.*dxIn*dqdx1[comp])/2.;
          Real Wjmm = -2.*Wjp - 3.*Wj + 6.*Wjm
            - inSign*6.*dxIn*dqdx1[comp];
          arrWcell[MD_OFFSETIX(i,+,inSign*o,comp)] = Wjm;
          arrWcell[MD_OFFSETIX(i,+,inSign*2*o,comp)] = Wjmm;
          arrInflow[MD_IX(i,comp)] = (7.*(Wjm + Wj) - (Wjmm + Wjp))/12.;
        }
      for(int comp = 0; comp != numWpcomp; ++comp)
        {
          Real Wj = arrWcell[MD_IX(i,comp)];
          Real Wjp = arrWcell[MD_OFFSETIX(i,-,outSign*q,comp)];
          Real Wjpp = arrWcell[MD_OFFSETIX(i,-,outSign*2*q,comp)];
          Real Wjm = -(Wjpp - 6.*Wjp + 3.*Wj
                       - outSign*6.*dxOut*dqdx2[comp])/2.;
          Real Wjmm =  -2.*Wjp - 3.*Wj + 6.*Wjm
            - outSign*6.*dxOut*dqdx2[comp];
          arrWcell[MD_OFFSETIX(i,+,outSign*q,comp)] = Wjm;
          arrWcell[MD_OFFSETIX(i,+,outSign*2*q,comp)] = Wjmm;
          arrOutflow[MD_IX(i,comp)] = (7.*(Wjm + Wj)- (Wjmm + Wjp))/12.;
        }
    }
}
#endif /* ! defined _BDRYCHARACTERISTICS_H_ */
