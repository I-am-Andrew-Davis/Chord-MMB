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
 * \file SA.cpp
 *
 * \brief Member functions for SA
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "FArrayBox.H"
#include "ParmParse.H"
#include "RootSolver.H"
#include "UnitNormalsF_F.H"

//----- Internal -----//

#include "SA.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "DataTemp.H"
#include "CNSIBC.H"
#include "PatchMappedFunc.H"
#include "ViscousTensor4thOrderOpF_F.H"

/*******************************************************************************
 *
 * Class SA: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

SA::SA()
  :
  m_Cb1(0.1355),
  m_Cb2(0.622),
  m_Cw1(0.),
  m_sigma(2./3.),
  m_k(0.41),
  m_Cw2(0.3),
  m_Cw3(2.),
  m_Cv1(7.1),
  m_Ct3(1.2),
  m_Ct4(0.5),
  m_Cprod(2.),
  m_Prt(0.71),
  m_initTurb(-1.),
  m_inletTurb(-1.),
  m_refLen(-1.),
  m_yPlusLam(11.0),
  m_E(9.)
{
  m_Cw1 = m_Cb1/(m_k*m_k) + (1. + m_Cb2)/m_sigma;
  // Read in turbulent input parameters
  ParmParse ppTURB("turb");
  ppTURB.get("init_turb_intensity", m_initTurb);
  if(m_initTurb < 0. || m_initTurb > 1.)
    {
      CRD::msg << "Input: 'init_turb_intensity' must be between 0 and 1!"
               << CRD::error;
    }
  m_inletTurb = m_initTurb;
  ppTURB.query("inlet_turb_intensity", m_inletTurb);
  if(m_inletTurb < 0. || m_inletTurb > 1.)
    {
      CRD::msg << "Input: 'inlet_turb_intensity' must be between 0 and 1!"
               << CRD::error;
    }
  ppTURB.get("ref_length", m_refLen);
  if(m_refLen < 0.)
    {
      CRD::msg << "Input: 'ref_length' must be > 0!" << CRD::error;
    }
  ppTURB.query("prandtl_num", m_Prt);
  if(m_Prt < 0.)
    {
      CRD::msg << "Input: 'prandtl_num' must be > 0!" << CRD::error;
    }
  // Calculate the laminar y+
  for (int i = 0; i != 10; ++i)
    {
      m_yPlusLam = std::log(std::max(m_E*m_yPlusLam, 1.))/m_k;
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

SA::~SA()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Solve for the turbulent thermal conductivity and dynamic viscosity
/** \param[in]  a_box  Cell box.  Flux needs to be computed on
 *                     surrounding faces.
 *  \param[out] a_mutFab
 *                     FAB of dynamic eddy viscosity
 *  \param[out] a_kappatFab
 *                     FAB of turbulent thermal conductivity
 *  \param[in]  a_muFab
 *                     FAB of dynamic viscosity
 *  \param[in]  a_kappaFab
 *                     FAB of thermal conductivity
 *  \param[in]  a_WfacePntFab
 *                     FAB containing the face-centered primitive variables
 *  \param[in]  a_strainRateFab
 *                     FAB of face-centered strain rate tensor
 *                     This has SpaceDim^2 components when being used and
 *                     only 1 component when serving as a placeholder
 *                     Retrieve the components using 
 *                     ViscousTensor4thOrderOp::tensorIdxRowOrder(row, col)
 *  \param[in]  a_gridMetrics
 *                     Level grid metrics
 *  \param[in]  a_dir  Current direction of solution
 *//*-----------------------------------------------------------------*/

void
SA::calcCoeffKappatMut(const Box&              a_box,
                       FArrayBox&              a_mutFab,
                       FArrayBox&              a_kappatFab,
                       const FArrayBox&        a_muFab,
                       const FArrayBox&        a_kappaFab,
                       const FArrayBox&        a_WfacePntFab,
                       const FArrayBox&        a_strainRateFab,
                       const LevelGridMetrics& a_gridMetrics,
                       const int               a_dir) const
{
  CH_TIME("SA::calcCoeffKappatMut");
  CRD::msg << CRD::fv4 << "SA::calcCoeffKappatMut" << CRD::end;

  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  Real Cv13 = std::pow(m_Cv1, 3);
  std::vector<Real> cn(CRDparam::g_numSpecies, 0.);
  MD_ARRAY_RESTRICT(arrW, a_WfacePntFab);
  MD_ARRAY_RESTRICT(arrMu, a_muFab);
  MD_ARRAY_RESTRICT(arrMut, a_mutFab);
  MD_ARRAY_RESTRICT(arrKappat, a_kappatFab);

  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real T = arrW[MD_IX(i, tempIndx)];
      Real nut = arrW[MD_IX(i, nutIndx)];
      Real mu = arrMu[MD_IX(i, 0)];
      Real chi = rho*nut/mu;
      Real chi3 = std::pow(chi, 3);
      Real fv1 = chi3/(chi3 + Cv13);
      for (int comp = 0; comp != CRDparam::g_numSpecies; ++comp)
        {
          const int wComp = comp + wCompStart;
          cn[comp] = arrW[MD_IX(i, wComp)];
        }
      Real Cp = CRDparam::g_CRDPhysics->cp(T, cn.data());
      Real mut = rho*nut*fv1;
      arrMut[MD_IX(i, 0)] = mut;
      arrKappat[MD_IX(i, 0)] = mut*Cp/m_Prt;
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Compute the intertial flux from primitive variable values on a face
/** \param[out] a_flux  Flux on the faces
 *  \param[in]  a_WFace Primitive state on the face
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_box   Box on which to compute flux
 *//*-----------------------------------------------------------------*/

void
SA::addTurbConvectiveFlux(FArrayBox&       a_flux,
                          const FArrayBox& a_WFace,
                          const int&       a_dir,
                          const Box&       a_box) const
{
  CH_TIME("SA::addTurbConvectiveFlux");
  CRD::msg << CRD::fv4 << "SA::addTurbConvectiveFlux" << CRD::end;

  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int rhoNutIndx = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  MD_ARRAY_RESTRICT(arrW, a_WFace);
  MD_ARRAY_RESTRICT(arrFlux, a_flux);

  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real uu = arrW[MD_IX(i, velIndx + a_dir)];
      Real nut = arrW[MD_IX(i, nutIndx)];
      arrFlux[MD_IX(i, rhoNutIndx)] = rho*uu*nut;
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Calculates the final diffusive flux for a single turbulent state
/** \param[in]  a_box  Box of faces to solve over
 *  \param[in]  a_NGradNufacePntFab
 *                     FAB of the face-centered mapped gradients 
 *                     of turbulent variables
 *  \param[out] a_NGradNufacePntFab
 *                     FAB of gradients modified accordingly
 *  \param[in]  a_muFab
 *                     FAB of the dynamic viscosity on each face
 *  \param[in]  a_kappaFab
 *                     FAB of the thermal conductivity values
 *  \param[in]  a_WfacePntFab
 *                     FAB containing face-centered primitive variables
 *  \param[in]  a_dir  Face direction
 *  \param[in]  a_wTurbComp
 *                     Primitive component of the current turbulent state
 *  \param[in]  a_tComp
 *                     Gradient component for use in a_NGradNufacePntFab
 *//*-----------------------------------------------------------------*/

void
SA::calcTurbDiffFlux(const Box&       a_box,
                     FArrayBox&       a_NGradNufacePntFab,
                     const FArrayBox& a_muFab,
                     const FArrayBox& a_kappaFab,
                     const FArrayBox& a_WfacePntFab,
                     const int        a_dir,
                     const int        a_wTurbComp,
                     const int        a_tComp) const
{
  CH_TIME("SA::calcTurbDiffFlux");
  CRD::msg << CRD::fv4 << "SA::calcTurbDiffFlux" << CRD::end;
  // Modify a_NGradNufacePntFab by (mu + mut)/sigma
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      const int curComp = a_tComp + dir;
      a_NGradNufacePntFab.mult(a_muFab, a_box, 0, curComp, 1);
    }
  a_NGradNufacePntFab.mult(1./m_sigma, a_box, a_tComp, SpaceDim);
}

/*--------------------------------------------------------------------*/
//  Calculates the turbulent source terms
/** \param[in]  a_box  Box of cells to solve cell-centered source terms
 *  \param[in]  a_disjointBox
 *                     Disjoint box
 *  \param[out] a_turbSourcePntFab
 *                     FAB of mapped cell-centered turbulent source terms
 *  \param[in]  a_GradWcellPntFab
 *                     FAB of cell-centered gradients in physical space
 *  \param[in]  a_WcellPntFab
 *                     FAB containing cell-centered primitive variables
 *  \param[in]  a_dataIndx
 *                     Current data index
 *  \param[in]  a_gridMetrics
 *                     Level grid metrics
 *//*-----------------------------------------------------------------*/

void
SA::calcTurbSourceTerms(const Box&              a_box,
                        const Box&              a_disjointBox,
                        FArrayBox&              a_turbSourcePntFab,
                        const FArrayBox&        a_GradWcellPntFab,
                        const FArrayBox&        a_WcellPntFab,
                        const DataIndex&        a_dataIndx,
                        const LevelGridMetrics& a_gridMetrics) const
{
  CH_TIME("SA::calcTurbSourceTerms");
  CRD::msg << CRD::fv4 << "SA::calcTurbSourceTerms" << CRD::end;
  // Source terms exist only for the turbulent viscosity governing equation
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int nuUIndx = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  // Constants used in box loop
  Real Cw36 = std::pow(m_Cw3, 6);
  Real Cv13 = std::pow(m_Cv1, 3);
  Real k2 = m_k*m_k;
  Real sixth = 1./6.;
  // Get physical coordinates
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  FABSTACKTEMP(XiFab, a_box, SpaceDim);  // Cartesian coordinates
  CRDparam::g_CNSIBC->getCellCompCoordinates(a_box, XiFab, blockCoordSys);
  // Find the cell-centered dynamic viscosity
  FABSTACKTEMP(muFab, a_box, 1);
  FABSTACKTEMP(kappaFab, a_box, 1);
  CRDparam::g_CRDPhysics->calcCoeffKappaMu(a_box,
                                           muFab,
                                           kappaFab,
                                           a_WcellPntFab);
  // make sure source starts with zero
  a_turbSourcePntFab.setVal(0., a_box, 0, a_turbSourcePntFab.nComp()-1);
  // get the wall distances
  CH_assert(a_gridMetrics.distanceDefined());
  const FArrayBox& distField = a_gridMetrics.m_distance[a_dataIndx];
  CH_assert(distField.box().contains(a_box));
  
  MD_ARRAY_RESTRICT(arrWCell, a_WcellPntFab);
  MD_ARRAY_RESTRICT(arrGradW, a_GradWcellPntFab);
  MD_ARRAY_RESTRICT(arrMu, muFab);
  MD_ARRAY_RESTRICT(arrXi, XiFab);
  MD_ARRAY_RESTRICT(arrDist, distField);
  MD_ARRAY_RESTRICT(arrTSource, a_turbSourcePntFab);
  MD_BOXLOOP(a_box, i)
    {
      RealVect mappedLoc(D_DECL(arrXi[MD_IX(i, 0)],
                                arrXi[MD_IX(i, 1)],
                                arrXi[MD_IX(i, 2)]));
      Real jVal = blockCoordSys.pointwiseJ(mappedLoc);
      // Distance from nearest wall
      Real dis = arrDist[MD_IX(i, 0)];
      Real rho = arrWCell[MD_IX(i, rhoIndx)];
      // mu
      Real mu = arrMu[MD_IX(i, 0)];
      // ~nu
      Real nut = arrWCell[MD_IX(i, nutIndx)];
      Real chi = rho*nut/mu;
      Real chi3 = std::pow(chi, 3);
      // fv1
      Real fv1 = chi3/(chi3 + Cv13);
      Real fv2 = 1. - chi/(1. + chi*fv1);
      // |S_{ij}|
      Real SMag = 0.;
      // |omega_{ij}|
      Real OMag = 0.;
      // grad nut . grad nut
      Real gradnut2 = 0.;
      // grad rho. grad nut
      Real gradrhodnut = 0.;
      for (int dir1 = 0; dir1 != SpaceDim; ++dir1)
        {
          const int vel1Comp = velIndx + dir1;
          for (int dir2 = 0; dir2 != SpaceDim; ++dir2)
            {
              const int vel2Comp = velIndx + dir2;
              const int vel1Grad2 =
                ViscousTensor4thOrderOp::gradIdx(vel1Comp, dir2);
              const int vel2Grad1 =
                ViscousTensor4thOrderOp::gradIdx(vel2Comp, dir1);
              Real du1dx2 = arrGradW[MD_IX(i, vel1Grad2)];
              Real du2dx1 = arrGradW[MD_IX(i, vel2Grad1)];
              // Solve 2(1/2(du1dx2 + du2dx1))^2
              SMag += 0.5*(du1dx2 + du2dx1)*(du1dx2 + du2dx1);
              if (dir1 != dir2)
                {
                  // Solve 2(1/2(du1dx2 - du2dx1))^2
                  OMag += 0.5*(du1dx2 - du2dx1)*(du1dx2 - du2dx1);
                }
            }
          const int gradnutComp =
            ViscousTensor4thOrderOp::gradIdx(nutIndx, dir1);
          const int gradrhoComp =
            ViscousTensor4thOrderOp::gradIdx(rhoIndx, dir1);
          Real gradnutdir1 = arrGradW[MD_IX(i, gradnutComp)];
          Real gradrhodir1 = arrGradW[MD_IX(i, gradrhoComp)];
          gradnut2 += gradnutdir1*gradnutdir1;
          gradrhodnut += gradnutdir1*gradrhodir1;
        }
      SMag = std::sqrt(SMag);
      OMag = std::sqrt(OMag);
      Real ft2 = m_Ct3*std::exp(-m_Ct4*chi*chi);
      Real SVal = OMag + m_Cprod*std::min(0., SMag - OMag);
      Real destR = nut/(k2*dis*dis);
      Real Stilde = SVal + destR*fv2;
      // Production
      Real prodG = m_Cb1*rho*(1. - ft2)*Stilde*nut;
      Real rval = destR/Stilde;
      Real fw = 2.00517;
      if (rval < 1.8 && Stilde != 0)
        {
          Real gval = rval + m_Cw2*(std::pow(rval, 6) - rval);
          Real tVar = (1. + Cw36)/(std::pow(gval, 6) + Cw36);
          fw = gval*std::pow(tVar, sixth);
        }
      // Destruction
      Real destY = -rho*(m_Cw1*fw - m_Cb1/k2*ft2)*std::pow(nut/dis, 2); //-rho*(m_Cw1*fw);
      // Source 1
      Real source1 = rho*m_Cb2/m_sigma*gradnut2;
      // Source 2
      Real source2 = -(mu/rho + nut)/m_sigma*gradrhodnut;
      // Final source values
      Real finalSource = prodG + destY + source1 + source2;
      // Multiply the source term to get J S
      arrTSource[MD_IX(i, nuUIndx)] = jVal*finalSource;
    }
}

/*--------------------------------------------------------------------*/
//  Initialize the turbulent variables in the flow field
/** \param[out] a_U    Cell-centered conservative variables
 *  \param[in]  a_W    Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                     Current disjointBox
 *  \param[in]  a_box  Box to initialize over
 *//*-----------------------------------------------------------------*/

void
SA::turbInitialize(FArrayBox&              a_U,
                   const FArrayBox&        a_W,
                   const LevelGridMetrics& a_gridMetrics,
                   const FluxBox&          a_unitNormals,
                   const DataIndex&        a_didx,
                   const Box&              a_disjointBox,
                   const Box&              a_box) const
{
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int rhoNutIndx = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  D_TERM(
    const int velComp1 = velIndx;,
    const int velComp2 = velIndx + 1;,
    const int velComp3 = velIndx + 2;);
  const Real turbL = m_refLen*0.07;
  const Real scale = std::sqrt(3./2.);
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_ARRAY_RESTRICT(arrU, a_U);

  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrU[MD_IX(i, rhoIndx)];
      D_TERM(Real uu = arrW[MD_IX(i, velComp1)];,
             Real vv = arrW[MD_IX(i, velComp2)];,
             Real ww = arrW[MD_IX(i, velComp3)];);
      Real speed = std::sqrt(D_TERM(uu*uu, + vv*vv, + ww*ww));
      Real mut = rho*(scale*speed*m_initTurb*turbL);
      arrU[MD_IX(i, rhoNutIndx)] = mut;
    }
  // /// Spalart 1994 paper recommends initial mu_t <= rhoNuT/2
  // a_U.setVal(CRDparam::g_mu/2., rhoNutindx);
  return;
}

/*--------------------------------------------------------------------*/
//  Find primitive turbulent variables from conservative variables
/** \param[out] a_W    Primitive variables
 *  \param[in]  a_U    Conservative variables
 *  \param[in]  a_box  Box to solve over
 *//*-----------------------------------------------------------------*/

void
SA::turbConsToPrim(FArrayBox&       a_W,
                   const FArrayBox& a_U,
                   const Box&       a_box) const
{
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int rhoNutIndx = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_ARRAY_RESTRICT(arrU, a_U);

  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrU[MD_IX(i, rhoIndx)];
      Real rhoNut = arrU[MD_IX(i, rhoNutIndx)];
      Real nut = rhoNut/rho;
      arrW[MD_IX(i, nutIndx)] = nut;
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Set the primitive state for turbulent variables, gets overwritten by
//  wall model
/** \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state corrected for presence of wall
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current Time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
SA::setTurbulentBC(
  const Box&                    a_boundaryFaceBox,
  FArrayBox&                    a_Wface,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const LevelGridMetrics&       a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int uNormIndx = velIndx + a_dir;
  // Set values to 0 at walls
  if (a_domT & CRDparam::DomainBCTypeSWall)
    {
      a_Wface.setVal(0., a_boundaryFaceBox, nutIndx, 1);
    }
  else if (a_domT & CRDparam::DomainBCTypeOutflow)
    {
      return;
    }
  else
    {
      const Real scale = std::sqrt(3./2.);
      const Real turbLen = m_refLen*0.07;
      MD_ARRAY_RESTRICT(arrWFace, a_Wface);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          Real uvel = arrWFace[MD_IX(i, uNormIndx)];
          arrWFace[MD_IX(i, nutIndx)] = scale*uvel*m_inletTurb*turbLen;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Applies the wall model to the faces and ghost cells
/** \param[in]  a_boundaryFaceBox
 *                      Box on which to apply the wall model
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Box on which to apply the wall model for ghost cells
 *  \param[in]  a_disjointBox
 *                      Disjoint box for grid metric functions
 *  \param[in]  a_WfaceBdryFab
 *                      Face-averaged primitive state including on boundaries
 *  \param[out] a_WfaceBdryFab
 *                      Face-averaged primitive state modified with wall model
 *  \param[in]  a_WfaceAvgDirFab
 *                      Face-averaged primitive state including on boundaries
 *  \param[out] a_WfaceAvgDirFab
 *                      Boundary face state out to a_boundaryFaceBox with
 *                      turbulent variables modified by wall model
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive with extrapolated ghost cells
 *  \param[out] a_WcellAvgFab
 *                      Cell-averaged primitive modified with wall model
 *  \param[in]  a_unitNormalFab
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_bcTypes
 *                      Vector of domain types
 *  \param[in]  a_bcBoxes
 *                      Vector of corresponding boxes for each domain type
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of the domain
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
SA::applyWallModel(
  const Box&                    a_boundaryFaceBox,
  const Box&                    a_boundaryFaceGhostBox,
  const Box&                    a_disjointBox,
  FArrayBox&                    a_WfaceBdryFab,
  FArrayBox&                    a_WfaceAvgDirFab,
  FArrayBox&                    a_WcellAvgFab,
  const FArrayBox&              a_unitNormalBasisFab,
  const CRDparam::DomainBCType& a_bcType,
  const Box&                    a_bcBox,
  const LevelGridMetrics&       a_gridMetrics,
  const DataIndex&              a_dataIndx,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Real                    a_time,
  const int                     a_level) const
{
  CH_assert(a_WfaceBdryFab.box().contains(a_boundaryFaceGhostBox));
  CH_assert(a_unitNormalBasisFab.box().contains(a_boundaryFaceGhostBox));
  CH_TIME("SA::applyWallModel");
  CRD::msg << CRD::fv4 << "SA::applyWallModel" << CRD::end;
  //**FIXME: ghost cells of velocity are not filled yet. A one-sided velocity
  //         gradient calculation is required in order to make this function
  //         work. The other option is to pull "applyWallModel" out of
  //         Inertial4thOrderOp::setBoundaryConditions and place it in its own
  //         new function.
  CRD::msg << "SA::applyWallModel: velocity ghost cells are not filled yet!"
           << CRD::error;
  // First, check to make sure a wall boundary exists in this direction
  //**FIXME: MixedBC type has been removed, so changes here may be required
  //         when it is reinstated in a different way
  bool wallExists = false;
  if (!a_bcBox.isEmpty() && (a_bcType & CRDparam::DomainBCTypeSWall))
    {
      wallExists = true;
    }
  if (!wallExists)
    {
      return;
    }
  // define some constants needed
  const BlockDomain& domain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);
  const RealVect dxVect = a_gridMetrics.dxVect();
  const int lohiSign = sign(a_side);
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int nutIndx = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  D_TERM( const int NormIndx = a_dir;,
          const int Tan1Indx = (a_dir+1)%SpaceDim;,
          // const int Tan2Indx = (a_dir+2)%SpaceDim;
    )
  IntVect wallNormVec(IntVect::Zero);
  wallNormVec[NormIndx] = 1;
  // IntVect wallTanVec(IntVect::Unit - wallNormVec);
  // Get interior box
  Box interiorBox(a_boundaryFaceGhostBox);
  interiorBox.shiftHalf(NormIndx, -lohiSign);
  // Get physical coordinates
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Shift the faces to align with the first layer of interior cells
  a_WfaceBdryFab.shiftHalf(a_dir, -lohiSign);
  {
    // Calculate the viscosity and conductivity at the wall
    FABSTACKTEMP(muFab, interiorBox, 1);
    FABSTACKTEMP(kappaFab, interiorBox, 1);
    CRDparam::g_CRDPhysics->calcCoeffKappaMu(interiorBox,
                                             muFab,
                                             kappaFab,
                                             a_WfaceBdryFab);
    // FAB of face-averaged gradients of transverse velocities
    // Ordered as 'du/dy, dw/dy' where u and w are wall tangent velocities
    // and y is the wall normal direction
    FABSTACKTEMP(GradUfaceBdry, a_boundaryFaceGhostBox, SpaceDim*SpaceDim);
    // In generalized mapped space this takes the form
    // 'du/dy' = $((\nabla_{\xi} \vec{u})(N^T/J) .  \hat{e}_t ) . \hat{e}_n$
    // Buckle up because this is going to be messy and expensive to solve!

    // solve $\grad_{xi} u$  on the face
    // a box to solve cell averaged gradients in
    const int gradStencilWidth = 2;
    Box gradBox(interiorBox);
    gradBox.grow(gradStencilWidth*wallNormVec);
    // Cell avg gradients, stencil width 2 needed for face averaged gradient
    FABSTACKTEMP(GradUcellAvg, gradBox, SpaceDim*SpaceDim);
    // alias the velocity
    FArrayBox WcellVel(CRDparam::g_CRDPhysics->velocityInterval(),
                       a_WcellAvgFab);
    // Compute the tangent components of the cell-averaged velocity gradients 
    CH_assert(a_WcellAvgFab.box().contains(grow(interiorBox,
                                                gradStencilWidth)));
    for(int velDir = 0; velDir != SpaceDim; ++velDir)
      {
        for(int gradDir = 0; gradDir != SpaceDim; ++gradDir)
          {
            if (gradDir != NormIndx)
              {
                const int gradVelComp =
                  ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
                FORT_CELLGRADDIR4THO(CHF_FRA1(GradUcellAvg, gradVelComp),
                                     CHF_CONST_FRA1(WcellVel, velDir),
                                     CHF_BOX(gradBox),
                                     CHF_CONST_INT(gradDir),
                                     CHF_CONST_REAL(dxVect[gradDir]));
              }
          }
      }
    // Solve cell-averaged gradients adjacent to non-periodic boundaries
    gradBox &= a_disjointBox;
    CRDparam::g_CNSIBC->cellAvgGradBC(WcellVel,
                                      GradUcellAvg,
                                      gradBox,
                                      domain,
                                      dxVect,
                                      false,
                                      a_level);
    // Compute face averaged gradient from the cell averaged gradient
    for(int velDir = 0; velDir != SpaceDim; ++velDir)
      {
        for(int gradDir = 0; gradDir != SpaceDim; ++gradDir)
          {
            const int gradVelComp =
              ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
            // Compute wall normal gradients on tangent faces
            // stencil of 2 cell-avg values in normal direction
            if(gradDir == NormIndx)
              {
                FORT_FACEGRADNORMDIR4THO(
                  CHF_FRA1(GradUfaceBdry, gradVelComp),
                  CHF_CONST_FRA1(WcellVel, velDir),
                  CHF_BOX(a_boundaryFaceGhostBox),
                  CHF_CONST_INT(NormIndx),
                  CHF_CONST_REAL(dxVect[gradDir]));
              }
            // Set tangential gradients on tangent faces
            // stencil of 2 cell-avg gradients in tangent direction 
            else
              {
                FORT_FACEGRADTANGDIR4THO(
                  CHF_FRA1(GradUfaceBdry, gradVelComp),
                  CHF_CONST_FRA1(GradUcellAvg, gradVelComp),
                  CHF_BOX(a_boundaryFaceGhostBox),
                  CHF_CONST_INT(NormIndx));
              }
          }
      }
    {
      // do transformation to get $\grad_{x} u$ as $(\grad_{\xi} u)(N^T/J)$
      FLUXBOXSTACKTEMP(AvgJinv, interiorBox, 1);
      blockCoordSys.getAvgJinverse(AvgJinv, interiorBox);
      FABSTACKTEMP(GradXUfaceBdry, a_boundaryFaceGhostBox, SpaceDim*SpaceDim);
      FABSTACKTEMP(NtJ, a_boundaryFaceGhostBox, SpaceDim*SpaceDim);
      MD_ARRAY_RESTRICT(arrNtJ, NtJ);
      MD_ARRAY_RESTRICT(arrJinv, AvgJinv[NormIndx]);
      MD_ARRAY_RESTRICT(arrN, a_gridMetrics.m_N[a_dataIndx][NormIndx]);
      MD_BOXLOOP(a_boundaryFaceGhostBox, i)
        {
          for (int c=0; c!=SpaceDim*SpaceDim; ++c)
            {
              arrNtJ[MD_IX(i, c)] = arrN[MD_IX(i, c)] * arrJinv[MD_IX(i, 0)];
            }
        }
      PatchMappedFunc::gradientCStoPS(a_boundaryFaceGhostBox,
                                      GradXUfaceBdry,
                                      GradUfaceBdry,
                                      NtJ);
      // get $\hat{e}_n$, the wall normal direction per cell
      FABSTACKTEMP(normalVect, a_boundaryFaceGhostBox, SpaceDim);
      normalVect.setVal(0.0);
      normalVect.setVal(1.0, NormIndx);
      // Transform normalVect into Cartesian physical space
      FORT_REVERSETRANSFORMF(CHF_FRA(normalVect),
                             CHF_CONST_FRA(a_unitNormalBasisFab),
                             CHF_BOX(a_boundaryFaceGhostBox));
      // get $\hat{e}_t$, the wall tangent direction per cell
      // FIXME: This is not tested for 3D!! In theory it should work
      // This should be the tangent in the streamwise direction on the face
      FABSTACKTEMP(tangentVect, a_boundaryFaceGhostBox, SpaceDim);
      // in 2D this is simple, since there is only 1 tangent direction
      if (SpaceDim == 2)
        {
          tangentVect.setVal(0.0);
          tangentVect.setVal(1.0, Tan1Indx);
          // Transform tangentVect into Cartesian physical space
          FORT_REVERSETRANSFORMF(CHF_FRA(tangentVect),
                                 CHF_CONST_FRA(a_unitNormalBasisFab),
                                 CHF_BOX(a_boundaryFaceGhostBox));
        }
      // For 3D the transformation gives 2 tangent vectors aligned with the
      // grid, but we need a single tangent aligned with the flow.
      // So transform the velocity vector, remove the wall normal component
      // and result is the tangent flow vector
      else if (SpaceDim == 3)
        {
          const int tmpVelIndx = 0;
          FORT_FORWARDTRANSFORMGENF(CHF_FRA(tangentVect),
                                    CHF_INT(tmpVelIndx),
                                    CHF_CONST_FRA(WcellVel),
                                    CHF_INT(tmpVelIndx),
                                    CHF_CONST_FRA(a_unitNormalBasisFab),
                                    CHF_BOX(a_boundaryFaceGhostBox));
          tangentVect.setVal(0.0, NormIndx);
          // Transform tangentVect back into Cartesian physical space
          FORT_REVERSETRANSFORMF(CHF_FRA(tangentVect),
                                 CHF_CONST_FRA(a_unitNormalBasisFab),
                                 CHF_BOX(a_boundaryFaceGhostBox));
          // Normalize the wall-tangent vector -- streamwise unit vector
          PatchMappedFunc::normalize(a_boundaryFaceGhostBox,
                                     tangentVect,
                                     Interval(tmpVelIndx, SpaceDim));
        }
      // solve $((\nabla_{\xi} \vec{u})(N^T/J) .  \hat{e}_t ) . \hat{e}_n$
      // and store it in the first component of GradUfaceBdry
      MD_ARRAY_RESTRICT(arrGradU, GradUfaceBdry);
      MD_ARRAY_RESTRICT(arrGradXU, GradXUfaceBdry);
      MD_ARRAY_RESTRICT(arrNorm, normalVect);
      MD_ARRAY_RESTRICT(arrTang, tangentVect);
      MD_BOXLOOP(a_boundaryFaceGhostBox, i)
        {
          arrGradU[MD_IX(i, 0)] = 0.0;
          for (int e=0; e!=SpaceDim; ++e)
            {
              arrGradU[MD_IX(i,0)] +=
                arrTang[MD_IX(i, e)] * (
                  D_TERM(
                    arrGradXU[MD_IX(i,
                                    ViscousTensor4thOrderOp::tensorIdxRowOrder(
                                      e, 0))] * arrNorm[MD_IX(i, 0)],
                    + arrGradXU[MD_IX(i,
                                    ViscousTensor4thOrderOp::tensorIdxRowOrder(
                                      e, 1))] * arrNorm[MD_IX(i, 1)],
                    + arrGradXU[MD_IX(i,
                                    ViscousTensor4thOrderOp::tensorIdxRowOrder(
                                      e, 2))] * arrNorm[MD_IX(i, 2)])
                  );
            }
        }
    }
  
    // shift from face to cell
    GradUfaceBdry.shiftHalf(NormIndx, -lohiSign);
    // get the wall distances
    CH_assert(a_gridMetrics.distanceDefined());
    const FArrayBox& distField = a_gridMetrics.m_distance[a_dataIndx];

    MD_ARRAY_RESTRICT(arrdUdY, GradUfaceBdry);
    MD_ARRAY_RESTRICT(arrMu, muFab);
    MD_ARRAY_RESTRICT(arrWCell, a_WcellAvgFab);
    MD_ARRAY_RESTRICT(arrWFace, a_WfaceBdryFab);
    MD_ARRAY_RESTRICT(arrDist, distField);
    const int MD_ID(o, NormIndx);
    //**FIXME: MixedBC type no longer is used. Changes may be required here
    //         when it is reimplemented in a different form.
    for(int vecL = 0; vecL != 1; ++vecL)
      {
        Box curBox = a_bcBox;
        curBox &= a_boundaryFaceGhostBox;
        if ((!curBox.isEmpty()) &&
            (a_bcType & CRDparam::DomainBCTypeSWall))
          {
            // cell box next to face
            curBox.shiftHalf(NormIndx, -lohiSign);
            MD_BOXLOOP(curBox, i)
              {
                Real rho = arrWCell[MD_IX(i, rhoIndx)];
                Real mu = arrMu[MD_IX(i, 0)];
                Real nu = mu/rho;
                Real yDist = arrDist[MD_IX(i, 0)];
                // note this is now streamwise velocity in normal direction
                Real tw = mu*std::abs(arrdUdY[MD_IX(i, 0)]);
                Real ut = std::sqrt(tw/rho);
                Real yplus = yDist*ut*rho/mu;
                Real nutWall = 0.;
                // If the first yplus value is greater than yPlusLam, then
                // use the wall treatment
                if(yplus > m_yPlusLam)
                  {
                    nutWall = nu*(m_k*yplus/(std::log(m_E*yplus)) - 1.);
                  }
                else
                  {
                    nutWall = 0.;
                  }
                arrWFace[MD_IX(i, nutIndx)] = nutWall;
                // Extrapolate ghost cell values
                Real Wj   = arrWCell[MD_IX(i, nutIndx)];
                Real Wjp  = arrWCell[MD_OFFSETIX(i,-,lohiSign*o,nutIndx)];
                Real Wjpp = arrWCell[MD_OFFSETIX(i,-,2*lohiSign*o,nutIndx)];
                Real Wjm = 4.*nutWall - (Wjpp - 5.*Wjp + 13.*Wj)/3.;
                Real Wjmm = 20.*nutWall - (8.*Wjpp - 37.*Wjp + 83.*Wj)/3. - Wjm;
                arrWCell[MD_OFFSETIX(i,+,lohiSign*o,nutIndx)] = Wjm;
                arrWCell[MD_OFFSETIX(i,+,2*lohiSign*o,nutIndx)] = Wjmm; 
              }
          }
      }
  }
  // Shift faces back to original locations
  a_WfaceBdryFab.shiftHalf(NormIndx, lohiSign);

  // Overwrite the face values of the turbulent variables
  a_WfaceAvgDirFab.copy(a_boundaryFaceBox,
                        CRDparam::g_CRDPhysics->turbPrimInterval(),
                        a_boundaryFaceBox, a_WfaceBdryFab,
                        CRDparam::g_CRDPhysics->turbPrimInterval());
}
