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
 * \file ViscousTensor4thOrderOp.cpp
 *
 * \brief Non-inline definitions for classes in ViscousTensor4thOrderOp.H
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "MOLUtilFunc.H"
#include "LevelGridMetrics.H"
#include "PointwiseDotProdF_F.H"
#include "EdgeToCellF_F.H"

//----- Internal -----//

#include "DebugOut.H"
#include "DataTemp.H"
#include "CRDPhysics.H"
#include "CRDutil.H"
#include "CNSIBC.H"
#include "PatchMappedFunc.H"
#include "ViscousTensor4thOrderOp.H"
#include "ViscousTensor4thOrderOpF_F.H"


/*******************************************************************************
 *
 * Class ViscousTensor4thOrderOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *  \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *//*-----------------------------------------------------------------*/

ViscousTensor4thOrderOp::ViscousTensor4thOrderOp(
  LevelGridMetrics& a_levelGridMetrics)
  :
  m_levelGridMetrics(a_levelGridMetrics)
{ }


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Add in the viscous flux
/** Version of flux that allocates arrays for intermediate results on
 *  the stack
 *  \param[in]  a_box   The flux is computed on the faces of 'a_box'
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[out] a_FluxfaceAvgFxb
 *                      The face-averaged flux
 *  \param[out] a_FluxfacePntFxb
 *                      The face-centered flux
 *  \param[out] a_turbSourceAvgFab
 *                      FAB containing mapped cell-averaged turbulent sources
 *  \param[out] a_invDtFab
 *                      Inverse time step values
 *  \param[in]  a_WcellAvgFab
 *                      The average primitive state in the cells set
 *                      on domain and 3 ghost at interior and 2 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[in]  a_WcellPntFab
 *                      The cell-centered primitive variables, for use in 
 *                      turbulent source terms
 *  \param[out] a_timeAvgFab
 *                      Time-averaged data for post-processing turbulent cases
 *  \param[in]  a_WfaceAvgFxb
 *                      The face-averaged primitive variables
 *  \param[in]  a_WfacePntFxb
 *                      The face-centered primitive variables
 *  \param[in]  a_facePntVelGradFxb
 *                      The face-centered velocity gradient
 *  \param[in]  a_facePntDeltaC
 *                      Face-centered cell-cutoff length for LES
 *  \param[in]  a_unitNormalsFxb
 *                      Face unit-normal vectors used for LES wall-model
 *  \param[out] a_invDtFab
 *                      Inverse time step values
 *  \param[in]  a_dt    Time-step size on current level
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_timeAfterFilterInit
 *                      Total time after initializing time-averaging filter
 *  \param[in]  a_level Index of the AMR level
 *  \param[in]  a_minDiffDt
 *                      Previous minimum diffusive time step
 *  \param[out] a_minDiffDt
 *                      Updated minimum diffusive time step
 *  \param[in] a_minDiffDtCell
 *                      Previous cell with minimum diffusive time step
 *  \param[out] a_minDiffDtCell
 *                      Updated cell with minimum diffusive time step
 *
 *  \note
 *  <ul>
 *    <li> On entry, ghost cells of a_vel must be filled by exchange
 *    <li> The flux is computed on all faces using interior data.
 *         Boundary fluxes may be corrected later.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::flux(const Box&         a_box,
                              const BlockDomain& a_domain,
                              FluxBox&           a_FluxfaceAvgFxb,
                              FluxBox&           a_FluxfacePntFxb,
                              FArrayBox&         a_invDtFab,
                              FArrayBox&         a_turbSourceAvgFab,
                              FArrayBox&         a_WcellAvgFab,
                              FArrayBox&         a_WcellPntFab,
                              FArrayBox&         a_timeAvgFab,
                              FluxBox&           a_faceAvgPlotFxb,
                              FluxBox&           a_WfaceAvgFxb,
                              FluxBox&           a_WfacePntFxb,
                              FluxBox&           a_facePntVelGradFxb,
                              FluxBox&           a_stressFluxFxb,
                              const FluxBox&     a_facePntDeltaC,
                              const FluxBox&     a_faceCoord,
                              const FluxBox&     a_unitNormalsFxb,
                              const DataIndex&   a_dataIndx,
                              const Real         a_dt,
                              const Real         a_time,
                              const Real         a_timeAfterFilterInit,
                              const int          a_level,
                              Real&              a_minDiffDt,
                              IntVect&           a_minDiffDtCell) const
{
  CH_TIME("ViscousTensor4thOrderOp::flux");
  CRD::msg << CRD::fv4 << "ViscousTensor4thOrderOp::flux" << CRD::end;
  // a_FluxfaceAvgFxb and a_FluxfacePntFxb are indexed as contiguous
  // in space, as shown
  // 0: rho x-flux
  // 1: rho y-flux
  // 2: rho u x-flux
  // 3: rho u y-flux
  // 4: rho v x-flux
  // 5: rho v y-flux etc.
  const int stressTensorComp = SpaceDim*SpaceDim;

  const int numFluxes = CRDparam::g_CRDPhysics->numFluxes();
  CH_assert(a_turbSourceAvgFab.nComp() == numFluxes);

//--Domain boxes

  // A domain box grown by 5 in connected directions and 2 in directions of
  // domain BC
  Box domain5Pp2BCbox = grow(a_domain.domainBox(), 5);
  // A domain box grown by 5 in connected directions and reduced by 1 in
  // directions of domain BC
  Box domain5Pm1BCbox = grow(a_domain.domainBox(), 5);
  for (const int dir : EachDir)
    {
      if (!a_domain.isConnected(dir, Side::Lo))
        {
          domain5Pp2BCbox.growLo(dir,  2 - 5);
          domain5Pm1BCbox.growLo(dir, -1 - 5);
        }
      if (!a_domain.isConnected(dir, Side::Hi))
        {
          domain5Pp2BCbox.growHi(dir,  2 - 5);
          domain5Pm1BCbox.growHi(dir, -1 - 5);
        }
    }
  
//--Other boxes

  Box box2Dom = grow(a_box,2);
  box2Dom &= a_domain;
  Box box2Dom2 = grow(a_box,2);
  box2Dom2 &= domain5Pp2BCbox;
  Box box1Dom = grow(a_box,1);
  box1Dom &= a_domain;

  Box GpntBox(box1Dom);
  Box GavgBox(a_box);
  Box box4Dom = grow(a_box,4);
  box4Dom &= a_domain;
  Box gradBox(box2Dom2);

  // Gradient of a tensor are ordered for a 2D example as
  // 0: du/dx
  // 1: du/dy
  // 2: dv/dx
  // 3: dv/dy
  // FluxBox of the face-centered viscous stress tensor
  FLUXBOXSTACKTEMP(VSTfacePntFxb, GpntBox, stressTensorComp);
  // FluxBox of the face-centered physical-space velocity-gradient tensor
  FLUXBOXSTACKTEMP(NGradUfacePntFxb, GpntBox, stressTensorComp);
  // FluxBox of the face-centered strain rate tensor
  FLUXBOXSTACKTEMP(strainRateTensorFxb, GpntBox, stressTensorComp);
  // Gradient of a scalar are ordered for a 2D example as
  // 0: dT/dx
  // 1: dT/dy
  // FluxBox of the point values of N^T/J Grad T
  FLUXBOXSTACKTEMP(NGradTfacePntFxb, GpntBox, SpaceDim);

  // FAB of the cell-averaged primitive variables
  // Gradients are ordered for a 2D example as
  // 0: drho/dx
  // 1: drho/dy
  // 2: du/dx
  // 3: du/dy etc.
  const int numGrad = CRDparam::g_CRDPhysics->numPrimitive()*SpaceDim;
  FABSTACKTEMP(GradWcellAvgFab, gradBox, numGrad);

  // Compute face-centered tau/mu, must multiply by mu
  computeVST(a_box,
             a_domain,
             domain5Pp2BCbox,
             domain5Pm1BCbox,
             VSTfacePntFxb,
             strainRateTensorFxb,
             NGradUfacePntFxb,
             GradWcellAvgFab,
             a_WcellAvgFab,
             a_WfaceAvgFxb,
             a_unitNormalsFxb,
             a_dataIndx,
             a_time,
             a_level);
  {
    NGradTfacePntFxb.setVal(0.);
    // Compute face-centered N^T/J grad T
    const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
    const bool rowMult = false;
    // Component for filling NGradTfacePntFxb is 0, only contains grad T
    const int tempVecStart = 0;
    computeNGradPhi(a_box,
                    a_domain,
                    domain5Pp2BCbox,
                    domain5Pm1BCbox,
                    NGradTfacePntFxb,
                    GradWcellAvgFab,
                    a_WcellAvgFab,
                    a_WfaceAvgFxb,
                    a_dataIndx,
                    tempIndx,
                    tempVecStart,
                    rowMult,
                    a_time,
                    a_level);
  }
  // FXB that contains the LES SGS momentum flux on each face
  FLUXBOXSTACKTEMP(sgsMomentumFxb, box1Dom, stressTensorComp);
  sgsMomentumFxb.setVal(0.);
  // FXB that contains the LES SGS energy flux on each face
  FLUXBOXSTACKTEMP(sgsEnergyFxb, box1Dom, SpaceDim);
  sgsEnergyFxb.setVal(0.);
  // Compute additions to fluxes from turbulence models
  if (CRDparam::g_turbModelType)
    {
      a_facePntVelGradFxb.copy(NGradUfacePntFxb);
      CRDparam::g_CRDPhysics->sgsModelFlux(sgsMomentumFxb,
                                           sgsEnergyFxb,
                                           a_faceAvgPlotFxb,
                                           a_WfacePntFxb,
                                           a_WcellAvgFab,
                                           NGradUfacePntFxb,
                                           strainRateTensorFxb,
                                           NGradTfacePntFxb,
                                           a_facePntDeltaC,
                                           a_faceCoord,
                                           a_unitNormalsFxb,
                                           m_levelGridMetrics,
                                           a_dataIndx,
                                           a_domain,
                                           a_box);
    }
  // Since the molecular diffusivity of each species, D_n, requires the
  // dynamic viscosity of the flow, we must save it so it can be used later
  // FXB that contains the dynamic viscosity on each face
  FLUXBOXSTACKTEMP(muFxb, box1Dom, 1);
  // FXB that contains the thermal conductivity
  FLUXBOXSTACKTEMP(kappaFxb, box1Dom, 1);
  // Calculate the dynamic viscosity and thermal conductivity on each face
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      IntVect growVect(IntVect::Unit);
      growVect[dir] = 0;
      Box muSolveBox = grow(a_box, growVect);
      muSolveBox &= a_domain;
      muSolveBox.surroundingNodes(dir);
      CRDparam::g_CRDPhysics->calcCoeffKappaMu(muSolveBox,
                                               muFxb[dir],
                                               kappaFxb[dir],
                                               a_WfacePntFxb[dir]);
      if (CRDparam::g_turbModelType)
        {
          FABSTACKTEMP(muTurb, muSolveBox, 1);
          muTurb.setVal(0.);
          FABSTACKTEMP(kappaTurb, muSolveBox, 1);
          kappaTurb.setVal(0.);
          CRDparam::g_CRDPhysics->calcCoeffKappatMut(muSolveBox,
                                                     muTurb,
                                                     kappaTurb,
                                                     muFxb[dir],
                                                     kappaFxb[dir],
                                                     a_WfacePntFxb[dir],
                                                     strainRateTensorFxb[dir],
                                                     m_levelGridMetrics,
                                                     dir);
          muFxb[dir].plus(muTurb);
          kappaFxb[dir].plus(kappaTurb);
        }
    }
  // Solve for the viscous time step
  {
    // FluxBox of pointwise N^T/J from LevelGridMetrics class
    const FluxBox& NtJ = m_levelGridMetrics.m_NtJ[a_dataIndx];
    const RealVect dxVect = m_levelGridMetrics.dxVect();
    FABSTACKTEMP(muCell, a_box, 1);
    muCell.setVal(0.);
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        FArrayBox& muFab = muFxb[dir];
        FORT_EDGETOINCREMENTCELL(CHF_CONST_FRA1(muFab, 0),
                                 CHF_FRA1(muCell, 0),
                                 CHF_BOX(a_box),
                                 CHF_CONST_INT(dir));
      }
    // Average by dividing by SpaceDim
    muCell.divide(SpaceDim);
    Real lambdaDMax = 16./3.;
    CRDparam::g_CRDPhysics->ellipticDtEvalRHS(a_domain,
                                              a_box,
                                              a_box,
                                              a_invDtFab,
                                              a_WcellPntFab,
                                              muCell,
                                              lambdaDMax,
                                              NtJ,
                                              m_levelGridMetrics,
                                              dxVect,
                                              a_minDiffDt,
                                              a_minDiffDtCell);
  }
  // Multiply the VST by mu
  for (int n = 0; n != stressTensorComp; ++n)
    {
      VSTfacePntFxb.mult(muFxb, box1Dom, 0, n, 1);
    }
  // Multiply the grad T by kappa
  for (int n = 0; n != SpaceDim; ++n)
    {
      NGradTfacePntFxb.mult(kappaFxb, box1Dom, 0, n, 1);
    }
  // Perform some time-averaging of turbulent quantities
  //**NOTE: Perform here because VST and SGS stresses are available
  if (CNSIBC::s_firstRKStage && (a_time >= CRDparam::g_startTimeAvgTime))
    {
      PatchMappedFunc::timeAvgData(a_timeAvgFab,
                                   a_WcellAvgFab,
                                   a_WfaceAvgFxb,
                                   a_WfacePntFxb,
                                   VSTfacePntFxb,
                                   sgsMomentumFxb,
                                   a_domain,
                                   a_box,
                                   a_dt,
                                   a_time,
                                   a_timeAfterFilterInit);
    }
  // Must set to zero because values are added to the fluxes
  a_FluxfacePntFxb.setVal(0.);
  a_FluxfaceAvgFxb.setVal(0.);
  {
    Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
    FluxBox UfacePntFxb(a_WfacePntFxb.box(), SpaceDim,
                        D_DECL(a_WfacePntFxb[0].dataPtr(velIntv.begin()),
                               a_WfacePntFxb[1].dataPtr(velIntv.begin()),
                               a_WfacePntFxb[2].dataPtr(velIntv.begin())));
    // FluxBox of point values of elliptic and heat flux in the energy equation
    FLUXBOXSTACKTEMP(energyFxb, box1Dom, SpaceDim);
    energyFxb.setVal(0.);

    // This FluxBox contains the \vec{J}_n = rho*D_n N/J grad c_n 
    // for each species
    // Example of ordering for 2 species in 2D
    // Comp 0: D(c_0)/D(x_0) 
    // Comp 1: D(c_0)/D(x_1)
    // Comp 2: D(c_1)/D(x_0)
    // Comp 3: D(c_1)/D(x_1)

    if ((CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf) &&
        (CRDparam::g_numSpecies > 1))
      {
        const int numSpecies = CRDparam::g_numSpecies;
        const int JnNumComp = SpaceDim*numSpecies;
        // First, we will solve for \nabla c_n in JnfacePntFxb
        FLUXBOXSTACKTEMP(JnfacePntFxb, GpntBox, JnNumComp);
        JnfacePntFxb.setVal(0.);
        const bool rowMult = false;
        for (int n = 0; n != numSpecies; ++n)
          {
            // The primitive component number for the nth species
            const int indexCN = 
              CRDparam::g_CRDPhysics->speciesPrimInterval().begin() + n;
            // The first component of JnfacePntFxb to modify
            const int cnVecStart = SpaceDim*n;
            computeNGradPhi(a_box,
                            a_domain,
                            domain5Pp2BCbox,
                            domain5Pm1BCbox,
                            JnfacePntFxb,
                            GradWcellAvgFab,
                            a_WcellAvgFab,
                            a_WfaceAvgFxb,
                            a_dataIndx,
                            indexCN,
                            cnVecStart,
                            rowMult,
                            a_time,
                            a_level);
          }
        // Next, we multiply \rho*D_n*\nabla c_n to get J_n in JnfacePntFxb
        // We also add h_n*J_n to the energy in this function
        for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
          {
            Box solvePntBox(GpntBox);
            solvePntBox.surroundingNodes(faceDir);
            CRDparam::g_CRDPhysics->speciesDiffHeatFlux(solvePntBox,
                                                        JnfacePntFxb[faceDir],
                                                        energyFxb[faceDir],
                                                        muFxb[faceDir],
                                                        kappaFxb[faceDir],
                                                        a_WfacePntFxb[faceDir]);
          }
        
        // Now add the species fluxes to the G vector
        for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
          {
            FArrayBox& GfacePntFab = a_FluxfacePntFxb[faceDir];
            const FArrayBox& JnfacePntFab = JnfacePntFxb[faceDir];
            for (int space = 0; space != SpaceDim; ++space)
              {
                for (int n = 0; n != numSpecies; ++n)
                  {
                    const int nFluxComp = 
                      CRDparam::g_CRDPhysics->speciesConsInterval().begin() + n;
                    // Copy the species flux to the species component of G
                    const int Gcomp = retrieveFluxDyadComp(nFluxComp,space);
                    const int JnComp = n*SpaceDim + space;
                    GfacePntFab.copy(JnfacePntFab, JnComp, Gcomp);
                  }
              }
          }
      }
    // Compute the face-averaged dot product of Tau and velocity
    // into the energy FluxBox
    const int UVecStart = 0;
    const bool rowMult = false;
    computeFOMatVecProd(box1Dom,
                        energyFxb,
                        VSTfacePntFxb,
                        UfacePntFxb,
                        UVecStart,
                        rowMult);
    // Add N/J Grad T to the energy flux
    energyFxb.plus(NGradTfacePntFxb,
                   box1Dom,
                   0,
                   0,
                   SpaceDim);
    // Subtract LES SGS flux from energyFxb and VSTfacePntFxb
    if (CRDparam::g_turbModelType)
      {
        VSTfacePntFxb.minus(sgsMomentumFxb,
                            box1Dom,
                            0,
                            0,
                            SpaceDim*SpaceDim);
        energyFxb.minus(sgsEnergyFxb,
                        box1Dom,
                        0,
                        0,
                        SpaceDim);
        a_stressFluxFxb.copy(VSTfacePntFxb, 0, 0, SpaceDim*SpaceDim);
      }
    // Add any turbulent diffusive fluxes
    const int numTurbVar = CRDparam::g_CRDPhysics->numTurbVar();
    if ((CRDparam::g_turbModelType && numTurbVar > 0) &&
        (!(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)))
      {
        const int turbVarStartW =
          CRDparam::g_CRDPhysics->turbPrimInterval().begin();
        // This FluxBox contains grad nu_t for each turbulent variable
        // Example of ordering for 2 variables in 2D
        // Comp 0: D(k)/D(x_0) 
        // Comp 1: D(k)/D(x_1)
        // Comp 2: D(omega)/D(x_0)
        // Comp 3: D(omega)/D(x_1)
        const int NutNumComp = SpaceDim*numTurbVar;
        FLUXBOXSTACKTEMP(NGradNutfacePntFxb, GpntBox, NutNumComp);
        NGradNutfacePntFxb.setVal(0.);
        const bool rowMult = false;
        for (int tComp = 0; tComp != numTurbVar; ++tComp)
          {
            // The primitive component related to the tComp variable
            const int wTurbComp = turbVarStartW + tComp;
            // The first component in NutfacePntFxb to modify
            const int tVecStart = SpaceDim*tComp;
            computeNGradPhi(a_box,
                            a_domain,
                            domain5Pp2BCbox,
                            domain5Pm1BCbox,
                            NGradNutfacePntFxb,
                            GradWcellAvgFab,
                            a_WcellAvgFab,
                            a_WfaceAvgFxb,
                            a_dataIndx,
                            wTurbComp,
                            tVecStart,
                            rowMult,
                            a_time,
                            a_level);
            for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
              {
                Box solvePntBox(GpntBox);
                solvePntBox.surroundingNodes(faceDir);
                CRDparam::g_CRDPhysics->calcTurbDiffFlux(
                  solvePntBox,
                  NGradNutfacePntFxb[faceDir],
                  muFxb[faceDir],
                  kappaFxb[faceDir],
                  a_WfacePntFxb[faceDir],
                  faceDir,
                  wTurbComp,
                  tVecStart);
              }
          }
        const int turbVarStartU =
          CRDparam::g_CRDPhysics->turbConsInterval().begin();
        // Now add the turbulent fluxes to the G vector
        for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
          {
            FArrayBox& GfacePntFab = a_FluxfacePntFxb[faceDir];
            const FArrayBox& NGradNutfacePntFab = NGradNutfacePntFxb[faceDir];
            for (int space = 0; space != SpaceDim; ++space)
              {
                for (int tComp = 0; tComp != numTurbVar; ++tComp)
                  {
                    const int tFluxComp = turbVarStartU + tComp;
                    // Copy the species flux to the species component of G
                    const int Gcomp = retrieveFluxDyadComp(tFluxComp,space);
                    const int NutComp = tComp*SpaceDim + space;
                    GfacePntFab.copy(NGradNutfacePntFab, NutComp, Gcomp);
                  }
              }
          }
      }
    // Solve turbulent source terms
    if ((CRDparam::g_turbModelType) &&
        (!(CRDparam::g_sgsModelType & CRDparam::SGSModelStretchedVortex)))
      {
        const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
        const int rhoVecStart = 0;
        // Calculate the gradient of density
        computeNGradPhi(a_box,
                        a_domain,
                        domain5Pp2BCbox,
                        domain5Pm1BCbox,
                        NGradTfacePntFxb, // Not used
                        GradWcellAvgFab,
                        a_WcellAvgFab,
                        a_WfaceAvgFxb,
                        a_dataIndx,
                        rhoIndx,
                        rhoVecStart,
                        rowMult,
                        a_time,
                        a_level,
                        true);
        // Find cell-centered gradients from cell-averaged gradients
        FABSTACKTEMP(GradWcellPntFab, box1Dom, numGrad);
        GradWcellPntFab.copy(GradWcellAvgFab);
        //**FIXME: problem here with MMB, MOLUtilFunc::deconvolve breaks here
        MOLUtilFunc::deconvolve(GradWcellPntFab,
                                GradWcellAvgFab,
                                gradBox,
                                a_domain,
                                -1);
        // FAB of mapped cell-centered source term
        FABSTACKTEMP(turbSourcePntFab, box1Dom, numFluxes);
        turbSourcePntFab.setVal(0.);
        // Do the mapping on GradW
        FABSTACKTEMP(GradWcellPntMapFab, box1Dom, numGrad);
        PatchMappedFunc::gradientCStoPS(
          box1Dom,
          GradWcellPntMapFab,
          GradWcellPntFab,
          m_levelGridMetrics.m_cellNtJ[a_dataIndx]);
        // compute the source term
        CRDparam::g_CRDPhysics->calcTurbSourceTerms(box1Dom,
                                                    a_box,
                                                    turbSourcePntFab,
                                                    GradWcellPntMapFab,
                                                    a_WcellPntFab,
                                                    a_dataIndx,
                                                    m_levelGridMetrics);
        // Convolve find <JS> from (JS)
        a_turbSourceAvgFab.copy(turbSourcePntFab);
        MOLUtilFunc::deconvolve(a_turbSourceAvgFab,
                                turbSourcePntFab,
                                a_box,
                                a_domain,
                                1);
      }

    // Transfer the components of the VST into G since the middle
    // components of G are equal to the columns of the VST
    for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
      {
        FArrayBox& GfacePntFab = a_FluxfacePntFxb[faceDir];
        const FArrayBox& VSTfacePntFab = VSTfacePntFxb[faceDir];
        const FArrayBox& energyFab = energyFxb[faceDir];
        for (int space = 0; space != SpaceDim; ++space)
          {
            int gVar = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();
            for (int var = 0; var != SpaceDim; ++var, ++gVar)
              {
                const int vstComp = tensorIdxRowOrder(var,space);
                const int gComp = retrieveFluxDyadComp(gVar,space);
                GfacePntFab.copy(VSTfacePntFab,vstComp,gComp);
              }
            const int engFluxComp = CRDparam::g_CRDPhysics->energyFluxIndex();
            // Copy the energy flux to the final component of G
            const int enGcomp = retrieveFluxDyadComp(engFluxComp,space);
            GfacePntFab.copy(energyFab,
                             space,
                             enGcomp);
          }
      }
  }

  // Convert the face-centered values for G into the face-averaged values of G
  a_FluxfaceAvgFxb.copy(a_FluxfacePntFxb);
  MOLUtilFunc::deconvolveFace(a_FluxfaceAvgFxb,
                              a_FluxfacePntFxb,
                              a_box,
                              a_domain,
                              1);
}

/*--------------------------------------------------------------------*/
//  Compute the viscous stress tensor on the faces
/** Version of flux that allocates arrays for intermediate results on
 *  the stack
 *  \param[in]  a_box   The flux is computed on the faces of 'a_box'
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_domain5Pp2BCbox
 *                      A domain box grown by 5 in connected
 *                      directions and 2 in directions of domain BC
 *  \param[in]  a_domain5Pm1BCbox
 *                      A domain box grown by 5 in connected
 *                      directions and reduced by 1 in directions of
 *                      domain BC
 *  \param[out] a_VSTfacePntFxb
 *                      The pointwise viscous stress tensor on the face
 *                      on domain and 2 ghost at interior and 0 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[out] a_strainRateTensorFxb
 *                      The face-centered strain rate tensor.
 *  \param[out] a_NGradUfacePntFxb
 *                      Face-centered physical-space gradients of velocity
 *  \param[out] a_GradWcellAvgFab
 *                      Cell-averaged gradients of primitive variables
 *  \param[in]  a_WcellAvgFab
 *                      The average primitive state in the cells set
 *                      on domain and 5 ghost at interior and 2 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[in]  a_WfaceAvgFab
 *                      The average primitive state on the faces set
 *                      on domain and 5 ghost at interior and 2 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[in]  a_unitNormalsFxb
 *                      Face unit-normal vectors used for LES wall-model
 *  \param[in]  a_dataIndx
 *                      The data index for use with the metric terms
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Current grid level
 *  \note
 *  <ul>
 *    <li> On entry, ghost cells of a_vel must be filled by exchange
 *    <li> The flux is computed on all faces using interior data.
 *         Boundary fluxes may be corrected later.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeVST(const Box&         a_box,
                                    const BlockDomain& a_domain,
                                    const Box&         a_domain5Pp2BCbox,
                                    const Box&         a_domain5Pm1BCbox,
                                    FluxBox&           a_VSTfacePntFxb,
                                    FluxBox&           a_strainRateTensorFxb,
                                    FluxBox&           a_NGradUfacePntFxb,
                                    FArrayBox&         a_GradWcellAvgFab,
                                    FArrayBox&         a_WcellAvgFab,
                                    FluxBox&           a_WfaceAvgFxb,
                                    const FluxBox&     a_unitNormalsFxb,
                                    const DataIndex&   a_dataIndx,
                                    const Real         a_time,
                                    const int          a_level) const
{
  // Tip: If debugging, initialize the gradients to BASEFAB_REAL_SETVAL so you
  // can easily see where they are defined.
  CH_TIME("ViscousTensor4thOrderOp::computeVST");
  CRD::msg << CRD::fv4 << "ViscousTensor4thOrderOp::computeVST" << CRD::end;
  const RealVect dxVect = m_levelGridMetrics.dxVect();
  const int stressTensorComp = SpaceDim*SpaceDim;

  Box box4Dom = grow(a_box,4);
  Box box4Dom2 = grow(a_box,4);
  box4Dom &= a_domain;
  box4Dom2 &= a_domain5Pp2BCbox;
  Box box2Dom = grow(a_box,2);
  box2Dom &= a_domain;
  Box box2Dom2 = grow(a_box,2);
  box2Dom2 &= a_domain5Pp2BCbox;

  Box box1Dom = grow(a_box,1);
  box1Dom &= a_domain;
  Box UfaceBox(box4Dom);
  Box UcellBox(box4Dom2);
  Box gradBox(box2Dom2);

  // Steps to extract the velocity from the primitive variables into U
  Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  // Alias FAB of the cell-averaged velocities
  FArrayBox UcellAvgFab(velIntv, a_WcellAvgFab);
  // Alias FluxBox of the face-averaged velocity
  FluxBox UfaceAvgFxb(a_WfaceAvgFxb.box(), SpaceDim,
                      D_DECL(a_WfaceAvgFxb[0].dataPtr(velIntv.begin()),
                             a_WfaceAvgFxb[1].dataPtr(velIntv.begin()),
                             a_WfaceAvgFxb[2].dataPtr(velIntv.begin())));

  // Gradient components of the alias are ordered for 2D as:
  // 0: du/dx
  // 1: du/dy
  // 2: dv/dx
  // 3: dv/dy
  // Alias of the cell-averaged velocity gradient tensor
  const int startGradComp = gradIdx(velIntv.begin(), 0);
  const int endGradComp = gradIdx(velIntv.end(), SpaceDim - 1);
  Interval gradInt(startGradComp, endGradComp);
  FArrayBox GradUcellAvgFab(gradInt, a_GradWcellAvgFab);
  // GradUcellAvgFab.setVal(BASEFAB_REAL_SETVAL);

  // FluxBox of the face-averaged velocity gradient tensor
  FLUXBOXSTACKTEMP(GradUfaceAvgFxb, gradBox, stressTensorComp);
  // GradUfaceAvgFxb.setVal(BASEFAB_REAL_SETVAL);

  // Currently solving for velocity tensor, not a single variable
  bool solvePhi = false;
  if (!a_domain5Pm1BCbox.contains(a_box))
    {
      // Solve the cell-averaged gradients at the cells adjacent to
      // non-periodic boundaries.  Note that this only sets tangential gradients
      // to a given face.
      CRDparam::g_CNSIBC->cellAvgGradBC(UcellAvgFab,
                                        GradUcellAvgFab,
                                        a_box,
                                        a_domain,
                                        dxVect,
                                        solvePhi,
                                        a_level);
    }

  // Compute the cell-averaged velocity gradients using cell-averaged velocity
  // ---------------------------------------------------
  //  I/O           Var            Stencil      Ghosts
  // =====  ====================  ==========  ==========
  // read   UcellAvgFab            2 per dir    5 in dom
  // write  GradUcellAvgFab                1    3 in dom
  // ---------------------------------------------------
  computeInteriorCellAvgVelGrad(GradUcellAvgFab,
                                UcellAvgFab,
                                UfaceAvgFxb,
                                a_box,
                                a_domain,
                                dxVect);

  // Compute the face-averaged velocity gradients using cell-averaged 
  // velocity and cell-averaged velocity gradients
  // ---------------------------------------------------
  //  I/O           Var            Stencil      Ghosts
  // =====  ====================  ==========  ==========
  // read   UcellAvgFab            2 per dir    5 in dom
  // read   GradUcellAvgFab        2 per dir    3 in dom
  // write  GradUfaceAvgFxb                1    3 in dom
  // ---------------------------------------------------
  computeInteriorFaceAvgVelGrad(GradUfaceAvgFxb,
                                GradUcellAvgFab,
                                UcellAvgFab,
                                a_box,
                                a_domain,
                                dxVect);

  // Fluxbox for the pointwise face-centered velocity gradient
  FLUXBOXSTACKTEMP(GradUfacePntFxb, box1Dom, stressTensorComp);
  GradUfacePntFxb.copy(GradUfaceAvgFxb);

  // We must tranform the face-averaged velocity gradient into the
  // pointwise face-centered velocity gradient
  MOLUtilFunc::deconvolveFace(GradUfacePntFxb,
                              GradUfaceAvgFxb,
                              box1Dom,
                              a_domain,
                              -1);

  // FluxBox of pointwise N^T/J from LevelGridMetrics class
  const FluxBox& NtJ = m_levelGridMetrics.m_NtJ[a_dataIndx];

  // Initialize the viscous stress tensor flux box
  a_VSTfacePntFxb.setVal(0.);
  a_strainRateTensorFxb.setVal(0.);
  a_NGradUfacePntFxb.setVal(0.);
  {
    // Compute face-centered values of N^T/J Grad u
    computeMatrixProduct(a_NGradUfacePntFxb,
                         GradUfacePntFxb,
                         NtJ,
                         0);

    // Modify wall-face velocity-gradients if specified by the LES wall-model
    //**NOTE: The SSV implementation demonstrates how to enforce zero values
    //        for tangential-derivatives of velocity on no-slip walls
    if (CRDparam::g_enforceModeledWallShearStress)
      {
        if (CRDparam::g_turbModelType)
          {
            CRDparam::g_CRDPhysics->modeledNoSlipWallVelocityGradient(
              a_NGradUfacePntFxb,
              a_unitNormalsFxb,
              a_WcellAvgFab,
              m_levelGridMetrics,
              a_domain,
              box1Dom,
              a_box);
          }
        else
          {
            // Model the wall shear-stress simply using the law-of-the-wall
            // This is only applicable to turbulent boundary layers
            CRDutil::enforceLawOfTheWallShearStress(a_NGradUfacePntFxb,
                                                    a_unitNormalsFxb,
                                                    a_WfaceAvgFxb,
                                                    m_levelGridMetrics,
                                                    a_domain,
                                                    box1Dom,
                                                    a_box);
          }
      }

    // FluxBox of DIV N^T/J u
    FLUXBOXSTACKTEMP(faceDIVNtU, box1Dom, 1);
    faceDIVNtU.setVal(0.);
    // Solve for DIV N^T/J u
    for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
      {
        const FArrayBox& NTUFab = a_NGradUfacePntFxb[faceDir];
        FArrayBox& DIVNtUFab = faceDIVNtU[faceDir];
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            const int gradVelComp =
              ViscousTensor4thOrderOp::tensorIdxRowOrder(dir, dir);
            DIVNtUFab.plus(NTUFab,gradVelComp,0,1);
          }
      }

    faceDIVNtU *= -2.*CRDparam::g_lambda;

    // Solve for the pointwise viscous stress tensor
    for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
      {
        // div N^T/J u
        const FArrayBox& DIVNTUFab = faceDIVNtU[faceDir];

        const FArrayBox& NTUFab = a_NGradUfacePntFxb[faceDir];

        FArrayBox& VSTFab = a_VSTfacePntFxb[faceDir];
        FArrayBox& strainRateFab = a_strainRateTensorFxb[faceDir];
        for (int velDir = 0; velDir != SpaceDim; ++velDir)
          {
            for (int gradDir = 0; gradDir != SpaceDim; ++gradDir)
              {
                const int gradVelComp =
                  ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
                const int gradInvComp = 
                  ViscousTensor4thOrderOp::tensorIdxRowOrder(gradDir, velDir);
                // Add ${\rm{N^T}}/J \nabla \vec{u}$ to $\vec{\vec{{S}}$
                strainRateFab.plus(NTUFab, 1./2., gradVelComp, gradVelComp, 1);
                // Add $({\rm{N^T}}/J \nabla \vec{u})^{\rm{T}}$
                // to $\vec{\vec{{S}}$
                strainRateFab.plus(NTUFab, 1./2., gradInvComp, gradVelComp, 1);
                // Add the divergence of ${\rm{N^T}}/J u$ to $\vec{\vec{\tau}}$
                if (gradDir == velDir)
                  {
                    VSTFab.plus(DIVNTUFab,0,gradVelComp,1);
                  }
              }
          }
      }
    // Add the strain rate tensor times 2 to the viscous stress tensor
    for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
      {
        a_VSTfacePntFxb[faceDir].plus(a_strainRateTensorFxb[faceDir],
                                      2., 0, 0, stressTensorComp);
      }
  }
}

/*--------------------------------------------------------------------*/
//  Computes the cell-averaged velocity gradients for the interior grid
/** 
 *  \param[out] a_GradUcellAvgFab
 *                      FAB of the cell-averaged velocity gradients
 *  \param[in]  a_UfaceAvgFxb   
 *                      Fxb of the face-averaged velocities
 *  \param[in]  a_box   Interior box
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_dxVect
 *                      Vector of the computational grid spacing
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeInteriorCellAvgVelGrad(
  FArrayBox&         a_GradUcellAvgFab,
  const FArrayBox&   a_UcellAvgFab,
  const FluxBox&     a_UfaceAvgFxb,
  const Box&         a_box,
  const BlockDomain& a_domain,
  const RealVect&    a_dxVect) const
{
  for (int velDir = 0; velDir != SpaceDim; ++velDir)
    {
      for (int gradDir = 0; gradDir != SpaceDim; ++gradDir)
        {
          Box centerBox(a_box);
          centerBox.grow(c_tanGhost);
          centerBox &= a_domain;
          const int gradVelComp =
            ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
          const FArrayBox& UfaceAvgFab = a_UfaceAvgFxb[gradDir];

          // FORT_CELLGRADDIR4THO(CHF_FRA1(a_GradUcellAvgFab,gradVelComp),
          //                      CHF_CONST_FRA1(a_UcellAvgFab,velDir),
          //                      CHF_BOX(centerBox),
          //                      CHF_CONST_INT(gradDir),
          //                      CHF_CONST_REAL(a_dxVect[gradDir]));
          FORT_CELLGRADFROMFACEAVG(
            CHF_FRA1(a_GradUcellAvgFab, gradVelComp),
            CHF_CONST_FRA1(UfaceAvgFab,velDir),
            CHF_BOX(centerBox),
            CHF_CONST_INT(gradDir),
            CHF_CONST_REAL(a_dxVect[gradDir]));
        }
    }
}

/*--------------------------------------------------------------------*/
//  Computes the face-averaged velocity gradients for the interior grid
/** 
 *  \param[out] a_GradUfaceAvgFxb
 *                      FluxBox of the face-averaged velocity gradients
 *  \param[in]  a_GradUcellAvgFab
 *                      FAB of the cell-averaged velocity gradients
 *  \param[in]  a_UcellAvgFab   
 *                      FAB of the cell-averaged velocities
 *  \param[in]  a_box   Interior box
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_dxVect
 *                      Vector of the computational grid spacing
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeInteriorFaceAvgVelGrad(
  FluxBox&           a_GradUfaceAvgFxb,
  const FArrayBox&   a_GradUcellAvgFab,
  const FArrayBox&   a_UcellAvgFab,
  const Box&         a_box,
  const BlockDomain& a_domain,
  const RealVect&    a_dxVect) const
{
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      IntVect growVect(c_tanGhost*IntVect::Unit);
      growVect[faceDir] = 0;
      FArrayBox& faceGVFab = a_GradUfaceAvgFxb[faceDir];
      // Normal gradients are solved up to and including boundaries
      Box solveBox(a_box);
      solveBox.grow(growVect);
      solveBox &= a_domain;
      solveBox.surroundingNodes(faceDir);
      switch (CRDparam::g_diffusiveDerivativeOrder)
        {
          // Second-order
        case 2:
          for (const int velDir : EachDir)
            {
              for (const int gradDir : EachDir)
                {
                  const Real factor = 1./a_dxVect[gradDir];
                  const int gradVelComp =
                    ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
                  // Compute wall normal gradients on tangent faces
                  if (gradDir == faceDir)
                    {
                      const int MD_ID(o, gradDir);
                      MD_BOXLOOP(solveBox, i)
                        {
                          faceGVFab[MD_IX(i,gradVelComp)] =
                            factor*(a_UcellAvgFab[MD_IX(i,velDir)] -
                                    a_UcellAvgFab[MD_OFFSETIX(i,-,o,velDir)]);
                        }
                    }
                  // Set tangential gradients on tangent faces
                  else
                    {
                      const int MD_ID(o, faceDir);
                      MD_BOXLOOP(solveBox, i)
                        {
                          faceGVFab[MD_IX(i,gradVelComp)] = 0.5*(
                            a_GradUcellAvgFab[MD_IX(i,gradVelComp)] + 
                            a_GradUcellAvgFab[MD_OFFSETIX(i,-,o,gradVelComp)]);
                        }
                    }
                }
            }
          break;
          // Fourth-order
        default:
          for (const int velDir : EachDir)
            {
              for (const int gradDir : EachDir)
                {
                  const int gradVelComp =
                    ViscousTensor4thOrderOp::tensorIdxRowOrder(velDir, gradDir);
                  // Compute wall normal gradients on tangent faces
                  if (gradDir == faceDir)
                    {
                      FORT_FACEGRADNORMDIR4THO(
                        CHF_FRA1(faceGVFab, gradVelComp),
                        CHF_CONST_FRA1(a_UcellAvgFab, velDir),
                        CHF_BOX(solveBox),
                        CHF_CONST_INT(gradDir),
                        CHF_CONST_REAL(a_dxVect[gradDir]));
                    }
                  // Set tangential gradients on tangent faces
                  else
                    {
                      FORT_FACEGRADTANGDIR4THO(
                        CHF_FRA1(faceGVFab, gradVelComp),
                        CHF_CONST_FRA1(a_GradUcellAvgFab, gradVelComp),
                        CHF_BOX(solveBox),
                        CHF_CONST_INT(faceDir));
                    }
                }
            }
          break;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the face-averaged product of N^T/J and the gradient of Phi
//  where Phi is either temperature or species concentration
/** Version of flux that allocates arrays for intermediate results on
 *  the stack
 *  \param[in]  a_box   The flux is computed on the faces of 'a_box'
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_domain5Pp2BCbox
 *                      A domain box grown by 5 in connected
 *                      directions and 2 in directions of domain BC
 *  \param[in]  a_domain5Pm1BCbox
 *                      A domain box grown by 5 in connected
 *                      directions and reduced by 1 in directions of
 *                      domain BC
 *  \param[out] a_NGradPhifaceAvgFxb
 *                      The face-averaged <N^T/J Grad Phi> on domain
 *                      and 1 ghost at interior and 0 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[out] a_GradWcellAvgFab
 *                      FAB containing the cell-averaged gradients of W
 *  \param[in]  a_WcellAvgFab
 *                      The cell-averaged primitive variables set
 *                      on domain and 5 ghost at interior and 2 ghosts
 *                      orthogonal to an exterior boundary.
 *  \param[in]  a_dataIndx
 *                      The data index for use with the metric terms
 *  \param[in]  a_indexPrim
 *                      The primitive component number within a_WfaceAvgFxb
 *                      and a_WcellAvgFxb to find gradients of 
 *  \param[in]  a_vecStart
 *                      First component to modify in a_NGradPhifacePntFxb, 
 *                      will modify component
 *                      Interval(a_vecStart, a_vecStart + SpaceDim - 1)
 *  \param[in]  a_outputStartComp
 *                      The first component in a_NGradPhifacePntFxb to 
 *                      solve for
 *  \param[in]  a_rowMult
 *                      If set to true, the rows of the matrix are 
 *                      multiplied by the vector, if false, the columns
 *                      are multiplied by the vector
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_onlyCellGrad
 *                      If set to true, only a_GradWcellAvgFab is filled
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeNGradPhi(
  const Box&         a_box,
  const BlockDomain& a_domain,
  const Box&         a_domain5Pp2BCbox,
  const Box&         a_domain5Pm1BCbox,
  FluxBox&           a_NGradPhifacePntFxb,
  FArrayBox&         a_GradWcellAvgFab,
  FArrayBox&         a_WcellAvgFab,
  FluxBox&           a_WfaceAvgFxb,
  const DataIndex&   a_dataIndx,
  const int          a_indexPrim,
  const int          a_vecStart,
  const bool         a_rowMult,
  const Real         a_time,
  const int          a_level,
  const bool         a_onlyCellGrad) const
{
  CH_TIME("ViscousTensor4thOrderOp::computeNGradPhi");
  CRD::msg << CRD::fv4 << "ViscousTensor4thOrderOp::computeNGradPhi" <<
    CRD::end;
  CH_assert(a_NGradPhifacePntFxb.nComp() > a_vecStart);
  CH_assert(a_WfaceAvgFxb.nComp() >= a_indexPrim);
  
  const RealVect dxVect = m_levelGridMetrics.dxVect();
  Box box4Dom2 = grow(a_box,4);
  box4Dom2 &= a_domain5Pp2BCbox;
  Box box4Dom = grow(a_box,4);
  box4Dom &= a_domain;

  Box box2Dom = grow(a_box,2);
  box2Dom &= a_domain;
  Box box2Dom2 = grow(a_box,2);
  box2Dom2 &= a_domain5Pp2BCbox;

  Box box1Dom = grow(a_box,1);
  box1Dom &= a_domain;

  Box phiFaceBox(box4Dom);
  Box phiCellBox(box4Dom2);
  Box gradBox(box2Dom2);

  // Alias FAB of the cell-averaged temperature
  Interval primInt(a_indexPrim, a_indexPrim);
  FArrayBox PhicellAvgFab(primInt, a_WcellAvgFab);
  // Alias FluxBox of the face-averaged temperature
  FluxBox PhifaceAvgFxb(a_WfaceAvgFxb.box(), 1,
                        D_DECL(a_WfaceAvgFxb[0].dataPtr(a_indexPrim),
                               a_WfaceAvgFxb[1].dataPtr(a_indexPrim),
                               a_WfaceAvgFxb[2].dataPtr(a_indexPrim)));

  // FAB of the cell-averaged scalar gradient tensor
  const int startGradComp = gradIdx(a_indexPrim, 0);
  const int endGradComp = gradIdx(a_indexPrim, SpaceDim - 1);
  Interval gradInt(startGradComp, endGradComp);
  FArrayBox GradPhicellAvgFab(gradInt, a_GradWcellAvgFab);

  // FluxBox of the face-averaged scalar gradient tensor
  FLUXBOXSTACKTEMP(GradPhifaceAvgFxb, gradBox, SpaceDim);
  // Solving phi now instead of the stress tensor
  bool solvePhi = true;
  if (!a_domain5Pm1BCbox.contains(a_box))
    {
      // Solve the cell-averaged gradients at the cells adjacent to
      // non-periodic boundaries
      CRDparam::g_CNSIBC->cellAvgGradBC(PhicellAvgFab,
                                        GradPhicellAvgFab,
                                        a_box,
                                        a_domain,
                                        dxVect,
                                        solvePhi,
                                        a_level);
    }

  // Compute the cell-averaged velocity gradients using cell-averaged phi
  // ---------------------------------------------------
  //  I/O           Var            Stencil      Ghosts
  // =====  ====================  ==========  ==========
  // read   T                      2 per dir    4 in dom
  // write  GradTcellAvgFab                1    2 in dom
  // ---------------------------------------------------
  computeInteriorCellAvgPhiGrad(GradPhicellAvgFab,
                                PhicellAvgFab,
                                PhifaceAvgFxb,
                                a_box,
                                a_domain,
                                dxVect);
  if (a_onlyCellGrad)
    {
      return;
    }
  // Compute the face-averaged temperature gradients using cell-averaged 
  // velocity and cell-averaged temperature gradients
  // ---------------------------------------------------
  //  I/O           Var            Stencil      Ghosts
  // =====  ====================  ==========  ==========
  // read   T                      2 per dir   4 i/e dom
  // read   GradTcellAvgFab        2 per dir   2 i/e dom
  // write  GradTfaceAvgFxb                1   2 i/e dom
  // ---------------------------------------------------
  computeInteriorFaceAvgPhiGrad(GradPhifaceAvgFxb,
                                GradPhicellAvgFab,
                                PhicellAvgFab,
                                a_box,
                                a_domain,
                                dxVect);

  // FluxBox of the face-centered temperature gradient
  FLUXBOXSTACKTEMP(GradPhifacePntFxb, box1Dom, SpaceDim);
  // Turn the face-averaged temperature gradients into face-centered
  // temperature gradients
  GradPhifacePntFxb.copy(GradPhifaceAvgFxb);
  MOLUtilFunc::deconvolveFace(GradPhifacePntFxb,
                              GradPhifaceAvgFxb,
                              box1Dom,
                              a_domain,
                              -1);

  // FluxBox of N^T/J from LevelGridMetrics class
  const FluxBox& NtJ = m_levelGridMetrics.m_NtJ[a_dataIndx];

  computeFOMatVecProd(box1Dom,
                      a_NGradPhifacePntFxb,
                      NtJ,
                      GradPhifacePntFxb,
                      a_vecStart,
                      a_rowMult);
}

/*--------------------------------------------------------------------*/
//  Computes the cell-averaged temperature gradients for the interior 
//  grid
/** 
 *  \param[out] a_GradPhicellAvgFab
 *                      FAB of the cell-averaged temperature gradients
 *  \param[in]  a_PhicellAvgFab
 *                      FAB of the cell-averaged temperature
 *  \param[in]  a_box   The gradients are computed in cells of 'a_box'
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_dxVect
 *                      Vector of the computational grid spacing
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeInteriorCellAvgPhiGrad(
  FArrayBox&         a_GradPhicellAvgFab,
  const FArrayBox&   a_PhicellAvgFab,
  const FluxBox&     a_PhifaceAvgFxb,
  const Box&         a_box,
  const BlockDomain& a_domain,
  const RealVect&    a_dxVect) const
{
  for (int gradDir = 0; gradDir != SpaceDim; ++gradDir)
    {
      Box centerBox(a_box);
      centerBox.grow(c_tanGhost);
      centerBox &= a_domain;
      const FArrayBox& PhifaceAvgFab = a_PhifaceAvgFxb[gradDir];
      // FORT_CELLGRADDIR4THO(CHF_FRA1(a_GradPhicellAvgFab,gradDir),
      //                      CHF_CONST_FRA1(a_PhicellAvgFab,0),
      //                      CHF_BOX(centerBox),
      //                      CHF_CONST_INT(gradDir),
      //                      CHF_CONST_REAL(a_dxVect[gradDir]));
      FORT_CELLGRADFROMFACEAVG(CHF_FRA1(a_GradPhicellAvgFab,gradDir),
                               CHF_CONST_FRA1(PhifaceAvgFab,0),
                               CHF_BOX(centerBox),
                               CHF_CONST_INT(gradDir),
                               CHF_CONST_REAL(a_dxVect[gradDir]));
    }
}

/*--------------------------------------------------------------------*/
//  Computes the face-averaged temperature gradients for the interior grid
/** 
 *  \param[out] a_GradPhifaceAvgFxb
 *                      FluxBox of the face-averaged temperature gradients
 *  \param[in]  a_GradPhicellAvgFab
 *                      FAB of the cell-averaged temperature gradients
 *  \param[in]  a_PhicellAvgFab
 *                      FAB of the cell-averaged temperature
 *  \param[in]  a_box   The gradients are computed on faces of 'a_box'
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_dxVect
 *                      Vector of the computational grid spacing
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeInteriorFaceAvgPhiGrad(
  FluxBox&           a_GradPhifaceAvgFxb,
  const FArrayBox&   a_GradPhicellAvgFab,
  const FArrayBox&   a_PhicellAvgFab,
  const Box&         a_box,
  const BlockDomain& a_domain,
  const RealVect&    a_dxVect) const
{
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      IntVect growVect(c_tanGhost*IntVect::Unit);
      // We will need data in the tangential directions of the faces later
      growVect[faceDir] = 0;
      FArrayBox& GradPhifaceFab = a_GradPhifaceAvgFxb[faceDir];
      // Normal gradients are solved up to and including boundaries
      Box solveBox(a_box);
      solveBox.grow(growVect);
      solveBox &= a_domain;
      solveBox.surroundingNodes(faceDir);
      switch (CRDparam::g_diffusiveDerivativeOrder)
        {
          // Second-order
        case 2:
          for (const int gradDir : EachDir)
            {
              // Compute wall normal gradients on tangent faces
              if (gradDir == faceDir)
                {
                  const int MD_ID(o, gradDir);
                  const Real factor = 1./a_dxVect[gradDir];
                  MD_BOXLOOP(solveBox, i)
                    {
                      GradPhifaceFab[MD_IX(i,gradDir)] =
                        factor*(a_PhicellAvgFab[MD_IX(i,0)] - 
                                a_PhicellAvgFab[MD_OFFSETIX(i,-,o,0)]);
                    }
                }
              // Set tangential gradients on tangent faces
              else
                {
                  const int MD_ID(o, faceDir);
                  MD_BOXLOOP(solveBox, i)
                    {
                      GradPhifaceFab[MD_IX(i,gradDir)] = 0.5*(
                        a_GradPhicellAvgFab[MD_IX(i,gradDir)] +
                        a_GradPhicellAvgFab[MD_OFFSETIX(i,-,o,gradDir)]);
                    }
                }
            }
          break;
          // Fourth-order
        default:
          for (const int gradDir : EachDir)
            {
              // Compute wall normal gradients on tangent faces
              if (gradDir == faceDir)
                {
                  FORT_FACEGRADNORMDIR4THO(
                    CHF_FRA1(GradPhifaceFab, gradDir),
                    CHF_CONST_FRA1(a_PhicellAvgFab, 0),
                    CHF_BOX(solveBox),
                    CHF_CONST_INT(gradDir),
                    CHF_CONST_REAL(a_dxVect[gradDir]));
                }
              // Set tangential gradients on tangent faces
              else
                {
                  FORT_FACEGRADTANGDIR4THO(
                    CHF_FRA1(GradPhifaceFab, gradDir),
                    CHF_CONST_FRA1(a_GradPhicellAvgFab, gradDir),
                    CHF_BOX(solveBox),
                    CHF_CONST_INT(faceDir));
                }
            }
          break;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Computes the face-averaged product of a matrix and a vector
/**
 *  \param[in]  a_box   Box over which to compute the product
 *  \param[out] a_outputFxb
 *                      FluxBox of the face-centered product
 *  \param[in]  a_MatfaceAvgFxb
 *                      FXB of the face-centered viscous stress tensor
 *  \param[in]  a_VecfaceAvgFxb
 *                      FXB of the face-centered velocity
 *  \param[in]  a_vecStart
 *                      First component to modify in a_outputFxb
 *  \param[in]  a_rowMult
 *                      If set to true, the rows of the matrix are 
 *                      multiplied by the vector, if false, the columns
 *                      are multiplied by the vector
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeFOMatVecProd(
  const Box&     a_box,
  FluxBox&       a_outputFxb,
  const FluxBox& a_MatfacePntFxb,
  const FluxBox& a_VecfacePntFxb,
  const int      a_vecStart,
  const bool     a_rowMult) const
{
  CH_TIME("ViscousTensor4thOrderOp::computeFOMatVecProd");
  CH_assert(a_MatfacePntFxb.nComp() > a_VecfacePntFxb.nComp());
  int startStride = 1;
  int stride = SpaceDim;
  if (a_rowMult)
    {
      startStride = SpaceDim;
      stride = 1;
    }
  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      FArrayBox& outputFab = a_outputFxb[faceDir];
      const FArrayBox& MatFab = a_MatfacePntFxb[faceDir];
      const FArrayBox& VecFab = a_VecfacePntFxb[faceDir];
      Box solveBox(a_box);
      solveBox.surroundingNodes(faceDir);

      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          int matStart = dir*startStride;
          int outputComp = dir + a_vecStart;
          FORT_MATVECPROD(CHF_FRA1(outputFab,outputComp),
                          CHF_CONST_FRA(MatFab),
                          CHF_CONST_FRA(VecFab),
                          CHF_BOX(solveBox),
                          CHF_CONST_INT(stride),
                          CHF_CONST_INT(matStart));
        }
    }
}

/*--------------------------------------------------------------------*/
//  Computes the face-averaged dot product between matrices over the
//  a_C.box()
/** 
 *  \param[out] a_C     FluxBox of the final dot product
 *                      Initalize this to zero outside the function
 *  \param[in]  a_A     FluxBox of the first matrix
 *  \param[in]  a_B     FluxBox of the first matrix
 *  \param[in]  a_varStart
 *                      Variable to start for A and B. If A and B are 
 *                      SpaceDim*SpaceDim, varStart should be 0.
 *                      If A and B are tangential derivatives of the flux dyads
 *                      varStart should be 0 or SpaceDim*SpaceDim depending on
 *                      if the first or second tangential gradient is being
 *                      accessed (only true for 3D)
 *//*-----------------------------------------------------------------*/

void
ViscousTensor4thOrderOp::computeMatrixProduct(FluxBox&       a_C,
                                              const FluxBox& a_A,
                                              const FluxBox& a_B,
                                              const int      a_varStart) const
{
  CH_assert(a_B.box().contains(a_C.box()));
  CH_assert(a_A.nComp() == a_B.nComp());
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      FArrayBox& cFab = a_C[dir];
      const FArrayBox& aFab = a_A[dir];
      const FArrayBox& bFab = a_B[dir];
      
      Box faceBox(a_C.box());
      
      faceBox.surroundingNodes(dir);
      FORT_POINTDOTMATPROD(CHF_FRA(cFab),
                           CHF_CONST_FRA(aFab),
                           CHF_CONST_FRA(bFab),
                           CHF_CONST_INT(a_varStart),
                           CHF_BOX(faceBox));
    }
}


/*--------------------------------------------------------------------*/
//  Retrieves the component of the space contiguous flux dyad
/** 
 *  \param[in]  a_space The space dimension to retrieve (starting from 0)
 *  \param[in]  a_var   The variable number to retrieve (starting from 0)
 *//*-----------------------------------------------------------------*/

int
ViscousTensor4thOrderOp::retrieveFluxDyadComp(const int a_var,
                                              const int a_space) const
{
  const int nFlux = CRDparam::g_CRDPhysics->numFluxes();
  CH_assert(a_space <= SpaceDim - 1);
  CH_assert(a_var <= nFlux - 1);
  
  int comp = a_var*SpaceDim + a_space;
  
  return comp;
}
