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
 * \file SSV.cpp
 *
 * \brief Member functions for SSV
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "CONSTANTS.H"
#include "MOLUtilFunc.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "FourthOrderUtil.H"
#include "PointwiseDotProdF_F.H"
#include "UnitNormalsF_F.H"
#include "RootSolver.H"
#include "AMRIO.H"

//----- Internal -----//

#include "SSV.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDutil.H"
#include "CRDPhysics.H"
#include "PatchCNSOp.H"
#include "CNSIBC.H"
#include "DataTemp.H"
#include "PatchMappedFunc.H"

/*******************************************************************************
 *
 * Class SSV: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

SSV::SSV()
  :
  m_vortexProportion(0),
  m_vortexModel(0),
  m_sgsKineticEnergyApprox(0),
  m_useWallModel(0),
  m_useOnlyWallModel(0),
  m_bndryVSTCorrection(0),
  m_virtualWallHeight(0.18),
  m_enforceWallShearStress(0),
  m_suppressWallModel(0),
  m_sgsEnergyFlux(1),
  m_constantK1(0),
  m_etaMaxBound(1e50)
{
  ParmParse ppTURB("turb");
  if (ppTURB.contains("sgs_kinetic_energy_approximation"))
    {
      ppTURB.get("sgs_kinetic_energy_approximation",
                 m_sgsKineticEnergyApprox);
      if ((m_sgsKineticEnergyApprox < 0) ||
          (m_sgsKineticEnergyApprox > 1))
        {
          m_sgsKineticEnergyApprox = 0;
          CRD::msg << "LES SV SGS model: sgs_kinetic_energy_approximation must" 
                   << " be >= 0 and < 2,"
                   << " default is 0, setting to 0" << CRD::error;
        }
    }
  Real vortexProportion = 1.;
  ppTURB.query("vortex_proportion", vortexProportion);
  if ((vortexProportion >= 0.) && (vortexProportion <= 1.))
    {
      m_vortexProportion = vortexProportion;
    }
  else
    {
      m_vortexProportion = 1.;
      CRD::msg << "LES SV SGS model: vortex_proportion must be >= 0 and "
               << "<= 1, default is 1, setting to 1" << CRD::error;
    }
  int vortexModel = 0;
  ppTURB.query("vortex_model", vortexModel);
  if ((vortexModel >= 0) && (vortexModel <= 1))
    {
      m_vortexModel = vortexModel;
    }
  else
    {
      CRD::msg << "LES SV SGS model: vortex_model must be 0 or 1" << CRD::error;
    }
  if (ppTURB.contains("use_wall_model"))
    {
      ppTURB.get("use_wall_model", m_useWallModel);
      if ((m_useWallModel < 0) || (m_useWallModel > 1))
        {
          m_useWallModel = 0;
          CRD::msg << "SV SGS LES model: use_wall_model must be >=0 and "
                   << "<= 1, default is 0, setting to 0 with no wall-model"
                   << CRD::var;
        }
    }
  else
    {
      CRD::msg << "SV SGS LES model: default does not use wall-model. "
               << "If a wall-model is desired, please restart and specify one."
               << CRD::var;
    }
  if (ppTURB.contains("only_use_wall_model"))
    {
      ppTURB.get("only_use_wall_model", m_useOnlyWallModel);
      if ((m_useOnlyWallModel < 0) || (m_useOnlyWallModel > 1))
        {
          m_useOnlyWallModel = 0;
          CRD::msg << "SV SGS LES model: only_use_wall_model must be >=0 and "
                   << "<= 1, default is 0, setting to 0"
                   << CRD::var;
        }
    }
  if (ppTURB.contains("use_wall_model_vst_bndry_correction"))
    {
      ppTURB.get("use_wall_model_vst_bndry_correction", m_bndryVSTCorrection);
      if ((m_bndryVSTCorrection < 0) || (m_bndryVSTCorrection > 1))
        {
          m_bndryVSTCorrection = 0;
          CRD::msg << "SV SGS LES model: use_wall_model_vst_bndry_correction "
                   << "must be >=0 and <= 1, default is 0, setting to 0"
                   << CRD::var;
        }
    }
  Real virtualWallHeight = 0.18;
  ppTURB.query("virtual_wall_height", virtualWallHeight);
  if ((virtualWallHeight >= 0) && (virtualWallHeight <= 1))
    {
      m_virtualWallHeight = virtualWallHeight;
    }
  ppTURB.query("enforce_wall_model_shear_stress", m_enforceWallShearStress);
  if ((m_enforceWallShearStress < 0) || (m_enforceWallShearStress > 1))
    {
      CRD::msg << "LES SV SGS wall-model: enforce_wall_model_shear_stress"
               << " must be 0 or 1 -- setting to 0" << CRD::var;
    }
  ppTURB.query("suppress_wall_model", m_suppressWallModel);
  if ((m_suppressWallModel < 0) || (m_suppressWallModel > 1))
    {
      CRD::msg << "LES SV SGS wall-model: suppress_wall_model"
	       << " must be 0 or 1 -- setting to 0" << CRD::var;
      m_suppressWallModel = 0;
    }
  ppTURB.query("use_sgs_energy_model", m_sgsEnergyFlux);
  if ((m_sgsEnergyFlux < 0) || (m_sgsEnergyFlux > 1))
    {
      CRD::msg << "LES SV SGS energy-model: use_sgs_energy_model"
               << " must be 0 or 1 -- setting to 1" << CRD::var;
      m_sgsEnergyFlux = 1;
    }
  ppTURB.query("use_constant_wall_model_K1", m_constantK1);
  if ((m_constantK1 < 0) || (m_constantK1 > 1))
    {
      CRD::msg << "LES SV SGS wall-model: use_constant_wall_model_K1"
               << " must be 0 or 1 -- setting to 0" << CRD::var;
      m_constantK1 = 0;
    }
  ppTURB.query("wall_model_eta_max_bound", m_etaMaxBound);
  if (m_etaMaxBound <= 0)
    {
      CRD::msg << "LES SV SGS wall-model: wall_model_eta_max_bound"
               << " must be > 0" << CRD::error;
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

SSV::~SSV()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Calculate the momentum flux subgrid-scale correction
/** \param[out] a_sgsMomentumFxb
 *                      Subgrid-scale momentum fluxes on faces
 *  \param[out] a_sgsEnergyFxb
 *                      Subgrid-scale energy fluxes on faces
 *  \param[in]  a_WfacePntFxb
 *                      Face-centered primitive variables
 *  \param[in]  a_NGradUfacePntFxb
 *                      Face-centered physical-space velocity gradients
 *  \param[in]  a_strainRateTensorFxb
 *                      Face-centered, physical-space 
 *                      strain-rate tensor
 *  \param[in]  a_NGradTfacePntFxb
 *                      Face-centered, physical-space
 *                      temperature gradient
 *  \param[in]  a_facePntDeltaC
 *                      Face-centered cell-cutoff length for LES
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_dataIndx
 *                      Data index for the grid information
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Cell-centered box with faces over which the
 *                      subgrid-scale flux must be computed
 *//*-----------------------------------------------------------------*/

void
SSV::sgsModelFlux(FluxBox&                a_sgsMomentumFxb,
                  FluxBox&                a_sgsEnergyFxb,
                  FluxBox&                a_faceAvgPlotFxb,
                  const FluxBox&          a_WfacePntFxb,
                  const FArrayBox&        a_WcellAvgFab,
                  const FluxBox&          a_NGradUfacePntFxb,
                  const FluxBox&          a_strainRateTensorFxb,
                  const FluxBox&          a_NGradTfacePntFxb,
                  const FluxBox&          a_facePntDeltaC,
                  const FluxBox&          a_faceCoord,
                  const FluxBox&          a_unitNormalsDirFxb,
                  const LevelGridMetrics& a_gridMetrics,
                  const DataIndex&        a_dataIndx,
                  const ProblemDomain&    a_domain,
                  const Box&              a_box) const
{
  CH_TIME("SSV::sgsModelFlux");

  // If only the wall-model is being used, there's nothing to do here
  if (m_useOnlyWallModel) { return; }

  Box disjointBox = a_gridMetrics.getBoxes()[a_dataIndx];
  Box box1Dom = grow(a_box,1);
  box1Dom    &= a_domain;

  // Global constants
  const int cRho  = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cTemp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cKE   = CRDparam::g_CRDPhysics->turbPrimInterval().begin();
  const int cSpecies = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numTComp = SpaceDim*SpaceDim;           // number of tensor comps

  // Memory allocation
  FLUXBOXSTACKTEMP(eigenVect, box1Dom, SpaceDim);      // eigenvector
  FLUXBOXSTACKTEMP(orientationFxb, box1Dom, numTComp); // SGS-vortex orientation
  FLUXBOXSTACKTEMP(sgsFxb, box1Dom, 1);

  // If using no coarsening, computing SGS KE on faces is 3 times as expensive
  // as in the cells (3 unique faces per cell), but has demonstrated higher
  // quality in regions of high SGS KE gradient (e.g. within boundary layers).
  // The following also tests the use of localized coarsening.
  if (false)
    {
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Box faceBox = grow(a_box, 1);
          faceBox.grow(dir, -1);
          faceBox &= a_domain;
          faceBox.surroundingNodes(dir);
          // Extend interpolated face state (necessary for lots of coarsening)
          Box interpBox = faceBox;
          interpBox.grow(CRDparam::g_sgskeCrsLevFilterRatio);
          interpBox &= a_domain;
          FABSTACKTEMP(WfacePntFabTemp,
                       interpBox,
                       a_WfacePntFxb[dir].nComp());
          WfacePntFabTemp.copy(a_WfacePntFxb[dir]);
          if (true)
            {
              for (int velComp = 0; velComp != SpaceDim; ++velComp)
                {
                  int faceInterpOrder = CRDparam::g_faceInterpolationOrder;
                  PatchMappedFunc::faceAvgValFromCellAvgCS(
                    WfacePntFabTemp, a_WcellAvgFab, cVel+velComp, cVel+velComp,
                    a_domain, interpBox, dir, faceInterpOrder, false);
                }
            }
          int localCoarsening = 1;
          // Compute subgrid-scale kinetic energy
          sgsKineticEnergy(sgsFxb[dir],
                           eigenVect[dir],
                           a_NGradUfacePntFxb[dir],
                           a_strainRateTensorFxb[dir],
                           WfacePntFabTemp,
                           a_facePntDeltaC[dir],
                           a_faceCoord[dir],
                           a_domain,
                           faceBox,
                           localCoarsening);
        }
    }

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box faceBox = grow(a_box, 1);
      faceBox.grow(dir, -1);
      faceBox &= a_domain;
      faceBox.surroundingNodes(dir);
      // Compute eigenvectors/eigenvalues here
      eigenDecomp(eigenVect[dir], a_strainRateTensorFxb[dir], faceBox);
      // Compute SGS-vortex orientation vector
      sgsVortexVector(orientationFxb[dir],
                      a_NGradUfacePntFxb[dir],
                      eigenVect[dir],
                      a_gridMetrics,
                      a_domain,
                      faceBox);
    }

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // If the wall-model is being used, we compute the model on all faces,
      // even boundaries. The rationale is that the wall-model
      // replaces the no-slip condition with a slip-condition and the wall
      // is modeled as being inside the flow region at a lifted virtual wall.
      // From a practical perspective, the wall-transpiration of momentum
      // through viscosity is removed with the slip-velocity condition and must
      // be replaced with something else (the interior model in this case).
      // If the wall-model is not used, walls do not have the SGS component
      // applied as there is still a no-slip velocity condition being used.
      Box faceBox = grow(a_box, 1);
      faceBox.grow(dir, -1);
      faceBox &= a_domain;
      faceBox.surroundingNodes(dir);

      // Check for no-slip wall boundaries here
      if (!m_useWallModel)
        {
          for (int side = -1; side != 2; ++side) // There are 2 sides: -1,+1
            {
              if (!side) { continue; } // Side "0" doesn't exist
              Box testBox = faceBox;
              testBox.shift(dir, side); // Shift testBox outside domain
              if (!a_domain.contains(testBox))
                {
                  int isWall = 0;
                  Side::LoHiSide whichSide = (side == 1) ? Side::Hi : Side::Lo;
                  hasNoSlipWall(isWall,a_gridMetrics,disjointBox,dir,whichSide);
                  if (isWall) { faceBox.growDir(dir, whichSide, -1); }
                }
            }
        }

      const FArrayBox& arrSGS = sgsFxb[dir];
      const FArrayBox& arrW  = a_WfacePntFxb[dir];
      const FArrayBox& arrO  = orientationFxb[dir];
      const FArrayBox& arrT  = a_NGradTfacePntFxb[dir];
      const FArrayBox& arrDC = a_facePntDeltaC[dir];
      FArrayBox& arrF        = a_sgsMomentumFxb[dir];
      FArrayBox& arrE        = a_sgsEnergyFxb[dir];

      // If not using face-computation of SGS KE, fill with interpolated values
      if (true)
        {
          MD_BOXLOOP(faceBox, i)
            {
              arrSGS[MD_IX(i, 0)] = arrW[MD_IX(i, cKE)];
            }
        }
      // Remove the (11) redundant max calls and (2) redundant sqrt calls
      FABSTACKTEMP(absSGSKE, faceBox, 1);
      FABSTACKTEMP(sqrtSGSKE, faceBox, 1);
      MD_BOXLOOP(faceBox, i)
        {
          absSGSKE[MD_IX(i, 0)] = std::max(0., arrSGS[MD_IX(i, 0)]);
          sqrtSGSKE[MD_IX(i, 0)] = std::sqrt(absSGSKE[MD_IX(i, 0)]);
        }
      // Compute the SGS momentum flux
      for (int comp = 0; comp != numTComp; ++comp)
        {
          MD_BOXLOOP(faceBox, i)
            {
              arrF[MD_IX(i,comp)] =
                arrW[MD_IX(i,cRho)]*absSGSKE[MD_IX(i,0)]*arrO[MD_IX(i,comp)];
            }
        }
      // Compute the SGS energy flux
      if (m_sgsEnergyFlux)
        {
          FABSTACKTEMP(gammaFab, faceBox, 1);
          FABSTACKTEMP(RgasFab, faceBox, 1);
          gammaFab.setVal(CRDparam::g_gamma);
          RgasFab.setVal(CRDparam::g_R);
          // Compute multispecies form of gamma and Rgas
          if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
            {
              std::vector<Real> cn(CRDparam::g_numSpecies, 0.);
              MD_BOXLOOP(faceBox, i)
                {
                  for (int comp = 0; comp != CRDparam::g_numSpecies; ++comp)
                    {
                      const int wComp = comp + cSpecies;
                      cn[comp] = arrW[MD_IX(i, wComp)];
                    }
                  gammaFab[MD_IX(i, 0)] = CRDparam::g_CRDPhysics->gamma(
                    arrW[MD_IX(i, cTemp)],cn.data());
                  RgasFab[MD_IX(i, 0)] =
                    CRDparam::g_CRDPhysics->Rgas(cn.data());
                }
            }
          // Put the SGS energy flux components together
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              MD_BOXLOOP(faceBox, i)
                {
                  const Real Rgas = RgasFab[MD_IX(i, 0)];
                  const Real gamma = gammaFab[MD_IX(i, 0)];
                  const Real c_1 = -0.5*Rgas*gamma/(gamma - 1.);
                  arrE[MD_IX(i,comp)] = c_1*arrW[MD_IX(i,cRho)]*
                    sqrtSGSKE[MD_IX(i,0)]*arrDC[MD_IX(i,0)]*
                    (D_TERM(arrO[MD_IX(i,srtIdx(comp,0))]*arrT[MD_IX(i,0)],
                          + arrO[MD_IX(i,srtIdx(comp,1))]*arrT[MD_IX(i,1)],
                          + arrO[MD_IX(i,srtIdx(comp,2))]*arrT[MD_IX(i,2)]));
                }
            }
        }
    }

  //**NOTE: Most of this should go in a plotting/data-processing utility
  // Only the upper triangle of the SGS stress tensor is needed
  int plotVariables = CRDparam::g_plotLoFaceAvgComps;
  plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
  if ((plotVariables & CRDparam::PlotTurbulentComps) &&
      (CNSIBC::s_firstRKStage))
    {
      int cStart = 0; // Yes, we'll use indexing within a plot class or
                      // another improved method --- as it is, this is annoying
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          cStart = CRDparam::g_CRDPhysics->numPrimitive();
        }
      if ((plotVariables & CRDparam::PlotFluxes) &&
          (CRDparam::g_physicsModels & CRDparam::PhysicsInertial))
        {
          cStart += CRDparam::g_CRDPhysics->numFluxes();
        }
      if ((plotVariables & CRDparam::PlotFluxes) &&
          (CRDparam::g_physicsModels & CRDparam::PhysicsViscous))
        {
          cStart += CRDparam::g_CRDPhysics->numFluxes();
        }
      cStart += SpaceDim*(SpaceDim - 1.)/2. + SpaceDim; // <u_i*u_j>
      cStart += 3; // <p^2>, <T^2>, <rho^2>

      int numUpperComps = SpaceDim*(SpaceDim - 1.)/2. + SpaceDim;
      FLUXBOXSTACKTEMP(sgsFacePntFxb, box1Dom, numUpperComps);
      FLUXBOXSTACKTEMP(sgsFaceAvgFxb, a_box, numUpperComps);
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          IntVect growVect = IntVect_unit;
          growVect[faceDir] = 0;
          Box faceBox = a_box;
          faceBox.grow(growVect); // Grow in tangential directions only
          faceBox &= a_domain;
          faceBox.surroundingNodes(faceDir);
          // Compute tau_sgs_ij
          int cLoc = 0;
          FArrayBox& sgsFacePntFab = sgsFacePntFxb[faceDir];
          FArrayBox& arrF = a_sgsMomentumFxb[faceDir];
          const FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  int tensorComp = dir + sComp*SpaceDim;
                  MD_BOXLOOP(faceBox, i)
                    {
                      Real rho = std::max(1.e-20, WfacePntFab[MD_IX(i, cRho)]);
                      Real sgsVal = arrF[MD_IX(i, tensorComp)];
                      sgsFacePntFab[MD_IX(i, cLoc)] = sgsVal/rho;
                    }
                  ++cLoc;
                }
            }
        }
      // Convolve tau_sgs_ij
      CRDutil::convolveFace(sgsFaceAvgFxb, sgsFacePntFxb, a_box, a_domain,
                            Interval(0, numUpperComps-1), 4, false, false);
      // Fill a_faceAvgPlotFxb with <tau_sgs_ij>
      a_faceAvgPlotFxb.copy(sgsFaceAvgFxb, 0, cStart, numUpperComps);

      // Compute terms related to tau_w if wall-modeling
      if (m_useWallModel)
        {
          FLUXBOXSTACKTEMP(facePntVelFxb, a_box, SpaceDim);
          facePntVelFxb.copy(a_WfacePntFxb, cVel, 0, SpaceDim);
          // faceTauFxb contains <mag_tau_w>, <tau_w_streamwise>, <tau_w^2>
          FLUXBOXSTACKTEMP(faceTauFxb, a_box, 3);
          for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
            {
              FArrayBox& facePntVelFab = facePntVelFxb[faceDir];
              FArrayBox& faceTauFab = faceTauFxb[faceDir];
              const FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
              const FArrayBox& unitNormalsDirFab = a_unitNormalsDirFxb[faceDir];
              Box faceBox = a_box;
              faceBox.surroundingNodes(faceDir);
              // Transform face velocity into face-normal space
              FORT_FORWARDTRANSFORMF(CHF_FRA(facePntVelFab),
                                     CHF_CONST_FRA(unitNormalsDirFab),
                                     CHF_BOX(faceBox));
              // Set the wall-normal component to zero
              facePntVelFab.setVal(0., faceDir);
              // Normalize the vector
              PatchMappedFunc::normalize(faceBox, facePntVelFab,
                                         Interval(0, SpaceDim-1));
              // Transform global streamwise direction into face-normal space
              FABSTACKTEMP(globStreamLocTanVect, faceBox, SpaceDim);
              // Set components to match user-defined vector
              const RealVect streamwisePhysVect = CRDparam::g_globStreamwiseDir;
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  globStreamLocTanVect.setVal(streamwisePhysVect[dir], dir);
                }
              // Transform globStreamLocTanVect into local normal-tangent space
              FORT_FORWARDTRANSFORMF(CHF_FRA(globStreamLocTanVect),
                                     CHF_CONST_FRA(unitNormalsDirFab),
                                     CHF_BOX(faceBox));
              // Zero the normal component of globStreamLocTanVect
              globStreamLocTanVect.setVal(0., faceDir);
              // Normalize globStreamLocTanVect
              PatchMappedFunc::normalize(faceBox, globStreamLocTanVect,
                                         Interval(0, SpaceDim-1));
              // Project velocity onto streamwise vector and multiply by tau_w
              const int cTurbStart =
                CRDparam::g_CRDPhysics->turbConsInterval().begin();
              const int numSGSKEComps = 1;
              const int cEta = cTurbStart + numSGSKEComps + faceDir;
              MD_BOXLOOP(faceBox, i)
                {
                  const Real eta = WfacePntFab[MD_IX(i, cEta)];
                  RealVect vVect(D_DECL(facePntVelFab[MD_IX(i, 0)],
                                        facePntVelFab[MD_IX(i, 1)],
                                        facePntVelFab[MD_IX(i, 2)]));
                  RealVect sVect(D_DECL(globStreamLocTanVect[MD_IX(i, 0)],
                                        globStreamLocTanVect[MD_IX(i, 1)],
                                        globStreamLocTanVect[MD_IX(i, 2)]));
                  faceTauFab[MD_IX(i, 0)] = eta;
                  faceTauFab[MD_IX(i, 1)] = eta*(vVect.dotProduct(sVect));
                  faceTauFab[MD_IX(i, 2)] = eta*eta;
                }
            }
          // Fill a_faceAvgPlotFxb with additional components
          cStart += SpaceDim*(SpaceDim - 1.)/2. + SpaceDim;
          a_faceAvgPlotFxb.copy(faceTauFxb, 0, cStart, 3);
        }
    }
}


/*--------------------------------------------------------------------*/
//  Calculate the strain-rate tensor eigenvalues and eigenvectors
/** \param[out] a_eigenVect
 *                      Strongest eigenvector of strain-rate tensor
 *  \param[in]  a_strainRateTensorFab
 *                      Face-centered strain-rate tensor
 *  \param[in]  a_box   Cells in which to perform eigendecomposition
 *  \return             Eigenvector stored as row vector
 *//*-----------------------------------------------------------------*/

void
SSV::eigenDecomp(FArrayBox&       a_eigenVect,
                 const FArrayBox& a_strainRateTensorFab,
                 const Box&       a_box) const
{
  CH_TIME("SSV::eigenDecomp");
#if CH_SPACEDIM == 3
  const Real epsilon = 1.e-20;
  const Real epsilon2 = 1.e-4;
  const Real cosTheta2 = std::cos(2.*PI/3.);
  const Real cosTheta4 = std::cos(4.*PI/3.);
  const Real sinTheta2 = std::sin(2.*PI/3.);
  const Real sinTheta4 = std::sin(4.*PI/3.);
  MD_ARRAY_RESTRICT(arrS,a_strainRateTensorFab);
  const Real hiTol = 1.e100;
  const Real loTol = -1.e100;
  MD_BOXLOOP(a_box, i)
    {
      D_TERM(const Real Sxx = arrS[MD_IX(i,srtIdx(0,0))];
             CH_assert(Sxx < hiTol && Sxx > loTol);
             CH_assert(!std::isnan(Sxx) && !std::isinf(Sxx));,
             const Real Syx = arrS[MD_IX(i,srtIdx(1,0))];
             CH_assert(Syx < hiTol && Syx > loTol);
             CH_assert(!std::isnan(Syx) && !std::isinf(Syx));
             const Real Syy = arrS[MD_IX(i,srtIdx(1,1))];
             CH_assert(Syy < hiTol && Syy > loTol);
             CH_assert(!std::isnan(Syy) && !std::isinf(Syy));,
             const Real Szx = arrS[MD_IX(i,srtIdx(2,0))];
             CH_assert(Szx < hiTol && Szx > loTol);
             CH_assert(!std::isnan(Szx) && !std::isinf(Szx));
             const Real Szy = arrS[MD_IX(i,srtIdx(2,1))];
             CH_assert(Szy < hiTol && Szy > loTol);
             CH_assert(!std::isnan(Szy) && !std::isinf(Szy));
             const Real Szz = arrS[MD_IX(i,srtIdx(2,2))];
             CH_assert(Szz < hiTol && Szz > loTol);
             CH_assert(!std::isnan(Szz) && !std::isinf(Szz)););
      // lambda^3 - a*lambda^2 - b*lambda - c = 0
      const Real a = (-1.)*(D_TERM(Sxx, + Syy, + Szz));
      const Real b = (-1.)*
        (D_TERM(,Syx*Syx - Sxx*Syy, + Szx*Szx + Szy*Szy - Sxx*Szz - Syy*Szz));
      const Real c = (-1.)*
        (D_SELECT(0.,0.,Szx*(Syx*Szy - Szx*Syy) + Szy*(Szx*Syx - Sxx*Szy)
                  + Szz*(Sxx*Syy - Syx*Syx)));
      const Real Q = (3.*b - a*a)/9.;
      const Real R = (9.*a*b - 27.*c - 2.*a*a*a)/54.;

      // NOTE: Q is <= 0 if S is symmetric (we enforce this algorithmically)
      const Real sqrtQ = std::sqrt(-Q + epsilon);
      const Real costheta = R/((-Q + epsilon)*sqrtQ);
      Real theta = 0.;
      if (costheta > 1.)
        {
          theta = 0.;
        }
      else if (costheta < -1.)
        {
          theta = PI;
        }
      else
        {
          theta = std::acos(costheta);
        }
      const Real cosTheta3 = std::cos(theta/3.);
      const Real sinTheta3 = std::sin(theta/3.);
      const Real eigenVal1 = 2.*sqrtQ*cosTheta3 - (a/3.);
      const Real eigenVal2 = 2.*sqrtQ*(cosTheta3*cosTheta2
                                     - sinTheta3*sinTheta2) - (a/3.);
      const Real eigenVal3 = 2.*sqrtQ*(cosTheta3*cosTheta4
                                     - sinTheta3*sinTheta4) - (a/3.);
#ifndef NDEBUG
      // Note: Some quick checks to make sure the eigenvalues are correct
      // (1) The trace of S is the sum of the eigenvalues
      // (2) The determinant of S is the product of the eigenvalues
      //     NOTE: Test (2) is sensitive to the magnitude of the eigenvalues
      //           -- don't use this unless it's updated
      // (3) The determinant of (S - lambda*I) = 0
      const Real trace = -a; // trace of S_ij
      Real norm = 0;
      for (int j = 0; j != SpaceDim; ++j)
        {
          for (int k = 0; k != SpaceDim; ++k)
            {
              norm += a_strainRateTensorFab[MD_IX(i,srtIdx(j,k))]*
                a_strainRateTensorFab[MD_IX(i,srtIdx(j,k))];
            }
        }
      norm = std::sqrt(norm);
      const Real determ_S = D_SELECT(0.,0.,
                            Sxx*(Syy*Szz - Szy*Szy) - Syx*(Syx*Szz - Szx*Szy)
                               + Szx*(Syx*Szy - Syy*Szx));
      const Real determ_Eigen1 = D_SELECT(0.,0.,
              (Sxx - eigenVal1)*((Syy - eigenVal1)*(Szz - eigenVal1) - Szy*Szy)
             - Syx*(Syx*(Szz - eigenVal1) - Szx*Szy)
             + Szx*(Syx*Szy - (Syy - eigenVal1)*Szx));
      const Real determ_Eigen2 = D_SELECT(0.,0.,
              (Sxx - eigenVal2)*((Syy - eigenVal2)*(Szz - eigenVal2) - Szy*Szy)
             - Syx*(Syx*(Szz - eigenVal2) - Szx*Szy)
             + Szx*(Syx*Szy - (Syy - eigenVal2)*Szx));
      const Real determ_Eigen3 = D_SELECT(0.,0.,
              (Sxx - eigenVal3)*((Syy - eigenVal3)*(Szz - eigenVal3) - Szy*Szy)
             - Syx*(Syx*(Szz - eigenVal3) - Szx*Szy)
             + Szx*(Syx*Szy - (Syy - eigenVal3)*Szx));

      const Real prod = eigenVal1*eigenVal2*eigenVal3;
      const Real sum = eigenVal1 + eigenVal2 + eigenVal3;
      const Real normed_determ1 = determ_Eigen1/(norm + epsilon);
      const Real normed_determ2 = determ_Eigen2/(norm + epsilon);
      const Real normed_determ3 = determ_Eigen3/(norm + epsilon);
      const Real diff_determProd = std::abs(determ_S - prod);
      const Real diff_traceSum = std::abs(trace - sum);
      if ((normed_determ1 > epsilon2) || (normed_determ2 > epsilon2) ||
          (normed_determ3 > epsilon2) || (diff_traceSum  > epsilon2))
        {
          const RealVect eigenVals(D_DECL(eigenVal1, eigenVal2, eigenVal3));
          CRD::msg << "SSV::eigenDecomp: eigenvalues incorrect!"
                   << " Eigenvalues of S = " << eigenVals
                   << " , Trace of S - eigenvalue sum = "
                   << diff_traceSum << " and must be < " << epsilon2
                   << " , Determ_S - eigenvalue product = "
                   << diff_determProd << " and should be < " << epsilon2
                   << " , Normed_determ1 = " << normed_determ1
                   << " and must be < " << epsilon2
                   << " , Normed_determ2 = " << normed_determ2
                   << " and must be < " << epsilon2
                   << " , Normed_determ3 = " << normed_determ3
                   << " and must be < " << epsilon2
                   << CRD::abort;
        }
#endif
      // Eigenvector computation
      const RealVect subDet(
        D_DECL((D_SELECT(0.,0.,
                         (Syy - eigenVal1)*(Szz - eigenVal1) - (Szy*Szy))),
               (D_SELECT(0.,0.,
                         (Szz - eigenVal1)*(Sxx - eigenVal1) - (Szx*Szx))),
               (D_SELECT(0.,0.,
                         (Sxx - eigenVal1)*(Syy - eigenVal1) - (Syx*Syx)))));
      const RealVect absSubDet(D_DECL(std::abs(subDet[0]),
                                      std::abs(subDet[1]),
                                      std::abs(subDet[2])));
      RealVect eVec(RealVect::Zero);
      if ((absSubDet[0] >= absSubDet[1]) &&
          (absSubDet[0] >= absSubDet[2]))
        {
          D_SELECT(,,
                   eVec[0] = 1.;
                   eVec[1] = (-Syx*(Szz - eigenVal1) + (Szx*Szy))/subDet[0];
                   eVec[2] = (-Szx*(Syy - eigenVal1) + (Syx*Szy))/subDet[0];);
        }
      else if ((absSubDet[1] >= absSubDet[0]) &&
               (absSubDet[1] >= absSubDet[2]))
        {
          D_SELECT(,,
                   eVec[0] = (-Syx*(Szz - eigenVal1) + (Szy*Szx))/subDet[1];
                   eVec[1] = 1.;
                   eVec[2] = (-Szy*(Sxx - eigenVal1) + (Syx*Szx))/subDet[1];);
        }
      else if ((absSubDet[2] >= absSubDet[0]) &&
               (absSubDet[2] >= absSubDet[1]))
        {
          D_SELECT(,,
                   eVec[0] = (-Szx*(Syy - eigenVal1) + (Szy*Syx))/subDet[2];
                   eVec[1] = (-Szy*(Sxx - eigenVal1) + (Szx*Syx))/subDet[2];
                   eVec[2] = 1.;);
        }
      else
        {
          pout() << "Face/Cell = " << MD_GETIV(i) << std::endl;
          CRD::msg << "SSV::eigenDecomp: eigenvector problem" << "\n"
                   << D_TERM(
                     "Sxx = " << Sxx << "\n", <<
                     "Syx = " << Syx << ", Syy = " << Syy << "\n", <<
                     "Szx = " << Szx << ", Szy = " << Szy << ", Szz = " << Szz)
                   << CRD::abort;
        }
      // Normalize final eigenvector
      const Real eigenSize = std::sqrt(D_TERM(eVec[0]*eVec[0],
                                            + eVec[1]*eVec[1],
                                            + eVec[2]*eVec[2]) + epsilon);
      D_TERM(a_eigenVect[MD_IX(i,0)] = eVec[0]/eigenSize;,
             a_eigenVect[MD_IX(i,1)] = eVec[1]/eigenSize;,
             a_eigenVect[MD_IX(i,2)] = eVec[2]/eigenSize;);
    }
#endif
#if CH_SPACEDIM == 2
  const Real epsilon = 1.e-20;
  const Real hiTol = 1.e100;
  const Real loTol = -1.e100;
  MD_BOXLOOP(a_box, i)
    {
      D_TERM(const Real Sxx = a_strainRateTensorFab[MD_IX(i,srtIdx(0,0))];
             CH_assert(Sxx < hiTol && Sxx > loTol);
             CH_assert(!std::isnan(Sxx) && !std::isinf(Sxx));,
             const Real Syx = a_strainRateTensorFab[MD_IX(i,srtIdx(1,0))];
             CH_assert(Syx < hiTol && Syx > loTol);
             CH_assert(!std::isnan(Syx) && !std::isinf(Syx));
             const Real Syy = a_strainRateTensorFab[MD_IX(i,srtIdx(1,1))];
             CH_assert(Syy < hiTol && Syy > loTol);
             CH_assert(!std::isnan(Syy) && !std::isinf(Syy));,);
      // Trace of the 2x2 matrix
      const Real traceS = Sxx + Syy;
      const Real traceSsqrd = traceS*traceS;
      // Determinant of 2x2 matrix
      const Real determS = Sxx*Syy - Syx*Syx;
      // Largest eigenvalue
      const Real eigVal1 = 0.5*traceS + std::sqrt(0.25*traceSsqrd - determS);
      // Eigenvector associated with largest eigenvalue
      RealVect eVec = RealVect_zero;
      eVec[0] = 1.;
      if (std::abs(Syx) > 0.)
        {
          eVec[0] = eigVal1 - Syy;
          eVec[1] = Syx;
        }
      // Normalize final eigenvector
      const Real eigenSize = std::sqrt(D_TERM(eVec[0]*eVec[0],
                                            + eVec[1]*eVec[1],
                                            + eVec[2]*eVec[2]) + epsilon);
      D_TERM(a_eigenVect[MD_IX(i,0)] = eVec[0]/eigenSize;,
             a_eigenVect[MD_IX(i,1)] = eVec[1]/eigenSize;,
             a_eigenVect[MD_IX(i,2)] = eVec[2]/eigenSize;);
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Compute the modeled SGS-vortex orientation
/** \param[out] a_orientationFab
 *                      Face-centered modeled SGS-vortex 
 *                      orientation tensor
 *  \param[in]  a_NGradUfacePntFab
 *                      Face-centered physical-space velocity gradients
 *  \param[in]  a_eigenVect
 *                      Eigenvectors of strain-rate tensor
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Cell-centered box over which orientation
 *                      tensor is calculated
 *//*-----------------------------------------------------------------*/

void
SSV::sgsVortexVector(FArrayBox&              a_orientationFab,
                     const FArrayBox&        a_NGradUfacePntFab,
                     const FArrayBox&        a_eigenVect,
                     const LevelGridMetrics& a_gridMetrics,
                     const ProblemDomain&    a_domain,
                     const Box&              a_box) const
{
  CH_TIME("SSV::sgsVortexVector");
  // Model based on strain-rate tensor eigenvector with
  // largest eigenvalue and resolved-scale-vorticity vector
  const int cEnd = PatchMappedFunc::m_numCurlComps;
  FABSTACKTEMP(vort, a_box, cEnd);
  vort.setVal(0.);

  if ((m_vortexProportion < 1.) || (m_vortexModel == 1))
    {
      const int useMagVortComp = 0;
      FABSTACKTEMP(magVort, a_box, 1);
      // Compute physical space vorticity
      PatchMappedFunc::curlPS(a_box, vort, a_NGradUfacePntFab);
      // Compute magnitude of vorticity
      PatchMappedFunc::magnitude(a_box, magVort, 0, vort, Interval(0, cEnd-1));
      // Add small value to magVort in order to avoid division by zero
      const Real eps = 1.e-40;
      magVort.plus(eps, 0, 1);
      // Divide vorticity vector by magnitude of vorticity
      PatchMappedFunc::divideVec(
        a_box, vort, Interval(0, cEnd-1), magVort, useMagVortComp);
    }
  if (m_vortexModel == 0)
    {
      for (int j = 0; j != SpaceDim; ++j)
        {
          for (int k = 0; k != SpaceDim; ++k)
            {
              Real kronDel = 0.; // Kronecker delta_jk
              if (j == k)
                {
                  kronDel = 1.;
                }
#if CH_SPACEDIM == 3
              MD_BOXLOOP(a_box, i)
                {
                  a_orientationFab[MD_IX(i, srtIdx(j,k))] =
                    (m_vortexProportion*(
                     kronDel - a_eigenVect[MD_IX(i,j)]*a_eigenVect[MD_IX(i,k)]))
                  + ((1. - m_vortexProportion)*
                     (kronDel - vort[MD_IX(i,j)]*vort[MD_IX(i,k)]));
                }
#endif
#if CH_SPACEDIM == 2
              MD_BOXLOOP(a_box, i)
                {
                  a_orientationFab[MD_IX(i, srtIdx(j,k))] =
                    (kronDel - a_eigenVect[MD_IX(i,j)]*a_eigenVect[MD_IX(i,k)]);
                }
#endif
            }
        }
    }
  else if (m_vortexModel == 1)
    {
      for (int j = 0; j != SpaceDim; ++j)
        {
          for (int k = 0; k != SpaceDim; ++k)
            {
              Real kronDel = 0.; // Kronecker delta_jk
              if (j == k)
                {
                  kronDel = 1.;
                }
              MD_BOXLOOP(a_box, i)
                {
                  // Cross-product of eigenVect and vorticity
                  RealVect e(D_DECL(a_eigenVect[MD_IX(i,0)],
                                    a_eigenVect[MD_IX(i,1)],
                                    a_eigenVect[MD_IX(i,2)]));
                  RealVect w(
                    D_DECL(vort[MD_IX(i,0)],vort[MD_IX(i,1)],vort[MD_IX(i,2)]));
                  RealVect dir3 = RealVect_zero;
                  D_TERM(,,dir3[0] = e[1]*w[2] - e[2]*w[1];
                           dir3[1] = e[2]*w[0] - e[0]*w[2];
                           dir3[2] = e[0]*w[1] - e[1]*w[0];);
                  Real magDir3 = dir3.vectorLength();
                  dir3 /= magDir3;
                  a_orientationFab[MD_IX(i, srtIdx(j,k))] =
                    (kronDel - dir3[j]*dir3[k]);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute kappa_C
/** \param[out] a_kappaC
 *                      Face-centered kappa_c
 *  \param[in]  a_WfacePntFab
 *                      Face-centered point primitive variables
 *  \param[in]  a_eigenVect
 *                      Eigenvector associated with largest eigenvalue
 *                      of strain-rate tensor in a cell
 *  \param[in]  a_strainRateTensorFab
 *                      Face-centered strain-rate tensor
 *  \param[in]  a_deltaC
 *                      Cell cutoff length
 *  \param[in]  a_box   Face-centered box where kappa_c is computed
 *//*-----------------------------------------------------------------*/

void
SSV::computeKappaC(FArrayBox&       a_kappaC,
                   const FArrayBox& a_WfacePntFab,
                   const FArrayBox& a_eigenVect,
                   const FArrayBox& a_strainRateTensorFab,
                   const FArrayBox& a_deltaC,
                   const Box&       a_box) const
{
  CH_TIME("SSV::computeKappaC");
  const int cRho  = CRDparam::g_CRDPhysics->densityIndex();
  const Real eps = 1.e-15; // Tolerance for division and sqrt
  // Compute viscosity for cutoff term in kappaC
  FABSTACKTEMP(muFab, a_box, 1);
  muFab.setVal(CRDparam::g_mu);
  FABSTACKTEMP(kFab, a_box, 1);
  kFab.setVal(CRDparam::g_K);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      CRDparam::g_CRDPhysics->calcCoeffKappaMu(a_box,muFab,kFab,a_WfacePntFab);
    }

  MD_ARRAY_RESTRICT(arrEV, a_eigenVect);
  MD_ARRAY_RESTRICT(arrSRT, a_strainRateTensorFab);
#if CH_SPACEDIM == 3
  MD_BOXLOOP(a_box, i)
    {
      const Real k_c = PI/a_deltaC[MD_IX(i,0)]; // Wavenumber cutoff
      Real nu = (muFab[MD_IX(i, 0)])/(a_WfacePntFab[MD_IX(i, cRho)]);

      // Stretching aligned with SGS-vortex axis
      Real a = D_TERM(0., + 0., +
             arrEV[MD_IX(i,0)]*arrEV[MD_IX(i,0)]*arrSRT[MD_IX(i,0)]
        + 2.*arrEV[MD_IX(i,0)]*arrEV[MD_IX(i,1)]*arrSRT[MD_IX(i,1)]
        + 2.*arrEV[MD_IX(i,0)]*arrEV[MD_IX(i,2)]*arrSRT[MD_IX(i,2)]
        +    arrEV[MD_IX(i,1)]*arrEV[MD_IX(i,1)]*arrSRT[MD_IX(i,4)]
        + 2.*arrEV[MD_IX(i,1)]*arrEV[MD_IX(i,2)]*arrSRT[MD_IX(i,5)]
        +    arrEV[MD_IX(i,2)]*arrEV[MD_IX(i,2)]*arrSRT[MD_IX(i,8)];);
      a_kappaC[MD_IX(i,0)] = k_c*std::sqrt(2.*nu/(3.*std::abs(a) + eps));
    }
#endif
#if CH_SPACEDIM == 2
  MD_BOXLOOP(a_box, i)
    {
      const Real k_c = PI/a_deltaC[MD_IX(i,0)]; // Wavenumber cutoff
      Real nu = (muFab[MD_IX(i, 0)])/(a_WfacePntFab[MD_IX(i, cRho)]);

      // Stretching aligned with SGS-vortex axis
      Real a = D_TERM(
             arrEV[MD_IX(i,0)]*arrEV[MD_IX(i,0)]*arrSRT[MD_IX(i,0)],
        + 2.*arrEV[MD_IX(i,0)]*arrEV[MD_IX(i,1)]*arrSRT[MD_IX(i,1)]
        +    arrEV[MD_IX(i,1)]*arrEV[MD_IX(i,1)]*arrSRT[MD_IX(i,3)];,);
      a_kappaC[MD_IX(i,0)] = k_c*std::sqrt(2.*nu/(3.*std::abs(a) + eps));
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Compute the Kolmogorov prefactor for the SGS KE calculation
/** \param[out] a_KP
 *                      Face-centered Kolmogorov prefactor
 *  \param[in]  a_WfacePntFab
 *                      Face-centered point primitive variables
 *  \param[out] a_eigenVect
 *                      Row-ordered eigenvectors of strain-rate tensor
 *  \param[out] a_kappaC
 *                      Face-centered kappa_c
 *  \param[in]  a_XFab
 *                      Physical space face coordinates
 *  \param[out] a_strainRateTensor
 *                      Face-centered, physical-space 
 *                      strain-rate tensor
 *  \param[in]  a_deltaC
 *                      Cell cutoff length
 *  \param[in]  a_domain
 *                      Problem domain for block
 *  \param[in]  a_box   Face-centered box where kappa_c is computed
 *//*-----------------------------------------------------------------*/

void
SSV::prefactorKolmogorov(FArrayBox&           a_KP,
                         const FArrayBox&     a_WfacePntFab,
                         const FArrayBox&     a_eigenVect,
                         const FArrayBox&     a_kappaC,
                         const FArrayBox&     a_XFab,
                         const FArrayBox&     a_strainRateTensor,
                         const FArrayBox&     a_deltaC,
                         const ProblemDomain& a_domain,
                         const Box&           a_box,
                         const int            a_localCoarsening) const
{
  CH_TIME("SSV::prefactorKolmogorov");

  // Global constants
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rad = 1; // Radius of ensemble-averaging
  const Real epsilon = 1.e-14; // Prevent r from becoming small negative
  int f_0 = 1;
  if (a_localCoarsening) // If using local coarsening, expand stencil
    {
      f_0 = CRDparam::g_sgskeCrsLevFilterRatio;
    }

  a_KP.setVal(0.); // Set to zero in order to fill with structure function
  FABSTACKTEMP(Q, a_box, 1);
  Q.setVal(0.);

  MD_ARRAY_RESTRICT(arrW, a_WfacePntFab);

  if (!m_sgsKineticEnergyApprox) // SGS KE from Voelkl, Chung, and Pullin
    {
      for (int j = -rad; j != (rad + 1); ++j)
        {
          for (int k = -rad; k != (rad + 1); ++k)
            {
#if CH_SPACEDIM == 3
              for (int l = -rad; l != (rad + 1); ++l)
                {
#endif
                  IntVect os(D_DECL(f_0*j,f_0*k,f_0*l));
                  IntVect osReverse = -os;
                  // Create the correct box for ensemble averaging over
                  Box ensembleBox(a_box);
                  ensembleBox.shift(os);        // Shift in offset vector dir
                  ensembleBox &= a_domain;      // Trim to match domain
                  ensembleBox.shift(osReverse); // Shift back
                  // Add another ensembleBox
                  Box ensembleCoarseBox(a_box);
                  ensembleCoarseBox.shift((f_0*os));
                  ensembleCoarseBox &= a_domain;
                  ensembleCoarseBox.shift((f_0*osReverse));
                  MD_BOXLOOP(ensembleBox, i)
                    {
                      // Compute the 2nd-order structure function
                      // We set the numerator of a_KP to the structure function
                      a_KP[MD_IX(i,0)] += D_TERM(
                       (arrW[MD_IX(i,cVel)]-arrW[MD_OFFSETIV(i,+,os,cVel)])*
                       (arrW[MD_IX(i,cVel)]-arrW[MD_OFFSETIV(i,+,os,cVel)]),
                     + (arrW[MD_IX(i,cVel+1)]-arrW[MD_OFFSETIV(i,+,os,cVel+1)])*
                       (arrW[MD_IX(i,cVel+1)]-arrW[MD_OFFSETIV(i,+,os,cVel+1)]),
                     + (arrW[MD_IX(i,cVel+2)]-arrW[MD_OFFSETIV(i,+,os,cVel+2)])*
                       (arrW[MD_IX(i,cVel+2)] -
                        arrW[MD_OFFSETIV(i,+,os,cVel+2)]););
                    }
#if CH_SPACEDIM == 3
                }
#endif
            }
        }

      const Real c_1 = 0.873469;
      const Real c_2 = 7.4022;
      const Real c_3 = 1.82642;
      const Real c_4 = 12.2946;
      const Real c_5 = 2./3.;
      const Real c_6 = 0.573159;
      const Real c_7 = -3./2.;
      const Real c_8 = PI/4.;

      // 2nd-order polynomial approximation with 0.2% max error
      // Valid for Cartesian meshes, anistropic grid has d > sqrt(3)
      //**FIXME(?): This polynomial interpolant does not match smoothly
      //**FIXME: Find better approximation that covers anisotropic meshes
      //**FIXME: Find approximation that covers entire range of d values
      const Real p_2 = -6.597682919;
      const Real p_3 = 14.912189960;
      const Real p_4 = -2.416021978;

      for (int j = -rad; j != (rad + 1); ++j)
        {
          for (int k = -rad; k != (rad + 1); ++k)
            {
#if CH_SPACEDIM == 3
              for (int l = -rad; l != (rad + 1); ++l)
                {
#endif
                  IntVect os(D_DECL(f_0*j,f_0*k,f_0*l)); // Offset vector
                  IntVect osReverse = -os;
                  // Create the correct box for ensemble averaging over
                  Box ensembleBox(a_box);
                  ensembleBox.shift(os);        // Shift in dir of offset vector
                  ensembleBox &= a_domain;      // Trim to match domain
                  ensembleBox.shift(osReverse); // Shift back
                  // Add another ensembleBox
                  Box ensembleCoarseBox(a_box);
                  ensembleCoarseBox.shift((f_0*os));
                  ensembleCoarseBox &= a_domain;
                  ensembleCoarseBox.shift((f_0*osReverse));
                  MD_BOXLOOP(ensembleBox, i)
                    {
                      const RealVect eVec(D_DECL(a_eigenVect[MD_IX(i,0)],
                                                 a_eigenVect[MD_IX(i,1)],
                                                 a_eigenVect[MD_IX(i,2)]));
                      const RealVect osr(
                        D_DECL(
                          a_XFab[MD_IX(i,0)] - a_XFab[MD_OFFSETIV(i,+,os,0)],
                          a_XFab[MD_IX(i,1)] - a_XFab[MD_OFFSETIV(i,+,os,1)],
                          a_XFab[MD_IX(i,2)] - a_XFab[MD_OFFSETIV(i,+,os,2)]));
                      const Real osrVecLenSqrd = osr.radSquared();
                      const Real eProjVecLen = osr.dotProduct(eVec);
                      const Real eProjVecLenSqrd = eProjVecLen*eProjVecLen;
                      const Real r = std::sqrt(
                        osrVecLenSqrd - eProjVecLenSqrd
                      + epsilon*0.5*(osrVecLenSqrd + eProjVecLenSqrd));
                      const Real d = r/a_deltaC[MD_IX(i,0)];
                      if (CRDparam::g_cartesian)
                        {
                          if (d < c_1)
                            {
                              Q[MD_IX(i,0)] += c_2*d*d - c_3*d*d*d*d;
                            }
                          else
                            {
                              Q[MD_IX(i,0)] += p_2 + p_3*d + p_4*d*d;
                            }
                        }
                      else
                        {
                          if (d < c_1)
                            {
                              Q[MD_IX(i,0)] += c_2*d*d - c_3*d*d*d*d;
                            }
                          else
                            {
                              Q[MD_IX(i,0)] += c_4*std::pow(d, c_5) - 6.
                                - c_6*std::pow(d, c_7)*std::sin(PI*d - c_8);
                            }
                        }
                    }
#if CH_SPACEDIM == 3
                }
#endif
            }
        }
      MD_BOXLOOP(a_box, i)
        {
          a_KP[MD_IX(i,0)] = a_KP[MD_IX(i,0)]/Q[MD_IX(i,0)];
        }
    }
  else // SGS KE calculation method from Mattner
    {
      const Real f_1 = 4./3.;
      const Real f_2 = 2./3.;
      for (int j = 0; j != SpaceDim; ++j)
        {
          for (int k = 0; k != SpaceDim; ++k)
            {
              MD_BOXLOOP(a_box, i)
                {
                  // Fill a_KP with the appropriate structural model
                  a_KP[MD_IX(i, 0)] +=
                    a_strainRateTensor[MD_IX(i,srtIdx(j,k))]*
                    a_strainRateTensor[MD_IX(i,srtIdx(j,k))];
                }
            }
        }
      MD_BOXLOOP(a_box, i)
        {
          Real kappa_c = PI/a_deltaC[MD_IX(i,0)]; // Wavenumber cutoff
          Real lambda_v = a_kappaC[MD_IX(i,0)]/kappa_c;
          Real f_3 = std::pow(lambda_v, f_2);
          Real f_4 = std::pow(kappa_c, f_1);
          a_KP[MD_IX(i,0)] = f_1*f_3*a_KP[MD_IX(i,0)]/f_4; // Scale a_KP
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute cell-avg SGS kinetic energy estimate with minimal input
/** \param[out] a_U     Cell-averaged conservative state
 *  \param[out] a_JU    Mapped, cell-averaged conservative state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for this level
 *  \param[out] a_dataIndx
 *                      Data index for box within this level
 *  \param[in]  a_fourthOrder
 *                      If true, SGS KE is computed to 4th-order
 *  \param[in]  a_box   Centered box over which subgrid-scale
 *                      kinetic energy is calculated
 *//*-----------------------------------------------------------------*/

void
SSV::cellAvgSGSKineticEnergy(FArrayBox&              a_U,
                             FArrayBox&              a_JU,
                             const LevelGridMetrics& a_gridMetrics,
                             const DataIndex&        a_dataIndx,
                             const bool              a_fourthOrder,
                             const Box&              a_box) const
{
  CH_TIME("SSV::cellAvgSGSKineticEnergy");

  const Interval turbInterval = CRDparam::g_CRDPhysics->turbConsInterval();
  const int sgsComp = turbInterval.begin();
  const int tensorComps = SpaceDim*SpaceDim;
  FArrayBox sgsFab(turbInterval, a_U);
  FABSTACKTEMP(eigVect, a_box, SpaceDim);
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box); //**FIXME: disjointBox

  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;

  const RealVect dxVect = a_gridMetrics.dxVect();

  // Second-order SGS KE computation
  const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  FABSTACKTEMP(WFab, box1Dom, numWcomp);
  // 2nd-ord prim from cell-avg cons state (disjointBox + 2 gc)
  //   <U>  -->  W|_2
  PatchCNSOp::computeWpntCell(
    WFab, a_U, WFab, box1Dom, blockDomain, true, false);

  // 2nd-ord phys gradients of 2nd-ord prim state (disjointBox + 2)
  //   W|_2  -->  N^tJ*Grad(W|_2)
  FABSTACKTEMP(GradVel2O, a_box, tensorComps);
  PatchMappedFunc::gradient2OPS(
    a_box,
    GradVel2O,
    WFab,
    a_gridMetrics.m_cellNtJ[a_dataIndx],
    blockDomain,
    CRDparam::g_CRDPhysics->velocityInterval(),
    dxVect);

  // Compute the strain-rate tensor
  const Real maxBound = 1.e50;
  const Real minBound = -1.e50;
  FABSTACKTEMP(srtFab, a_box, tensorComps);
  for (int j = 0; j != SpaceDim; ++j)
    {
      for (int k = 0; k != SpaceDim; ++k)
        {
          MD_BOXLOOP(a_box, i)
            {
              Real srtVal = 0.5*(
                GradVel2O[MD_IX(i,srtIdx(j,k))]
              + GradVel2O[MD_IX(i,srtIdx(k,j))]);
              if ((srtVal > maxBound) || (srtVal < minBound) ||
                  (std::isnan(srtVal)) || (std::isinf(srtVal)))
                {
                  srtVal = 0.;
                }
              srtFab[MD_IX(i,srtIdx(j,k))] = srtVal;
            }
        }
    }

  // Compute delta_c
  FABSTACKTEMP(deltaC, a_box, 1); // Cell cutoff length
  const Real cellVol = a_gridMetrics.dxVect().product();
  const Real c_0  = 1./SpaceDim;
  const FArrayBox& JFab = a_gridMetrics.m_J[a_dataIndx];
  if (CRDparam::g_cartesian)
    {
      const Real deltaCvalCartesian = std::pow(cellVol, c_0);
      MD_BOXLOOP(a_box, i)
        {
          deltaC[MD_IX(i,0)] = deltaCvalCartesian;
        }
    }
  else
    {
      MD_BOXLOOP(a_box, i)
        {
          deltaC[MD_IX(i,0)] = std::pow(JFab[MD_IX(i, 0)]*cellVol, c_0);
        }
    }

  const FArrayBox& XFab = a_gridMetrics.m_cellCoord[a_dataIndx];

  // Compute subgrid-scale kinetic energy
  sgsKineticEnergy(sgsFab, eigVect, GradVel2O, srtFab, WFab, deltaC, XFab,
                   blockDomain, a_box);

  // Compute JU version of subgrid-scale kinetic energy estimate
  MD_BOXLOOP(a_box, i)
    {
      a_JU[MD_IX(i, sgsComp)] = JFab[MD_IX(i, 0)]*sgsFab[MD_IX(i, 0)];
    }
}

/*--------------------------------------------------------------------*/
//  Compute cell-avg SGS kinetic energy estimate with minimal input
/** \param[out] a_SGSKE Cell-averaged SGS KE estimate
 *  \param[in]  a_JW    Mapped, cell-averaged conservative state
 *  \param[in]  a_JGradW
 *                      Mapped, cell-averaged velocity gradient
 *  \param[in]  a_cellPntJ
 *                      Cell-centered physical cell volume
 *  \param[in]  a_crsXFab
 *                      Cell-centered physical space coordinates
 *  \param[in]  a_crsDeltaC
 *                      Cell-centered grid-scale cutoff-length
 *  \param[out] a_probDom
 *                      Problem domain
 *//*-----------------------------------------------------------------*/

void
SSV::crsLevCellAvgSGSKineticEnergy(
  LevelData<FArrayBox>&       a_SGSKE,
  const LevelData<FArrayBox>& a_JW,
  const LevelData<FArrayBox>& a_JGradW,
  const LevelData<FArrayBox>& a_cellPntJ,
  const LevelData<FArrayBox>& a_crsXFab,
  const LevelData<FArrayBox>& a_crsDeltaC,
  const ProblemDomain&        a_probDom) const
{
  CH_TIME("SSV::crsLevCellAvgSGSKineticEnergy");

  const DisjointBoxLayout crsBoxes = a_SGSKE.getBoxes();
  const int tensorComps = SpaceDim*SpaceDim;
  for (DataIterator dit = crsBoxes.dataIterator(); dit.ok(); ++dit)
    {
      Box disjointBox = crsBoxes[dit];
      Box box1Dom = grow(disjointBox, 1);
      box1Dom &= a_probDom;
      Box box2Dom = grow(disjointBox, 2);
      box2Dom &= a_probDom;
      const FArrayBox& JW = a_JW[dit];
      const FArrayBox& JFab = a_cellPntJ[dit];
      FArrayBox& sgsFab = a_SGSKE[dit];
      const FArrayBox& JGradW = a_JGradW[dit];
      const FArrayBox& XFab = a_crsXFab[dit];
      const FArrayBox& deltaC = a_crsDeltaC[dit];

      // 1) Deconvolve the gradients
      FABSTACKTEMP(cellPntVelGrad, box1Dom, tensorComps);
      CRDutil::deconvolve(cellPntVelGrad, JGradW, box1Dom, a_probDom,
                          cellPntVelGrad.interval(), 4, 1, false);

      // 2) Bring the gradients back to physical space (divide by J)
      for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
        {
          MD_BOXLOOP(box1Dom, i)
            {
              cellPntVelGrad[MD_IX(i, comp)] /= JFab[MD_IX(i, 0)];
            }
        }

      // 3) Deconvolve the density and velocity
      int cellPntWComps = 1 + SpaceDim;
      FABSTACKTEMP(cellPntW, box2Dom, cellPntWComps);
      CRDutil::deconvolve(cellPntW, JW, box2Dom, a_probDom,
                          Interval(0,cellPntWComps-1), 4, 1, false);

      // 4) Bring the density and velocity back to physical space (divide by J)
      for (int comp = 0; comp != (SpaceDim+1); ++comp)
        {
          MD_BOXLOOP(box2Dom, i)
            {
              cellPntW[MD_IX(i, comp)] /= JFab[MD_IX(i, 0)];
            }
        }

      // 5) Compute the strain-rate tensor
      FABSTACKTEMP(srtFab, box1Dom, tensorComps);
      for (int j = 0; j != SpaceDim; ++j)
        {
          for (int k = 0; k != SpaceDim; ++k)
            {
              MD_BOXLOOP(box1Dom, i)
                {
                  srtFab[MD_IX(i,srtIdx(j,k))] = 0.5*(
                    cellPntVelGrad[MD_IX(i,srtIdx(j,k))]
                  + cellPntVelGrad[MD_IX(i,srtIdx(k,j))]);
                }
            }
        }

      // 6) compute subgrid-scale kinetic energy
      FABSTACKTEMP(eigVect, box1Dom, SpaceDim);
      sgsKineticEnergy(sgsFab, eigVect, cellPntVelGrad, srtFab, cellPntW,
                       deltaC, XFab, a_probDom, box1Dom);

      // 7)  Compute cell-averaged J*SGSKE
      FABSTACKTEMP(cellAvgSGSFab, disjointBox, 1);
      // 7a) Multiply sgsFab by pntJ
      sgsFab.mult(JFab);
      // 7b) Convolve sgsFab to get <J*SGSKE>
      CRDutil::deconvolve(cellAvgSGSFab, sgsFab, disjointBox, a_probDom,
                          cellAvgSGSFab.interval(), 4, -1, false);
      // 7c) Copy <J*SGSKE> into the correct Fab
      sgsFab.copy(cellAvgSGSFab);
    }
}

/*--------------------------------------------------------------------*/
//  Compute SGS kinetic energy estimate
/** \param[out] a_sgsKE Subgrid scale kinetic energy estimate
 *  \param[out] a_eigVect
 *                      Row-ordered eigenvectors of strain-rate tensor
 *  \param[in]  a_NGradUfacePntFab
 *                      Face-centered physical-space velocity gradients
 *  \param[out] a_srtFab
 *                      Face-centered, physical-space 
 *                      strain-rate tensor
 *  \param[in]  a_WfacePntFab
 *                      Face-centered point primitive variables
 *  \param[in]  a_deltaC
 *                      Cell cutoff length
 *  \param[in]  a_XFab
 *                      Physical space face coordinates
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Face-centered box over which subgrid-scale
 *                      kinetic energy is calculated
 *//*-----------------------------------------------------------------*/

void
SSV::sgsKineticEnergy(FArrayBox&           a_sgsKE,
                      FArrayBox&           a_eigVect,
                      const FArrayBox&     a_NGradUfacePntFab,
                      const FArrayBox&     a_srtFab,
                      const FArrayBox&     a_WfacePntFab,
                      const FArrayBox&     a_deltaC,
                      const FArrayBox&     a_XFab,
                      const ProblemDomain& a_domain,
                      const Box&           a_box,
                      const int            a_localCoarsening) const
{
  CH_TIME("SSV::sgsKineticEnergy");

  FABSTACKTEMP(kappaC, a_box, 1);
  // NOTE: We now compute the model on all faces, even boundaries
  eigenDecomp(a_eigVect, a_srtFab, a_box);
  computeKappaC(kappaC, a_WfacePntFab, a_eigVect, a_srtFab, a_deltaC, a_box);
  // Use a_sgsKE as temporary space for computing Kolmogorov prefactor
  prefactorKolmogorov(a_sgsKE, a_WfacePntFab, a_eigVect, kappaC, a_XFab,
                      a_srtFab, a_deltaC, a_domain, a_box, a_localCoarsening);
  const Real f_1 = 1./3.;
  MD_BOXLOOP(a_box, i)
    {
      // Compute the incomplete-gamma-function factor
      const Real k_c2 = kappaC[MD_IX(i,0)]*kappaC[MD_IX(i,0)];
      const Real k_c4 = k_c2*k_c2;
      const Real k_c6 = k_c2*k_c4;
      Real gammaFactor = 0.;
      if (k_c2 < 2.42806)
        {
          gammaFactor = ((3. + 2.5107*k_c2 + 0.330357*k_c4 + 0.0295481*k_c6)/
             (1. + 0.336901*k_c2 + 0.0416684*k_c4 + 0.00187191*k_c6))
            - 4.06235*std::pow(k_c2, f_1);
        }
      else
        {
          gammaFactor = ((13.111151 + 8.666661*k_c2 + k_c4)/
             (10.370367 + 23.333326*k_c2 + 10.*k_c4 + k_c6))*std::exp(-k_c2);
        }
      // SGS kinetic energy
      a_sgsKE[MD_IX(i,0)] = 0.5*a_sgsKE[MD_IX(i,0)]*gammaFactor;
    }
}

/*--------------------------------------------------------------------*/
//  Compute SGS kinetic energy correction estimate for consToPrim
/** \param[out] a_WcellFab
 *                      Cell primitive variables with pressure
 *                      corrected by SGS KE estimate
 *  \param[in]  a_box   Cell-centered box over which pressure is
 *                      corrected after consToPrim calculation
 *//*-----------------------------------------------------------------*/

void
SSV::consToPrimCorrectionLES(FArrayBox& a_WcellFab,
                             const Box& a_box) const
{
  CH_TIME("SSV::consToPrimCorrectionLES");
  if (!CRDparam::g_useConsToPrimCorrection)
    {
      return;
    }

  // 1) FIXME: Does not work with thermally-perfect physics. To fix,
  //           (a) Update the calculations of gamma and R
  //           (b) Add multi-species SGS corrections (if necessary)

  // Global constants
  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const int cPress = CRDparam::g_CRDPhysics->pressureIndex();
  const int cSGSKE = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const Real gamma = CRDparam::g_gamma;    // Ratio of specific heats

  CH_assert(a_WcellFab.box().contains(a_box));

  MD_BOXLOOP(a_box, i)
    {
      const Real rho = a_WcellFab[MD_IX(i,cRho)];
      const Real sgsKE = a_WcellFab[MD_IX(i, cSGSKE)];
      a_WcellFab[MD_IX(i, cPress)] -= (gamma - 1.)*rho*sgsKE;
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
SSV::applyWallModel(const Box&              a_boundaryFaceBox,
                    const Box&              a_boundaryFaceGhostBox,
                    const Box&              a_disjointBox,
                    FArrayBox&              a_WfaceBdryFab,
                    FArrayBox&              a_WfaceAvgDirFab,
                    FArrayBox&              a_WcellAvgFab,
                    const FArrayBox&        a_unitNormalBasisFab,
                    const BCInfo&           a_bcTypes,
                    const Box&              a_bcBoxes,
                    const LevelGridMetrics& a_gridMetrics,
                    const DataIndex&        a_dataIndx,
                    const int               a_dir,
                    const Side::LoHiSide&   a_side,
                    const Real              a_time,
                    const int               a_level) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::applyWallModel");
}

/*--------------------------------------------------------------------*/
//  Applies the wall model to the faces and ghost cells
/** \param[out] a_bndryVelFab
 *                      Slip velocities on the boundary
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Domain of the current block
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_bndryNtJFab
 *                      Face-centered NtJ
 *  \param[in]  a_cellFaceAvgBdryBox
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_cellBoxWcellBdry
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_disjointBox
 *                      Disjoint box for grid metric functions
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *//*-----------------------------------------------------------------*/

void
SSV::turbModelSlipVel(FArrayBox&              a_bndryVelFab,
                      FArrayBox&              a_WcellAvgFab,
                      const FluxBox&          a_WfaceAvgFxb,
                      const LevelGridMetrics& a_gridMetrics,
                      const ProblemDomain&    a_domain,
                      const DataIndex&        a_dataIndx,
                      const FArrayBox&        a_unitNormalsDirFab,
                      const FArrayBox&        a_bndryNtJFab,
                      const Box&              a_boundaryFaceBox,
                      const Box&              a_boundaryFaceGhostBox,
                      const Box&              a_disjointBox,
                      const int               a_dir,
                      const Side::LoHiSide&   a_side) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::turbModelSlipVel");

  computeSlipVel(a_bndryVelFab,
                 a_WcellAvgFab,
                 a_WfaceAvgFxb,
                 a_gridMetrics,
                 a_domain,
                 a_dataIndx,
                 a_unitNormalsDirFab,
                 a_bndryNtJFab,
                 a_boundaryFaceBox,
                 a_boundaryFaceGhostBox,
                 a_disjointBox,
                 a_dir,
                 a_side);
}

/*--------------------------------------------------------------------*/
//  Applies the wall model to the faces and ghost cells
/** \param[out] a_bndryVelFab
 *                      Slip velocities on the boundary
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Domain of the current block
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_bndryNtJFab
 *                      Face-centered NtJ Fab
 *  \param[in]  a_cellFaceAvgBdryBox
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_cellBoxWcellBdry
 *                      Box of boundary faces requiring cell-averages
 *  \param[in]  a_disjointBox
 *                      Disjoint box for grid metric functions
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of box on which to work (low or high)
 *//*-----------------------------------------------------------------*/

void
SSV::computeSlipVel(FArrayBox&              a_bndryVelFab,
                    FArrayBox&              a_WcellAvgFab,
                    const FluxBox&          a_WfaceAvgFxb,
                    const LevelGridMetrics& a_gridMetrics,
                    const ProblemDomain&    a_domain,
                    const DataIndex&        a_dataIndx,
                    const FArrayBox&        a_unitNormalsDirFab,
                    const FArrayBox&        a_bndryNtJFab,
                    const Box&              a_boundaryFaceBox,
                    const Box&              a_boundaryFaceGhostBox,
                    const Box&              a_disjointBox,
                    const int               a_dir,
                    const Side::LoHiSide&   a_side) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::computeSlipVel");

  // We need to figure out if this is a wall or not
  BoundaryIndex bcIdx;
  bcIdx.define(
    a_gridMetrics.getCoordSys().whichBlock(a_disjointBox), a_dir, a_side);
  auto domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
  if (!(CRDparam::DomainBCTypeSWall & domBC.m_type)) { return; }

  // Choose a global, physical-space streamwise vector
  const RealVect streamwisePhysVect = CRDparam::g_globStreamwiseDir;

  // Global indices
  const int cRho        = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel        = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cTurbStart  = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int cSGSKE      = cTurbStart;
  const int cEta        = cTurbStart + a_dir + 1;
  const int tensorComps = SpaceDim*SpaceDim;
  // Constants
  const Real mu = CRDparam::g_mu;
  const Real h_nu_plus = 11.; // Empirical constant
  const Real gamma_11 = 0.45; // Empirical constant
  const Real eps = 1.e-40; // Small value (invalid divide by zero, etc.)

  // Opposite side of current side (use for opposite shift direction)
  Side::LoHiSide oppSide = flip(a_side);
  // Sign for shifting into the domain: +1 if lo, -1 if hi
  const int offsetSign = sign(oppSide);
  // Sign for shifting out of the domain: -1 if lo, +1 if hi
  const int normSign = sign(a_side);

  // Create box for first layer of interior cells
  Box firstInteriorCellBox = a_boundaryFaceGhostBox;
  // Grow the face box by one on the interior side of the box
  firstInteriorCellBox.growDir(a_dir, oppSide, 1);
  // Turn the face box into a cell box
  firstInteriorCellBox.enclosedCells(a_dir);

  // Create box for first interior layer of wall-tangential faces
  Box firstInteriorFaceBox = a_boundaryFaceGhostBox;
  // Shift the wall-faces box into the domain by 1
  firstInteriorFaceBox.shift(a_dir, offsetSign);

  // Create box for wall-faces and first layer of interior wall-tangential faces
  Box wallDistanceBox = a_boundaryFaceGhostBox;
  // Grow the wall-faces box to include the first interior cell
  wallDistanceBox.growDir(a_dir, oppSide, 1);

  // We need the coordinate system for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));

  // Compute the physical distance vector btwn the wall and first interior faces
  FABSTACKTEMP(XiFab, wallDistanceBox, SpaceDim); // computational space coords
  FABSTACKTEMP(XFab, wallDistanceBox, SpaceDim);  // physical space coordinates
  // Get the coordinates on the wall and first interior faces
  CRDparam::g_CNSIBC->getFaceCoordinates(
    wallDistanceBox, XiFab, XFab, a_dir, blockCoordSys);
  // Compute the difference between XFab on the top and bottom faces of the cell
  FABSTACKTEMP(wallDistanceVect, firstInteriorFaceBox, SpaceDim);
  const int MD_ID(ii, a_dir);
  MD_BOXLOOP(firstInteriorFaceBox, i)
    {
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          Real xCompInterior = XFab[MD_IX(i, comp)];
          Real xCompWall = XFab[MD_OFFSETIX(i, +, normSign*ii, comp)];
          // Compute (interior - wall)*(side compensation) (1 if lo, -1 if hi)
          wallDistanceVect[MD_IX(i, comp)] =
            offsetSign*(xCompInterior - xCompWall);
        }
    }
  // Shift the wallDistanceVect fab to the wall
  wallDistanceVect.shift(a_dir, normSign);
  // Transform the wall-distance vector into normal-tangent space
  FORT_FORWARDTRANSFORMF(CHF_FRA(wallDistanceVect),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Shift the wallDistanceVect fab to the first interior cell, not face
  wallDistanceVect.shiftHalf(a_dir, offsetSign);

  // Compute u_tau and h_0_plus
  const FArrayBox& WfaceAvgFab = a_WfaceAvgFxb[a_dir];
  FABSTACKTEMP(wallRhoFab, wallDistanceBox, 1);
  // Interpolate the cell-averaged density and SGS KE to the faces
  const int orderOfAccuracy = 4;
  PatchMappedFunc::faceAvgValFromCellAvgCS(
    wallRhoFab, a_WcellAvgFab, cRho, 0, a_domain, wallDistanceBox,
    a_dir, orderOfAccuracy, false);
  // Temporarily copy the density from the wall-face and interior face
  wallRhoFab.copy(WfaceAvgFab, cRho, 0, 1);
  wallRhoFab.shiftHalf(a_dir, offsetSign); // Shift wall face to the cell
  FABSTACKTEMP(uTau, firstInteriorCellBox, 1);
  FABSTACKTEMP(h0Plus, firstInteriorCellBox, 1);
  FABSTACKTEMP(h0, firstInteriorCellBox, 1);
  MD_BOXLOOP(firstInteriorCellBox, i)
    {
      Real wallRho = wallRhoFab[MD_IX(i, 0)];
      Real nu = mu/wallRho;
      Real invNu = wallRho/mu;
      // Note that eta_0 must be >= 0
      const Real eta_0 = std::max(0., a_WcellAvgFab[MD_IX(i, cEta)]);
      // Compute h_0 -- use the normal component of wallDistanceVect
      const Real dNormal = std::abs(wallDistanceVect[MD_IX(i, a_dir)]);
      // We set the virtual wall to be 18% of the cell normal-direction height
      const Real h_0 = m_virtualWallHeight*dNormal;
      h0[MD_IX(i, 0)] = h_0;
      uTau[MD_IX(i, 0)] = std::sqrt(nu*eta_0); // Compute u_tau
      h0Plus[MD_IX(i, 0)] = h_0*uTau[MD_IX(i, 0)]*invNu; // Compute h_0_plus
    }

  wallRhoFab.shiftHalf(a_dir, normSign); // Shift wall face back to wall
  h0.shiftHalf(a_dir, normSign); // Shift h0 to the wall
  // Compute the spatial gradient of spatially filtered eta_0
  // (a) In the wall plane, spatially average eta_0
  IntVect growTanVect = IntVect_unit;
  growTanVect[a_dir] = 0;
  Box cellBoxTan1 = firstInteriorCellBox;
  cellBoxTan1.grow(growTanVect);
  cellBoxTan1 &= a_domain;
  Box cellBoxTan2 = cellBoxTan1;
  cellBoxTan2.grow(growTanVect);
  cellBoxTan2 &= a_domain;
  FABSTACKTEMP(eta, cellBoxTan2, 1);
  MD_BOXLOOP(cellBoxTan2, i)
    {
      const Real eta_0 = std::max(0., a_WcellAvgFab[MD_IX(i, cEta)]);
      eta[MD_IX(i, 0)] = eta_0;
    }
  FABSTACKTEMP(filteredEta, cellBoxTan1, 1);
  MD_BOXLOOP(cellBoxTan1, i)
    {
      // Loop over a neighbor box and sum the values
      Box neighborBox(MD_GETIV(i), MD_GETIV(i));
      neighborBox.grow(growTanVect);
      neighborBox &= a_domain;
      Real sumEta = 0.;
      int numPnts = 0;
      MD_BOXLOOP(neighborBox, j)
        {
          sumEta += eta[MD_IX(j, 0)];
          ++numPnts;
        }
      filteredEta[MD_IX(i, 0)] = sumEta/numPnts;
    }
  // (b) In the wall plane, compute the derivatives of filteredEta
  FABSTACKTEMP(etaGrad, firstInteriorCellBox, SpaceDim);
  etaGrad.setVal(0.);
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      if (comp != a_dir)
        {
          PatchMappedFunc::cellAvgDerivFromCellAvgCS(
            etaGrad, filteredEta, 0, comp, a_domain, firstInteriorCellBox,
            a_gridMetrics.dxVect(), comp, false, false);
        }
    }
  // (c) Map the derivatives to physical space
  // Shift etaGrad to the first interior face
  etaGrad.shiftHalf(a_dir, offsetSign);
  FABSTACKTEMP(physGradEta, firstInteriorFaceBox, SpaceDim);
  PatchMappedFunc::gradientCStoPS(
    firstInteriorFaceBox, physGradEta, etaGrad, a_bndryNtJFab);
  // (d) Get the wall-tangential magnitude
  // Shift physGradEta to the wall face
  physGradEta.shift(a_dir, normSign);
  FORT_FORWARDTRANSFORMF(CHF_FRA(physGradEta),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Set the wall-normal component to zero
  physGradEta.setVal(0.0, a_dir);
  // Bring physGradEta back to physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(physGradEta),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Interpolate the cell-averaged velocity
  FABSTACKTEMP(intFaceVel, firstInteriorFaceBox, SpaceDim);
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      PatchMappedFunc::faceAvgValFromCellAvgCS(
        intFaceVel, a_WcellAvgFab, cVel+comp, comp, a_domain,
        firstInteriorFaceBox, a_dir, orderOfAccuracy, false);
    }
  // Now fill most of intFaceVel with a_WfaceAvgFxb
  intFaceVel.copy(a_WfaceAvgFxb[a_dir], cVel, 0, SpaceDim);
  // Shift intFaceVel to the wall
  intFaceVel.shift(a_dir, normSign);
  // Transform intFaceVel into wall-normal coordinates
  FORT_FORWARDTRANSFORMF(CHF_FRA(intFaceVel),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Set the wall-normal component to zero
  intFaceVel.setVal(0.0, a_dir);
  // Compute the magnitude of intFaceVel
  FABSTACKTEMP(faceVelMag, a_boundaryFaceGhostBox, 1);
  PatchMappedFunc::magnitudeSpace(
    a_boundaryFaceGhostBox, faceVelMag, 0, intFaceVel, 0);
  // Shift eta to the wall
  eta.shiftHalf(a_dir, normSign);

  // Compute tangential velocity magnitude in first interior cell
  FABSTACKTEMP(interiorCellVel, firstInteriorCellBox, SpaceDim);
  interiorCellVel.setVal(0.);
  // Now copy a_WfaceAvgFxb into WfaceVel to cover as much as possible
  interiorCellVel.copy(a_WcellAvgFab, cVel, 0, SpaceDim);
  // Shift interiorCellVel to the wall for transformation
  interiorCellVel.shiftHalf(a_dir, normSign);
  // Transform the velocity vector
  FORT_FORWARDTRANSFORMF(CHF_FRA(interiorCellVel),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Set the normal component to zero
  interiorCellVel.setVal(0.0, a_dir);
  // Reverse transform
  FORT_REVERSETRANSFORMF(CHF_FRA(interiorCellVel),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));

  // Compute strain-rate tensor
  // a) compute the face-normal gradient of velocity
  // b) compute the face-tangential gradient of velocity
  // c) interpolate the face-tangential gradient of velocity to the face
  Box cellAvgGradVelBox = firstInteriorCellBox;
  // Grow the cellAvgGradVelBox three cells into the domain
  cellAvgGradVelBox.growDir(a_dir, oppSide, 3);
  FABSTACKTEMP(cellAvgGradVel, cellAvgGradVelBox, tensorComps);
  FABSTACKTEMP(faceAvgGradVel, firstInteriorFaceBox, tensorComps);
  for (int i = 0; i != SpaceDim; ++i)
    {
      const int velComp = cVel + i;
      const int gradComp = i*SpaceDim + a_dir;
      PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(
        faceAvgGradVel, a_WcellAvgFab, velComp, gradComp,
        a_domain, firstInteriorFaceBox, a_gridMetrics.dxVect(), a_dir,
        true,   // 4th-order if true, 2nd-order if false
        false); // only interior stencil if true
    }
  // Compute the cell-averaged gradient of velocity
  for (int row = 0; row != SpaceDim; ++row)
    {
      const int velComp = cVel + row;
      for (int col = 0; col != SpaceDim; ++col)
        {
          if (col != a_dir)
            {
              const int gradComp = row*SpaceDim + col;
              PatchMappedFunc::cellAvgDerivFromCellAvgCS(
                cellAvgGradVel, a_WcellAvgFab, velComp, gradComp, a_domain,
                cellAvgGradVelBox, a_gridMetrics.dxVect(), col, true, false);
            }
        }
    }
  // Interpolate the cell-averaged gradients to the faces
  for (int i = 0; i != SpaceDim; ++i)
    {
      for (int j = 0; j != SpaceDim; ++j)
        {
          int comp = i*SpaceDim + j;
          if (j != a_dir)
            {
              PatchMappedFunc::faceAvgValFromCellAvgCS(
                faceAvgGradVel, cellAvgGradVel, comp, comp, a_domain,
                firstInteriorFaceBox, a_dir, orderOfAccuracy, false);
            }
        }
    }

  // Transform the computational space gradient into physical space
  FABSTACKTEMP(faceAvgPhysGradVel, firstInteriorFaceBox, tensorComps);
  PatchMappedFunc::gradientCStoPS(firstInteriorFaceBox,
                                  faceAvgPhysGradVel,
                                  faceAvgGradVel,
                                  a_bndryNtJFab);

  // Compute the strain-rate tensor
  FABSTACKTEMP(srtFab, firstInteriorFaceBox, tensorComps);
  for (int j = 0; j != SpaceDim; ++j)
    {
      for (int k = 0; k != SpaceDim; ++k)
        {
          MD_BOXLOOP(firstInteriorFaceBox, i)
            {
              srtFab[MD_IX(i,srtIdx(j,k))] = 0.5*(
                faceAvgPhysGradVel[MD_IX(i,srtIdx(j,k))]
              + faceAvgPhysGradVel[MD_IX(i,srtIdx(k,j))]);
            }
        }
    }

  // Compute the SGS momentum flux on the face
  // Memory allocation
  FABSTACKTEMP(eigenVect, firstInteriorFaceBox, SpaceDim);
  FABSTACKTEMP(orientationFab, firstInteriorFaceBox, tensorComps);
  FABSTACKTEMP(sgsMomentumFlux, firstInteriorFaceBox, tensorComps);

  // Compute eigenvectors/eigenvalues here
  eigenDecomp(eigenVect, srtFab, firstInteriorFaceBox);
  // Compute SGS-vortex orientation vector
  sgsVortexVector(orientationFab,
                  faceAvgGradVel,
                  eigenVect,
                  a_gridMetrics,
                  a_domain,
                  firstInteriorFaceBox);

  // Interpolate the cell-averaged density and SGS KE to the faces
  FABSTACKTEMP(faceDensity, firstInteriorFaceBox, 1);
  FABSTACKTEMP(faceSGSKE, firstInteriorFaceBox, 1);
  PatchMappedFunc::faceAvgValFromCellAvgCS(
    faceDensity, a_WcellAvgFab, cRho, 0, a_domain, firstInteriorFaceBox,
    a_dir, orderOfAccuracy, false);
  PatchMappedFunc::faceAvgValFromCellAvgCS(
    faceSGSKE, a_WcellAvgFab, cSGSKE, 0, a_domain, firstInteriorFaceBox,
    a_dir, orderOfAccuracy, false);

  // Now copy a_WfaceAvgFxb to cover most of the faces -- in the end, only
  // one face should be not covered by original a_WfaceAvgFxb
  faceDensity.copy(WfaceAvgFab, cRho, 0, 1);
  faceSGSKE.copy(WfaceAvgFab, cSGSKE, 0, 1);
  // Set the faceSGSKE fab to always be zero or more
  MD_BOXLOOP(firstInteriorFaceBox, i)
    {
      Real sgsKE = faceSGSKE[MD_IX(i, 0)];
      faceSGSKE[MD_IX(i, 0)] = std::max(0., sgsKE);
    }
  // Now finish computing the SGS momentum flux
  for (int comp = 0; comp != (SpaceDim*SpaceDim); ++comp)
    {
      MD_BOXLOOP(firstInteriorFaceBox, i)
        {
          const Real sgsKE = faceSGSKE[MD_IX(i, 0)];
          sgsMomentumFlux[MD_IX(i, comp)] =
            faceDensity[MD_IX(i,0)]*sgsKE*orientationFab[MD_IX(i,comp)];
        }
    }

  // Compute Kolmogorov-like parameter
  FABSTACKTEMP(K1, firstInteriorFaceBox, 1);
  // Get the face-normal vector, $\hat{e}_n$, on the wall face
  //**NOTE: the vector begins in normal-tangent space
  FABSTACKTEMP(normalVect, a_boundaryFaceGhostBox, SpaceDim);
  normalVect.setVal(0.0);
  normalVect.setVal(1.0, a_dir);
  // Transform normalVect into physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(normalVect),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  normalVect.shift(a_dir, offsetSign); // Shift normalVect to interior face
  // Compute the local wall-tangent components of the global streamwise vector
  FABSTACKTEMP(globStreamLocTanVect, a_boundaryFaceGhostBox, SpaceDim);
  // Set the components to match what is given by the user-defined vector
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      globStreamLocTanVect.setVal(streamwisePhysVect[dir], dir);
    }
  // Transform globStreamLocTanVect into the local normal-tangent space
  FORT_FORWARDTRANSFORMF(CHF_FRA(globStreamLocTanVect),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  // Zero the normal component of globStreamLocTanVect
  globStreamLocTanVect.setVal(0., a_dir);
  // Transform globStreamLocTanVect back to Cartesian physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(globStreamLocTanVect),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  globStreamLocTanVect.shift(a_dir, offsetSign); // Shift to interior face
  // Normalize the wall-tangent vector -- this is the streamwise unit vector
  PatchMappedFunc::normalize(firstInteriorFaceBox,
                             globStreamLocTanVect,
                             Interval(0, SpaceDim-1));
  // Compute the local wall-tangent components of the local velocity vector
  FABSTACKTEMP(locStreamLocTanVect, a_boundaryFaceGhostBox, SpaceDim);
  // Get the velocity components from the first interior face
  FABSTACKTEMP(WfaceVel, firstInteriorFaceBox, SpaceDim);
  // Interpolate the cell-averaged velocity to the faces
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      PatchMappedFunc::faceAvgValFromCellAvgCS(
        WfaceVel, a_WcellAvgFab, cVel + dir, dir, a_domain,
        firstInteriorFaceBox, a_dir, orderOfAccuracy, false);
    }
  // Now copy a_WfaceAvgFxb into WfaceVel to cover as much as possible
  WfaceVel.copy(WfaceAvgFab, cVel, 0, SpaceDim);
  // Transform the velocity vector
  const int transformIndx = 0;
  WfaceVel.shift(a_dir, normSign); // Shift to the wall
  FORT_FORWARDTRANSFORMGENF(CHF_FRA(locStreamLocTanVect),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(WfaceVel),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(a_unitNormalsDirFab),
                            CHF_BOX(a_boundaryFaceGhostBox));
  // Set the normal component to zero to obtain the wall-tangent velocity vector
  locStreamLocTanVect.setVal(0.0, a_dir);
  // Transform locStreamLocTanVect back into Cartesian physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(locStreamLocTanVect),
                         CHF_CONST_FRA(a_unitNormalsDirFab),
                         CHF_BOX(a_boundaryFaceGhostBox));
  locStreamLocTanVect.shift(a_dir, offsetSign); // Shift to interior face
  WfaceVel.shift(a_dir, offsetSign); // Shift to interior face
  // Get the local velocity magnitude for comparison with q before proceeding
  FABSTACKTEMP(velMagFab, firstInteriorFaceBox, 1);
  MD_BOXLOOP(firstInteriorFaceBox, i)
    {
      RealVect velMagVect(D_DECL(locStreamLocTanVect[MD_IX(i, 0)],
                                 locStreamLocTanVect[MD_IX(i, 1)],
                                 locStreamLocTanVect[MD_IX(i, 2)]));
      velMagFab[MD_IX(i, 0)] = velMagVect.vectorLength();
    }
  // Normalize the wall-tangent vector -- this is the streamwise unit vector
  PatchMappedFunc::normalize(firstInteriorFaceBox,
                             locStreamLocTanVect,
                             Interval(0, SpaceDim-1));
  MD_BOXLOOP(firstInteriorFaceBox, i)
    {
      // Compute the wall-normal, streamwise component of tau_sgs
      Real tau_ns = 0.;
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          tau_ns += globStreamLocTanVect[MD_IX(i, dir)]*(
            D_TERM(
              sgsMomentumFlux[MD_IX(i,srtIdx(dir,0))]*normalVect[MD_IX(i,0)],
            + sgsMomentumFlux[MD_IX(i,srtIdx(dir,1))]*normalVect[MD_IX(i,1)],
            + sgsMomentumFlux[MD_IX(i,srtIdx(dir,2))]*normalVect[MD_IX(i,2)]));
        }
      const Real sgsKE = faceSGSKE[MD_IX(i, 0)];
      // Make sure that tau_planar_z is always >= 0
      const Real tau_planar_z = std::abs(tau_ns) + eps;
      K1[MD_IX(i, 0)] = gamma_11*std::sqrt(sgsKE)/(2.*std::sqrt(tau_planar_z));
    }

  if (m_constantK1) { K1.setVal(0.41); }

  // Compute the slip-velocity for the wall
  K1.shift(a_dir, normSign);
  velMagFab.shift(a_dir, normSign);
  globStreamLocTanVect.shift(a_dir, normSign);
  locStreamLocTanVect.shift(a_dir, normSign);
  uTau.shiftHalf(a_dir, normSign);
  h0Plus.shiftHalf(a_dir, normSign);
  normalVect.shift(a_dir, normSign);
  MD_BOXLOOP(a_boundaryFaceGhostBox, i)
    {
      // Tangential components of the global streamwise vector
      const RealVect globalStreamDir(D_DECL(globStreamLocTanVect[MD_IX(i, 0)],
                                            globStreamLocTanVect[MD_IX(i, 1)],
                                            globStreamLocTanVect[MD_IX(i, 2)]));
      // Tangential components of the local streamwise vector
      const RealVect localStreamDir(D_DECL(locStreamLocTanVect[MD_IX(i, 0)],
                                           locStreamLocTanVect[MD_IX(i, 1)],
                                           locStreamLocTanVect[MD_IX(i, 2)]));
      // Local wall-normal unit vector
      const RealVect normVect(D_DECL(normalVect[MD_IX(i, 0)],
                                     normalVect[MD_IX(i, 1)],
                                     normalVect[MD_IX(i, 2)]));
      // Determine if local backflow exists
      const Real velDotProd = globalStreamDir.dotProduct(localStreamDir);
      bool separatedFlow = false;
      if (velDotProd <= 0.)
        {
          separatedFlow = true;
        }
      // If streamwisePhysVect is over 80 degrees, assume separation
      //**NOTE: The following corresponds roughly to theta = 81.89 degrees
      //**NOTE: This is probably very generous, separation should occur earlier
      //**NOTE: Move the following to documentation after initial development
      // Angle btwn vectors a,b: $\theta = \text{acos}\left(\frac{
      // \vec{a}\dot\vec{b}}{\lVert\vec{a}\rVert\lVert\vec{b}\rVert}\right)$
      // Since a,b are unit vectors, we know
      // $\theta = \text{acos}\left(\vec{a}\dot\vec{b}\right)$
      // Since we choose $\theta$ = 81.89, we get $\vec{a}\dot\vec{b} = 0.99$
      //**NOTE: We also assume the same thing for flow impinging on a surface
      const Real normDotProd = std::abs(globalStreamDir.dotProduct(normVect));
      if (normDotProd >= 0.99)
        {
          separatedFlow = true;
        }

      const Real kappa_1 = std::max(eps, K1[MD_IX(i, 0)]);
      const Real u_tau = uTau[MD_IX(i, 0)];
      const Real h_0_plus = h0Plus[MD_IX(i, 0)];
      Real q_0 = 0.;
      if (separatedFlow || (h_0_plus < h_nu_plus))
        {
          q_0 = u_tau*h_0_plus;
        }
      else // This increases robustness of code -- no issues with log
        {
          q_0 = u_tau*((1./kappa_1)*std::log(h_0_plus/h_nu_plus) + h_nu_plus);
        }
      const RealVect velBound(D_DECL(interiorCellVel[MD_IX(i, 0)],
                                     interiorCellVel[MD_IX(i, 1)],
                                     interiorCellVel[MD_IX(i, 2)]));
      const Real localVelMag = velBound.vectorLength();
      q_0 = std::min(q_0, localVelMag);
      // Assign the velocities in Cartesian physical space
      D_TERM(a_bndryVelFab[MD_IX(i, 0)] = q_0*localStreamDir[0];,
             a_bndryVelFab[MD_IX(i, 1)] = q_0*localStreamDir[1];,
             a_bndryVelFab[MD_IX(i, 2)] = q_0*localStreamDir[2];);

      // Add a wall-normal velocity component
      const RealVect etaGradVect(D_DECL(physGradEta[MD_IX(i,0)],
                                        physGradEta[MD_IX(i,1)],
                                        physGradEta[MD_IX(i,2)]));
      const Real etaGradMag = globalStreamDir.dotProduct(etaGradVect);
      Real y_vel =
        -faceVelMag[MD_IX(i,0)]*h0[MD_IX(i,0)]*etaGradMag/(2.*eta[MD_IX(i,0)]);
      y_vel = std::min(0.1*q_0, std::max(-0.1*q_0, y_vel));
      D_TERM(a_bndryVelFab[MD_IX(i, 0)] += y_vel*normVect[0];,
             a_bndryVelFab[MD_IX(i, 1)] += y_vel*normVect[1];,
             a_bndryVelFab[MD_IX(i, 2)] += y_vel*normVect[2];);
    }
}

/*--------------------------------------------------------------------*/
//  Compute source term for eta_0 in order to advance eta_0 in time
/** \param[out] a_turbSourceAvgFab
 *                      Cell-averaged turbulent source terms
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_facePntVelGradFxb
 *                      The face-centered velocity gradient
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_RHSfab
 *                      Cell-averaged conservative state update
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_domain
 *                      Domain of this block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics on this level
 *  \param[in]  a_dataIndx
 *                      Index of current box in DisjointBoxLayout
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_box   Cell-centered box over which eta_0 source-term
 *                      is computed
 *//*-----------------------------------------------------------------*/

void
SSV::updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                   FArrayBox&              a_WcellAvgFab,
                   const FArrayBox&        a_JUFab,
                   const FluxBox&          a_facePntVelGradFxb,
                   const FluxBox&          a_WfaceAvgFxb,
                   const FluxBox&          a_fluxFxb,
                   FluxBox&                a_stressFluxFxb,
                   const FArrayBox&        a_RHSfab,
                   const FluxBox&          a_unitNormals,
                   const ProblemDomain&    a_domain,
                   const LevelGridMetrics& a_gridMetrics,
                   const DataIndex&        a_dataIndx,
                   const Box&              a_disjointBox,
                   const Box&              a_box,
                   const Real              a_dt) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::updateEtaZero");

  // Fuse a dir-loop and a side-loop
  for( int dirSide = 0; dirSide != 2*SpaceDim; ++dirSide)
    {
      const int dir = dirSide/2;    // Translate dirSide index into dir
      const int side = dirSide % 2; // Translate dirSide index into side
      Side::LoHiSide whichSide = (side == 1) ? Side::Hi : Side::Lo;
      Side::LoHiSide oppSide = flip(whichSide);

      Box bndryFaceBox;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        bndryFaceBox, a_box, a_domain, dir, whichSide);

      if (!bndryFaceBox.isEmpty())
        {
          CH_assert(bndryFaceBox.ixType().test(dir)); // Check if face box

          Box firstInteriorCellBox = bndryFaceBox; // Wall faces
          firstInteriorCellBox.growDir(dir, oppSide, 1); // +1 interior faces
          firstInteriorCellBox.enclosedCells(dir); // Now, first layer of cells

          int isWall = 0; // Check boundary condition
          hasNoSlipWall(isWall,a_gridMetrics,a_disjointBox,dir,whichSide);

          // We only support the use of walls in this function
          if (isWall)
            {
              const int etaUpdateMethod = 2;
              if (etaUpdateMethod == 1)
                {
                  newUpdateEtaZero(a_turbSourceAvgFab, a_WcellAvgFab, a_JUFab,
                                   a_facePntVelGradFxb[dir], a_WfaceAvgFxb[dir],
                                   a_WfaceAvgFxb, a_stressFluxFxb, a_RHSfab,
                                   a_unitNormals[dir], a_unitNormals, a_domain,
                                   a_gridMetrics, a_dataIndx, whichSide, dir,
                                   a_disjointBox, firstInteriorCellBox,
                                   bndryFaceBox, a_dt);
                }
              else if (etaUpdateMethod == 2)
                {
                  updateEtaZeroLocal(a_turbSourceAvgFab, a_WcellAvgFab, a_JUFab,
                                     a_facePntVelGradFxb[dir],
                                     a_WfaceAvgFxb[dir], a_WfaceAvgFxb,
                                     a_fluxFxb, a_stressFluxFxb, a_RHSfab,
                                     a_unitNormals[dir], a_unitNormals,
                                     a_domain, a_gridMetrics, a_dataIndx,
                                     whichSide, dir, a_disjointBox,
                                     firstInteriorCellBox, bndryFaceBox, a_dt);
                }
              else
                {
                  updateEtaZero(a_turbSourceAvgFab, a_WcellAvgFab,
                                a_facePntVelGradFxb[dir], a_WfaceAvgFxb[dir],
                                a_RHSfab, a_unitNormals[dir], a_domain,
                                a_gridMetrics, a_dataIndx, whichSide, dir,
                                a_disjointBox, firstInteriorCellBox,
                                bndryFaceBox, a_dt);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute source term for eta_0 in order to advance eta_0 in time
/** \param[out] a_turbSourceAvgFab
 *                      Cell-averaged turbulent source terms
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_facePntVelGradFab
 *                      The face-centered velocity gradient
 *  \param[in]  a_WfaceAvgFab
 *                      Face-averaged primitive state
 *  \param[in]  a_RHSfab
 *                      Cell-averaged conservative state update
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_domain
 *                      Domain of this block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics on this level
 *  \param[in]  a_dataIndx
 *                      Index of current box in DisjointBoxLayout
 *  \param[in]  a_side  Current side of box to work on (low or high)
 *  \param[in]  a_dir   Direction of face on which to update eta_0
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_box   Cell-centered box over which eta_0 source-term
 *                      is computed
 *  \param[in] a_boundaryFaceBox
 *                      Box of boundary faces next to which eta_0
 *                      source term is computed
 *//*-----------------------------------------------------------------*/

void
SSV::updateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                   FArrayBox&              a_WcellAvgFab,
                   const FArrayBox&        a_facePntVelGradFab,
                   const FArrayBox&        a_WfaceAvgFab,
                   const FArrayBox&        a_RHSfab,
                   const FArrayBox&        a_unitNormals,
                   const ProblemDomain&    a_domain,
                   const LevelGridMetrics& a_gridMetrics,
                   const DataIndex&        a_dataIndx,
                   const Side::LoHiSide&   a_side,
                   const int               a_dir,
                   const Box&              a_disjointBox,
                   const Box&              a_box,
                   const Box&              a_boundaryFaceBox,
                   const Real              a_dt) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::updateEtaZero");

  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int turbCompBegin = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int cEta = turbCompBegin + a_dir + 1;
  const Real eps = 1.e-50;

  Side::LoHiSide oppSide = flip(a_side);
  const int normSign = sign(a_side);
  const int offsetSign = sign(oppSide);

  // Create the interior-face-box
  Box interiorFaceBox = a_boundaryFaceBox;
  // Shift the box to the interior face
  interiorFaceBox.shift(a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);
  // Get rho*u at the interior face and shift it to the wall
  FABSTACKTEMP(interiorFaceRhoU, interiorFaceBox, 1);
  MD_BOXLOOP(interiorFaceBox, i)
    {
      interiorFaceRhoU[MD_IX(i, 0)] =
        a_WfaceAvgFab[MD_IX(i, cRho)]*a_WfaceAvgFab[MD_IX(i, cVel)];
    }
  interiorFaceRhoU.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // We need to transform the velocity at the interior face
  FABSTACKTEMP(interiorFaceVel, interiorFaceBox, SpaceDim);
  // Copy the velocity from a_WfaceAvgFab into interiorFaceVel
  interiorFaceVel.copy(
    a_WfaceAvgFab, interiorFaceBox, cVel, interiorFaceBox, 0, SpaceDim);
  // Shift the fab to the wall so we can transform it using unit normals
  interiorFaceVel.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform the fab using the unit normals on the wall face
  FORT_FORWARDTRANSFORMF(CHF_FRA(interiorFaceVel),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));

  // Expand a_boundaryFaceBox to include ghost cells
  Box bndryBox = a_boundaryFaceBox;
  IntVect growVect = IntVect_unit;
  growVect[a_dir] = 0;
  bndryBox.grow(growVect);
  bndryBox &= a_domain;
  // Make a box to cover the first interior cells adjacent to the bndry face
  Box firstInteriorCellBox = bndryBox; // Set to cover the ghost cells
  firstInteriorCellBox.growDir(a_dir,oppSide,1); // Include first interior face
  firstInteriorCellBox.enclosedCells(a_dir); // Turn this into cell box
  // Store the velocity for shifting purposes (a_WcellAvgFab stays const)
  FABSTACKTEMP(velFab, firstInteriorCellBox, SpaceDim);
  velFab.copy(a_WcellAvgFab, cVel, 0, SpaceDim); // Fill with a_WcellAvgFab
  velFab.shiftHalf(a_dir, normSign); // Shift velFab to the wall face
  // Store the eta_0 value for shifting purposes
  FABSTACKTEMP(eta, firstInteriorCellBox, 1);
  eta.copy(a_WcellAvgFab, cEta, 0, 1); // Fill with a_WcellAvgFab
  eta.shiftHalf(a_dir, normSign); // Shift eta to the wall face
  // Transform the velocity vector
  const int transformIndx = 0;
  FABSTACKTEMP(locStreamLocTanVect, bndryBox, SpaceDim);
  FORT_FORWARDTRANSFORMGENF(CHF_FRA(locStreamLocTanVect),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(velFab),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(a_unitNormals),
                            CHF_BOX(bndryBox));
  // Set the normal component to zero to obtain the wall-tangent velocity vector
  locStreamLocTanVect.setVal(0.0, a_dir);
  // Normalize the wall-tangent vector -- this is the streamwise unit vector
  PatchMappedFunc::normalize(bndryBox,
                             locStreamLocTanVect,
                             Interval(0, SpaceDim-1));

  // Create a normal vector in normal-tangent space
  FABSTACKTEMP(normVect, bndryBox, SpaceDim);
  normVect.setVal(0.);
  normVect.setVal(1., a_dir);
  // Transform normalVect into Cartesian physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(normVect),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(bndryBox));
  // Compute the wall-normal derivative of the wall-normal velocity
  FABSTACKTEMP(modVelGrad, bndryBox, SpaceDim*SpaceDim);
  modVelGrad.setVal(0.); // Most values need to be zero
  MD_BOXLOOP(bndryBox, i)
    {
      Real normDeriv = 0.;
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          normDeriv += normVect[MD_IX(i, comp)]*(
            D_TERM(
              a_facePntVelGradFab[MD_IX(i,srtIdx(comp,0))]*
                normVect[MD_IX(i,0)],
            + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,1))]*
                normVect[MD_IX(i,1)],
            + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,2))]*
                normVect[MD_IX(i,2)]));
        }
      modVelGrad[MD_IX(i, srtIdx(a_dir,a_dir))] = normDeriv;
    }
  // Fill the normal-derivs of the tangential velocities using model
  for (int velComp = 0; velComp != SpaceDim; ++velComp)
    {
      if (velComp != a_dir)
        {
          const int cGrad = a_dir + SpaceDim*velComp;
          MD_BOXLOOP(bndryBox, i)
            {
              modVelGrad[MD_IX(i, cGrad)] =
                locStreamLocTanVect[MD_IX(i,velComp)]*eta[MD_IX(i,0)];
            }
        }
    }
  // Transform back to Cartesian space
  FABSTACKTEMP(tempState, bndryBox, SpaceDim*SpaceDim);
  // Right multiply grad(u) by a_unitNormals
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int cGrad = col + SpaceDim*row;
          MD_BOXLOOP(bndryBox, i)
            {
              tempState[MD_IX(i, cGrad)] =
                D_TERM(
                  modVelGrad[MD_IX(i,srtIdx(row,0))]*
                    a_unitNormals[MD_IX(i,srtIdx(0,col))],
                + modVelGrad[MD_IX(i,srtIdx(row,1))]*
                    a_unitNormals[MD_IX(i,srtIdx(1,col))],
                + modVelGrad[MD_IX(i,srtIdx(row,2))]*
                    a_unitNormals[MD_IX(i,srtIdx(2,col))]);
            }
        }
    }
  // Left multiply temp state by a_unitNormals^T
  Interval tensorIntv(0, SpaceDim*SpaceDim - 1);
  FArrayBox cartGradVel(tensorIntv, modVelGrad); // Alias to avoid allocation
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int cGrad = col + SpaceDim*row;
          MD_BOXLOOP(bndryBox, i)
            {
              cartGradVel[MD_IX(i, cGrad)] =
                D_TERM(
                  a_unitNormals[MD_IX(i,srtIdx(0,row))]*
                    tempState[MD_IX(i,srtIdx(0,col))],
                + a_unitNormals[MD_IX(i,srtIdx(1,row))]*
                    tempState[MD_IX(i,srtIdx(1,col))],
                + a_unitNormals[MD_IX(i,srtIdx(2,row))]*
                    tempState[MD_IX(i,srtIdx(2,col))]);
            }
        }
    }
  // Calculate the viscosity (and thermal conductivity unfortunately)
  FABSTACKTEMP(muFab, bndryBox, 1);
  FABSTACKTEMP(kappaFab, bndryBox, 1);
  CRDparam::g_CRDPhysics->calcCoeffKappaMu(
    bndryBox, muFab, kappaFab, a_WfaceAvgFab);
  // Fill the new viscous stress tensor with the strain-rate tensor
  FArrayBox vstFab(tensorIntv, tempState); // Alias to avoid allocation
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int vstComp = row*SpaceDim + col;
          MD_BOXLOOP(bndryBox, i)
            {
              vstFab[MD_IX(i, vstComp)] = 
                cartGradVel[MD_IX(i, srtIdx(row, col))]
              + cartGradVel[MD_IX(i, srtIdx(col, row))];
            }
        }
    }
  // Now subtract the divergence of velocity
  MD_BOXLOOP(bndryBox, i)
    {
      Real divVel = 2.*CRDparam::g_lambda*(D_TERM(
        cartGradVel[MD_IX(i, srtIdx(0,0))],
      + cartGradVel[MD_IX(i, srtIdx(1,1))],
      + cartGradVel[MD_IX(i, srtIdx(2,2))]));
      for (int row = 0; row != SpaceDim; ++row)
        {
          vstFab[MD_IX(i, srtIdx(row,row))] -= divVel;
        }
    }
  // Now multiply by the viscosity
  for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
    {
      MD_BOXLOOP(bndryBox, i)
        {
          vstFab[MD_IX(i, comp)] *= muFab[MD_IX(i, 0)];
        }
    }
  // Store vstFab in vstDiffFab for difference calculation later
  FABSTACKTEMP(vstDiffFab, bndryBox, SpaceDim*SpaceDim);
  vstDiffFab.copy(vstFab);
  // Now convolve the new viscous stress tensor
  FArrayBox vstAvgFab(tensorIntv, cartGradVel); // Alias to avoid allocation
  MOLUtilFunc::deconvolveCenterFace(
    vstAvgFab,vstFab,a_boundaryFaceBox,a_domain,a_dir);
  // Compute the original viscous stress tensor for comparison
  FArrayBox origVSTFab(tensorIntv, vstFab); // Alias to avoid allocation
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int vstComp = row*SpaceDim + col;
          MD_BOXLOOP(bndryBox, i)
            {
              origVSTFab[MD_IX(i, vstComp)] = 
                a_facePntVelGradFab[MD_IX(i, srtIdx(row, col))]
              + a_facePntVelGradFab[MD_IX(i, srtIdx(col, row))];
            }
        }
    }
  // Now subtract the divergence of velocity
  MD_BOXLOOP(bndryBox, i)
    {
      Real divVel = 2.*CRDparam::g_lambda*(D_TERM(
        a_facePntVelGradFab[MD_IX(i, srtIdx(0,0))],
      + a_facePntVelGradFab[MD_IX(i, srtIdx(1,1))],
      + a_facePntVelGradFab[MD_IX(i, srtIdx(2,2))]));
      for (int row = 0; row != SpaceDim; ++row)
        {
          origVSTFab[MD_IX(i, srtIdx(row,row))] -= divVel;
        }
    }
  // Now multiply by the viscosity
  for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
    {
      MD_BOXLOOP(bndryBox, i)
        {
          origVSTFab[MD_IX(i, comp)] *= muFab[MD_IX(i, 0)];
        }
    }
  // Subtract the two face-centered viscous stress tensors
  for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
    {
      MD_BOXLOOP(bndryBox, i)
        {
          vstDiffFab[MD_IX(i, comp)] -= origVSTFab[MD_IX(i, comp)];
        }
    }
  // Store vstAvgFab in vstDiffAvgFab for difference calculation later
  FABSTACKTEMP(vstDiffAvgFab, bndryBox, SpaceDim*SpaceDim);
  vstDiffAvgFab.copy(vstAvgFab);
  // Now convolve the original viscous stress tensor
  FArrayBox origVSTAvgFab(tensorIntv, vstAvgFab); // Alias to avoid allocation
  MOLUtilFunc::deconvolveCenterFace(
    origVSTAvgFab, origVSTFab, a_boundaryFaceBox, a_domain, a_dir);
  // Subtract the two viscous stress tensors
  for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
    {
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          vstDiffAvgFab[MD_IX(i, comp)] -= origVSTAvgFab[MD_IX(i, comp)];
        }
    }
  // Map the viscous stress tensors
  Box mappedBox = a_boundaryFaceBox;
  mappedBox.growDir(a_dir, oppSide, 1); // Include first interior face
  mappedBox.enclosedCells(a_dir); // Turn this into a cell-box
  Box dataBox = mappedBox;
  dataBox.grow(1); // Grow all directions by 1
  const FluxBox& Nfxb = a_gridMetrics.m_N[a_dataIndx];
  FLUXBOXSTACKTEMP(mappedFlux, mappedBox, SpaceDim);
  FLUXBOXSTACKTEMP(pntFlux, dataBox, SpaceDim*SpaceDim);
  FLUXBOXSTACKTEMP(avgFlux, dataBox, SpaceDim*SpaceDim);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      mappedFlux[dir].setVal(0.);
      pntFlux[dir].setVal(0.);
      avgFlux[dir].setVal(0.);
      if (dir == a_dir)
        {
          pntFlux[dir].copy(vstDiffFab);
          avgFlux[dir].copy(vstDiffAvgFab);
        }
    }
  Interval fluxIntv(0, SpaceDim - 1);
  bool fourthOrder = false;
  if (CRDparam::g_faceConvolveFlatten < 2 ||
      CRDparam::g_faceDeconvolveFlatten < 2)
    {
      fourthOrder = true;
    }
  const BlockCoordSys* blockCoordSys =
    a_gridMetrics.getCoordSys(a_disjointBox);
  blockCoordSys->computeMetricTermProductAverage(
    mappedFlux, avgFlux, Nfxb, SpaceDim, pntFlux, mappedBox, fourthOrder,
    fluxIntv, fluxIntv, 1, &a_domain);

  // Divide this by the discretization size
  FABSTACKTEMP(divTau, a_boundaryFaceBox, SpaceDim);
  const Real invDx = 1./(a_gridMetrics.dxVect()[0]); // Assuming isotropic dx
  const FArrayBox& mappedFluxFab = mappedFlux[a_dir];
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          divTau[MD_IX(i, comp)] = normSign*invDx*mappedFluxFab[MD_IX(i,comp)];
        }
    }

  // Bound eta_0 update with reasonable values
  // We assume that a local solution to the law-of-the-wall provides
  // a reasonable estimate for the instantaneous, local eta value
  FABSTACKTEMP(etaEstFab, bndryBox, 1);
  etaEst(etaEstFab, a_WcellAvgFab, a_gridMetrics, a_unitNormals,
         bndryBox, a_disjointBox, a_dir, a_side);
  FABSTACKTEMP(etaBound, bndryBox, 2);
  etaBound.setVal(0.);
  MD_BOXLOOP(bndryBox, i)
    {
      etaBound[MD_IX(i, 0)] = 0.1*(etaEstFab[MD_IX(i, 0)]);
      etaBound[MD_IX(i, 1)] = 10.*(etaEstFab[MD_IX(i, 0)]);
    }

  // Switch between original method and temporary method
  const int method = 0;

  // Get cell J from gridMetrics
  const FArrayBox& cellJ = a_gridMetrics.m_J[a_dataIndx];
  // Now we loop over the wall-parallel interior-cell planes
  FABSTACKTEMP(etaUpdateFab, a_box, 2); // One comp for eta, one for J
  etaUpdateFab.setVal(0.);
  for (int j = 0; j != CRDparam::g_wallModelCrsRatio; ++j)
    {
      // Get the wall-parallel interior-cell box
      Box cellBox = a_box;
      cellBox.shift(a_dir, offsetSign*j);
      // We transform the RHS(momentum) into normal, streamwise, spanwise
      FABSTACKTEMP(flowMomRHS, cellBox, SpaceDim);
      // Fill the fab with the RHS data for momentum
      flowMomRHS.copy(a_RHSfab, cellBox, cVel, cellBox, 0, SpaceDim);
      // Shift the fab to the wall so we can transform it using unit normals
      flowMomRHS.shiftHalf(a_dir, normSign*(2*j + 1));
      // Before transforming, add the difference of the VST
      // Only do this if this is the wall-adjacent layer of cells
      if (j == 0 && m_bndryVSTCorrection)
        {
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              MD_BOXLOOP(a_boundaryFaceBox, i)
                {
                  flowMomRHS[MD_IX(i, comp)] += divTau[MD_IX(i, comp)];
                }
            }
        }
      // Transform the fab using the unit normals on the wall face
      FORT_FORWARDTRANSFORMF(CHF_FRA(flowMomRHS),
                             CHF_CONST_FRA(a_unitNormals),
                             CHF_BOX(a_boundaryFaceBox));
      // Shift the fab back to the first interior cell
      flowMomRHS.shiftHalf(a_dir, offsetSign*(2*j + 1));
      // We also need to transform velocity into normal, streamwise, spanwise
      FABSTACKTEMP(flowVel, cellBox, SpaceDim);
      // Fill the fab with velocity data
      flowVel.copy(a_WcellAvgFab, cellBox, cVel, cellBox, 0, SpaceDim);
      // Shift the fab to the wall so we can transform it using unit normals
      flowVel.shiftHalf(a_dir, normSign*(2*j + 1));
      // Transform the fab using the unit normals on the wall face
      FORT_FORWARDTRANSFORMF(CHF_FRA(flowVel),
                             CHF_CONST_FRA(a_unitNormals),
                             CHF_BOX(a_boundaryFaceBox));
      // Shift the fab back to the first interior cell
      flowVel.shiftHalf(a_dir, offsetSign*(2*j + 1));
      const int wallAdjCellIndx = a_box.smallEnd(a_dir);
      // Compute the eta_0 update
      if (method == 0) // Previous method
        {
          MD_BOXLOOP(cellBox, i)
            {
              const Real rho = a_WcellAvgFab[MD_IX(i, cRho)];
              const Real invRho = 1./rho;
              const Real rhoRHS = a_RHSfab[MD_IX(i, cRho)];
              // The stream-coordinate oriented velocity in first interior cell
              RealVect vel(D_DECL(flowVel[MD_IX(i, 0)],
                                  flowVel[MD_IX(i, 1)],
                                  flowVel[MD_IX(i, 2)]));
              // Get stream-coord oriented RHS for just velocity (not momentum)
              RealVect streamVelRHS(
                D_DECL(invRho*(flowMomRHS[MD_IX(i, 0)] - rhoRHS*vel[0]),
                       invRho*(flowMomRHS[MD_IX(i, 1)] - rhoRHS*vel[1]),
                       invRho*(flowMomRHS[MD_IX(i, 2)] - rhoRHS*vel[2])));
              // Now, we zero the normal component of velocity
              vel[a_dir] = 0.;
              // Similarly, we zero the normal component of streamVelRHS
              streamVelRHS[a_dir] = 0.;
              // Compute the wall-tangential velocity magnitude
              const Real q = std::sqrt(vel.radSquared()) + eps;
              // Now we essentially normalize the velocity by q
              const RealVect normedVel = vel/q;
              // Here's time-derivative of the wall-tangent velocity magnitude
              const Real q_timeDeriv = normedVel.dotProduct(streamVelRHS);
              // Get eta_0 in the first cell (this is where we store it for now)
              IntVect bndryIV = MD_GETIV(i);
              bndryIV[a_dir] = wallAdjCellIndx;
              Real eta_0 = a_WcellAvgFab[MD_IV(bndryIV, cEta)];
              // Finally, the eta_0 update
              Real etaUpdate = 2.*eta_0*q_timeDeriv;
              if ((eta_0 <= 0.) && (etaUpdate < 0.))
                {
                  etaUpdate = std::abs(2.*q_timeDeriv);
                }
              etaUpdateFab[MD_IV(bndryIV, 0)] += etaUpdate;
              // Sum the J values to update interior cells appropriately
              etaUpdateFab[MD_IV(bndryIV, 1)] += cellJ[MD_IX(i, 0)];
            }
        }
      else if (method == 1) // Preserving as much conservation as possible
        {
          MD_BOXLOOP(cellBox, i)
            {
              const Real rhoU_RHS = a_RHSfab[MD_IX(i, cVel)];
              // Get eta_0 in the first cell (this is where we store it for now)
              IntVect bndryIV = MD_GETIV(i);
              bndryIV[a_dir] = wallAdjCellIndx;
              const Real eta_0 = a_WcellAvgFab[MD_IV(bndryIV, cEta)];
              // Finally, the eta_0 update
              Real etaUpdate = 2.*eta_0*rhoU_RHS;
              if ((eta_0 <= 0.) && (etaUpdate < 0.))
                {
                  etaUpdate = std::abs(2.*rhoU_RHS);
                }
              etaUpdateFab[MD_IV(bndryIV, 0)] += etaUpdate;
              // Sum the J values to update interior cells appropriately
              etaUpdateFab[MD_IV(bndryIV, 1)] += cellJ[MD_IX(i, 0)];
            }
        }
    }
  if (method == 0)
    {
      // Now, loop over the near-wall box and set eta updates
      MD_BOXLOOP(a_box, i)
        {
          // The stream-coordinate oriented velocity at the interior face
          RealVect faceVel(D_DECL(interiorFaceVel[MD_IX(i, 0)],
                                  interiorFaceVel[MD_IX(i, 1)],
                                  interiorFaceVel[MD_IX(i, 2)]));
          // Now, we zero the normal component of velocity
          faceVel[a_dir] = 0.;
          // Compute the wall-tangential vel magnitude at interior face
          const Real q_h = faceVel.vectorLength() + eps;
          // Create a box of wall-normal cells to fill with etaUpdate
          Box wallNormalBox = a_box;
          // Get the current IntVect
          IntVect currVect = MD_GETIV(i);
          if (a_side == Side::Lo)
            {
              // If bndry is low side, set low side of wallNormalBox to currVect
              wallNormalBox.setSmall(currVect);
              // Also, get the wall-normal component of a_disjointBox high end
              int highSide = a_disjointBox.bigEnd()[a_dir];
              // Set the wall-normal component of currVect to highSide
              currVect.setVal(a_dir, highSide);
              // Set the high side of wallNormalBox to the new currVect
              wallNormalBox.setBig(currVect);
            }
          else
            {
              // If bndry is high, set high side of wallNormalBox to currVect
              wallNormalBox.setBig(currVect);
              // Also, get the wall-normal component of a_disjointBox low end
              int lowSide = a_disjointBox.smallEnd()[a_dir];
              // Set the wall-normal component of currVect to lowSide
              currVect.setVal(a_dir, lowSide);
              // Set the low side of wallNormalBox to the new currVect
              wallNormalBox.setSmall(currVect);
            }
          Real localJ = cellJ[MD_IX(i, 0)];
          Real Jsum = etaUpdateFab[MD_IX(i, 1)];
          Real invQh = 1./q_h;
          Real etaUpdate = invQh*etaUpdateFab[MD_IX(i, 0)]*localJ/Jsum;
          // We now enforce eta boundedness
          const Real etaMaxBound = etaBound[MD_IX(i, 1)];
          const Real etaMinBound = etaBound[MD_IX(i, 0)];
          const Real eta_0 = eta[MD_IX(i, 0)];
          const Real etaUpdateMaxBound = localJ*(etaMaxBound - eta_0)/a_dt;
          const Real etaUpdateMinBound = localJ*(etaMinBound - eta_0)/a_dt;
          etaUpdate =
            std::max(etaUpdateMinBound, std::min(etaUpdate, etaUpdateMaxBound));
          // Now loop over wallNormalBox and assign etaUpdate to all cells
          MD_BOXLOOP(wallNormalBox, j)
            {
              a_turbSourceAvgFab[MD_IX(j, cEta)] = etaUpdate;
            }
        }
    }
  else if (method == 1)
    {
      // Now, loop over the near-wall box and set eta updates
      MD_BOXLOOP(a_box, i)
        {
          // Interior face streamwise momentum
          const Real rhoU = interiorFaceRhoU[MD_IX(i, 0)];
          // Create a box of wall-normal cells to fill with etaUpdate
          Box wallNormalBox = a_box;
          // Get the current IntVect
          IntVect currVect = MD_GETIV(i);
          if (a_side == Side::Lo)
            {
              // If bndry is low side, set low side of wallNormalBox to currVect
              wallNormalBox.setSmall(currVect);
              // Also, get the wall-normal component of a_disjointBox high end
              int highSide = a_disjointBox.bigEnd()[a_dir];
              // Set the wall-normal component of currVect to highSide
              currVect.setVal(a_dir, highSide);
              // Set the high side of wallNormalBox to the new currVect
              wallNormalBox.setBig(currVect);
            }
          else
            {
              // If bndry is high, set high side of wallNormalBox to currVect
              wallNormalBox.setBig(currVect);
              // Also, get the wall-normal component of a_disjointBox low end
              int lowSide = a_disjointBox.smallEnd()[a_dir];
              // Set the wall-normal component of currVect to lowSide
              currVect.setVal(a_dir, lowSide);
              // Set the low side of wallNormalBox to the new currVect
              wallNormalBox.setSmall(currVect);
            }
          Real localJ = cellJ[MD_IX(i, 0)];
          Real Jsum = etaUpdateFab[MD_IX(i, 1)];
          Real invQh = 1./rhoU;
          // Real c_0 = 0.1; // Bound how rapidly eta_0 responds to
          //                 // local flow variation
          Real etaUpdate = invQh*etaUpdateFab[MD_IX(i, 0)]*localJ/Jsum;
          // Now we check that the update
          // 1) doesn't push eta_0 below the physically resolved gradient
          // 2) doesn't push eta_0 above the estimated gradient -- we can get
          //    this value from the Musker or Spalding profiles
          // Now loop over wallNormalBox and assign etaUpdate to all cells
          MD_BOXLOOP(wallNormalBox, j)
            {
              a_turbSourceAvgFab[MD_IX(j, cEta)] = etaUpdate;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute source term for eta_0 in order to advance eta_0 in time
/** \param[out] a_turbSourceAvgFab
 *                      Cell-averaged turbulent source terms
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_facePntVelGradFab
 *                      The face-centered velocity gradient
 *  \param[in]  a_WfaceAvgFab
 *                      Face-averaged primitive state
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state (FluxBox)
 *  \param[in]  a_stressFluxFxb
 *                      Combined stresses (SGS and viscous) for LES wall-model
 *  \param[in]  a_RHSfab
 *                      Cell-averaged conservative state update
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_domain
 *                      Domain of this block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics on this level
 *  \param[in]  a_dataIndx
 *                      Index of current box in DisjointBoxLayout
 *  \param[in]  a_side  Current side of box to work on (low or high)
 *  \param[in]  a_dir   Direction of face on which to update eta_0
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_box   Cell-centered box over which eta_0 source-term
 *                      is computed
 *  \param[in] a_boundaryFaceBox
 *                      Box of boundary faces next to which eta_0
 *                      source term is computed
 *//*-----------------------------------------------------------------*/
void
SSV::newUpdateEtaZero(FArrayBox&              a_turbSourceAvgFab,
                      FArrayBox&              a_WcellAvgFab,
                      const FArrayBox&        a_JUFab,
                      const FArrayBox&        a_facePntVelGradFab,
                      const FArrayBox&        a_WfaceAvgFab,
                      const FluxBox&          a_WfaceAvgFxb,
                      FluxBox&                a_stressFluxFxb,
                      const FArrayBox&        a_RHSfab,
                      const FArrayBox&        a_unitNormals,
                      const FluxBox&          a_unitNormalsFxb,
                      const ProblemDomain&    a_domain,
                      const LevelGridMetrics& a_gridMetrics,
                      const DataIndex&        a_dataIndx,
                      const Side::LoHiSide&   a_side,
                      const int               a_dir,
                      const Box&              a_disjointBox,
                      const Box&              a_box,
                      const Box&              a_boundaryFaceBox,
                      const Real              a_dt) const
{
  if (!m_useWallModel) return;
  CH_TIME("SSV::newUpdateEtaZero");

  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cPres = CRDparam::g_CRDPhysics->pressureIndex();
  const int cTurbBegin = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  const int cEta = cTurbBegin + a_dir + 1;
  const Real invDx = 1./(a_gridMetrics.dxVect()[a_dir]);
  const Real eps = 1e-50;

  Side::LoHiSide oppSide = flip(a_side);
  const int normSign = sign(a_side);
  const int offsetSign = sign(oppSide);

  // Get the face-plane-box where all outer-flow variables will be evaluated
  Box innerFaceBox = a_boundaryFaceBox;
  // Shift the box to the interior face
  innerFaceBox.shift(a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);
  // Box of cells surrounding innerFaceBox
  Box innerCellBox = innerFaceBox;
  innerCellBox.grow(a_dir, 1);
  innerCellBox.enclosedCells(a_dir);

  // As a hack for all tangential derivatives on faces, we're going to
  // temporarily shift all the face values into the lower adjacent cells
  // so we can maintain consistency with the cell-derivative function
  Box adjacentCellBox = innerFaceBox;
  adjacentCellBox.shiftHalf(a_dir, -1);

  // Get cell J from gridMetrics
  const FArrayBox& cellJ = a_gridMetrics.m_J[a_dataIndx];
  // Get face N^T/J from gridMetrics
  const FluxBox& NtJ = a_gridMetrics.m_NtJ[a_dataIndx];

  // Compute the pressure gradient
  // a) using only differences of face values on opposite sides of cells,
  //    (this requires face-normal interpolation),
  // b) or, by using cell-average values and face-tangent values,
  //    (evaluating all derivatives at face, no interpolation necessary)
  // // ----------------------------------------------------------------- //
  // // Method (a)
  // FABSTACKTEMP(pressGradAFace, innerFaceBox, SpaceDim);
  // FABSTACKTEMP(pressGradACell, innerCellBox, SpaceDim);
  // // Compute cell-averaged pressure gradients
  // for (int deriv = 0; deriv != SpaceDim; ++deriv)
  //   {
  //     const int MD_ID(ii, deriv);
  //     const FArrayBox& WfaceAvgFabDir = a_WfaceAvgFxb[deriv];
  //     MD_BOXLOOP(innerCellBox, i)
  //       {
  //         pressGradACell[MD_IX(i, deriv)] = invDx*(
  //           WfaceAvgFabDir[MD_OFFSETIX(i, +, ii, cPres)] -
  //           WfaceAvgFabDir[MD_IX(i, cPres)]);
  //       }
  //   }
  // // Interpolate cell-averaged pressure gradients to the faces
  // for (int comp = 0; comp != SpaceDim; ++comp)
  //   {
  //     int orderOfAccuracy = 2;
  //     PatchMappedFunc::faceAvgValFromCellAvgCS(
  //       pressGradAFace, pressGradACell, comp, comp, a_domain, innerFaceBox,
  //       a_dir, orderOfAccuracy, false);
  //   }
  // // ----------------------------------------------------------------- //
  // Method (b)
  FABSTACKTEMP(pressGradB, innerFaceBox, SpaceDim);
  for (int deriv = 0; deriv != SpaceDim; ++deriv)
    {
      if (deriv == a_dir)
        {
          PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(
            pressGradB, a_WcellAvgFab, cPres, deriv, a_domain, innerFaceBox,
            a_gridMetrics.dxVect(), a_dir, false, false);
        }
      else
        {
          Box copyFaceBox = innerFaceBox;
          copyFaceBox.grow(deriv, 1);
          copyFaceBox &= a_domain;
          FABSTACKTEMP(facePres, copyFaceBox, 1);
          const FArrayBox& WfaceAvgFabDir = a_WfaceAvgFxb[a_dir];
          facePres.copy(WfaceAvgFabDir,copyFaceBox,cPres,copyFaceBox,0,1);
          // Just shift into the next lower cell (really doesn't matter where
          facePres.shiftHalf(a_dir, -1);
          pressGradB.shiftHalf(a_dir, -1);
          PatchMappedFunc::cellAvgDerivFromCellAvgCS(
            pressGradB, facePres, 0, deriv, a_domain, adjacentCellBox,
            a_gridMetrics.dxVect(), deriv, false, false);
          // Now make sure to shift things back as necessary
          pressGradB.shiftHalf(a_dir, 1);
        }
    }
  // It seems that method (b) works better, so use it for future derivatives

  // Transform the pressure gradient into physical space using N^T/J
  FABSTACKTEMP(physGradPres, innerFaceBox, SpaceDim);
  const FArrayBox& NtJFab = NtJ[a_dir];
  for (int col = 0; col != SpaceDim; ++col)
    {
      MD_BOXLOOP(innerFaceBox, i)
        {
          physGradPres[MD_IX(i, col)] = D_TERM(
            pressGradB[MD_IX(i,0)]*NtJFab[MD_IX(i,srtIdx(0,col))],
            + pressGradB[MD_IX(i,1)]*NtJFab[MD_IX(i,srtIdx(1,col))],
            + pressGradB[MD_IX(i,2)]*NtJFab[MD_IX(i,srtIdx(2,col))]);
        }
    }
  // Shift pressure gradient to the wall
  physGradPres.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform the pressure gradient into wall-normal space using a_unitNormals
  FORT_FORWARDTRANSFORMF(CHF_FRA(physGradPres),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Zero the normal derivative component
  physGradPres.setVal(0., a_dir);
  // Transform the pressure gradient back into physical space with unitNormals
  FORT_REVERSETRANSFORMF(CHF_FRA(physGradPres),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Shift pressure gradient away from wall
  physGradPres.shift(a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);

  // We need two pieces for all further fluxes
  // (a) tangential derivatives of fluxes
  // (b) normal derivatives of fluxes (evaluated from wall to interior face)

  // Combine the velocity, SGS, and molecular stress tensors together
  Box grownFaceBox = innerFaceBox;
  grownFaceBox.grow(1); // Grow by one face in all directions, even normal
  grownFaceBox.grow(a_dir, -1);
  grownFaceBox &= a_domain;
  // Compute 2nd-order rho*u*u and subtract the viscous + SGS fluxes
  FArrayBox& viscStress = a_stressFluxFxb[a_dir];
  FABSTACKTEMP(combinedStressTensor, grownFaceBox, SpaceDim*SpaceDim);
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int comp = srtIdx(row, col);
          MD_BOXLOOP(grownFaceBox, i)
            {
              const Real rho = a_WfaceAvgFab[MD_IX(i, cRho)];
              const Real vel1 = a_WfaceAvgFab[MD_IX(i, cVel+row)];
              const Real vel2 = a_WfaceAvgFab[MD_IX(i, cVel+col)];
              const Real viscFlux = viscStress[MD_IX(i, comp)];
              combinedStressTensor[MD_IX(i, comp)] = rho*vel1*vel2 - viscFlux;
            }
        }
    }

  // Expand a_boundaryFaceBox to include ghost cells
  Box bndryBox = a_boundaryFaceBox;
  IntVect growVect = IntVect_unit;
  growVect[a_dir] = 0;
  bndryBox.grow(growVect);
  bndryBox &= a_domain;
  // Make a box to cover the first interior cells adjacent to the bndry face
  Box firstInteriorCellBox = bndryBox; // Set to cover the ghost cells
  firstInteriorCellBox.growDir(a_dir,oppSide,1); // Include 1st interior face
  firstInteriorCellBox.enclosedCells(a_dir); // Turn this into cell box
  // Store the eta_0 value for shifting purposes
  FABSTACKTEMP(eta, firstInteriorCellBox, 1);
  // Fill using a_WcellAvgFab
  MD_BOXLOOP(firstInteriorCellBox, i)
    {
      eta[MD_IX(i, 0)] = a_WcellAvgFab[MD_IX(i, cEta)];
    }
  eta.shiftHalf(a_dir, normSign); // Shift eta to the wall face
  FABSTACKTEMP(tempState, bndryBox, SpaceDim*SpaceDim);
  // We need to transform the velocity at the interior face
  FABSTACKTEMP(interiorFaceVel, innerFaceBox, SpaceDim);
  // Copy the velocity from a_WfaceAvgFab into interiorFaceVel
  interiorFaceVel.copy(
    a_WfaceAvgFab, innerFaceBox, cVel, innerFaceBox, 0, SpaceDim);
  // Shift the fab to the wall so we can transform it using unit normals
  interiorFaceVel.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform the fab using the unit normals on the wall face
  FORT_FORWARDTRANSFORMF(CHF_FRA(interiorFaceVel),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Compute the wall viscous flux
  {
    // Store the velocity for shifting purposes (a_WcellAvgFab stays const)
    FABSTACKTEMP(velFab, firstInteriorCellBox, SpaceDim);
    velFab.copy(a_WcellAvgFab, cVel, 0, SpaceDim); // Fill w/ a_WcellAvgFab
    velFab.shiftHalf(a_dir, normSign); // Shift velFab to the wall face
    // Transform the velocity vector
    const int transformIndx = 0;
    FABSTACKTEMP(locStreamLocTanVect, bndryBox, SpaceDim);
    FORT_FORWARDTRANSFORMGENF(CHF_FRA(locStreamLocTanVect),
                              CHF_INT(transformIndx),
                              CHF_CONST_FRA(velFab),
                              CHF_INT(transformIndx),
                              CHF_CONST_FRA(a_unitNormals),
                              CHF_BOX(bndryBox));
    // Set normal component to zero to obtain the wall-tangent velocity vector
    locStreamLocTanVect.setVal(0.0, a_dir);
    // Normalize the wall-tangent vector -- this is the streamwise unit vector
    PatchMappedFunc::normalize(bndryBox,
                               locStreamLocTanVect,
                               Interval(0, SpaceDim-1));

    // Create a normal vector in normal-tangent space
    FABSTACKTEMP(normVect, bndryBox, SpaceDim);
    normVect.setVal(0.);
    normVect.setVal(1., a_dir);
    // Transform normalVect into Cartesian physical space
    FORT_REVERSETRANSFORMF(CHF_FRA(normVect),
                           CHF_CONST_FRA(a_unitNormals),
                           CHF_BOX(bndryBox));
    // Compute the wall-normal derivative of the wall-normal velocity
    FABSTACKTEMP(modVelGrad, bndryBox, SpaceDim*SpaceDim);
    modVelGrad.setVal(0.); // Most values need to be zero
    MD_BOXLOOP(bndryBox, i)
      {
        Real normDeriv = 0.;
        for (int comp = 0; comp != SpaceDim; ++comp)
          {
            normDeriv += normVect[MD_IX(i, comp)]*(
              D_TERM(
                a_facePntVelGradFab[MD_IX(i,srtIdx(comp,0))]*
                normVect[MD_IX(i,0)],
                + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,1))]*
                normVect[MD_IX(i,1)],
                + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,2))]*
                normVect[MD_IX(i,2)]));
          }
        modVelGrad[MD_IX(i, srtIdx(a_dir,a_dir))] = normDeriv;
      }
    // Fill the normal-derivs of the tangential velocities using model
    for (int velComp = 0; velComp != SpaceDim; ++velComp)
      {
        if (velComp != a_dir)
          {
            const int cGrad = a_dir + SpaceDim*velComp;
            MD_BOXLOOP(bndryBox, i)
              {
                modVelGrad[MD_IX(i, cGrad)] =
                  locStreamLocTanVect[MD_IX(i,velComp)]*eta[MD_IX(i,0)];
              }
          }
      }
    // Transform back to Cartesian space
    // Right multiply grad(u) by a_unitNormals
    for (int row = 0; row != SpaceDim; ++row)
      {
        for (int col = 0; col != SpaceDim; ++col)
          {
            const int cGrad = col + SpaceDim*row;
            MD_BOXLOOP(bndryBox, i)
              {
                tempState[MD_IX(i, cGrad)] =
                  D_TERM(
                    modVelGrad[MD_IX(i,srtIdx(row,0))]*
                    a_unitNormals[MD_IX(i,srtIdx(0,col))],
                    + modVelGrad[MD_IX(i,srtIdx(row,1))]*
                    a_unitNormals[MD_IX(i,srtIdx(1,col))],
                    + modVelGrad[MD_IX(i,srtIdx(row,2))]*
                    a_unitNormals[MD_IX(i,srtIdx(2,col))]);
              }
          }
      }
    // Left multiply temp state by a_unitNormals^T
    Interval tensorIntv(0, SpaceDim*SpaceDim - 1);
    FArrayBox cartGradVel(tensorIntv, modVelGrad); // Alias to avoid allocation
    for (int row = 0; row != SpaceDim; ++row)
      {
        for (int col = 0; col != SpaceDim; ++col)
          {
            const int cGrad = col + SpaceDim*row;
            MD_BOXLOOP(bndryBox, i)
              {
                cartGradVel[MD_IX(i, cGrad)] =
                  D_TERM(
                    a_unitNormals[MD_IX(i,srtIdx(0,row))]*
                    tempState[MD_IX(i,srtIdx(0,col))],
                    + a_unitNormals[MD_IX(i,srtIdx(1,row))]*
                    tempState[MD_IX(i,srtIdx(1,col))],
                    + a_unitNormals[MD_IX(i,srtIdx(2,row))]*
                    tempState[MD_IX(i,srtIdx(2,col))]);
              }
          }
      }
    // Calculate the viscosity (and thermal conductivity unfortunately)
    FABSTACKTEMP(muFab, bndryBox, 1);
    FABSTACKTEMP(kappaFab, bndryBox, 1);
    CRDparam::g_CRDPhysics->calcCoeffKappaMu(
      bndryBox, muFab, kappaFab, a_WfaceAvgFab);
    // Fill the new viscous stress tensor with the strain-rate tensor
    FArrayBox vstFab(tensorIntv, tempState); // Alias to avoid allocation
    for (int row = 0; row != SpaceDim; ++row)
      {
        for (int col = 0; col != SpaceDim; ++col)
          {
            const int vstComp = row*SpaceDim + col;
            MD_BOXLOOP(bndryBox, i)
              {
                vstFab[MD_IX(i, vstComp)] = 
                  cartGradVel[MD_IX(i, srtIdx(row, col))]
                  + cartGradVel[MD_IX(i, srtIdx(col, row))];
              }
          }
      }
    // Now subtract the divergence of velocity
    MD_BOXLOOP(bndryBox, i)
      {
        Real divVel = 2.*CRDparam::g_lambda*(
          D_TERM(cartGradVel[MD_IX(i, srtIdx(0,0))],
                 + cartGradVel[MD_IX(i, srtIdx(1,1))],
                 + cartGradVel[MD_IX(i, srtIdx(2,2))]));
        for (int row = 0; row != SpaceDim; ++row)
          {
            vstFab[MD_IX(i, srtIdx(row,row))] -= divVel;
          }
      }
    // Now multiply by the viscosity
    for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
      {
        MD_BOXLOOP(bndryBox, i)
          {
            vstFab[MD_IX(i, comp)] *= muFab[MD_IX(i, 0)];
          }
      }
  }

  // Compute the velocity stress tensor update
  FABSTACKTEMP(velUpdateDerivs, innerFaceBox, SpaceDim*SpaceDim);
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          // Note that for each stress, the full gradient is necessary
          // for the physical space transformation
          const int comp = srtIdx(row, col);
          FABSTACKTEMP(physGradFaceStress, innerFaceBox, SpaceDim);
          FABSTACKTEMP(gradFaceStress, innerFaceBox, SpaceDim);
          for (int deriv = 0; deriv != SpaceDim; ++deriv)
            {
              // Compute the face-normal derivative
              if (deriv == a_dir)
                {
                  // Evaluate the difference between the current face and
                  // the wall face. Then, using the computational space
                  // length between the current face and the wall, compute
                  // the derivative

                  // First, shift to the wall
                  gradFaceStress.shift(
                    a_dir, normSign*CRDparam::g_wallModelCrsRatio);
                  combinedStressTensor.shift(
                    a_dir, normSign*CRDparam::g_wallModelCrsRatio);
                  const Real numC = CRDparam::g_wallModelCrsRatio;
                  MD_BOXLOOP(a_boundaryFaceBox, i)
                    {
                      const Real interiorFaceVal =
                        combinedStressTensor[MD_IX(i, comp)];
                      const Real wallFaceVal = -tempState[MD_IX(i, comp)];
                      gradFaceStress[MD_IX(i, deriv)] =
                        offsetSign*numC*invDx*(interiorFaceVal - wallFaceVal);
                    }
                  // Finally, shift away from the wall
                  gradFaceStress.shift(
                    a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);
                  combinedStressTensor.shift(
                    a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);
                }
              else // Compute the tangential derivative
                {
                  Box copyFaceBox = innerFaceBox;
                  copyFaceBox.grow(deriv, 1);
                  copyFaceBox &= a_domain;
                  FABSTACKTEMP(faceStress, copyFaceBox, 1);
                  faceStress.copy(
                    combinedStressTensor,copyFaceBox,comp,copyFaceBox,0,1);
                  faceStress.shiftHalf(a_dir, -1);
                  gradFaceStress.shiftHalf(a_dir, -1);
                  PatchMappedFunc::cellAvgDerivFromCellAvgCS(
                    gradFaceStress, faceStress, 0, deriv, a_domain,
                    adjacentCellBox, a_gridMetrics.dxVect(), deriv, false,
                    false);
                  // Now make sure to shift things back as necessary
                  gradFaceStress.shiftHalf(a_dir, 1);
                }
            }
          // Transform the gradient into physical space
          for (int indx = 0; indx != SpaceDim; ++indx)
            {
              MD_BOXLOOP(innerFaceBox, i)
                {
                  physGradFaceStress[MD_IX(i,indx)] = D_TERM(
                    gradFaceStress[MD_IX(i,0)]*NtJFab[MD_IX(i,srtIdx(0,indx))],
                  + gradFaceStress[MD_IX(i,1)]*NtJFab[MD_IX(i,srtIdx(1,indx))],
                  + gradFaceStress[MD_IX(i,2)]*NtJFab[MD_IX(i,srtIdx(2,indx))]);
                }
            }

          // Place the necessary component in the data container
          MD_BOXLOOP(innerFaceBox, i)
            {
              velUpdateDerivs[MD_IX(i, comp)] =
                physGradFaceStress[MD_IX(i, col)];
            }
        }
    }

  FABSTACKTEMP(uDt, innerFaceBox, SpaceDim);
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(innerFaceBox, i)
        {
          uDt[MD_IX(i, comp)] = -physGradPres[MD_IX(i, comp)]
            - (D_TERM(velUpdateDerivs[MD_IX(i, srtIdx(comp, 0))],
                    + velUpdateDerivs[MD_IX(i, srtIdx(comp, 1))],
                    + velUpdateDerivs[MD_IX(i, srtIdx(comp, 2))]));
        }
    }
  // Shift to the wall
  uDt.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform uDt using unitNormals
  FORT_FORWARDTRANSFORMF(CHF_FRA(uDt),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Shift to the first interior cell
  uDt.shiftHalf(a_dir, offsetSign);
  // Shift interior face velocity to the first interior cell
  interiorFaceVel.shiftHalf(a_dir, offsetSign);
  // Bound eta_0 update with reasonable values
  // We assume that a local solution to the law-of-the-wall provides
  // a reasonable estimate for the instantaneous, local eta value
  FABSTACKTEMP(etaEstFab, bndryBox, 1);
  etaEst(etaEstFab, a_WcellAvgFab, a_gridMetrics, a_unitNormals,
         bndryBox, a_disjointBox, a_dir, a_side);
  FABSTACKTEMP(etaBound, bndryBox, 2);
  etaBound.setVal(0.);
  MD_BOXLOOP(bndryBox, i)
    {
      etaBound[MD_IX(i, 0)] = 1.e-5; // Essentially adding a low bound of 0
      etaBound[MD_IX(i, 1)] = 100.*(etaEstFab[MD_IX(i, 0)]);
    }
  // Now compute the eta_0 update
  MD_BOXLOOP(a_box, i)
    {
      Real rho = a_WcellAvgFab[MD_IX(i, cRho)];
      RealVect vel(D_DECL(interiorFaceVel[MD_IX(i, 0)],
                          interiorFaceVel[MD_IX(i, 1)],
                          interiorFaceVel[MD_IX(i, 2)]));
      RealVect uDtVec(D_DECL(uDt[MD_IX(i,0)],uDt[MD_IX(i,1)],uDt[MD_IX(i,2)]));
      // Zero the normal components of both vel and uDtVec
      vel[a_dir] = 0.;
      uDtVec[a_dir] = 0.;
      // Compute the wall-tangential velocity magnitude
      const Real q = rho*(std::sqrt(vel.radSquared()) + eps);
      // Now we essentially normalize the velocity by q and divide by q again
      const RealVect normedVel = vel/(q*q);
      // Here's time-derivative of the wall-tangent velocity magnitude
      const Real q_timeDeriv = normedVel.dotProduct(uDtVec);
      // Get eta_0 in the first cell
      const Real eta_0 = eta[MD_IX(i, 0)];
      // We need some value that incorporates the ideas of outer flow influence
      // on the inner flow profile and how fast that happens. This should be
      // incorporated as some term that makes the etaUpdate smaller as the
      // region of integration (i.e. the near wall cell volume) increases.
      // For now, we're just going to pick something.
      Real damping_0 = 1.;
      // Finally, the eta_0 update
      Real etaUpdate = damping_0*2.*eta_0*q_timeDeriv;
      // Get the current IntVect
      IntVect currVect = MD_GETIV(i);
      // Create a box of wall-normal cells to fill with etaUpdate
      Box wallNormalBox = a_box;
      if (a_side == Side::Lo)
        {
          // If bndry is low side, set low side of wallNormalBox to currVect
          wallNormalBox.setSmall(currVect);
          // Also, get the wall-normal component of a_disjointBox high end
          int highSide = a_disjointBox.bigEnd()[a_dir];
          // Set the wall-normal component of currVect to highSide
          currVect.setVal(a_dir, highSide);
          // Set the high side of wallNormalBox to the new currVect
          wallNormalBox.setBig(currVect);
        }
      else
        {
          // If bndry is high, set high side of wallNormalBox to currVect
          wallNormalBox.setBig(currVect);
          // Also, get the wall-normal component of a_disjointBox low end
          int lowSide = a_disjointBox.smallEnd()[a_dir];
          // Set the wall-normal component of currVect to lowSide
          currVect.setVal(a_dir, lowSide);
          // Set the low side of wallNormalBox to the new currVect
          wallNormalBox.setSmall(currVect);
        }
      // We now enforce eta boundedness
      Real wallLocalJ = cellJ[MD_IX(i, 0)];
      const Real etaMaxBound = etaBound[MD_IX(i, 1)];
      const Real etaMinBound = etaBound[MD_IX(i, 0)];
      const Real etaUpdateMaxBound = wallLocalJ*(etaMaxBound - eta_0)/a_dt;
      const Real etaUpdateMinBound = wallLocalJ*(etaMinBound - eta_0)/a_dt;
      etaUpdate =
        std::max(etaUpdateMinBound, std::min(etaUpdate, etaUpdateMaxBound));
      // Now loop over wallNormalBox and assign etaUpdate to all cells
      MD_BOXLOOP(wallNormalBox, j)
        {
          // Real localJ = cellJ[MD_IX(j, 0)];
          a_turbSourceAvgFab[MD_IX(j, cEta)] = etaUpdate;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute source term for eta_0 in order to advance eta_0 in time
/** \param[out] a_turbSourceAvgFab
 *                      Cell-averaged turbulent source terms
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_facePntVelGradFab
 *                      The face-centered velocity gradient
 *  \param[in]  a_WfaceAvgFab
 *                      Face-averaged primitive state
 *  \param[in]  a_RHSfab
 *                      Cell-averaged conservative state update
 *  \param[in]  a_unitNormals
 *                      Face-averaged unit-normal basis
 *  \param[in]  a_domain
 *                      Domain of this block
 *  \param[in]  a_gridMetrics
 *                      Grid metrics on this level
 *  \param[in]  a_dataIndx
 *                      Index of current box in DisjointBoxLayout
 *  \param[in]  a_side  Current side of box to work on (low or high)
 *  \param[in]  a_dir   Direction of face on which to update eta_0
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_box   Cell-centered box over which eta_0 source-term
 *                      is computed
 *  \param[in] a_boundaryFaceBox
 *                      Box of boundary faces next to which eta_0
 *                      source term is computed
 *//*-----------------------------------------------------------------*/

void
SSV::updateEtaZeroLocal(FArrayBox&              a_turbSourceAvgFab,
                        FArrayBox&              a_WcellAvgFab,
                        const FArrayBox&        a_JUFab,
                        const FArrayBox&        a_facePntVelGradFab,
                        const FArrayBox&        a_WfaceAvgFab,
                        const FluxBox&          a_WfaceAvgFxb,
                        const FluxBox&          a_fluxFxb,
                        FluxBox&                a_stressFluxFxb,
                        const FArrayBox&        a_RHSfab,
                        const FArrayBox&        a_unitNormals,
                        const FluxBox&          a_unitNormalsFxb,
                        const ProblemDomain&    a_domain,
                        const LevelGridMetrics& a_gridMetrics,
                        const DataIndex&        a_dataIndx,
                        const Side::LoHiSide&   a_side,
                        const int               a_dir,
                        const Box&              a_disjointBox,
                        const Box&              a_box,
                        const Box&              a_boundaryFaceBox,
                        const Real              a_dt) const
{
  if (!m_useWallModel) return;
  CH_TIME("SSV::updateEtaZeroLocal");

  // Basically, we try to preserve the velocity dynamics as much as possible
  // and just look at a single cell which is not too close to the wall

  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int turbCompBegin = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  int etaIndx = turbCompBegin + a_dir + 1;
  const Real eps = 1e-50;

  Side::LoHiSide oppSide = flip(a_side);
  const int normSign = sign(a_side);
  const int offsetSign = sign(oppSide);

  // Create the interior-face-box
  Box interiorFaceBox = a_boundaryFaceBox;
  // Shift the box to the interior face
  interiorFaceBox.shift(a_dir, offsetSign*CRDparam::g_wallModelCrsRatio);
  // Get rho*u at the interior face and shift it to the wall
  FABSTACKTEMP(interiorFaceRhoU, interiorFaceBox, SpaceDim);
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(interiorFaceBox, i)
        {
          interiorFaceRhoU[MD_IX(i, comp)] = a_WfaceAvgFab[MD_IX(i, rhoIndx)]*
            a_WfaceAvgFab[MD_IX(i, WvelIndx + comp)];
        }
    }
  interiorFaceRhoU.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform the fab using the unit normals on the wall face
  FORT_FORWARDTRANSFORMF(CHF_FRA(interiorFaceRhoU),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Shift interiorFaceRhoU to the wall-adjacent cell
  interiorFaceRhoU.shiftHalf(a_dir, offsetSign);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  // ------------------------ //
  //     Do we need this?     //
  // ------------------------ //

  // We need to transform the velocity at the interior face
  FABSTACKTEMP(interiorFaceVel, interiorFaceBox, SpaceDim);
  // Copy the velocity from a_WfaceAvgFab into interiorFaceVel
  interiorFaceVel.copy(
    a_WfaceAvgFab, interiorFaceBox, WvelIndx, interiorFaceBox, 0, SpaceDim);
  // Shift the fab to the wall so we can transform it using unit normals
  interiorFaceVel.shift(a_dir, normSign*CRDparam::g_wallModelCrsRatio);
  // Transform the fab using the unit normals on the wall face
  FORT_FORWARDTRANSFORMF(CHF_FRA(interiorFaceVel),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  // Step 1: Create the wall-face momentum flux
  // Expand a_boundaryFaceBox to include ghost cells
  Box bndryBox = a_boundaryFaceBox;
  IntVect growVect = IntVect_unit;
  growVect[a_dir] = 0;
  bndryBox.grow(growVect);
  bndryBox &= a_domain;
  // Make a box to cover the first interior cells adjacent to the bndry face
  Box firstInteriorCellBox = bndryBox; // Set to cover the ghost cells
  firstInteriorCellBox.growDir(a_dir,oppSide,1); // Include first interior face
  firstInteriorCellBox.enclosedCells(a_dir); // Turn this into cell box
  // Store the velocity for shifting purposes (a_WcellAvgFab stays const)
  FABSTACKTEMP(velFab, firstInteriorCellBox, SpaceDim);
  velFab.copy(a_WcellAvgFab, WvelIndx, 0, SpaceDim); // Fill with a_WcellAvgFab
  velFab.shiftHalf(a_dir, normSign); // Shift velFab to the wall face
  // Store the eta_0 value for shifting purposes
  FABSTACKTEMP(eta, firstInteriorCellBox, 1);
  eta.copy(a_WcellAvgFab, etaIndx, 0, 1); // Fill with a_WcellAvgFab

  // -------------------------------------------------------------------- //

  // Make sure eta_0 is >= 0
  MD_BOXLOOP(firstInteriorCellBox, i)
    {
      eta[MD_IX(i, 0)] = std::max(0., eta[MD_IX(i, 0)]);
    }

  // Create box for first interior layer of wall-tangential faces
  Box firstInteriorFaceBox = bndryBox;
  // Shift the wall-faces box into the domain by 1
  firstInteriorFaceBox.shift(a_dir, offsetSign);

  // Create box for wall-faces + first layer of interior wall-tangential faces
  Box wallDistanceBox = bndryBox;
  // Grow the wall-faces box to include the first interior cell
  wallDistanceBox.growDir(a_dir, oppSide, 1);

  // We need the coordinate system for the block
  const BlockCoordSys& blockCoordSysCoords =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Compute physical distance vector btwn the wall and first interior faces
  FABSTACKTEMP(XiFab, wallDistanceBox, SpaceDim); // computational space coords
  FABSTACKTEMP(XFab, wallDistanceBox, SpaceDim); // physical space coordinates
  // Get the coordinates on the wall and first interior faces
  CRDparam::g_CNSIBC->getFaceCoordinates(
    wallDistanceBox, XiFab, XFab, a_dir, blockCoordSysCoords);
  // Compute difference between XFab on the top and bottom faces of the cell
  FABSTACKTEMP(wallDistanceVect, firstInteriorFaceBox, SpaceDim);
  const int MD_ID(o, a_dir);
  MD_BOXLOOP(firstInteriorFaceBox, i)
    {
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          Real xCompInterior = XFab[MD_IX(i, comp)];
          Real xCompWall = XFab[MD_OFFSETIX(i, +, normSign*o, comp)];
          // Compute (interior - wall)*(side compensation) (1 if lo, -1 if hi)
          wallDistanceVect[MD_IX(i, comp)] =
            offsetSign*(xCompInterior - xCompWall);
        }
    }
  // Shift the wallDistanceVect fab to the wall
  wallDistanceVect.shift(a_dir, normSign);
  // Transform the wall-distance vector into normal-tangent space
  FORT_FORWARDTRANSFORMF(CHF_FRA(wallDistanceVect),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(bndryBox));
  // Shift the wallDistanceVect fab to the first interior cell, not face
  wallDistanceVect.shiftHalf(a_dir, offsetSign);

  // Assuming that the wall is at a raised virtual wall, the virtual eta
  // should be different from the actual eta. Compute the virtual eta and
  // test the eta_0 update with this
  FABSTACKTEMP(etaVirtual, firstInteriorCellBox, 1);
  const Real mu = CRDparam::g_mu; //**NOTE: This won't work with multispecies
  const Real virtualWallHeight = m_virtualWallHeight;
  MD_BOXLOOP(firstInteriorCellBox, i)
    {
      const Real k_0 = 0.41;
      const Real C_0 = 0.001093; // From Musker 1979
      Real nu = mu/(a_WcellAvgFab[MD_IX(i, rhoIndx)]);
      // Compute h_0 -- use the normal component of wallDistanceVect
      const Real dNormal = wallDistanceVect[MD_IX(i, a_dir)];
      // Compute u_tau
      const Real uTau = std::sqrt(nu*eta[MD_IX(i, 0)]);
      // Compute y_plus at the virtual wall
      const Real y_p = (virtualWallHeight*dNormal)*uTau/nu;
      const Real y_p_sq = y_p*y_p;
      // Compute the fraction of eta_0 that eta is at the virtual wall
      const Real eta_fraction =
        (k_0 + C_0*y_p_sq)/(k_0 + C_0*y_p_sq + C_0*k_0*y_p*y_p_sq);
      etaVirtual[MD_IX(i, 0)] = eta_fraction*eta[MD_IX(i, 0)];
    }

  eta.shiftHalf(a_dir, normSign); // Shift eta to the wall face
  etaVirtual.shiftHalf(a_dir, normSign); // Shift etaVirtual to the wall face

  // ----------------------------------------------------------------------- //

  // Transform the velocity vector
  const int transformIndx = 0;
  FABSTACKTEMP(locStreamLocTanVect, bndryBox, SpaceDim);
  FORT_FORWARDTRANSFORMGENF(CHF_FRA(locStreamLocTanVect),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(velFab),
                            CHF_INT(transformIndx),
                            CHF_CONST_FRA(a_unitNormals),
                            CHF_BOX(bndryBox));
  // Set the normal component to zero to obtain the wall-tangent velocity vector
  locStreamLocTanVect.setVal(0.0, a_dir);
  // Normalize the wall-tangent vector -- this is the streamwise unit vector
  PatchMappedFunc::normalize(bndryBox,
                             locStreamLocTanVect,
                             Interval(0, SpaceDim-1));

  // Create a normal vector in normal-tangent space
  FABSTACKTEMP(normVect, bndryBox, SpaceDim);
  normVect.setVal(0.);
  normVect.setVal(1., a_dir);
  // Transform normalVect into Cartesian physical space
  FORT_REVERSETRANSFORMF(CHF_FRA(normVect),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(bndryBox));
  // Compute the wall-normal derivative of the wall-normal velocity
  FABSTACKTEMP(modVelGrad, bndryBox, SpaceDim*SpaceDim);
  modVelGrad.setVal(0.); // Most values need to be zero
  MD_BOXLOOP(bndryBox, i)
    {
      Real normDeriv = 0.;
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          normDeriv += normVect[MD_IX(i, comp)]*(
            D_TERM(
              a_facePntVelGradFab[MD_IX(i,srtIdx(comp,0))]*
                normVect[MD_IX(i,0)],
            + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,1))]*
                normVect[MD_IX(i,1)],
            + a_facePntVelGradFab[MD_IX(i,srtIdx(comp,2))]*
                normVect[MD_IX(i,2)]));
        }
      modVelGrad[MD_IX(i, srtIdx(a_dir,a_dir))] = normDeriv;
    }
  // Fill the normal-derivs of the tangential velocities using model
  for (int velComp = 0; velComp != SpaceDim; ++velComp)
    {
      if (velComp != a_dir)
        {
          const int cGrad = a_dir + SpaceDim*velComp;
          MD_BOXLOOP(bndryBox, i)
            {
              modVelGrad[MD_IX(i, cGrad)] =
                locStreamLocTanVect[MD_IX(i,velComp)]*eta[MD_IX(i,0)];
            }
        }
    }
  // Transform back to Cartesian space
  FABSTACKTEMP(tempState, bndryBox, SpaceDim*SpaceDim);
  // Right multiply grad(u) by a_unitNormals
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int cGrad = col + SpaceDim*row;
          MD_BOXLOOP(bndryBox, i)
            {
              tempState[MD_IX(i, cGrad)] =
                D_TERM(
                  modVelGrad[MD_IX(i,srtIdx(row,0))]*
                    a_unitNormals[MD_IX(i,srtIdx(0,col))],
                  + modVelGrad[MD_IX(i,srtIdx(row,1))]*
                    a_unitNormals[MD_IX(i,srtIdx(1,col))],
                  + modVelGrad[MD_IX(i,srtIdx(row,2))]*
                    a_unitNormals[MD_IX(i,srtIdx(2,col))]);
            }
        }
    }
  // Left multiply temp state by a_unitNormals^T
  Interval tensorIntv(0, SpaceDim*SpaceDim - 1);
  FArrayBox cartGradVel(tensorIntv, modVelGrad); // Alias to avoid allocation
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int cGrad = col + SpaceDim*row;
          MD_BOXLOOP(bndryBox, i)
            {
              cartGradVel[MD_IX(i, cGrad)] =
                D_TERM(
                  a_unitNormals[MD_IX(i,srtIdx(0,row))]*
                  tempState[MD_IX(i,srtIdx(0,col))],
                + a_unitNormals[MD_IX(i,srtIdx(1,row))]*
                  tempState[MD_IX(i,srtIdx(1,col))],
                + a_unitNormals[MD_IX(i,srtIdx(2,row))]*
                  tempState[MD_IX(i,srtIdx(2,col))]);
            }
        }
    }
  // Calculate the viscosity (and thermal conductivity unfortunately)
  FABSTACKTEMP(muFab, bndryBox, 1);
  FABSTACKTEMP(kappaFab, bndryBox, 1);
  CRDparam::g_CRDPhysics->calcCoeffKappaMu(
    bndryBox, muFab, kappaFab, a_WfaceAvgFab);
  // Fill the new viscous stress tensor with the strain-rate tensor
  FArrayBox vstFab(tensorIntv, tempState); // Alias to avoid allocation
  for (int row = 0; row != SpaceDim; ++row)
    {
      for (int col = 0; col != SpaceDim; ++col)
        {
          const int vstComp = row*SpaceDim + col;
          MD_BOXLOOP(bndryBox, i)
            {
              vstFab[MD_IX(i, vstComp)] = 
                cartGradVel[MD_IX(i, srtIdx(row, col))]
              + cartGradVel[MD_IX(i, srtIdx(col, row))];
            }
        }
    }
  // Now subtract the divergence of velocity
  MD_BOXLOOP(bndryBox, i)
    {
      Real divVel = 2.*CRDparam::g_lambda*(D_TERM(
        cartGradVel[MD_IX(i, srtIdx(0,0))],
      + cartGradVel[MD_IX(i, srtIdx(1,1))],
      + cartGradVel[MD_IX(i, srtIdx(2,2))]));
      for (int row = 0; row != SpaceDim; ++row)
        {
          vstFab[MD_IX(i, srtIdx(row,row))] -= divVel;
        }
    }
  // Now multiply by the viscosity
  for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
    {
      MD_BOXLOOP(bndryBox, i)
        {
          vstFab[MD_IX(i, comp)] *= muFab[MD_IX(i, 0)];
        }
    }
  // Store vstFab in vstDiffFab for difference calculation later
  FABSTACKTEMP(vstDiffFab, bndryBox, SpaceDim*SpaceDim);
  vstDiffFab.copy(vstFab);
  // Now convolve the new viscous stress tensor
  FArrayBox vstAvgFab(tensorIntv, cartGradVel); // Alias to avoid allocation
  vstAvgFab.copy(vstFab);
  MOLUtilFunc::deconvolveCenterFace(
    vstAvgFab,vstFab,a_boundaryFaceBox,a_domain,a_dir);
  // DEBUG ONLY
  // // Compute the original viscous stress tensor for comparison
  // FArrayBox origVSTFab(tensorIntv, vstFab); // Alias to avoid allocation
  // for (int row = 0; row != SpaceDim; ++row)
  //   {
  //     for (int col = 0; col != SpaceDim; ++col)
  //       {
  //         const int vstComp = row*SpaceDim + col;
  //         MD_BOXLOOP(bndryBox, i)
  //           {
  //             origVSTFab[MD_IX(i, vstComp)] = 
  //               a_facePntVelGradFab[MD_IX(i, srtIdx(row, col))]
  //             + a_facePntVelGradFab[MD_IX(i, srtIdx(col, row))];
  //           }
  //       }
  //   }
  // // Now subtract the divergence of velocity
  // MD_BOXLOOP(bndryBox, i)
  //   {
  //     Real divVel = 2.*CRDparam::g_lambda*(D_TERM(
  //       a_facePntVelGradFab[MD_IX(i, srtIdx(0,0))],
  //     + a_facePntVelGradFab[MD_IX(i, srtIdx(1,1))],
  //     + a_facePntVelGradFab[MD_IX(i, srtIdx(2,2))]));
  //     for (int row = 0; row != SpaceDim; ++row)
  //       {
  //         origVSTFab[MD_IX(i, srtIdx(row,row))] -= divVel;
  //       }
  //   }
  // // Now multiply by the viscosity
  // for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
  //   {
  //     MD_BOXLOOP(bndryBox, i)
  //       {
  //         origVSTFab[MD_IX(i, comp)] *= muFab[MD_IX(i, 0)];
  //       }
  //   }
  // END DEBUG ONLY
  // DEBUG ONLY -- actually, this is likely to be always used soon
  // FArrayBox& viscStress = a_stressFluxFxb[a_dir];
  // for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
  //   {
  //     MD_BOXLOOP(bndryBox, i)
  //       {
  //         vstDiffFab[MD_IX(i, comp)] -= viscStress[MD_IX(i, comp)];
  //       }
  //   }
  // END DEBUG ONLY
  // DEBUG ONLY
  // // Subtract the two face-centered viscous stress tensors
  // for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
  //   {
  //     MD_BOXLOOP(bndryBox, i)
  //       {
  //         vstDiffFab[MD_IX(i, comp)] -= origVSTFab[MD_IX(i, comp)];
  //       }
  //   }
  // END DEBUG ONLY
  // Store vstAvgFab in vstDiffAvgFab for difference calculation later
  // FABSTACKTEMP(vstDiffAvgFab, bndryBox, SpaceDim*SpaceDim);
  // vstDiffAvgFab.copy(vstAvgFab);
  // DEBUG ONLY
  // // Now convolve the original viscous stress tensor
  // FArrayBox origVSTAvgFab(tensorIntv, vstAvgFab); // Avoiding allocation
  // origVSTAvgFab.copy(origVSTFab);
  // MOLUtilFunc::deconvolveCenterFace(
  //   origVSTAvgFab, origVSTFab, a_boundaryFaceBox, a_domain, a_dir);
  // END DEBUG ONLY
  // DEBUG ONLY
  // Now convolve the original viscous stress tensor
  // FArrayBox origVSTAvgFab(tensorIntv, vstAvgFab); // Avoid allocation
  // origVSTAvgFab.copy(viscStress);
  // MOLUtilFunc::deconvolveCenterFace(
  //   origVSTAvgFab, viscStress, a_boundaryFaceBox, a_domain, a_dir);
  // END DEBUG ONLY
  // Subtract the two viscous stress tensors
  // for (int comp = 0; comp != SpaceDim*SpaceDim; ++comp)
  //   {
  //     MD_BOXLOOP(a_boundaryFaceBox, i)
  //       {
  //         vstDiffAvgFab[MD_IX(i, comp)] -= origVSTAvgFab[MD_IX(i, comp)];
  //       }
  //   }
  // Map the viscous stress tensors
  Box mappedBox = a_boundaryFaceBox;
  // Enclose all faces up to the evaluation face
  mappedBox.growDir(a_dir, oppSide, CRDparam::g_wallModelCrsRatio);
  mappedBox.enclosedCells(a_dir); // Turn this into a cell-box
  Box dataBox = mappedBox;
  dataBox.grow(1); // Grow all directions by 1
  const FluxBox& Nfxb = a_gridMetrics.m_N[a_dataIndx];
  FLUXBOXSTACKTEMP(mappedFlux, mappedBox, SpaceDim);
  FLUXBOXSTACKTEMP(pntFlux, dataBox, SpaceDim*SpaceDim);
  FLUXBOXSTACKTEMP(avgFlux, dataBox, SpaceDim*SpaceDim);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      mappedFlux[dir].setVal(0.);
      pntFlux[dir].setVal(0.);
      avgFlux[dir].setVal(0.);
      if (dir == a_dir)
        {
          pntFlux[dir].copy(vstFab);
          avgFlux[dir].copy(vstAvgFab);
        }
    }
  Interval fluxIntv(0, SpaceDim - 1);
  bool fourthOrder = false;
  if (CRDparam::g_faceConvolveFlatten < 2 ||
      CRDparam::g_faceDeconvolveFlatten < 2)
    {
      fourthOrder = true;
    }
  const BlockCoordSys* blockCoordSys = a_gridMetrics.getCoordSys(a_disjointBox);
  blockCoordSys->computeMetricTermProductAverage(
    mappedFlux, avgFlux, Nfxb, SpaceDim, pntFlux, mappedBox, fourthOrder,
    fluxIntv, fluxIntv, 1, &a_domain);

  // Step 2: Evaluate the wall-normal and wall-tangential RHS(momentum)
  // RHS(momentum) comes from difference of face fluxes divided by dx
  Box rhsBox = a_box;
  rhsBox.growDir(a_dir, oppSide, (CRDparam::g_wallModelCrsRatio-1));
  // This holds RHS contribution derived from wall-parallel faces
  FABSTACKTEMP(wallNormalRHS, rhsBox, SpaceDim);
  wallNormalRHS.setVal(0.);
  // This holds RHS contribution derived from wall-perpendicular faces
  FABSTACKTEMP(wallTanRHS, rhsBox, SpaceDim);
  wallTanRHS.setVal(0.);
  // Create the wall-normal fluxes
  FABSTACKTEMP(vertFluxFab, mappedFlux[a_dir].box(), a_fluxFxb[a_dir].nComp());
  vertFluxFab.copy(a_fluxFxb[a_dir]);
  // Zero the wall-normal momentum fluxes right on the wall face
  vertFluxFab.setVal(0., a_boundaryFaceBox, WvelIndx, SpaceDim);
  // Set the wall-normal momentum fluxes on the wall face
  vertFluxFab.minus(mappedFlux[a_dir], 0, WvelIndx, SpaceDim);
  // Add the pressure flux at the wall? -- but it gets zeroed out no matter what
  // anyway, so it doesn't have to be added
  // For each velocity component, loop over each face direction
  for (int velComp = 0; velComp != SpaceDim; ++velComp)
    {
      // For each face dir, compute the contribution of the faces to RHS
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          const FArrayBox& fluxFab = a_fluxFxb[faceDir];
          // const FArrayBox& mappedFluxFab = mappedFlux[faceDir];
          const int MD_ID(ii, faceDir);
          // Remember that RHS is -(1/dx)*(F_R - F_L) -- so, we build in the
          // negative sign into the invDx term to make things easier
          const Real invDx = -1./(a_gridMetrics.dxVect()[faceDir]);
          // Also, we account for the virtual raised wall here
          const Real virtualWallHeight = m_virtualWallHeight;
          const Real w_0 = 1. + virtualWallHeight;
          const Real invWallDx = -1./(w_0*(a_gridMetrics.dxVect()[faceDir]));
          // If faces are perpendicular to wall, sum both face contributions
          if (faceDir != a_dir)
            {
              MD_BOXLOOP(rhsBox, i)
                {
                  const Real hiFlux =
                    fluxFab[MD_OFFSETIX(i, +, ii, WvelIndx+velComp)];
                  const Real loFlux = fluxFab[MD_IX(i, WvelIndx+velComp)];
                  wallTanRHS[MD_IX(i, velComp)] += invDx*(hiFlux - loFlux);
                }
            }
          else // If faces are wall-parallel, assign the face contributions
            {
              //**NOTE: FluxFix should be subtracted in this case
              MD_BOXLOOP(rhsBox, i)
                {
                  // DEBUG ONLY
                  // const Real loFluxFix = mappedFluxFab[MD_IX(i, velComp)];
                  // const Real hiFluxFix = 
                  //   mappedFluxFab[MD_OFFSETIX(i, +, ii, velComp)];
                  // const Real loFlux = fluxFab[MD_IX(i, WvelIndx+velComp)]
                  //   - loFluxFix;
                  // const Real hiFlux =
                  //   fluxFab[MD_OFFSETIX(i, +, ii, WvelIndx+velComp)]
                  //   - hiFluxFix;
                  // END DEBUG ONLY
                  //**FIXME: when we switch back to coarsening, the following
                  //         will be wrong
                  const Real loFlux = vertFluxFab[MD_IX(i, WvelIndx+velComp)];
                  const Real hiFlux =
                    vertFluxFab[MD_OFFSETIX(i, +, ii, WvelIndx+velComp)];
                  Real localInvDx = invDx;
                  const int wallAdjCellIndx =
                    firstInteriorCellBox.smallEnd()[a_dir];
                  if (MD_GETIV(i)[a_dir] == wallAdjCellIndx)
                    {
                      localInvDx = invWallDx;
                    }
                  wallNormalRHS[MD_IX(i, velComp)] =
                    localInvDx*(hiFlux - loFlux);
                }
            }
        }
    }

  // Step 3: Combine the two

  // First, combine the vertical fluxes
  const FArrayBox& cellJ = a_gridMetrics.m_J[a_dataIndx];
  FABSTACKTEMP(vertFlux, a_box, SpaceDim);
  FABSTACKTEMP(sumCellJ, a_box, 1);
  vertFlux.setVal(0.);
  sumCellJ.setVal(0.);
  for (int cell = 0; cell != CRDparam::g_wallModelCrsRatio; ++cell)
    {
      // Get the wall-parallel interior-cell box
      Box cellBox = a_box;
      cellBox.shift(a_dir, offsetSign*cell);
      // Shift the summed vertFlux to the cell of interest
      vertFlux.shift(a_dir, offsetSign*cell);
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          MD_BOXLOOP(cellBox, i)
            {
              vertFlux[MD_IX(i, comp)] += wallNormalRHS[MD_IX(i, comp)];
            }
        }
      // Shift the summed vertFlux back to the wall-adjacent cell
      vertFlux.shift(a_dir, normSign*cell);
      // Shift the summed sumCellJ to the cell of interest
      sumCellJ.shift(a_dir, offsetSign*cell);
      MD_BOXLOOP(cellBox, i)
        {
          sumCellJ[MD_IX(i, 0)] += cellJ[MD_IX(i, 0)];
        }
      // Shift the summed sumCellJ back to the wall-adjacent cell
      sumCellJ.shift(a_dir, normSign*cell);
    }
  // Now, weight the vertical fluxes
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(a_box, i)
        {
          const Real sumJ = sumCellJ[MD_IX(i,0)];
          vertFlux[MD_IX(i,comp)] = (vertFlux[MD_IX(i,comp)])/sumJ;
        }
    }
  // Get the box to work over
  Box cellBoxInterior = a_box;
  cellBoxInterior.shift(a_dir, offsetSign*(CRDparam::g_wallModelCrsRatio-1));
  // If the coarsening uses wall-adjacent data for the wall-perpendicular
  // fluxes, we need to sum these and weight them as well
  const bool sumWallTanFlux = true;
  FABSTACKTEMP(tanFlux, a_box, SpaceDim);
  if (sumWallTanFlux)
    {
      tanFlux.setVal(0.);
      for (int cell = 0; cell != CRDparam::g_wallModelCrsRatio; ++cell)
        {
          // Get the wall-parallel interior-cell box
          Box cellBox = a_box;
          cellBox.shift(a_dir, offsetSign*cell);
          // Shift the summed tanFlux to the cell of interest
          tanFlux.shift(a_dir, offsetSign*cell);
          for (int comp = 0; comp != SpaceDim; ++comp)
            {
              MD_BOXLOOP(cellBox, i)
                {
                  tanFlux[MD_IX(i, comp)] += wallTanRHS[MD_IX(i, comp)];
                }
            }
          // Shift the summed tanFlux back to the wall-adjacent cell
          tanFlux.shift(a_dir, normSign*cell);
        }
      // Now, weight the vertical fluxes
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          MD_BOXLOOP(a_box, i)
            {
              const Real sumJ = sumCellJ[MD_IX(i,0)];
              tanFlux[MD_IX(i,comp)] = (tanFlux[MD_IX(i,comp)])/sumJ;
            }
        }
      // Shift tanFlux to the cell of interest away from the wall
      tanFlux.shift(a_dir, offsetSign*(CRDparam::g_wallModelCrsRatio-1));
    }
  else
    {
      // Shift tanFlux to the cell of interest away from the wall
      tanFlux.shift(a_dir, offsetSign*(CRDparam::g_wallModelCrsRatio-1));
      // Fill tanFlux with wallTanRHS
      tanFlux.copy(wallTanRHS, 0, 0, SpaceDim);
      // Weight tanFlux
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          MD_BOXLOOP(cellBoxInterior, i)
            {
              const Real invJval = 1./(cellJ[MD_IX(i, 0)]);
              tanFlux[MD_IX(i, comp)] = invJval*tanFlux[MD_IX(i, comp)];
            }
        }
    }
  FABSTACKTEMP(totalFlux, a_box, SpaceDim);
  // Shift totalFlux to the cell of interest away from the wall
  totalFlux.shift(a_dir, offsetSign*(CRDparam::g_wallModelCrsRatio-1));
  // Shift vertFlux to the cell of interest away from the wall
  vertFlux.shift(a_dir, offsetSign*(CRDparam::g_wallModelCrsRatio-1));
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(cellBoxInterior, i)
        {
          totalFlux[MD_IX(i, comp)] = vertFlux[MD_IX(i, comp)]
            + tanFlux[MD_IX(i, comp)];
        }
    }
  // Shift totalFlux back to the wall-adjacent cell
  totalFlux.shift(a_dir, normSign*(CRDparam::g_wallModelCrsRatio-1));
  // Multiply totalFlux by the wall-adjacent cell J value
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(a_box, i)
        {
          totalFlux[MD_IX(i, comp)] *= cellJ[MD_IX(i, 0)];
        }
    }

  // Step 4: Update eta_0

  // Bound eta_0 update with reasonable values
  // We assume that a local solution to the law-of-the-wall provides
  // a reasonable estimate for the instantaneous, local eta value
  FABSTACKTEMP(etaEstFab, bndryBox, 1);
  etaEst(etaEstFab, a_WcellAvgFab, a_gridMetrics, a_unitNormals,
         bndryBox, a_disjointBox, a_dir, a_side);
  FABSTACKTEMP(etaBound, bndryBox, 2);
  etaBound.setVal(0.);
  MD_BOXLOOP(bndryBox, i)
    {
      etaBound[MD_IX(i, 0)] = 1.e-5; // Essentially adding a low bound of 0
      etaBound[MD_IX(i, 1)] = std::min(2.*etaEstFab[MD_IX(i,0)], m_etaMaxBound);
    }

  // DEBUG ONLY
  FABSTACKTEMP(flowMomRHS, a_box, SpaceDim);
  // Fill the fab with the RHS data for momentum
  flowMomRHS.copy(a_RHSfab, a_box, WvelIndx, a_box, 0, SpaceDim);
  // END DEBUG ONLY

  // DEBUG ONLY
  FABSTACKTEMP(compareTotalRHS, rhsBox, SpaceDim);
  for (int comp = 0; comp != SpaceDim; ++comp)
    {
      MD_BOXLOOP(rhsBox, i)
        {
          compareTotalRHS[MD_IX(i, comp)] = wallTanRHS[MD_IX(i, comp)]
            + wallNormalRHS[MD_IX(i, comp)];
        }
    }
  // END DEBUG ONLY

  // Shift the fab to the wall so we can transform it using unit normals
  totalFlux.shiftHalf(a_dir, normSign);
  // Transform the fab using the unit normals on the wall face
  FORT_FORWARDTRANSFORMF(CHF_FRA(totalFlux),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_boundaryFaceBox));
  // Shift the fab back to the first interior cell
  totalFlux.shiftHalf(a_dir, offsetSign);
  // Compute the eta_0 update
  MD_BOXLOOP(a_box, i)
    {
      const Real eta_0 = eta[MD_IX(i, 0)];
      // const Real eta_0_virtual = etaVirtual[MD_IX(i, 0)];
      // Momentum at evaluation face in the interior
      RealVect rhoU(D_DECL(interiorFaceRhoU[MD_IX(i, 0)],
                           interiorFaceRhoU[MD_IX(i, 1)],
                           interiorFaceRhoU[MD_IX(i, 2)]));
      // RHS(momentum)
      RealVect momRHS(D_DECL(totalFlux[MD_IX(i, 0)],
                             totalFlux[MD_IX(i, 1)],
                             totalFlux[MD_IX(i, 2)]));
      // Now, we zero the normal component of momentum
      rhoU[a_dir] = 0.;
      // Similarly, we zero the normal component of momRHS
      momRHS[a_dir] = 0.;
      // Compute the wall-tangential velocity magnitude
      const Real q = std::sqrt(rhoU.radSquared()) + eps;
      // Now we essentially normalize the velocity by q
      const RealVect normedVel = rhoU/(q*q);
      // Here's time-derivative of the wall-tangent velocity magnitude
      const Real q_timeDeriv = normedVel.dotProduct(momRHS);
      // The eta_0 update
      // const Real alpha_0 = 1.0; // DEBUG ONLY -- representing time filtering
      Real etaUpdate = 2.*eta_0*q_timeDeriv;

      // Create a box of wall-normal cells to fill with etaUpdate
      Box wallNormalBox = a_box;
      // Get the current IntVect
      IntVect currVect = MD_GETIV(i);
      if (a_side == Side::Lo)
        {
          // If bndry is low side, set low side of wallNormalBox to currVect
          wallNormalBox.setSmall(currVect);
          // Also, get the wall-normal component of a_disjointBox high end
          int highSide = a_disjointBox.bigEnd()[a_dir];
          // Set the wall-normal component of currVect to highSide
          currVect.setVal(a_dir, highSide);
          // Set the high side of wallNormalBox to the new currVect
          wallNormalBox.setBig(currVect);
        }
      else
        {
          // If bndry is high, set high side of wallNormalBox to currVect
          wallNormalBox.setBig(currVect);
          // Also, get the wall-normal component of a_disjointBox low end
          int lowSide = a_disjointBox.smallEnd()[a_dir];
          // Set the wall-normal component of currVect to lowSide
          currVect.setVal(a_dir, lowSide);
          // Set the low side of wallNormalBox to the new currVect
          wallNormalBox.setSmall(currVect);
        }
      Real localJ = cellJ[MD_IX(i, 0)];
      // We now enforce eta boundedness
      const Real etaMaxBound = etaBound[MD_IX(i, 1)];
      const Real etaMinBound = etaBound[MD_IX(i, 0)];
      const Real etaUpdateMaxBound = localJ*(etaMaxBound - eta_0)/a_dt;
      const Real etaUpdateMinBound = localJ*(etaMinBound - eta_0)/a_dt;
      etaUpdate =
        std::max(etaUpdateMinBound, std::min(etaUpdate, etaUpdateMaxBound));
      // Now loop over wallNormalBox and assign etaUpdate to all cells
      MD_BOXLOOP(wallNormalBox, j)
        {
          a_turbSourceAvgFab[MD_IX(j, etaIndx)] = etaUpdate;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the inertial flux from primitive variable values on a face
/** \param[out] a_flux  Flux on the faces
 *  \param[in]  a_WFace Primitive state on the face
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_box   Box on which to compute flux
 *//*-----------------------------------------------------------------*/

void
SSV::addTurbConvectiveFlux(FArrayBox&       a_flux,
                           const FArrayBox& a_WFace,
                           const int&       a_dir,
                           const Box&       a_box) const
{
  CH_TIME("SSV::addTurbConvectiveFlux");
  CRD::msg << CRD::fv4 << "SSV::addTurbConvectiveFlux" << CRD::end;
  const Interval turbIntv = CRDparam::g_CRDPhysics->turbConsInterval();
  const int turbBeginIndx = turbIntv.begin();
  const int turbEndIndx = turbBeginIndx + turbIntv.size();

  for (int comp = turbBeginIndx; comp != turbEndIndx; ++comp)
    {
      MD_BOXLOOP(a_box, i)
        {
          a_flux[MD_IX(i, comp)] = 0.;
        }
    }
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
SSV::calcTurbSourceTerms(const Box&              a_box,
                         const Box&              a_disjointBox,
                         FArrayBox&              a_turbSourcePntFab,
                         const FArrayBox&        a_GradWcellPntFab,
                         const FArrayBox&        a_WcellPntFab,
                         const DataIndex&        a_dataIndx,
                         const LevelGridMetrics& a_gridMetrics) const
{
  CH_TIME("SSV::calcTurbSourceTerms");
  CRD::msg << CRD::fv4 << "SSV::calcTurbSourceTerms" << CRD::end;
  const Interval turbIntv = CRDparam::g_CRDPhysics->turbConsInterval();
  // Set all turbulent source terms to zero
  a_turbSourcePntFab.setVal(0., a_box, turbIntv.begin(), turbIntv.size());
}

/*--------------------------------------------------------------------*/
//  Initialize the turbulent variables in the flow field
/** \param[out] a_U     Cell-centered conservative variables
 *  \param[in]  a_W     Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_didx  Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                      Current disjointBox
 *  \param[in]  a_box   Box to initialize over
 *//*-----------------------------------------------------------------*/

void
SSV::turbInitialize(FArrayBox&              a_U,
                    const FArrayBox&        a_W,
                    const LevelGridMetrics& a_gridMetrics,
                    const FluxBox&          a_unitNormals,
                    const DataIndex&        a_didx,
                    const Box&              a_disjointBox,
                    const Box&              a_box) const
{
  if (!m_useWallModel || m_suppressWallModel) return;
  CH_TIME("SSV::turbInitialize");

  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  // We need to set the values for a_U in the eta components to zero
  const int cTurb = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  for (int d = 0; d != SpaceDim; ++d) { a_U.setVal(0., a_box, (cTurb+d+1)); }

  // Initialize the eta_0 component for the wall-model
  IntVect boundedBox = IntVect_zero;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      const int cEta = cTurb + dir + 1;
      // We need to loop over every possible boundary face direction
      Box boundaryFaceBoxLo;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxLo, a_box, blockDomain, dir, Side::Lo);
      Box boundaryFaceBoxHi;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxHi, a_box, blockDomain, dir, Side::Hi);

      // Check if a single box has both low and high wall boundaries
      if (!boundaryFaceBoxLo.isEmpty() && !boundaryFaceBoxHi.isEmpty())
        {
          // Check low boundary condition
          BoundaryIndex bcIdxLo;
          bcIdxLo.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
                         dir, Side::Lo);
          BCInfo domBCLo = CRDparam::g_CNSIBC->getDomainBC(bcIdxLo);
          // Check high boundary condition
          BoundaryIndex bcIdxHi;
          bcIdxHi.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
                         dir, Side::Hi);
          BCInfo domBCHi = CRDparam::g_CNSIBC->getDomainBC(bcIdxHi);
          // If both boundaries are walls, we need to restrict the
          // initialization to half the box or less
          if ((CRDparam::DomainBCTypeSWall & domBCHi.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBCHi.m_type) &&
              (CRDparam::DomainBCTypeSWall & domBCLo.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBCLo.m_type))
            {
              boundedBox[dir] = 1;
            }
        }

      if (!boundaryFaceBoxLo.isEmpty())
        {
          // Check boundary condition
          BoundaryIndex bcIdx;
          bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
                       dir, Side::Lo);
          BCInfo domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
          // We only support the use of walls in this function
          if ((CRDparam::DomainBCTypeSWall & domBC.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBC.m_type))
            {
              FABSTACKTEMP(etaFab, boundaryFaceBoxLo, 1);
              etaEst(etaFab, a_W, a_gridMetrics, a_unitNormals[dir],
                     boundaryFaceBoxLo, a_disjointBox, dir, Side::Lo);
              etaFab.shiftHalf(dir, 1); // Shift to first interior cell
              // Set the number of cells to be filled with etaFab
              int numCellsIter = a_disjointBox.size(dir);
              if (boundedBox[dir]) // If hi + lo walls, fill half the box here
                {
                  numCellsIter = 0.5*a_disjointBox.size(dir);
                }
              Box etaBox = boundaryFaceBoxLo; // Create box to fill with etaFab
              etaBox.shiftHalf(dir, 1); // First interior cell box
              // Now fill the cells with etaFab
              for (int cell = 0; cell != numCellsIter; ++cell)
                {
                  MD_BOXLOOP(etaBox, i)
                    {
                      a_U[MD_IX(i, cEta)] = etaFab[MD_IX(i, 0)];
                    }
                  etaBox.shift(dir, 1); // Shift box to next layer of cells
                  etaFab.shift(dir, 1); // Shift fab to next layer of cells
                }
            }
        }

      if (!boundaryFaceBoxHi.isEmpty())
        {
          // Check boundary condition
          BoundaryIndex bcIdx;
          bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
                       dir, Side::Hi);
          BCInfo domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
          // We only support the use of walls in this function
          if ((CRDparam::DomainBCTypeSWall & domBC.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBC.m_type))
            {
              FABSTACKTEMP(etaFab, boundaryFaceBoxHi, 1);
              etaEst(etaFab, a_W, a_gridMetrics, a_unitNormals[dir],
                     boundaryFaceBoxHi, a_disjointBox, dir, Side::Hi);
              etaFab.shiftHalf(dir, -1); // Shift to first interior cell
              // Set the number of cells to be filled with etaFab
              int numCellsIter = a_disjointBox.size(dir);
              if (boundedBox[dir]) // If hi + lo walls, fill half the box here
                {
                  numCellsIter = 0.5*a_disjointBox.size(dir);
                }
              Box etaBox = boundaryFaceBoxLo; // Create box to fill with etaFab
              etaBox.shiftHalf(dir, -1); // First interior cell box
              // Now fill the cells with etaFab
              for (int cell = 0; cell != numCellsIter; ++cell)
                {
                  MD_BOXLOOP(etaBox, i)
                    {
                      a_U[MD_IX(i, cEta)] = etaFab[MD_IX(i, 0)];
                    }
                  etaBox.shift(dir, -1); // Shift box to next layer of cells
                  etaFab.shift(dir, -1); // Shift fab to next layer of cells
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Find primitive turbulent variables from conservative variables
/** \param[out] a_W    Primitive variables
 *  \param[in]  a_U    Conservative variables
 *  \param[in]  a_box  Box to solve over
 *//*-----------------------------------------------------------------*/

void
SSV::turbConsToPrim(FArrayBox&       a_W,
                    const FArrayBox& a_U,
                    const Box&       a_box) const
{
  CH_TIME("SSV::turbConsToPrim");
  CRD::msg << CRD::fv4 << "SSV::turbConsToPrim" << CRD::end;
  const Interval turbIntv = CRDparam::g_CRDPhysics->turbConsInterval();
  // In this case, turbulent primitive and conservative variables are identical
  a_W.copy(a_U,a_box,turbIntv.begin(),a_box,turbIntv.begin(),turbIntv.size());
}

/*--------------------------------------------------------------------*/
//  Store SGS kinetic energy from current time to act as source term
/** \param[out] a_sgsKE
 *                     Cell-averaged, mapped SGS kinetic energy
 *  \param[in]  a_JU   Mapped conservative variables
 *  \param[in]  a_store
 *                     Store the SGS kinetic energy; otherwise, update
 *//*-----------------------------------------------------------------*/

void
SSV::storeSGSKineticEnergy(LevelData<FArrayBox>& a_sgsKE,
                           LevelData<FArrayBox>& a_JU,
                           const bool            a_store) const
{
  CH_TIME("SSV::storeSGSKineticEnergy");
  CRD::msg << CRD::fv4 << "SSV::storeSGSKineticEnergy" << CRD::end;

  const int sgsKEcomp = CRDparam::g_CRDPhysics->turbConsInterval().begin();
  if (a_store)
    {
      a_JU.copyTo(Interval(sgsKEcomp, sgsKEcomp), a_sgsKE, Interval(0,0));
    }
  else
    {
      for (DataIterator dit = a_JU.getBoxes().dataIterator(); dit.ok(); ++dit)
        {
          const Box disjointBox = a_JU.getBoxes()[dit];
          FArrayBox& keSrc = a_sgsKE[dit];
          FArrayBox& JU = a_JU[dit];
          MD_BOXLOOP(disjointBox, i)
            {
              keSrc[MD_IX(i, 0)] = JU[MD_IX(i, sgsKEcomp)] - keSrc[MD_IX(i, 0)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set wall-model data from finer level
/** \param[in]  a_JU    \<\JU\> state on this level
 *  \param[in]  a_blockDomain
 *                      Domain for this block on this level
 *  \param[in]  a_box   Current disjoint box
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the current level
 *  \param[in]  a_dataIndex
 *                      DataIndex for the current box
 *//*-----------------------------------------------------------------*/

void
SSV::setWallModelFromFiner(FArrayBox&              a_JU,
                           const ProblemDomain&    a_blockDomain,
                           const Box&              a_box,
                           const LevelGridMetrics& a_gridMetrics,
                           const DataIndex&        a_dataIndex) const
{
  // Check for boxes with both high and low wall boundaries
  IntVect boundedBox = IntVect_zero;
  const FArrayBox& JFab = a_gridMetrics.m_J[a_dataIndex];
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // Global constants
      const int cTurb = CRDparam::g_CRDPhysics->turbConsInterval().begin();
      const int cEta = cTurb + dir + 1;
      int hasHiWall = 0; // Check for high wall boundary
      int hasLoWall = 0; // Check for low wall boundary
      Box boundaryFaceBoxLo; // Low-side bndry-face box
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxLo, a_box, a_blockDomain, dir, Side::Lo);
      Box boundaryFaceBoxHi; // High-side bndry-face box
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxHi, a_box, a_blockDomain, dir, Side::Hi);
      // If boundaryFaceBoxLo is not empty, check for wall-boundary
      if (!boundaryFaceBoxLo.isEmpty())
        {
          // Check boundary condition
          BoundaryIndex bcIdx;
          bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_box),
                       dir, Side::Lo);
          BCInfo domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
          // We only support the use of walls in this function
          if ((CRDparam::DomainBCTypeSWall & domBC.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBC.m_type))
            {
              hasLoWall = 1;
            }
        }
      // If boundaryFaceBoxHi is not empty, check for wall-boundary
      if (!boundaryFaceBoxHi.isEmpty())
        {
          // Check boundary condition
          BoundaryIndex bcIdx;
          bcIdx.define(a_gridMetrics.getCoordSys().whichBlock(a_box),
                       dir, Side::Hi);
          BCInfo domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
          // We only support the use of walls in this function
          if ((CRDparam::DomainBCTypeSWall & domBC.m_type) &&
              !(CRDparam::DomainBCTypeMixed & domBC.m_type))
            {
              hasHiWall = 1;
            }
        }
      // If a_box has both high and low wall-boundaries, set boundedBox
      if (hasLoWall && hasHiWall)
        {
          boundedBox[dir] = 1;
        }
      // If low wall boundary, set the data
      if (hasLoWall)
        {
          Box etaBox = boundaryFaceBoxLo;
          etaBox.growHi(dir, 1); // Add one face on the high side of the box
          etaBox.enclosedCells(dir); // Turn this box into a cell-box
          CH_assert(a_JU.box().contains(etaBox));
          MD_BOXLOOP(etaBox, i)
            {
              Real eta = a_JU[MD_IX(i, cEta)]/JFab[MD_IX(i, 0)];
              Box wallNormalBox = etaBox; // Box of wall-normal cells for eta
              IntVect currVect = MD_GETIV(i); // Get the current IntVect
              wallNormalBox.setSmall(currVect); // Set low side to currVect
              int highSide = a_box.bigEnd()[dir]; // Wall-norm comp of high end
              if (boundedBox[dir]) // If both sides walls, only set half the box
                {
                  int halfBoxLength = 0.5*a_box.size(dir);
                  highSide = a_box.smallEnd()[dir] + halfBoxLength - 1;
                }
              // Set the wall-normal component of currVect to highSide
              currVect.setVal(dir, highSide);
              // Set the high side of wallNormalBox to the new currVect
              wallNormalBox.setBig(currVect);
              // Now loop over wallNormalBox and assign etaUpdate
              MD_BOXLOOP(wallNormalBox, j)
                {
                  Real J = JFab[MD_IX(j, 0)];
                  a_JU[MD_IX(j, cEta)] = J*eta;
                }
            }
        }
      // If high wall boundary, set the data
      if (hasHiWall)
        {
          Box etaBox = boundaryFaceBoxHi;
          etaBox.growLo(dir, 1); // Add one face on the low side of the box
          etaBox.enclosedCells(dir); // Turn this box into a cell-box
          CH_assert(a_JU.box().contains(etaBox));
          MD_BOXLOOP(etaBox, i)
            {
              Real eta = a_JU[MD_IX(i, cEta)]/JFab[MD_IX(i, 0)];
              Box wallNormalBox = etaBox; // Box of wall-normal cells for eta
              IntVect currVect = MD_GETIV(i); // Get the current IntVect
              wallNormalBox.setBig(currVect); // Set high side to currVect
              int lowSide = a_box.smallEnd()[dir]; // Wall-norm comp of low end
              if (boundedBox[dir]) // If both sides walls, only set half the box
                {
                  int halfBoxLength = 0.5*a_box.size(dir);
                  lowSide = a_box.bigEnd()[dir] - halfBoxLength + 1;
                }
              // Set the wall-normal component of currVect to lowSide
              currVect.setVal(dir, lowSide);
              // Set the low side of wallNormalBox to the new currVect
              wallNormalBox.setSmall(currVect);
              // Now loop over wallNormalBox and assign etaUpdate
              MD_BOXLOOP(wallNormalBox, j)
                {
                  Real J = JFab[MD_IX(j, 0)];
                  a_JU[MD_IX(j, cEta)] = J*eta;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Estimate eta_0 based on near-wall state
/** \param[out] a_etaFab
 *                      Eta_0 estimate on the wall faces (face box)
 *  \param[in]  a_W     Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_box   Box to initialize over
 *  \param[in]  a_disjointBox
 *                      Disjoint box of this dataIndex
 *  \param[in]  a_dir   The face-direction for initialization
 *  \param[in]  a_side  The side of the boundary for initialization
 *//*-----------------------------------------------------------------*/

void
SSV::etaEst(FArrayBox&              a_etaFab,
            const FArrayBox&        a_WcellAvgFab,
            const LevelGridMetrics& a_gridMetrics,
            const FArrayBox&        a_unitNormals,
            const Box&              a_box,
            const Box&              a_disjointBox,
            const int               a_dir,
            const Side::LoHiSide&   a_side) const
{
  CH_TIME("SSV::etaEst");

  // Constant indices
  const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int cRho = CRDparam::g_CRDPhysics->densityIndex();

  // Opposite side of current side (use for opposite shift direction)
  Side::LoHiSide oppSide = flip(a_side);
  // Sign for shifting into the domain: +1 if lo, -1 if hi
  const int offsetSign = sign(oppSide);
  // Sign for shifting out of the domain: -1 if lo, +1 if hi
  const int normSign = sign(a_side);

  // Get the block coordinate system
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));

  // Get the physical space coordinates on the wall face
  FABSTACKTEMP(faceXFab, a_box, SpaceDim);
  FABSTACKTEMP(faceXiFab, a_box, SpaceDim);
  CRDparam::g_CNSIBC->getFaceCoordinates(
    a_box, faceXiFab, faceXFab, a_dir, blockCoordSys);

  // Set up a box covering the second layer of interior cells
  Box secondInteriorCellBox = a_box;
  // Grow the face box by one on the interior side of the box
  secondInteriorCellBox.growDir(a_dir, oppSide, 1);
  // Turn the face box into a cell box
  secondInteriorCellBox.enclosedCells(a_dir);
  // Shift the cell box to the second interior cell (currently at first)
  secondInteriorCellBox.shift(a_dir, offsetSign);

  // Get the physical space coordinates in the second interior cell
  FABSTACKTEMP(cellXFab, secondInteriorCellBox, SpaceDim);
  FABSTACKTEMP(cellXiFab, secondInteriorCellBox, SpaceDim);
  CRDparam::g_CNSIBC->getCellCoordinates(
    secondInteriorCellBox, cellXiFab, cellXFab, blockCoordSys);
  // Shift cellXFab to the wall
  cellXFab.shiftHalf(a_dir, normSign); // First interior face
  cellXFab.shift(a_dir, normSign); // Wall face

  // Compute the distance between the wall and the cell
  MD_BOXLOOP(a_box, i)
    {
      for (int comp = 0; comp != SpaceDim; ++comp)
        {
          faceXFab[MD_IX(i, comp)] =
            offsetSign*(cellXFab[MD_IX(i, comp)] - faceXFab[MD_IX(i, comp)]);
        }
    }
  // Transform the wall-distance vector into normal-tangent space
  FORT_FORWARDTRANSFORMF(CHF_FRA(faceXFab),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_box));
  // Compute the wall-tangential component of the velocity
  FABSTACKTEMP(streamwiseTanVel, secondInteriorCellBox, SpaceDim);
  // Fill with values from the domain interior
  streamwiseTanVel.copy(a_WcellAvgFab, secondInteriorCellBox, cVel,
                        secondInteriorCellBox, 0, SpaceDim);
  // Shift streamwiseTanVel to the wall
  streamwiseTanVel.shiftHalf(a_dir, normSign);
  streamwiseTanVel.shift(a_dir, normSign);
  // Transform streamwiseTanVel int normal-tangent space
  FORT_FORWARDTRANSFORMF(CHF_FRA(streamwiseTanVel),
                         CHF_CONST_FRA(a_unitNormals),
                         CHF_BOX(a_box));
  // Set the normal component to zero
  streamwiseTanVel.setVal(0., a_dir);
  // Calculate the viscosity (and thermal conductivity unfortunately)
  FABSTACKTEMP(muCell, secondInteriorCellBox, 1);
  FABSTACKTEMP(kappaCell, secondInteriorCellBox, 1);
  CRDparam::g_CRDPhysics->calcCoeffKappaMu(
    secondInteriorCellBox,muCell,kappaCell,a_WcellAvgFab);
  // Turn mu into nu
  MD_BOXLOOP(secondInteriorCellBox, i)
    {
      muCell[MD_IX(i, 0)] /= a_WcellAvgFab[MD_IX(i, cRho)];
    }
  // Shift mu to the wall face
  muCell.shiftHalf(a_dir, normSign); // First interior face
  muCell.shift(a_dir, normSign); // Wall face
  // Compute eta_0
  MD_BOXLOOP(a_box, i)
    {
      Real nu = muCell[MD_IX(i, 0)];
      Real yVal = faceXFab[MD_IX(i, a_dir)];
      RealVect streamwiseVel(D_DECL(streamwiseTanVel[MD_IX(i, 0)],
                                    streamwiseTanVel[MD_IX(i, 1)],
                                    streamwiseTanVel[MD_IX(i, 2)]));
      Real u_tau = 1.e-10;
      Real uVel = streamwiseVel.vectorLength();
      if (uVel >= 0.) // If uVel = 0, u_tau is indeterminant
        {
          // Friction velocity (based on Spalding)
          Real yPlusMaxGuess = 1000000.;
          Real uTauMin = std::sqrt(uVel*nu/yVal);
          Real uTauMax = yPlusMaxGuess*nu/yVal;
          int iterBrent  = 0;
          int errorBrent = 0;

          const UTauSpaldingFunc& f = UTauSpaldingFunc(yVal,nu,uVel);
          u_tau = RootSolver::BrentER(iterBrent,errorBrent,f,uTauMin,uTauMax);
          if (errorBrent != 0 || u_tau != u_tau)
            {
              CRD::msg << "SSV: Bad uTau value: " << u_tau << CRD::error;
            }
        }
      a_etaFab[MD_IX(i, 0)] = u_tau*u_tau/nu;
    }
}

/*--------------------------------------------------------------------*/
//  Check if a box has a no-slip wall boundary
/** \param[out] a_etaFab
 *                      Eta_0 estimate on the wall faces (face box)
 *  \param[in]  a_W     Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_box   Box to initialize over
 *  \param[in]  a_disjointBox
 *                      Disjoint box of this dataIndex
 *  \param[in]  a_dir   The face-direction for initialization
 *  \param[in]  a_side  The side of the boundary for initialization
 *//*-----------------------------------------------------------------*/

void
SSV::hasNoSlipWall(int&                    a_hasNoSlipWall,
                   const LevelGridMetrics& a_gridMetrics,
                   const Box&              a_disjointBox,
                   const int               a_dir,
                   const Side::LoHiSide&   a_side) const
{
  CH_TIME("SSV::hasNoSlipWall");
  a_hasNoSlipWall = 0; // Default to not being a no-slip wall boundary
  BoundaryIndex bcIdx;
  bcIdx.define(
    a_gridMetrics.getCoordSys().whichBlock(a_disjointBox), a_dir, a_side);
  BCInfo domBC = CRDparam::g_CNSIBC->getDomainBC(bcIdx);
  if ((CRDparam::DomainBCTypeSWall & domBC.m_type) &&
      !(CRDparam::DomainBCTypeMixed & domBC.m_type))
    {
      a_hasNoSlipWall = 1;
    }
}

/*--------------------------------------------------------------------*/
//  Modify velocity gradient to incorporate wall-shear-stress from model
/** \param[out] a_NGradUfacePntFxb
 *                      Physical-space velocity-gradient
 *  \param[in]  a_unitNormalsFxb
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Box over which wall velocity-gradients
 *                      are modified based on wall-model
 *  \param[in]  a_disjointBox
 *                      Disjoint box of this dataIndex
 *//*-----------------------------------------------------------------*/

void
SSV::modeledNoSlipWallVelocityGradient(
  FluxBox&                a_NGradUfacePntFxb,
  const FluxBox&          a_unitNormalsFxb,
  const FArrayBox&        a_WcellAvgFab,
  const LevelGridMetrics& a_gridMetrics,
  const ProblemDomain&    a_domain,
  const Box&              a_box,
  const Box&              a_disjointBox) const
{
  // if (!m_enforceWallShearStress) return;
  // CH_TIME("SSV::modeledNoSlipWallVelocityGradient");

  // // Loop over all face directions
  // for (int dirSide = 0; dirSide != 2*SpaceDim; ++dirSide)
  //   {
  //     // The following few lines save us an indent by merging a "side loop"
  //     // with a "dir loop"
  //     // Translate the dirSide index into dir
  //     const int dir = dirSide/2;
  //     // Translate the dirSide index into side
  //     const int side = dirSide % 2;
  //     // Current side
  //     Side::LoHiSide whichSide = (side == 1) ? Side::Hi : Side::Lo;
  //     // Opposite side of current side (use for opposite shift direction)
  //     Side::LoHiSide oppSide = flip(whichSide);
  //     // Sign for shifting out of the domain: -1 if lo, +1 if hi
  //     const int normSign = sign(whichSide);
  //     // Constant indices
  //     const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
  //     const int cEta = CRDparam::g_CRDPhysics->turbConsInterval().begin()+dir+1;
  //     // Check if there is a boundary box
  //     Box bndryFaceBox;
  //     CRDparam::g_CNSIBC->getBoundaryFaces(
  //       bndryFaceBox, a_box, a_domain, dir, whichSide);

  //     int isWall = 0;
  //     Side::LoHiSide whichSide = (side == 1) ? Side::Hi : Side::Lo;
  //     hasNoSlipWall(isWall,a_gridMetrics,disjointBox,dir,whichSide);

  //     // Check if this is a no-slip wall
  //     if (isWall)
  //       {
  //         // Face-centered physical-space velocity-gradients
  //         FArrayBox& facePntVelGrad = a_NGradUfacePntFxb[dir];
  //         const FArrayBox& unitNormalsFab = a_unitNormalsFxb[dir];
  //         // Velocity gradients modified by LES wall-model
  //         FABSTACKTEMP(modelVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
  //         // Everything except for normal derivs should be zero
  //         modelVelGrad.setVal(0.);
  //         // Get the face-normal vector at the boundary
  //         FABSTACKTEMP(normVect, bndryFaceBox, SpaceDim);
  //         normVect.setVal(0.);
  //         normVect.setVal(1., dir);
  //         // Transform normalVect into Cartesian physical space
  //         FORT_REVERSETRANSFORMF(CHF_FRA(normVect),
  //                                CHF_CONST_FRA(unitNormalsFab),
  //                                CHF_BOX(bndryFaceBox));
  //         // Compute wall-normal derivative of wall-normal velocity
  //         // Just extract it from a_NGradUfacePntFxb
  //         MD_BOXLOOP(bndryFaceBox, i)
  //           {
  //             Real normDeriv = 0.;
  //             for (int comp = 0; comp != SpaceDim; ++comp)
  //               {
  //                 normDeriv += normVect[MD_IX(i, comp)]*(
  //                   D_TERM(
  //                     facePntVelGrad[MD_IX(i,srtIdx(comp,0))]*
  //                     normVect[MD_IX(i,0)],
  //                     + facePntVelGrad[MD_IX(i,srtIdx(comp,1))]*
  //                     normVect[MD_IX(i,1)],
  //                     + facePntVelGrad[MD_IX(i,srtIdx(comp,2))]*
  //                     normVect[MD_IX(i,2)]));
  //               }
  //             modelVelGrad[MD_IX(i,srtIdx(dir, dir))] = normDeriv;
  //           }
  //         // Create box covering first layer of wall-adjacent cells
  //         Box firstInteriorCellBox = bndryFaceBox;
  //         // Grow the face box by one on the interior side of the box
  //         firstInteriorCellBox.growDir(dir, oppSide, 1);
  //         // Turn the face box into a cell box
  //         firstInteriorCellBox.enclosedCells(dir);
  //         // Store velocity for shifting (a_WcellAvgFab stays const)
  //         FABSTACKTEMP(streamTanVect, firstInteriorCellBox, SpaceDim);
  //         // Fill with a_WcellAvgFab
  //         streamTanVect.copy(a_WcellAvgFab, cVel, 0, SpaceDim);
  //         streamTanVect.shiftHalf(dir, normSign); // Shift to wall face
  //         // Store the eta_0 value for shifting purposes
  //         FABSTACKTEMP(eta, firstInteriorCellBox, 1);
  //         // Fill with a_WcellAvgFab
  //         eta.copy(a_WcellAvgFab, cEta, 0, 1);
  //         eta.shiftHalf(dir, normSign); // Shift to wall face
  //         // Transform the velocity vector (in-place transformation)
  //         FORT_FORWARDTRANSFORMF(CHF_FRA(streamTanVect),
  //                                CHF_CONST_FRA(unitNormalsFab),
  //                                CHF_BOX(bndryFaceBox));
  //         // Zero the normal comp to get wall-tangent velocity vector
  //         streamTanVect.setVal(0.0, dir);
  //         // Normalize wall-tangent vector -- streamwise unit vector
  //         PatchMappedFunc::normalize(
  //           bndryFaceBox, streamTanVect, Interval(0,SpaceDim-1));
  //         // Fill normal-derivs of tangential velocities using model
  //         for (int velComp = 0; velComp != SpaceDim; ++velComp)
  //           {
  //             if (velComp != dir)
  //               {
  //                 const int cGrad = dir + SpaceDim*velComp;
  //                 MD_BOXLOOP(bndryFaceBox, i)
  //                   {
  //                     modelVelGrad[MD_IX(i, cGrad)] =
  //                       streamTanVect[MD_IX(i,velComp)]*eta[MD_IX(i,0)];
  //                   }
  //               }
  //           }
  //         // Transform back to Cartesian space (a lot here)
  //         FABSTACKTEMP(tempWallVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
  //         // First, right multiply grad(u) by a_unitNormals
  //         for (int row = 0; row != SpaceDim; ++row)
  //           {
  //             for (int col = 0; col != SpaceDim; ++col)
  //               {
  //                 const int cGrad = col + SpaceDim*row;
  //                 MD_BOXLOOP(bndryFaceBox, i)
  //                   {
  //                     tempWallVelGrad[MD_IX(i, cGrad)] =
  //                       D_TERM(
  //                         modelVelGrad[MD_IX(i,srtIdx(row,0))]*
  //                         unitNormalsFab[MD_IX(i,srtIdx(0,col))],
  //                         + modelVelGrad[MD_IX(i,srtIdx(row,1))]*
  //                         unitNormalsFab[MD_IX(i,srtIdx(1,col))],
  //                         + modelVelGrad[MD_IX(i,srtIdx(row,2))]*
  //                         unitNormalsFab[MD_IX(i,srtIdx(2,col))]);
  //                   }
  //               }
  //           }
  //         // Next, left multiply by a_unitNormals^T
  //         FABSTACKTEMP(wallVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
  //         for (int row = 0; row != SpaceDim; ++row)
  //           {
  //             for (int col = 0; col != SpaceDim; ++col)
  //               {
  //                 const int cGrad = col + SpaceDim*row;
  //                 MD_BOXLOOP(bndryFaceBox, i)
  //                   {
  //                     wallVelGrad[MD_IX(i, cGrad)] =
  //                       D_TERM(
  //                         unitNormalsFab[MD_IX(i,srtIdx(0,row))]*
  //                         tempWallVelGrad[MD_IX(i,srtIdx(0,col))],
  //                         + unitNormalsFab[MD_IX(i,srtIdx(1,row))]*
  //                         tempWallVelGrad[MD_IX(i,srtIdx(1,col))],
  //                         + unitNormalsFab[MD_IX(i,srtIdx(2,row))]*
  //                         tempWallVelGrad[MD_IX(i,srtIdx(2,col))]);
  //                   }
  //               }
  //           }
  //         // Finally, copy the result into a_NGradUfacePntFxb
  //         for (int c = 0; c != SpaceDim*SpaceDim; ++c)
  //           {
  //             MD_BOXLOOP(bndryFaceBox, i)
  //               {
  //                 facePntVelGrad[MD_IX(i, c)] = wallVelGrad[MD_IX(i, c)];
  //               }
  //           }
  //       }
  //   }
}
