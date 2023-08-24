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
 * \file Smagorinsky.cpp
 *
 * \brief Member functions for Smagorinksy LES SGS model
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "FArrayBox.H"

//----- Internal -----//

#include "Smagorinsky.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"

/*******************************************************************************
 *
 * Class Smagorinsky: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

Smagorinsky::Smagorinsky()
  :
  m_Cs(0.1),
  m_Prt(0.71)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

Smagorinsky::~Smagorinsky()
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
 *  \param[out] a_kappaFab
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
Smagorinsky::calcCoeffKappatMut(const Box&              a_box,
                                FArrayBox&              a_mutFab,
                                FArrayBox&              a_kappatFab,
                                const FArrayBox&        a_muFab,
                                const FArrayBox&        a_kappaFab,
                                const FArrayBox&        a_WfacePntFab,
                                const FArrayBox&        a_strainRateFab,
                                const LevelGridMetrics& a_gridMetrics,
                                const int               a_dir) const
{
  CH_TIME("Smagorinsky::calcCoeffKappatMut");
  CRD::msg << CRD::fv4 << "Smagorinsky::calcCoeffKappatMut" << CRD::end;

  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int numStrainComp = SpaceDim*SpaceDim;
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  std::vector<Real> cn(CRDparam::g_numSpecies, 0.);
  Real dx = a_gridMetrics.dxVect()[a_dir];
  CH_assert(numStrainComp == a_strainRateFab.nComp());
  MD_ARRAY_RESTRICT(arrW, a_WfacePntFab);
  MD_ARRAY_RESTRICT(arrMut, a_mutFab);
  MD_ARRAY_RESTRICT(arrKappat, a_kappatFab);
  MD_ARRAY_RESTRICT(arrSRT, a_strainRateFab);

  MD_BOXLOOP(a_box, i)
    {
      const Real rho = arrW[MD_IX(i, rhoIndx)];
      const Real T = arrW[MD_IX(i, tempIndx)];
      for(int comp = 0; comp != CRDparam::g_numSpecies; ++comp)
        {
          const int wComp = comp + wCompStart;
          cn[comp] = arrW[MD_IX(i, wComp)];
        }
      Real Cp = CRDparam::g_CRDPhysics->cp(T, cn.data());
      Real MagSRT = 0.0;
      for(int comp = 0; comp != numStrainComp; comp++)
        {
          MagSRT += std::pow(arrSRT[MD_IX(i, comp)], 2);
        }
      MagSRT = std::sqrt(2.*MagSRT);
      Real mut = rho*m_Cs*m_Cs*dx*dx*MagSRT;
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
Smagorinsky::addTurbulentConvectiveFlux(FArrayBox&       a_flux,
                                        const FArrayBox& a_WFace,
                                        const int&       a_dir,
                                        const Box&       a_box) const
{
  return;
}
