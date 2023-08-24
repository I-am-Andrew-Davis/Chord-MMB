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
 * \file CNSIBCReferenceCubeBC.cpp
 *
 * \brief Member functions for CNSIBCReferenceCubeBC
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "SetCentersF_F.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"
#include "CNSIBCReferenceCubeBCF_F.H"
#include "CRDPhysics.H"

/*******************************************************************************
 *
 * Class CNSIBC: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor set BC at domain extents
/** This will probably move to the CS for multiblock
 *  \param[in]  a_domainBC
 *                      An array of BC at domain extents.  The array
 *                      is ordered [SpaceDim][Side] in C ordering.
 *  \param[in]  a_wnum  Number of primitive states
 *  \param[in]  a_rhoRef
 *                      Reference (freestream) density
 *  \param[in]  a_uRef  Reference (freestream) fluid velocity.  This
 *                      must be greater than zero.  Use sound speed if
 *                      in doubt.
 *//*-----------------------------------------------------------------*/

CNSIBCReferenceCubeBC::CNSIBCReferenceCubeBC(
  const int  a_wnum,
  const Real a_rhoRef,
  const Real a_uRef)
:
  CNSIBC()
{
  CRD::msg << "!! Inheriting from CNSIBCReferenceCube is deprecated !!\n Update your ibc file to inherit directly from CNSIBC or CNSIBCGeneralized instead." << CRD::error;

  
  m_refPrimStateBC.define(SpaceDim, 2, a_wnum + CRDparam::g_numSpecies);
  // Set m_primStateBC to default values
  Real *p;
  const int densityComp = URHO;
  p = &m_refPrimStateBC(0, 0, densityComp);
  for (int n = 2*SpaceDim; n--;) *p++ = CRDparam::g_rho;
  D_TERM(
    p = &m_refPrimStateBC(0, 0, WVELX);
    for (int n = 2*SpaceDim; n--;) *p++ = 0.;,
    p = &m_refPrimStateBC(0, 0, WVELY);
    for (int n = 2*SpaceDim; n--;) *p++ = 0.;,
    p = &m_refPrimStateBC(0, 0, WVELZ);
    for (int n = 2*SpaceDim; n--;) *p++ = 0.;)
  const Real pres = CRDparam::g_rho*CRDparam::g_R*CRDparam::g_T;
  const int presComp = WPRES;
  p = &m_refPrimStateBC(0, 0, presComp);
  for (int n = 2*SpaceDim; n--;) *p++ = pres;
  // p = &m_refPrimStateBC(0, 0, WTMPT);
  // for (int n = 2*SpaceDim; n--;) *p++ = CRDparam::g_T;
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCReferenceCubeBC::~CNSIBCReferenceCubeBC()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/

const char *const
CNSIBCReferenceCubeBC::IBCName() const
{
  return "reference cube";
}

/*--------------------------------------------------------------------*/
//  Set primitive state on a domain boundary
/** \param[in]  a_dir   Direction normal to boundary
 *  \param[in]  a_side  Side of boundary (-1 or 0 = low, 1 = high)
 *  \param[in]  a_comp  Component to set
 *  \param[in]  a_val   New value for component
 *//*-----------------------------------------------------------------*/

void
CNSIBCReferenceCubeBC::setReferenceBCState(const int  a_dir,
                                           const int  a_side,
                                           const int  a_comp,
                                           const Real a_val)
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side >= -1 && a_side <= 1);
  CH_assert(a_comp >= 0);
  m_refPrimStateBC(a_dir, (a_side+1)/2, a_comp) = a_val;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the imposed (exterior or farfield) primitive state at flow BCC
/** State is set according to reference conditions
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
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
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCReferenceCubeBC::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const int                     a_dir,
  const Side::LoHiSide&         a_side,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const CRDparam::DomainBCType& a_domT) const
{
  const int lohiSign = sign(a_side);
  const int lohiIdx  = (lohiSign + 1)/2;
  for (int comp = 0; comp != CRDparam::g_CRDPhysics->numPrimitive(); 
       ++comp)
    {
      a_Wface.setVal(m_refPrimStateBC(a_dir, lohiIdx, comp),
                     a_boundaryFaceBox,
                     comp);
    }
}

  

