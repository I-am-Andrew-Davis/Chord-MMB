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
 * \file CNSIBCReferenceCubeBC.H
 *
 * \brief IBC implementing BC methods for regular BC in a (hyper-)cube
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCREFERENCECUBEBC_H_
#define _CNSIBCREFERENCECUBEBC_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "CHArray.H"

//----- Internal -----//

#include "CNSIBC.H"

//---- Debug ----//
#include "CRDhardconfig.H"

/*******************************************************************************
 */
///  BC and initial conditions on the lattice.
/**  This class serves as reference for setting initial and boundary conditions
 *   on the lattice.  It is expected that more complicated setups derive from
 *   this class.  Note that pure function 'initialize' is not implemented so
 *   this class is not sufficient as an IBC.
 *
 *//*+*************************************************************************/

class CNSIBCReferenceCubeBC : public CNSIBC
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor sets BC at domain extents
  CNSIBCReferenceCubeBC(
    const int  a_wnum = WNUM,
    const Real a_rhoRef = CRDparam::g_gamma,
    const Real a_uRef = 1.);

  /// Destructor
  virtual ~CNSIBCReferenceCubeBC();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCReferenceCubeBC(const CNSIBCReferenceCubeBC&);

//--Assignment not permitted

  CNSIBCReferenceCubeBC &operator=(const CNSIBCReferenceCubeBC&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Set primitive state on a domain boundary
  virtual void setReferenceBCState(const int  a_dir,
                                   const int  a_side,
                                   const int  a_comp,
                                   const Real a_val);

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the imposed (exterior or farfield) primitive state at flow BC
  virtual void setImposedBCprimState(
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
    const CRDparam::DomainBCType& a_domT) const;
  
/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  CHArray<Real, 3, ArZeroCol, ArSp::NewArrayAlloc<Real> > m_refPrimStateBC;
                                      ///< The primitive state on the boundary
                                      ///< (set and use only as required)
                                      ///< (dir, side, comp)
};

#endif /* ! defined _CNSIBCREFERENCECUBEBC_H_ */
