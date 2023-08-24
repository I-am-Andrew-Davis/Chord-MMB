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
 * \file Smagorinsky.H
 *
 * \brief Class for Large Eddy Simulation
 *
 *//*+*************************************************************************/

#ifndef _SMAGORINSKY_H_
#define _SMAGORINSKY_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "LES.H"

//--Forward declarations

/*******************************************************************************
 */
///  Smagorinsky derived from LES base class
/**
 *   This class is an implementation of Large Eddy Simulation
 *
 *//*+*************************************************************************/

class Smagorinsky : public LES
{
  
/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  Smagorinsky();
  

  /// Destructor
  ~Smagorinsky();

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  Smagorinsky(const Smagorinsky&);

//--Assignment not permitted

  Smagorinsky &operator=(const Smagorinsky&);

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--Variables

  /// Return the number of turbulent variables
  virtual int numTurbVar() const;

//--Misc

  /// Output turbulence model info
  virtual void writeTurbModelInfo() const;

  /// Output turbulence model name
  virtual const char *const modelName() const;

//--Algorithms (on a patch)

  /// Calculate SGS thermal conductivity and dynamic viscosity
  virtual void calcCoeffKappatMut(const Box&              a_box,
                                  FArrayBox&              a_mutFab,
                                  FArrayBox&              a_kappatFab,
                                  const FArrayBox&        a_muFab,
                                  const FArrayBox&        a_kappaFab,
                                  const FArrayBox&        a_WfacePntFab,
                                  const FArrayBox&        a_strainRateFab,
                                  const LevelGridMetrics& a_gridMetrics,
                                  const int               a_dir) const;

  /// Calculate the turbulent convective fluxes
  virtual void addTurbulentConvectiveFlux(FArrayBox&       a_flux,
                                          const FArrayBox& a_WFace,
                                          const int&       a_dir,
                                          const Box&       a_box) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  Real m_Cs;                          ///< Some constant
  Real m_Prt;                         ///< Turbulent Prandtl number
};

/*******************************************************************************
 *
 * Class Smagorinsky: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of turbulence variables
/** 
 *//*-----------------------------------------------------------------*/

inline int
Smagorinsky::numTurbVar() const
{
  return 0;
}

/*--------------------------------------------------------------------*/
//  Write any information related to the turbulence model
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

inline void
Smagorinsky::writeTurbModelInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Cs\n" << m_Cs << CRD::var;
  CRD::msg << "Turbulent Prandtl number\n" << m_Prt << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Write the name of the turbulence model
/**
 *//*-----------------------------------------------------------------*/

inline const char *const
Smagorinsky::modelName() const
{
  return "LES: Standard Smagorinsky SGS model";
}

#endif  /* ! defined _SMAGORINSKY_H_ */
