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
 * \file CRDState.H
 *
 * \brief State information, often used for IBC
 *
 *//*+*************************************************************************/

#ifndef _CRDSTATE_H_
#define _CRDSTATE_H_

//----- Chombo Library -----//

#include "RealVect.H"
#include "CH_Hash.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "TurbModeling.H"
#include "LES.H"
#include "CRDparam.H"
#include "CNSIBC.H"
#include "CRDPhysics.H"


/*******************************************************************************
 */
///  State information
/**
 *   This class stores state information in a global sense.
 *
 *   WARNING: Adding new states can invalidate pointers to existing states.  For
 *   quick lookup, it is best to save state indicies and only use pointers or
 *   references in blocks where new states are not added.
 *
 *//*+*************************************************************************/

class CRDState
{


/*==============================================================================
 * Types
 *============================================================================*/

  enum ThermoDef : unsigned
  {
    ThermoDefUndefined   = 0,
    ThermoDefDensity     = (1<<0),
    ThermoDefPressure    = (1<<1),
    ThermoDefTemperature = (1<<2)
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  CRDState();

  /// Copy
  CRDState(const CRDState& a_other);

  /// Assignment
  CRDState &operator=(const CRDState& a_other);

  /// Destructor
  ~CRDState() = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// Read all state information
  static void readStateInfo();

  /// Write all state information
  static void writeStateInfo();

  /// Get the unique integer index for a given state name
  static int nameIndex(const std::string a_name);

  /// Get a state for a given name
  static CRDState& get(const std::string a_name);

  /// Get a state for a given index
  static CRDState& get(const int a_idx);

  /// Get a particular component
  const Real& operator()(const int a_comp, const int a_idxPhys = 0) const;

  /// Get a modifiable particular component
  Real& operator()(const int a_comp, const int a_idxPhys = 0);
  
  /// Get the density
  const Real& density(const int a_idxPhys = 0) const;

  /// Get the modifiable density
  Real& density(const int a_idxPhys = 0);

  /// Get the velocity vector (does not return a RealVect, use auto)
  const RealVectConstAlias velocity(const int a_idxPhys = 0) const;

  /// Get the modifiable velocity vector (does not return a RealVect, use auto)
  RealVectAlias velocity(const int a_idxPhys = 0);

  /// Get the pressure
  const Real& pressure(const int a_idxPhys = 0) const;

  /// Get the modifiable pressure
  Real& pressure(const int a_idxPhys = 0);

  /// Get the temperature
  const Real& temperature(const int a_idxPhys = 0) const;

  /// Get the modifiable temperature
  Real& temperature(const int a_idxPhys = 0);

  /// Set extra thermodynamic states (if 2 are defined)
  void setExtraThermo(const int a_idxPhys = 0);


/*==============================================================================
 * Data members
 *============================================================================*/

public:

  std::string m_name;                 ///< Name of the state
  CHArray<Real, 2, ArZeroRow, ArSp::NewArrayAlloc<Real>> m_state;
                                      ///< Vector of state information
                                      ///< (physics_set, comp)
  RealVect m_frameVelocity;           ///< Reference frame velocity for working
                                      ///< with stagnation conditions
  unsigned m_thermoUser = ThermoDefUndefined;
                                      ///< Thermodynamics states defined by the
                                      ///< user
  unsigned m_thermoDefs = ThermoDefUndefined;
                                      ///< Thermodyanic states that are defined

private:

  static std::vector<CRDState> s_states;
  static std::unordered_map<std::string, int,
                            CH_Hash::google_CityHash<std::string>> s_names;
                                      ///< A map of state names to integer keys
};


/*******************************************************************************
 *
 * Class CRDState: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get a state for a given index
/** \param[in]  a_idx   Index
 *  \return             State
 *//*-----------------------------------------------------------------*/

inline CRDState&
CRDState::get(const int a_idx)
{
  CH_assert (a_idx < s_states.size());
  return s_states[a_idx];
}

/*--------------------------------------------------------------------*/
//  Get a particular component
/** \param[in]  a_comp  Component index
 *  \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Density
 *//*-----------------------------------------------------------------*/

inline const Real&
CRDState::operator()(const int a_comp, const int a_idxPhys) const
{
  CH_assert(a_comp < CRDparam::CRDPhysicsSet(a_idxPhys)->numPrimitive());
  return m_state(a_idxPhys, a_comp);
}

/*--------------------------------------------------------------------*/
//  Get a modifiable particular component
/** \param[in]  a_comp  Component index
 *  \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Modifiable density
 *//*-----------------------------------------------------------------*/

inline Real&
CRDState::operator()(const int a_comp, const int a_idxPhys)
{
  CH_assert(a_comp < CRDparam::CRDPhysicsSet(a_idxPhys)->numPrimitive());
  return m_state(a_idxPhys, a_comp);
}

/*--------------------------------------------------------------------*/
//  Get the density
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Density
 *//*-----------------------------------------------------------------*/

inline const Real&
CRDState::density(const int a_idxPhys) const
{
  return m_state(a_idxPhys, CRDparam::CRDPhysicsSet(a_idxPhys)->densityIndex());
}

/*--------------------------------------------------------------------*/
//  Get the modifiable density
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Modifiable density
 *//*-----------------------------------------------------------------*/

inline Real&
CRDState::density(const int a_idxPhys)
{
  m_thermoUser |= ThermoDefDensity;
  m_thermoDefs |= ThermoDefDensity;
  return m_state(a_idxPhys, CRDparam::CRDPhysicsSet(a_idxPhys)->densityIndex());
}

/*--------------------------------------------------------------------*/
//  Get the velocity vector (does not return a RealVect, use auto)
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Velocity vector
 *
 *  The returned type will behave like a RealVect but is actually an
 *  alias
 *//*-----------------------------------------------------------------*/

inline const RealVectConstAlias
CRDState::velocity(const int a_idxPhys) const
{
  return RealVectConstAlias(
    &m_state(a_idxPhys,
             CRDparam::CRDPhysicsSet(a_idxPhys)->velocityInterval().begin()));
}

/*--------------------------------------------------------------------*/
//  Get the modifiable velocity vector (does not return a RealVect,
//  use auto)
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Modifiable velocity vector
 *
 *  The returned type will behave like a RealVect but is actually an
 *  alias
 *//*-----------------------------------------------------------------*/

inline RealVectAlias
CRDState::velocity(const int a_idxPhys)
{
  return RealVectAlias(
    &m_state(a_idxPhys,
             CRDparam::CRDPhysicsSet(a_idxPhys)->velocityInterval().begin()));
}

/*--------------------------------------------------------------------*/
//  Get the pressure
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Pressure
 *//*-----------------------------------------------------------------*/

inline const Real&
CRDState::pressure(const int a_idxPhys) const
{
  return m_state(a_idxPhys,
                 CRDparam::CRDPhysicsSet(a_idxPhys)->pressureIndex());
}

/*--------------------------------------------------------------------*/
//  Get the modifiable pressure
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Modifiable pressure
 *//*-----------------------------------------------------------------*/

inline Real&
CRDState::pressure(const int a_idxPhys)
{
  m_thermoUser |= ThermoDefPressure;
  m_thermoDefs |= ThermoDefPressure;
  return m_state(a_idxPhys,
                 CRDparam::CRDPhysicsSet(a_idxPhys)->pressureIndex());
}

/*--------------------------------------------------------------------*/
//  Get the temperature
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Temperature
 *//*-----------------------------------------------------------------*/

inline const Real&
CRDState::temperature(const int a_idxPhys) const
{
  return m_state(a_idxPhys,
                 CRDparam::CRDPhysicsSet(a_idxPhys)->temperatureIndex());
}

/*--------------------------------------------------------------------*/
//  Get the modifiable temperature
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *  \return             Modifiable temperature
 *//*-----------------------------------------------------------------*/

inline Real&
CRDState::temperature(const int a_idxPhys)
{
  m_thermoUser |= ThermoDefTemperature;
  m_thermoDefs |= ThermoDefTemperature;
  return m_state(a_idxPhys,
                 CRDparam::CRDPhysicsSet(a_idxPhys)->temperatureIndex());
}

#endif  /* ! defined _CRDSTATE_H_ */
