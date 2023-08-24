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
 * \file CRDState.cpp
 *
 * \brief Member functions for CRDState
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "ParmParse.H"

//----- Internal -----//

#include "CRDState.H"
#include "CRDmsg.H"


/*******************************************************************************
 *
 * Class CRDState: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Static data member initialization
 *============================================================================*/

std::vector<CRDState> CRDState::s_states;
std::unordered_map<std::string, int,
                   CH_Hash::google_CityHash<std::string>> CRDState::s_names(64);


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

CRDState::CRDState()
  :
  m_frameVelocity(RealVect_zero)
{
  const int numPhysics = CRDparam::numCRDPhysics();
  int maxPrim = 0;
  for (int idxPhys = 0; idxPhys != numPhysics; ++idxPhys)
    {
      maxPrim = std::max(maxPrim,
                         CRDparam::CRDPhysicsSet(idxPhys)->numPrimitive());
    }
  m_state.define(numPhysics, maxPrim);
  m_state = 0.;
}

/*--------------------------------------------------------------------*/
//  Copy
/**
 *//*-----------------------------------------------------------------*/

CRDState::CRDState(const CRDState& a_other)
  :
  m_name(a_other.m_name),
  m_state(a_other.m_state.size(1), a_other.m_state.size(0)),
  m_frameVelocity(a_other.m_frameVelocity),
  m_thermoUser(a_other.m_thermoUser),
  m_thermoDefs(a_other.m_thermoDefs)
{
  const Real *q = a_other.m_state.begin();
  for (Real* p = m_state.begin(), *p_end = m_state.end(); p != p_end;)
    {
      *p++ = *q++;
    }
}

/*--------------------------------------------------------------------*/
//  Assignment
/**
 *//*-----------------------------------------------------------------*/

CRDState& CRDState::operator=(const CRDState& a_other)
{
  if (&a_other != this)
    {
      m_name = a_other.m_name;
      m_state.define(a_other.m_state.size(1), a_other.m_state.size(0));
      const Real *q = a_other.m_state.begin();
      for (Real* p = m_state.begin(), *p_end = m_state.end(); p != p_end;)
        {
          *p++ = *q++;
        }
      m_frameVelocity = a_other.m_frameVelocity;
      m_thermoUser = a_other.m_thermoUser;
      m_thermoDefs = a_other.m_thermoDefs;
    }
  return *this;
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Read all state information
/**
 *//*-----------------------------------------------------------------*/

void
CRDState::readStateInfo()
{
  ParmParse ppSTATE("state");
  char physSetName[64] = { '\0' };
  std::string stateName;
  const int numState = ppSTATE.countval("names");
  for (int idxState = 0; idxState != numState; ++idxState)
    {
      ppSTATE.get("names", stateName, idxState);
      if (CRDState::nameIndex(stateName) != -1)
        {
          CRD::msg << "Duplicate state names in input file.  Previous state "
            "with name \"" << stateName << "\" will be overwritten."
                   << CRD::warn;
        }
      CRDState& state = CRDState::get(stateName);
      {
        
        std::string varName(stateName + "_frame-velocity");
        std::vector<Real> value(SpaceDim, 0.0);
        ppSTATE.queryarr(varName.c_str(), value, 0, SpaceDim);
        for (const int idxVel : EachDir)
          {
            state.m_frameVelocity[idxVel] = value[idxVel];
          }
      }
      for (int idxPhys = 0; idxPhys != CRDparam::numCRDPhysics(); ++idxPhys)
        {
          physSetName[0] = '\0';
          if (idxPhys > 0)
            {
              sprintf(physSetName, "_phys%d", idxPhys);
            }
          const CRDPhysics* physics = CRDparam::CRDPhysicsSet(idxPhys);
          for (int comp = 0; comp != physics->numPrimitive(); ++comp)
            {
              if (physics->velocityInterval().contains(comp))
                {
                  const int cVel = physics->velocityInterval().begin();
                  if (comp == cVel)
                    {
                      std::string varName(stateName + physSetName +
                                          "_velocity");
                      if (ppSTATE.contains(varName))
                        {
                          std::vector<Real> value(SpaceDim);
                          ppSTATE.getarr(varName.c_str(), value, 0,
                                         SpaceDim);
                          for (const int idxVel : EachDir)
                            {
                              state(comp + idxVel, idxPhys) = value[idxVel];
                            }
                        }
                    }
                }
              else
                {
                  std::string varName(stateName + physSetName + '_' +
                                      physics->primStateName(comp));
                  if (ppSTATE.contains(varName))
                    {
                      Real value = 0.;
                      ppSTATE.get(varName.c_str(), value);
                      state(comp, idxPhys) = value;
		      const int densityIdx = physics->densityIndex();
		      const int pressIdx = physics->pressureIndex();
                      if (comp == densityIdx)
			{
			  state.m_thermoUser |= ThermoDefDensity;
                          state.m_thermoDefs |= ThermoDefDensity;
			}
		      else if (comp == pressIdx)
			{
                          state.m_thermoUser |= ThermoDefPressure;
                          state.m_thermoDefs |= ThermoDefPressure;
			}
		      else
			{
			  state.m_thermoUser |= ThermoDefTemperature;
			  state.m_thermoDefs |= ThermoDefTemperature;
			}
                    }
                }
            }  // Loop over components
          state.setExtraThermo();
        }  // Loop over physics sets
    }
}

/*--------------------------------------------------------------------*/
//  Write all state information
/**
 *//*-----------------------------------------------------------------*/

void
CRDState::writeStateInfo()
{
  for (const auto names : s_names)
    {
      CRD::msg << "State: " << names.first << " (" << nameIndex(names.first)
               << ')' << CRD::end;
      CRDState& state = CRDState::get(names.second);
      CRD::msg << "  frame-velocity\n(";
      for (const int idxVel : EachDir)
        {
          if (idxVel != 0) CRD::msg << ',';
          CRD::msg << state.m_frameVelocity[idxVel];
        }
      CRD::msg << ')' << CRD::var;
      for (int idxPhys = 0; idxPhys != CRDparam::numCRDPhysics(); ++idxPhys)
        {
          const CRDPhysics* physics = CRDparam::CRDPhysicsSet(idxPhys);
          for (int comp = 0; comp != CRDparam::g_CRDPhysics->numPrimitive();
               ++comp)
            {
              if (physics->velocityInterval().contains(comp))
                {
                  const int cVel = physics->velocityInterval().begin();
                  if (comp == cVel)
                    {
                      CRD::msg << "  velocity\n(";
                      for (const int idxVel : EachDir)
                        {
                          if (idxVel != 0) CRD::msg << ',';
                          CRD::msg << state(comp + idxVel, idxPhys);
                        }
                      CRD::msg << ')' << CRD::var;
                    }
                }
              else
                {
                  CRD::msg << "  " << physics->primStateName(comp)
                           << '\n' << state(comp, idxPhys)
                           << CRD::var;
                }
            }  // Loop over components
        }  // Loop over physics
    }
}

/*--------------------------------------------------------------------*/
//  Get the unique integer index for a given state name
/** \param[in]  a_name  Name of the state
 *  \return             Integer index
 *//*-----------------------------------------------------------------*/

int
CRDState::nameIndex(const std::string a_name)
{
  auto iter = s_names.find(a_name);
  if (iter == s_names.end())
    {
      return -1;
    }
  else
    {
      return iter->second;
    }
}

/*--------------------------------------------------------------------*/
//  Get a state for a given name
/** This will allocate a new state if it does not yet exist
 *  \param[in]  a_name  Name of the state
 *  \return             State
 *//*-----------------------------------------------------------------*/

CRDState&
CRDState::get(const std::string a_name)
{
  auto ret = s_names.insert(std::make_pair(a_name, s_states.size()));
  if (ret.second)  // New insertion
    {
      s_states.push_back(CRDState{});
    }
  return s_states[ret.first->second];
}

/*--------------------------------------------------------------------*/
//  Set extra thermodynamic states (if 2 are defined)
/** \param[in]  a_idxPhys
 *                      Physics set (default 0)
 *//*-----------------------------------------------------------------*/

void
CRDState::setExtraThermo(const int a_idxPhys)
{
  const CRDPhysics* physics = CRDparam::CRDPhysicsSet(a_idxPhys);
  const Real R = physics->Rgas(
    &(this->operator()(std::max(0, physics->speciesPrimInterval().begin()),
                       a_idxPhys)),
    1);
  const unsigned setThermo =
    m_thermoDefs & (  ThermoDefDensity
                    | ThermoDefPressure
                    | ThermoDefTemperature);
  // All using ideal gas law
  switch (setThermo)
    {
    case (ThermoDefDensity | ThermoDefPressure):
      temperature(a_idxPhys) = pressure(a_idxPhys)/(R*density(a_idxPhys));
      m_thermoDefs |= ThermoDefTemperature;
      break;
    case (ThermoDefPressure | ThermoDefTemperature):
      density(a_idxPhys) = pressure(a_idxPhys)/(R*temperature(a_idxPhys));
      m_thermoDefs |= ThermoDefDensity;
      break;
    case (ThermoDefTemperature | ThermoDefDensity):
      pressure(a_idxPhys) = density(a_idxPhys)*R*temperature(a_idxPhys);
      m_thermoDefs |= ThermoDefPressure;
      break;
    }
}
