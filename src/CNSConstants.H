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
 * \file CNSConstants.H
 *
 * \brief Relavant constants for compressible Navier-Stokes
 *
 *//*+*************************************************************************/

#ifndef _CNSCONSTANTS_H_
#define _CNSCONSTANTS_H_

//----- Chombo Library -----//

#include "REAL.H"
#include "UsingNamespace.H"

//**FIXME Probably want to check and reference ISBN 92-9194-004-6
/// International standard atmosphere
namespace ISA
{
  const Real g_R              = 287.1;                      // (J/(kg K))
  const Real g_gamma          = 1.4;                        //
  const Real g_rho            = 1.225;                      // (kg/m^3)
  const Real g_pressure       = 101325.;                    // (Pa)
  const Real g_temperature    = 288.15;                     // (K)
  const Real g_mu             = 1.789E-5;                   // (kg/m s)
  const Real g_kappa          = 2.5326E-2;                  // (W/(m K))
  const Real g_gravity        = 9.80665;                    // (m/s^2)
}

#endif /* ! defined _CNSCONSTANTS_H_ */
