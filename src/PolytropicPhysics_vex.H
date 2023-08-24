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
 * \file CNSPhysics.H
 *
 * \brief Vector implementations of kernels in PolytropicPhysics
 *
 *//*+*************************************************************************/

#ifndef _POLYTROPICPHYSICS_VEX_H_
#define _POLYTROPICPHYSICS_VEX_H_

#include "UsingNamespace.H"
#include "NamespaceVar.H"

//-- Forward declarations from Chombo

#include "NamespaceHeader.H"
class Box;
class FArrayBox;
#include "NamespaceFooter.H"

#ifdef CRD_USE_VEX

/// Compute the solution to the Riemann problem
void riemann_vex(FArrayBox&       a_Wstar,
                 const FArrayBox& a_Wleft,
                 const FArrayBox& a_Wright,
                 const int        a_dir,
                 const Box&       a_box);

#endif  /* defined CRD_USE_VEX */
#endif  /* ! defined _POLYTROPICPHYSICS_VEX_H_ */
