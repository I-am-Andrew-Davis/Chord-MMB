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
 * \file VEx.H
 *
 * \brief Support for vector extensions
 *
 *//*+*************************************************************************/

#ifndef _VEX_H_
#define _VEX_H_

//----- Chombo Library -----//

#include "CH_config.H"
#include "UsingNamespace.H"

//--Check if we should use vector operations

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
#define CRD_USE_VEX
#endif

#endif  /* ! defined _VEX_H_ */
