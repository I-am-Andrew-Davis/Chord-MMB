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
 * \file chord.H
 *
 * \brief Prototype if Chord is built as a library
 *
 *//*+*************************************************************************/

/* Normally the interface to Chord is the main function and this prototype is
   meaningless.  However Chord can also be built as a library which is often
   done for tests.  If so, this is the prototype to the function.  The test
   or caller must define main.
*/

#ifndef _CHORD_H_
#define _CHORD_H_

//----- Standard Library -----//

#include <functional>

//--Forward declarations from Chombo

#include "BaseNamespaceHeader.H"
template <class T>
class Vector;
#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"
class Box;
class AMRLevel;
#include "NamespaceFooter.H"

//--Forward declarations

class CNSIBCFactory;
class CRDPhysics;

/// Interface to the Chord program
int chord(int a_argc, const char* a_argv[],
          CNSIBCFactory *const a_ibcFactory = nullptr,
          CRDPhysics *const    a_physics = nullptr,
          std::function<Vector<Vector<Box>>(const AMR&)> *const
                               a_getFixedHierarchy = nullptr);

#endif  /* ! defined _CHORD_H_ */
