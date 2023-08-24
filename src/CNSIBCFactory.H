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
 * \file CNSIBCFactory.H
 *
 * \brief Abstract base class for constructing IBC objects
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCFACTORY_H_
#define _CNSIBCFACTORY_H_

//--Forward declarations

class CNSIBC;


/*******************************************************************************
 */
/// Factory for a new CNSIBC object
/** This class serves as abstract base for a factory that constructs a
 *  CNSIBC object.  It is only used when the Chord program is configured as
 *  a library and a custom CNSIBC object is required.  The construction must
 *  be defered until ParmParse has been initialized (in chord.cpp) which is the
 *  reason for this class.
 *
 *//*+*************************************************************************/

class CNSIBCFactory
{
public:

  virtual ~CNSIBCFactory() = default;

  /// Make a CNSIBC object
  virtual CNSIBC* new_CNSIBC() const = 0;
};
  
#endif /* ! defined _CNSIBCFACTORY_H_ */
