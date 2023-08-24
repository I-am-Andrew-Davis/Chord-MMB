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
 * \file CRDhardconf.H
 *
 * \brief Hard-coded configurations options
 *
 * These are only indended for development and research purposes.  Users of
 * Chord should not have to adjust these parameters.  As well, they should be
 * temporary in the sense that there should be some plan to select a preferred
 * approach and remove the configuration parameter.
 *
 *//*+*************************************************************************/

#ifndef _CRDHARDCONF_H_
#define _CRDHARDCONF_H_

// Defined for testing mapped viscous operations.  The mapped viscous terms
// require a massive number of ghost cells in state and metrics.  Instead, we
// break up the operations and use an exchange.  Whereas the Cartesian
// case is local to a patch, the mapped approach requires data stored in
// LevelData and exchanges.  We are keeping the Cartesian approach for
// validation until much more is know about how to implement the mapped approach
// 0 - Cartesian approach to computing viscous flux
// 1 - Mapped approach to computing viscous flux
#define TEST_USEMAPPEDVISC 1

#endif /* ! defined _CRDHARDCONF_H_ */
