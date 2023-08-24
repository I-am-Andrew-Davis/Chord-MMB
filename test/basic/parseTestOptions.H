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
 * \file CRDutil.H
 *
 * \brief Utilities for Chord
 *
 *//*+*************************************************************************/

#ifndef _PARSETESTOPTIONS_H_
#define _PARSETESTOPTIONS_H_

#include <cstring>

inline void
parseTestOptions(int argc, const char* argv[], bool& a_verbose)
{
  for (int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp(argv[i], "-v", 3) == 0)
            {
              a_verbose = true ;
            }
          else if (strncmp(argv[i], "-q", 3) == 0)
            {
              a_verbose = false ;
            }
          else
            {
              break ;
            }
        }
    }
}

#endif /* ! defined _PARSETESTOPTIONS_H_ */
