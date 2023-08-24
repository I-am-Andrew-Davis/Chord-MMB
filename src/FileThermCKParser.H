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
 * \file FileThermCKParser.H
 *
 * \brief Parser for thermodynamic, transport, and reaction data
 *
 *//*+*************************************************************************/

#ifndef _FILETHERMCKPARSER_H_
#define _FILETHERMCKPARSER_H_

//----- Chombo Library -----//

#include "CH_HDF5.H"
#include "RefCountedPtr.H"

//----- Internal -----//

#include "CRDparam.H"
#include "FileThermParser.H"

/*******************************************************************************
 */
///  Parser for thermodynamic data, transport data, and reaction mechanism
/**
 *   This class allows Chemkin thermodynamic data to be read in
 *
 *//*+*************************************************************************/

class FileThermCKParser : public FileThermParser
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  FileThermCKParser(ThermPhysics& a_thPhys)
    :
    FileThermParser(a_thPhys)
    { }

  /// Destructor
  virtual ~FileThermCKParser()
    { }

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  FileThermCKParser(const FileThermCKParser&);

//--Assignment not permitted

  FileThermCKParser &operator=(const FileThermCKParser&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--Functions for thermodynamic data

  /// Read the thermodynamic file
  virtual void readThermFile();

  /// Fill thermodynamic lookup tables
  virtual void fillThermoTables();

//--Functions for reaction data

  /// Assign lookup tables for K_fwd and K_bkwd
  virtual void fillReactTables(std::vector<Real>& a_preAF,
                               std::vector<Real>& a_betai,
                               std::vector<Real>& a_EAR,
                               std::vector<Real>& a_REVpreAB,
                               std::vector<Real>& a_REVbetai,
                               std::vector<Real>& a_REVEAR);

};

#endif  /* ! defined _FILETHERMCKPARSER_H_ */
