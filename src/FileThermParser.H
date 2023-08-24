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
 * \file FileThermParser.H
 *
 * \brief Parser for thermodynamic, transport, and reaction data
 *
 *//*+*************************************************************************/

#ifndef _FILETHERMPARSER_H_
#define _FILETHERMPARSER_H_

//----- Chombo Library -----//

#include "CH_HDF5.H"
#include "RefCountedPtr.H"

//----- Internal -----//

#include "CRDparam.H"
#include "ThermPhysics.H"

/*******************************************************************************
 */
///  Parser for thermodynamic data, transport data, and reaction mechanism
/**
 *   This class defaults to reading NASA specified thermodynamic and transport
 *   data and Chemkin format reaction mechanism. This also fills lookup tables
 *   for relevant data
 *
 *//*+*************************************************************************/

class FileThermParser
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  FileThermParser(ThermPhysics& a_thPhys)
    :
    m_thPhys(a_thPhys)
    { }

  /// Destructor
  virtual ~FileThermParser()
    { }

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  FileThermParser(const FileThermParser&);

//--Assignment not permitted

  FileThermParser &operator=(const FileThermParser&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--General lookup table functionality

  /// Set lookup table values
  void setLookupTableVals();

//--Functions for thermodynamic data

  /// Read the thermodynamic file
  virtual void readThermFile();

  /// Calculate (\Delta_f H)^0 - (\Delta_f H)^298.15 for a given species
  Real solveDelHF(const std::string&              a_speciesName,
                  const std::vector<std::string>& a_elementList,
                  const std::vector<std::string>& a_refStates,
                  const std::vector<int>&         a_checkEl,
                  const std::vector<Real>&        a_numEinR,
                  const std::vector<Real>&        a_refHVals,
                  const int&                      a_curSpec) const;

  /// Fill thermodynamic lookup tables
  virtual void fillThermoTables();

//--Functions for reaction data

  /// Read the reaction file
  void readReactionFile();

  /// Count the different types and number of each type of reactions
  void countNumReactions();

  /// Size the vectors for reaction values
  void sizeReactionVectors();

  /// Extract all coefficients from the reaction file(s)
  void extractAllReactCoeffs(std::vector<Real>& a_preAF,
                             std::vector<Real>& a_betai,
                             std::vector<Real>& a_EAR,
                             std::vector<Real>& a_REVpreAB,
                             std::vector<Real>& a_REVbetai,
                             std::vector<Real>& a_REVEAR);

  /// Extract the species from the string and assign coeff values
  void extractSpecies(std::string& a_str,
                      const int&   a_ProdReact,
                      const int&   a_cRN);

  /// Extract the species and values for the third body coeffs
  void extractTBCoeff(std::string&  a_str,
                      const int&    a_tbRN);

  /// Extract arbitrary order reaction values
  void extractARBCoeff(std::string&  a_str,
                       const int&    a_cRN,
                       const int&    a_arbRN,
                       const bool    a_FORD,
                       const bool    a_firstArb);

  /// Extract the PR species, ie (+H2O)
  void extractPRSpeciesAlpha(const std::string a_str,
                             const int&        a_tbRN);

  /// Extract the coeff for Troe methodology
  void extractPRTroeValues(std::string&  a_str,
                           const int&    a_prRN);

  /// Assign lookup tables for K_fwd and K_bkwd
  virtual void fillReactTables(std::vector<Real>& a_preAF,
                               std::vector<Real>& a_betai,
                               std::vector<Real>& a_EAR,
                               std::vector<Real>& a_REVpreAB,
                               std::vector<Real>& a_REVbetai,
                               std::vector<Real>& a_REVEAR);

//--Functions for transport data

  /// Read transport file
  void readTransportFile(std::vector<std::array<Real, 4> >& a_muVecL,
                         std::vector<std::array<Real, 4> >& a_muVecH,
                         std::vector<std::array<Real, 4> >& a_kVecL,
                         std::vector<std::array<Real, 4> >& a_kVecH);
  
  /// Fill transport lookup tables
  void fillTransportTables(std::vector<std::array<Real, 4> >& a_muVecL,
                           std::vector<std::array<Real, 4> >& a_muVecH,
                           std::vector<std::array<Real, 4> >& a_kVecL,
                           std::vector<std::array<Real, 4> >& a_kVecH);

//--General extraction and string modifying functions

  /// Determines what separator is used
  size_t findSeparator(const std::string& a_str) const;

  /// Finds and replaces in a string
  void findReplace(std::string&       a_str,
                   const std::string& a_Rstr,
                   const std::string& a_Nstr) const;

  /// Extracts values from a string delimited by a_dl
  Real extractValue(std::string&      a_str,
                    const std::string a_dl) const;

  /// Extracts values from a string where values are of length a_nc
  Real extractValue(std::string& a_str,
                    const int    a_nc) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  ThermPhysics& m_thPhys;             ///< Reference to ThermPhysics

};

/*--------------------------------------------------------------------*/
//  Determines what the separator is used in reaction mechanism files
/** \param[out] loc     Size of the separator
 *  \param[in]  a_str   String to be read
 *//*-----------------------------------------------------------------*/

inline size_t
FileThermParser::findSeparator(const std::string& a_str) const
{
  size_t loc = 0;
  int count = 0;
  char c(' ');
  for (int i = 0; i != a_str.length(); ++i)
    {
      if (a_str[i] == c)
        {
          ++count;
        }
      if (count >= 3)
        {
          loc = i;
          break;
        }
    }
  return loc;
}

/*--------------------------------------------------------------------*/
//  Finds and replaces in a string
/** \param[out] a_str   Modified string
 *  \param[in]  a_str   String to be read
 *  \param[in]  a_Rstr  String of what is to be replaced
 *  \param[in]  a_Nstr  String to replace it with
 *//*-----------------------------------------------------------------*/

inline void
FileThermParser::findReplace(std::string&       a_str,
                             const std::string& a_Rstr,
                             const std::string& a_Nstr) const
{
  std::string::size_type pos = 0u;
  while ((pos = a_str.find(a_Rstr, pos)) != std::string::npos)
    {
      a_str.replace(pos, a_Rstr.length(), a_Nstr);
      pos += a_Nstr.length();
    }
}

/*--------------------------------------------------------------------*/
//  Extracts the values from the string delimited by a_dl
/** \param[in]  a_str   String of values
 *  \param[in]  a_dl    The delimiter
 *//*-----------------------------------------------------------------*/

inline Real
FileThermParser::extractValue(std::string&      a_str,
                              const std::string a_dl) const
{
  size_t locaif = a_str.find_first_of("0123456789-");
  size_t locail = a_str.find(a_dl,locaif);
  if (locail == std::string::npos)
    {
      locail = a_str.length();
    }
  std::string freqFact = a_str.substr(locaif, locail);
  Real outputVal = std::atof(freqFact.c_str());
  // Remove the value from the a_str
  a_str = a_str.substr(locail, a_str.length());
  return outputVal;
}

/*--------------------------------------------------------------------*/
//  Extracts the values from the string where values are of length a_nc
/** \param[in]  a_str   String of values
 *  \param[in]  a_nc    Number of chars the value is
 *//*-----------------------------------------------------------------*/

inline Real
FileThermParser::extractValue(std::string& a_str,
                              const int    a_nc) const
{
  std::string freqFact = a_str.substr(0, a_nc);
  Real outputVal = std::atof(freqFact.c_str());
  // Remove the value from the a_str
  a_str = a_str.substr(a_nc, a_str.length());
  return outputVal;
}

#endif  /* ! defined _FILETHERMPARSER_H_ */
