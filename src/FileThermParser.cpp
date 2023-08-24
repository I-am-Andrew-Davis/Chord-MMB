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
 * \file FileThermParser.cpp
 *
 * \brief Member functions for FileThermParser
 *
 *//*+*************************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <array>

//----- Chombo -----//

#include "CH_System.H"
#include "ParmParse.H"

//----- Internal -----//

#include "FileThermParser.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "DataTemp.H"

/*******************************************************************************
 *
 * Class FileThermParser: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set values related to lookup tables
/** 
 *//*-----------------------------------------------------------------*/

void
FileThermParser::setLookupTableVals()
{
  m_thPhys.m_lookupLoT = 250.;
  m_thPhys.m_lookupHiT = 5600.;
  m_thPhys.m_midLookupT = 1000.;
  ParmParse ppTHERM("therm");
  // Set stencil for using Lagrange polynomial interpolation
  if (ppTHERM.contains("interp_points"))
    {
      // Number of points to use during interpolation
      m_thPhys.m_interpN = ppTHERM.countval("interp_points");
      m_thPhys.m_diffCL.assign(m_thPhys.m_interpN,0);
      ppTHERM.getarr("interp_points",m_thPhys.m_diffCL, 0, m_thPhys.m_interpN);
    }
  else
    {
      m_thPhys.m_interpN = 5;
      m_thPhys.m_lookupDelT = 1.;
      m_thPhys.m_diffCL.assign(m_thPhys.m_interpN,0);
      // Set the points to use for interpolation, +/- the closest point
      m_thPhys.m_diffCL[0] = 0; // Use the closest point, point 0
      m_thPhys.m_diffCL[1] = 1; // Use next point ahead of point 0
      m_thPhys.m_diffCL[2] = 2; // Use point two ahead of point 0
      m_thPhys.m_diffCL[3] = -1; // Use point behind point 0
      m_thPhys.m_diffCL[4] = -2; // Use point two behind point 0
    }
  ppTHERM.query("interp_lookup_delt",m_thPhys.m_lookupDelT);
  m_thPhys.m_interpDenoms.assign(m_thPhys.m_interpN,0.);
  m_thPhys.m_lookupSize = int((m_thPhys.m_lookupHiT - m_thPhys.m_lookupLoT)/
                              m_thPhys.m_lookupDelT) + 1; 
  // Set the denominators for the Lagragian points
  for (int n = 0; n != m_thPhys.m_interpN; ++n)
    {
      m_thPhys.m_interpDenoms[n] = 1.;
      for (int k = 0; k != m_thPhys.m_interpN; ++k)
        {
          if (n != k)
            {
              m_thPhys.m_interpDenoms[n] *= m_thPhys.m_lookupDelT*
                (m_thPhys.m_diffCL[n] - m_thPhys.m_diffCL[k]);
            }
        }
      m_thPhys.m_interpDenoms[n] = 1./m_thPhys.m_interpDenoms[n];
    }
}

/*--------------------------------------------------------------------*/
//  Read in the therm data from a file
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::readThermFile()
{
  if (!CH_System::fileExists(m_thPhys.m_thermoFile.c_str()))
    {
      CRD::msg << "Input: 'thermo_file' " << m_thPhys.m_thermoFile
               << " does not exist!" << CRD::error;
    }
  const int numSpecies = CRDparam::g_numSpecies;
  std::vector<std::string> speciesNames(numSpecies, "U");
  std::vector<std::string> altSpecNames(numSpecies, "U");
  // FIXME: The element and reference lists should be comprehensive
  // List of elements that could exist in the species
  std::vector<std::string> elementList = {"C", "H", "O", "N", "Ar", "He"};
  // List of the corresponding reference states for each element by
  // name to be looked up in reference file
  std::vector<std::string> refStates = {"C(gr)", "H2", "O2", "N2", "Ar", "He"};
  // Vector that contains the # elements in ref state
  std::vector<Real> numEinR = {1., 2., 2., 2., 1., 1.};
  CH_assert(refStates.size() == elementList.size());
  CH_assert(numEinR.size() == refStates.size());
  const int numRefElements = refStates.size();
  // Vector to show if the reference element is in the list of species
  // If -1, element is not in species,
  // -2, if element is in species list but reference state is not
  // #, the species # that corresponds to the reference state
  std::vector<int> checkEl(numRefElements, -1);
  for (int i = 0; i != numSpecies; ++i)
    {
      speciesNames[i] = CRDparam::g_speciesNames[i];
      altSpecNames[i] = m_thPhys.m_altSpecNames[i];
      for (int j = 0; j != numRefElements; ++j)
        {
          // Check if the species is a reference state
          if (speciesNames[i] == refStates[j] ||
              altSpecNames[i] == refStates[j])
            {
              checkEl[j] = i;
            }
          // Check if the species contains one of the elements
          else if (speciesNames[i].find(elementList[j]) != std::string::npos &&
                   altSpecNames[i].find(elementList[j]) != std::string::npos &&
                   checkEl[j] == -1)
            {
              checkEl[j] = -2;
            }
        }
    }
  // Add the reference states to the species list
  int k = 0;
  for (int j = 0; j != numRefElements; ++j)
    {
      if (checkEl[j] == -2)
        {
          speciesNames.push_back(refStates[j]);
          altSpecNames.push_back(refStates[j]);
          checkEl[j] = numSpecies + k++;
        }
    }
  // Retrieve the electron information 
  speciesNames.push_back("e-");
  altSpecNames.push_back("e-");
  // Number of species and number of reference states
  const int numTotalSpecies = speciesNames.size();
  // Reference enthalpies for all species
  std::vector<Real> refEnthalpy(numTotalSpecies, 0.);
  ifstream readFile(m_thPhys.m_thermoFile);
  // Molar mass location and number of chars
  const int mmLoc = 55;
  const int mmLength = 10;
  // The number of chars in each value we are extracting is 16
  const int valLength = 16;
  // Set up here to prevent compiler warning
  Real dummyspot = 0.;
  if (readFile.is_open())
    {
      // The current number corresponding to the species
      int curSpec = -1;
      int curLine = -1;
      // Bool for if the species is not a reference species
      bool solSpec = false;
      std::string readLine;
      while (std::getline(readFile, readLine))
        {
          // Skip commented out lines
          if (readLine.find("!") <= 3)
            {
              continue;
            }
          if (curSpec == -1 && readLine.find(" ") == 0)
            {
              continue;
            }
          std::string specstring = readLine.substr(0,readLine.find(" ",1));
          for (int i = 0; i != numTotalSpecies; ++i)
            {
              if (specstring == speciesNames[i] ||
                  specstring == altSpecNames[i])
                {
                  curSpec = i;
                  curLine = 0;
                  solSpec = false;
                  if (curSpec < numSpecies)
                    {
                      solSpec = true;
                    }
                  break;
                }
            }
          if (curSpec == -1)
            {
              continue;
            }
          // Extract the molar mass for non-reference species
          if (curLine == 1 && solSpec)
            {
              readLine = readLine.substr(mmLoc, readLine.length());
              m_thPhys.m_molMass[curSpec] = extractValue(readLine,
                                                         mmLength)/1000.;
            }
          // Extract H(298.15) - H(0) for all species
          else if (curLine == 2)
            {
              std::string hform = readLine.substr(
                readLine.find_last_of(" ")+1,readLine.size());
              refEnthalpy[curSpec] = std::atof(hform.c_str());
              if (!solSpec)
                {
                  curSpec = -1;
                  continue;
                }
            }
          // Replace any D notation with E
          findReplace(readLine, "D", "E");
          // Extract the first portion of the high enthalpy coefficients
          if (curLine == 3)
            {
              m_thPhys.m_hnA1L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA2L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA3L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA4L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA5L[curSpec] = extractValue(readLine,valLength);
            }
          else if (curLine == 4)
            {
              m_thPhys.m_hnA6L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA7L[curSpec] = extractValue(readLine,valLength);
              dummyspot += extractValue(readLine,valLength);
              m_thPhys.m_hnB1L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnB2L[curSpec] = extractValue(readLine,valLength);
            }
          else if (curLine == 6)
            {
              m_thPhys.m_hnA1H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA2H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA3H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA4H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA5H[curSpec] = extractValue(readLine,valLength);
            }
          else if (curLine == 7)
            {
              m_thPhys.m_hnA6H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA7H[curSpec] = extractValue(readLine,valLength);
              dummyspot += extractValue(readLine,valLength);
              m_thPhys.m_hnB1H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnB2H[curSpec] = extractValue(readLine,valLength);
              curSpec = -1;
            }
          curLine++;
        }
    }
  readFile.close();
  // Check if any species were not found
  for (int i = 0; i != numSpecies; ++i)
    {
      // Solve for \Delta_f H(0) - \Delta_f H(298) + (H(0) - H(298))
      m_thPhys.m_H0[i] = 0.;
      m_thPhys.m_Rn[i] = m_thPhys.m_Rmol/m_thPhys.m_molMass[i];
      if (m_thPhys.m_hnA1L[i] == -1.)
        {
          CRD::msg << "FileThermParser::readThermFile: No inputs found for "
                   << CRDparam::g_speciesNames[i] << "!" << CRD::error;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Calculate the difference in heat of formation from 298.15 to 0 using
//  the elements in each species relative to the reference species.
/** \param[in]  a_speciesName
 *                      String containing the species name
 *  \param[in]  a_elementList
 *                      Vector of strings of element names
 *  \param[in]  a_refStates
 *                      Vector of strings of names of reference states 
 *                      for elements. IE H2 is ref state for H
 *  \param[in]  a_checkEl
 *                      Vector of component numbers of ref states for elements
 *  \param[in]  a_numEinR
 *                      Vector of number of elements in ref state
 *                      IE N2 has 2 N's while C(gr) has 1 C
 *  \param[in]  a_refHVals
 *                      Extracted values of H^0_n = (H(0) - H(298.15))_n
 *  \param[in]  a_curSpec
 *                      Component number of current species
 *  \return     hval    Returns \Delta_f^0 H(0) - \Delta_f^0 H(298.15) + 
 *                      (H(0) - H(298.15)) = \sum \nu_i H_i = H_n
 *  A few examples:
 *  Species - HCO, Elements - [H,C,O], Refs States - [H2,C(gr),O2]
 *  H^0_{HCO} = 1/2 H_{H2} + 1 H_{C(gr)} + 1/2 H_{O2}
 *  Species - CO2, Elements - [C,O], Refs States - [C(gr),O2]
 *  H^0_{CO2} = 1 H_{C(gr)} + 1 H_{O2}
 *  Species - OH+, Elements - [H,O], Refs States - [H2,O2,e-]
 *  H^0_{OH+} = 1/2 H_{H2} + 1/2 H_{O2} - H_{e-}
 *//*-----------------------------------------------------------------*/

Real
FileThermParser::solveDelHF(const std::string&              a_speciesName,
                            const std::vector<std::string>& a_elementList,
                            const std::vector<std::string>& a_refStates,
                            const std::vector<int>&         a_checkEl,
                            const std::vector<Real>&        a_numEinR,
                            const std::vector<Real>&        a_refHVals,
                            const int&                      a_curSpec) const
{
  const int numRefSpecies = a_elementList.size();
  const int numTotalSpecies = a_refHVals.size();
  // Final component of a_refHVals is the enthalpy for e-
  const int compEM = numTotalSpecies - 1;
  const Real EMHVal = a_refHVals[compEM];
  Real hval = 0.;
  // Loop over all reference species
  for (int j = 0; j != numRefSpecies; ++j)
    {
      // Component in a_refHVals
      int specRef = a_checkEl[j];
      if (specRef == -1)
        {
          continue;
        }
      // Number of elements in species
      int numEl = 0;
      int sizeEl = a_elementList[j].length();
      for (int i = 0; i != a_speciesName.length(); ++i)
        {
          std::string cutst = a_speciesName.substr(i,sizeEl);
          if (cutst == a_elementList[j])
            {
              numEl++;
              int diffLen = a_speciesName.length() - (i+sizeEl);
              if (diffLen > 0)
                {
                  char nc = a_speciesName.at(i+sizeEl);
                  // Do not count if the next character is lowercase
                  if (std::islower(nc))
                    {
                      numEl--;
                      continue;
                    }
                  else if (std::isupper(nc))
                    {
                      continue;
                    }
                  else if (std::isdigit(nc))
                    {
                      numEl = nc - '0';
                      // Check if the # of elements is double digits
                      if (diffLen > 1)
                        {
                          char nc2 = a_speciesName.at(i+sizeEl+1);
                          if (std::isdigit(nc2))
                            {
                              numEl = numEl*10 + (nc2 - '0');
                            }
                        }
                    }
                }
            }
        }
      Real ratio = numEl/a_numEinR[j];
      hval += a_refHVals[specRef]*ratio;
    }
  if (a_speciesName.find("+") != std::string::npos &&
      a_speciesName.find(",") == std::string::npos)
    {
      hval -= EMHVal;
    }
  else if (a_speciesName.find("-") != std::string::npos &&
           a_speciesName.find(",") == std::string::npos)
    {
      hval += EMHVal;
    }
  return hval;
}

/*--------------------------------------------------------------------*/
//  Set the values in the lookup tables for thermodynamic properties
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::fillThermoTables()
{
  const int nS = CRDparam::g_numSpecies;
  const Real delT = m_thPhys.m_lookupDelT;
  const Real loT = m_thPhys.m_lookupLoT;
  std::vector<Real> initVec(nS,0.);
  m_thPhys.m_lookupH.assign(nS,initVec);
  for (int sp = 0; sp != nS; ++sp)
    {
      m_thPhys.m_lookupH[sp].assign(m_thPhys.m_lookupSize, 0.);
      for (int i = 0; i != m_thPhys.m_lookupSize; ++i)
        {
          const Real T = loT + (Real)(i)*delT;
          Real logT = std::log(T);
          if (T >= m_thPhys.m_midLookupT)
            {
              m_thPhys.m_lookupH[sp][i] = m_thPhys.m_Rn[sp]*
                (m_thPhys.m_hnB1H[sp] - m_thPhys.m_hnA1H[sp]/T +
                 m_thPhys.m_hnA2H[sp]*logT +
                 T*(m_thPhys.m_hnA3H[sp] +
                    T*(m_thPhys.m_hnA4H[sp]/2. + T*
                       (m_thPhys.m_hnA5H[sp]/3. +
                        T*(m_thPhys.m_hnA6H[sp]/4.
                           + T*m_thPhys.m_hnA7H[sp]/5.)))));
            }
          else
            {
              m_thPhys.m_lookupH[sp][i] = m_thPhys.m_Rn[sp]*
                (m_thPhys.m_hnB1L[sp] - m_thPhys.m_hnA1L[sp]/T +
                 m_thPhys.m_hnA2L[sp]*logT +
                 T*(m_thPhys.m_hnA3L[sp] +
                    T*(m_thPhys.m_hnA4L[sp]/2. + T*
                       (m_thPhys.m_hnA5L[sp]/3. +
                        T*(m_thPhys.m_hnA6L[sp]/4.
                           + T*m_thPhys.m_hnA7L[sp]/5.)))));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Read the reaction input file
//  When reading the reaction file, be sure the reaction section begins
//  with the word REACTIONS. In the line containing REACTIONS, be sure
//  to include any unit modifications for A_r and E_a
/** 
 *//*-----------------------------------------------------------------*/

void
FileThermParser::readReactionFile()
{
  CRD::msg << CRD::fv3 << "FileThermParser::readReactionFile" << CRD::end;
  // Number of reactions
  m_thPhys.m_numReactions = 0;
  // Number of third body reactions
  m_thPhys.m_numTBReact = 0;
  // Number of pressure dependent reactions
  m_thPhys.m_numPRReact = 0;
  countNumReactions();
  // Size all the members vectors
  sizeReactionVectors();
  std::vector<Real> preAF(m_thPhys.m_numReactions);
  std::vector<Real> betai(m_thPhys.m_numReactions);
  std::vector<Real> EAR(m_thPhys.m_numReactions);
  std::vector<Real> REVpreAB(m_thPhys.m_numREVReact);
  std::vector<Real> REVbetai(m_thPhys.m_numREVReact);
  std::vector<Real> REVEAR(m_thPhys.m_numREVReact);
  extractAllReactCoeffs(preAF, betai, EAR, REVpreAB, REVbetai, REVEAR);
  fillReactTables(preAF, betai, EAR, REVpreAB, REVbetai, REVEAR);
}

/*--------------------------------------------------------------------*/
//  Count the number of reactions from the file(s)
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::countNumReactions()
{
  std::string reactLine("REACTIONS");
  std::string reactLine2("REAC");
  int rfLen = m_thPhys.m_reactionFile.length();
  char *filename = new char[rfLen + 1];
  std::strcpy(filename, m_thPhys.m_reactionFile.c_str());
  ifstream readFile(filename);
  delete [] filename;
  if (readFile.is_open())
    {
      std::string readLine;
      // Look for the word REACTIONS
      while (std::getline(readFile, readLine))
        {
          size_t loc, loc2;
          // Look for the word REACTIONS or REAC
          loc = readLine.find(reactLine);
          loc2 = readLine.find(reactLine2);
          if (loc != std::string::npos || loc2 != std::string::npos)
            {
              break;
            }
        }
      while (std::getline(readFile, readLine))
        {
          size_t fexl;
          bool arbReact = false; // True when reading an arbitrary order react
          fexl = readLine.find("!");
          // Remove any part that follows an !
          readLine = readLine.substr(0,fexl);
          // Check if line starts with a ! and skip said line
          if (fexl != 0 && readLine.length() > 4)
            {
              // Check that Chord can handle all reaction mechanisms present
              if (readLine.find("LT") != std::string::npos ||
                  readLine.find("JAN") != std::string::npos ||
                  readLine.find("HV") != std::string::npos ||
                  readLine.find("TDEP") != std::string::npos ||
                  readLine.find("XSMI") != std::string::npos ||
                  readLine.find("EXCI") != std::string::npos ||
                  readLine.find("MOME") != std::string::npos)
                {
                  CRD::msg << "Reaction #" << m_thPhys.m_numReactions
                           << " is not acceptable for the current "
                           << "implementation of Chord" << CRD::error;
                }
              // Ensure line is not TB or PR coeffs
              if (readLine.find("/") == std::string::npos &&
                  readLine.find("=") !=  std::string::npos)
                {
                  ++m_thPhys.m_numReactions;
                  arbReact = false;
                  size_t prc = readLine.find("(+");
                  size_t tbc = readLine.find("M");
                  // Check if eq is pressure-dependent
                  if (prc != std::string::npos)
                    {
                      ++m_thPhys.m_numPRReact;
                      ++m_thPhys.m_numTBReact;
                    }
                  // If the reaction has an M, then it is third body
                  else if (tbc != std::string::npos)
                    {
                      ++m_thPhys.m_numTBReact;
                    }
                }
              if (readLine.find("REV") != std::string::npos)
                {
                  ++m_thPhys.m_numREVReact;
                }
              if ((readLine.find("FORD") != std::string::npos ||
                   readLine.find("RORD") != std::string::npos) && !arbReact)
                {
                  ++m_thPhys.m_numArbReact;
                  arbReact = true;
                }
            }
        }
      readFile.close();
    }
  if (m_thPhys.m_numReactions != CRDparam::g_numReactions)
    {
      CRD::msg << "FileThermParser: Number of reactions "
               << "specified in input file is incorrect. There are "
               << m_thPhys.m_numReactions << " reactions in reaction file."
               << CRD::error;
    }
}

/*--------------------------------------------------------------------*/
// Size all vectors for reaction data
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::sizeReactionVectors()
{
  // Number of species
  const int nS = CRDparam::g_numSpecies;
  // Reaction variables
  // Number of reactions
  const int nRCT = m_thPhys.m_numReactions;
  // Number of reactions with third body
  const int nRTB = m_thPhys.m_numTBReact;
  CH_assert(nRTB <= nRCT);
  // Number of pressure dependent reactions
  const int nRPR = m_thPhys.m_numPRReact;
  // Number of reactants for each reaction
  m_thPhys.m_numRects.assign(nRCT,1);
  // Number of products for each reaction
  m_thPhys.m_numProds.assign(nRCT,1);
  // Reversible or irreversible reaction
  m_thPhys.m_revReact.assign(nRCT,1);
  // The reference vectors tell us which species are associated with which
  // reaction.
  // Index #: reaction #, 1st, 2nd, or 3rd species
  // 0: 1st reaction, 1st species
  // 1: 1st reaction, 2nd species
  // 2: 1st reaction, 3rd species
  // 3: 2nd reaction, 1st species
  // 4: 2nd reaction, 2nd species
  // 5: 2nd reaction, 3rd species
  // etc.
  // IE if we had the reactions
  // H2+H2+O2=>H2O+H2O
  // OH+OH=>O+H2O
  // and the species were ordered
  // as 0-H2, 1-O2, 2-O, 3-OH, and 4-H2O; the reactants reference vector is:
  // m_refNup[0] = 0;
  // m_refNup[1] = 1;
  // m_refNup[2] = -1; This means there are only 2 reactant species
  // m_refNup[3] = 3;
  // m_refNup[4] = -1; This means there is only 1 reactant species
  // m_refNup[5] = -1;
  // The product reference vector is:
  // m_refNupp[0] = 4;
  // m_refNupp[1] = -1;
  // m_refNupp[2] = -1;
  // m_refNupp[3] = 2;
  // m_refNupp[4] = 4;
  // m_refNupp[5] = -1;
  // Number for reference vectors
  const int nRefRS = m_thPhys.m_maxSpecPerReact*nRCT;
  m_thPhys.m_refNup.assign(nRefRS,-1);
  m_thPhys.m_refNupp.assign(nRefRS,-1);
  m_thPhys.m_nup.assign(nRefRS,0.);
  m_thPhys.m_nupp.assign(nRefRS,0.);
  // Vector to tell if the reaction is third body or not
  m_thPhys.m_thirdBody.assign(nRCT,0);
  // Size of vector for coeffs that vary based on TB reactions and species
  // Similarly indexed as reaction coeffs
  // Index #: species #, third body reaction #
  const int nRSTB = nRTB*nS;
  // Vector containing coefficients for third body reactions
  m_thPhys.m_tbAlpha.assign(nRSTB,1.);
  // Check if there are pressure-dependent reactions
  if (nRPR > 0)
    {
      // Vector to tell if reaction is pressure-dependent
      m_thPhys.m_presDep.assign(nRCT,0);
      // Vector of high or low reaction values EAR, beta and A
      m_thPhys.m_PREAR.assign(nRPR,0.);
      m_thPhys.m_PRbetai.assign(nRPR,0.);
      m_thPhys.m_PRpreAF.assign(nRPR,0.);

      m_thPhys.m_TroeAlpha.assign(nRPR,0.);
      m_thPhys.m_TroeTS3.assign(nRPR,0.);
      m_thPhys.m_TroeTS2.assign(nRPR,0.);
      m_thPhys.m_TroeTS.assign(nRPR,0.);

      m_thPhys.m_SRIa.assign(nRPR,0.);
      m_thPhys.m_SRIb.assign(nRPR,0.);
      m_thPhys.m_SRIc.assign(nRPR,0.);
      m_thPhys.m_SRId.assign(nRPR,1.);
      m_thPhys.m_SRIe.assign(nRPR,0.);
    }
  // Arbitrary order reactions
  const int nARB = m_thPhys.m_numArbReact;
  if (nARB > 0)
    {
      const int nRefArb = m_thPhys.m_maxSpecPerReact*nARB;
      m_thPhys.m_ford.assign(nRefArb, 0.);
      m_thPhys.m_rord.assign(nRefArb, 0.);
    }
}

/*--------------------------------------------------------------------*/
//  Extract all reaction coefficients
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::extractAllReactCoeffs(std::vector<Real>& a_preAF,
                                       std::vector<Real>& a_betai,
                                       std::vector<Real>& a_EAR,
                                       std::vector<Real>& a_REVpreAB,
                                       std::vector<Real>& a_REVbetai,
                                       std::vector<Real>& a_REVEAR)
{
  // Universal gas constant in units cal/(mol K) to match E_a
  Real avoNum = 6.0221409E23;
  Real R_ea = 1.9872036;
  const int numEaTypes = 4;
  std::vector<std::string> EaUnitTypes =
    {"KCAL/MOLE", "JOULES/MOLE", "KELVINS", "EVOLTS"};
  std::vector<Real> EaUnits = {R_ea/1000., m_thPhys.m_Rmol, 1.,
                               R_ea*3.8293E-20};
  // Units for A (Default is cm-s-K-mol)
  m_thPhys.m_preAmodUnits = 1.E-6;
  Real modA1 = avoNum;   // For when units are MOLECULES
  Real modA2 = 1./1000.; // For when units are KMOLES
  std::string reactLine("REACTIONS");
  std::string reactLine2("REAC");
  int rfLen = m_thPhys.m_reactionFile.length();
  char *filename = new char[rfLen + 1];
  std::strcpy(filename, m_thPhys.m_reactionFile.c_str());
  ifstream readFile(filename);
  delete [] filename;
  if (readFile.is_open())
    { 
      std::string readLine;
      while (std::getline(readFile, readLine))
        {
          size_t loc, loc2;
          // Look for the word REACTIONS or REAC
          loc = readLine.find(reactLine);
          loc2 = readLine.find(reactLine2);
          if (loc != std::string::npos || loc2 != std::string::npos)
            {
              // Check if units are given in m^3 instead of cm^3
              if (readLine.find("METER") != std::string::npos)
                {
                  m_thPhys.m_preAmodUnits = 1.;
                }
              // Check to see if units for A are specified
              if (readLine.find("MOLECULES") != std::string::npos)
                {
                  m_thPhys.m_preAmodUnits *= modA1;
                }
              else if (readLine.find("KMOLES") != std::string::npos)
                {
                  m_thPhys.m_preAmodUnits *= modA2;
                }
              // Check to see if units for E_a are specified
              for (int eaType = 0; eaType != numEaTypes; ++eaType)
                {
                  if (readLine.find(EaUnitTypes[eaType]) != std::string::npos)
                    {
                      R_ea = EaUnits[eaType];
                    }
                }
              break;
            }
        }
      // Current reaction number
      int cRN = 0;
      // Current third body reaction number
      // Starts at a negative index because it is increased before
      // being used
      int tbRN = -1;
      // Current pressure-dependent reaction number
      int prRN = -1;
      // Current REV reaction number
      int revRN = -1;
      // Current arbitrary order reaction number
      int arbRN = -1;
      // Check if still reading an arbitrary order reaction
      bool curArbReact = false;
      // Pressure-dependent line,
      // 0 - line to look for LOW or HIGH
      // 1 - line possibly containing type, either TROE or SRI
      // 2,3,4... - Lines containing alpha coefficients
      int prLine = -1;
      while (std::getline(readFile, readLine))
        {
          size_t fexl;
          fexl = readLine.find("!");
          // Remove any part that follows an !
          readLine = readLine.substr(0,fexl);
          // Check if line starts with a ! and skip said line
          if (fexl != 0 && readLine.length() > 4)
            {
              // Check to ensure it isn't a TB or PR coeff line
              if (readLine.find("/") == std::string::npos &&
                  readLine.find("=") != std::string::npos)
                {
                  // Split into two strings, one with equation, one with
                  // coeffs
                  std::string coeffstring;
                  size_t loc = readLine.find('\t');
                  // If tabs aren't used, look for 3 consecutive spaces
                  if (loc == std::string::npos)
                    {
                      loc = findSeparator(readLine);
                    }
                  coeffstring = readLine.substr(loc+1, readLine.length());
                  // Make a string of the just the stoichiometric eqn
                  std::string eqstring = readLine.substr(0,loc);
                  // Remove white space from string
                  eqstring.erase(
                    std::remove_if(eqstring.begin(),
                                   eqstring.end(),
                                   [](char c)
                                   {
                                     return std::isspace(
                                       static_cast<unsigned char>(c));
                                   }),
                    eqstring.end());
                  // Separate the reactants and products
                  // Location of the equal sign, either <=>, =>, or =
                  size_t loceq = 0;
                  // Length of the equal sign
                  int eqlength = 0;
                  if (eqstring.find("<") != std::string::npos &&
                      eqstring.find(">") != std::string::npos)
                    {
                      loceq = eqstring.find("<");
                      eqlength = 3;
                      m_thPhys.m_revReact[cRN] = 1;
                    }
                  else if (eqstring.find(">") != std::string::npos)
                    {
                      loceq = eqstring.find("=");
                      eqlength = 2;
                      // Set this reaction as irreversible
                      m_thPhys.m_revReact[cRN] = -1;
                    }
                  else if (eqstring.find("=") != std::string::npos)
                    {
                      loceq = eqstring.find("=");
                      eqlength = 1;
                      m_thPhys.m_revReact[cRN] = 1;
                    }
                  // String of reactants
                  std::string rside = eqstring.substr(0,loceq);
                  // String of products
                  std::string pside =
                    eqstring.substr(loceq+eqlength, eqstring.length());
                  size_t prcR = rside.find("(+");
                  size_t tbc = rside.find("M");
                  // Check if eq is pressure-dependent and third body
                  if (prcR != std::string::npos && tbc != std::string::npos)
                    {
                      m_thPhys.m_presDep[cRN] = 1;
                      ++prRN;
                      // Reset prLine
                      prLine = 0;
                      m_thPhys.m_thirdBody[cRN] = 1;
                      // Find the length of the species in parenthesis
                      size_t prdR = rside.find(")") - (prcR+2);
                      // Extract the pressure-dependent species
                      std::string prstring = rside.substr(prcR+2,prdR);
                      size_t prcP = pside.find("(");
                      // Remove the pressure-dependent species from the
                      // products and reactants
                      rside = rside.substr(0, prcR);
                      pside = pside.substr(0, prcP);
                      ++tbRN;
                    }
                  // Check if eq is just third body
                  else if (tbc != std::string::npos)
                    {
                      // Remove the third-body species from the
                      // products and reactants
                      rside = rside.substr(0, rside.size()-2);
                      pside = pside.substr(0, pside.size()-2);
                      m_thPhys.m_thirdBody[cRN] = 1;
                      ++tbRN;
                    }
                  // Check if eq is pressure-dependent but not third body
                  else if (prcR != std::string::npos)
                    {
                      m_thPhys.m_presDep[cRN] = 1;
                      m_thPhys.m_thirdBody[cRN] = -1;
                      // Find the length of the species in parenthesis
                      size_t prdR = rside.find(")") - (prcR+2);
                      // Extract the pressure-dependent species
                      std::string prstring = rside.substr(prcR+2,prdR);
                      size_t prcP = pside.find("(");
                      // Remove the pressure-dependent species from the
                      // products and reactants
                      rside = rside.substr(0, prcR);
                      pside = pside.substr(0, prcP);
                      ++prRN;
                      ++tbRN;
                      // Extract the species from the string and assign
                      // the alpha TB value to one for that species
                      extractPRSpeciesAlpha(prstring, tbRN);
                      // Reset prLine
                      prLine = 0;
                    }
                  // Extract the reaction species
                  extractSpecies(rside, -1, cRN);
                  // Extract the product species
                  extractSpecies(pside, 1, cRN);
                  // Extract the reaction coefficients
                  a_preAF[cRN] = extractValue(coeffstring," ");
                  Real beta = extractValue(coeffstring," ");
                  a_betai[cRN] = beta;
                  m_thPhys.m_betai[cRN] = beta;
                  Real EAR = extractValue(coeffstring," ")/R_ea;
                  a_EAR[cRN] = EAR;
                  m_thPhys.m_EAR[cRN] = EAR;
                  // Advance the current reaction number
                  ++cRN;
                  curArbReact = false;
                }
              // If slash is found, extract the PR coeffs from line
              else if (readLine.find("/") != std::string::npos)
                {
                  bool fordReact = false;
                  bool rordReact = false;
                  if (readLine.find("FORD") != std::string::npos)
                    {
                      fordReact = true;
                    }
                  else if (readLine.find("RORD") != std::string::npos)
                    {
                      rordReact = true;
                    }
                  // Check if line contains REV
                  if (readLine.find("REV") != std::string::npos)
                    {
                      CH_assert(m_thPhys.m_numREVReact > 0);
                      ++revRN;
                      // Set to -1 for REV reactions
                      m_thPhys.m_revReact[cRN-1] = 0;
                      std::string revcoefstr = readLine;
                      size_t locslash = revcoefstr.find("/");
                      size_t lastslash = revcoefstr.find_last_of("/")
                        - locslash;
                      // Make string only the coeffs
                      revcoefstr = revcoefstr.substr(locslash+1,
                                                     lastslash-1);
                      // Extract values into a_REVEAR, a_REVpreAB, a_REVbetai
                      a_REVpreAB[revRN] = extractValue(revcoefstr," ");
                      a_REVbetai[revRN] = extractValue(revcoefstr," ");
                      a_REVEAR[revRN] = extractValue(revcoefstr," ")/R_ea;
                    }
                  // Check if arbitrary reaction order reaction
                  else if (fordReact || rordReact)
                    {
                      // Since multiple lines following the reaction can be
                      // specified for FORD or RORD, care must be taken
                      bool firstArb = false;
                      if (!curArbReact)
                        {
                          curArbReact = true;
                          firstArb = true;
                          ++arbRN;
                          m_thPhys.m_revReact[cRN-1] *= 2;
                        }
                      std::string fordcoefstr = readLine;
                      size_t locslash = fordcoefstr.find("/");
                      size_t lastslash = fordcoefstr.find_last_of("/")
                        - locslash;
                      // Make string only the specs and coeffs
                      fordcoefstr = fordcoefstr.substr(locslash+1,
                                                       lastslash-1);
                      // Extract and assign coeff values
                      extractARBCoeff(fordcoefstr, cRN-1, arbRN, fordReact,
                                      firstArb);
                    }
                  // If the last equation was pressure dependent
                  else if (prLine != -1)
                    {
                      CH_assert(m_thPhys.m_presDep[cRN-1] != 0);
                      // Check for LOW or HIGH
                      // (unimolecular or bimolecular)
                      // and extract A, beta, and EAR
                      if (prLine == 0)
                        {
                          size_t lowPR;
                          lowPR = readLine.find("LOW");
                          // If low is found, presDep should be negative
                          if (lowPR != std::string::npos)
                            {
                              m_thPhys.m_presDep[cRN-1] = -1;
                            }
                          // Extract the PR A, beta, and EAR values
                          std::string prcoefstr = readLine;
                          size_t locslash = prcoefstr.find("/");
                          size_t lastslash = prcoefstr.find_last_of("/")
                            - locslash;
                          // Make string only the coeffs
                          prcoefstr = prcoefstr.substr(locslash+1,
                                                       lastslash-1);
                          // Extract values from string into m_PREAR,
                          // m_PRpreAF, and m_PRbetai
                          m_thPhys.m_PRpreAF[prRN] =
                            extractValue(prcoefstr," ");
                          m_thPhys.m_PRbetai[prRN] =
                            extractValue(prcoefstr," ");
                          m_thPhys.m_PREAR[prRN] =
                            extractValue(prcoefstr," ")/R_ea;
                          prLine++;
                        }
                      // Check for Troe method
                      else if (prLine == 1 &&
                               readLine.find("TROE") != std::string::npos)
                        {
                          // Assign value of 2 or -2 for Troe
                          m_thPhys.m_presDep[cRN-1] *= 2;
                          std::string prcoefstr = readLine;
                          size_t locslash = prcoefstr.find("/");
                          size_t lastslash = prcoefstr.find_last_of("/")
                            - locslash;
                          // Make string only the coeffs
                          prcoefstr = prcoefstr.substr(locslash+1, lastslash-1);
                          // Extract Troe values alpha, T*, T**, T***
                          m_thPhys.m_TroeAlpha[prRN] =
                            extractValue(prcoefstr," ");
                          m_thPhys.m_TroeTS3[prRN] =
                            extractValue(prcoefstr," ");
                          m_thPhys.m_TroeTS[prRN] = extractValue(prcoefstr," ");
                          // If T** is listed, set it to a negative value
                          size_t localif =
                            prcoefstr.find_first_of("0123456789-");
                          if (localif == std::string::npos)
                            {
                              m_thPhys.m_TroeTS2[prRN]= -100.;
                            }
                          else
                            {
                              m_thPhys.m_TroeTS2[prRN] =
                                extractValue(prcoefstr," ");
                            }
                          // Set to -1 to find third body values
                          // on next line
                          prLine = -1;
                        }
                      // Check for SRI method
                      else if (prLine == 1 &&
                               readLine.find("SRI") != std::string::npos)
                        {
                          m_thPhys.m_presDep[cRN-1] *= 3;
                          std::string prcoefstr = readLine;
                          size_t locslash = prcoefstr.find("/");
                          size_t lastslash = prcoefstr.find_last_of("/")
                            - locslash;
                          // Make string only the coeffs
                          prcoefstr = prcoefstr.substr(locslash+1,
                                                       lastslash-1);
                          // Extract SRI values for a, b, and c
                          m_thPhys.m_SRIa[prRN] = extractValue(prcoefstr," ");
                          m_thPhys.m_SRIb[prRN] = extractValue(prcoefstr," ");
                          m_thPhys.m_SRIc[prRN] = extractValue(prcoefstr," ");
                          // Check if d and e variables are listed
                          size_t localif =
                            prcoefstr.find_first_of("0123456789-");
                          if (localif != std::string::npos)
                            {
                              m_thPhys.m_SRId[prRN] =
                                extractValue(prcoefstr," ");
                              localif = prcoefstr.find_first_of("0123456789-");
                              if (localif != std::string::npos)
                                {
                                  m_thPhys.m_SRIe[prRN] =
                                    extractValue(prcoefstr," ");
                                }
                            }
                          prLine = -1;
                        }
                      // Otherwise it is the Lindemann method
                      else if (m_thPhys.m_thirdBody[cRN-1] == 1)
                        {
                          std::string tbcoeffstring = readLine;
                          // Remove white space from string
                          tbcoeffstring.erase(
                            std::remove_if(
                              tbcoeffstring.begin(),
                              tbcoeffstring.end(),
                              [](char c)
                              {
                                return std::isspace(
                                  static_cast<unsigned char>(c));
                              }),
                            tbcoeffstring.end());
                          // Extract and assign third body coeff values
                          extractTBCoeff(tbcoeffstring, tbRN);
                          prLine = -1;
                        }
                      else
                        {
                          prLine = -1;
                        }
                    }
                  else if (m_thPhys.m_numTBReact > 0)
                    {
                      // Check if last equation was third body
                      if (m_thPhys.m_thirdBody[cRN-1] == 1)
                        {
                          std::string tbcoeffstring = readLine;
                          // Remove white space from string
                          tbcoeffstring.erase(
                            std::remove_if(
                              tbcoeffstring.begin(),
                              tbcoeffstring.end(),
                              [](char c)
                              {
                                return std::isspace(
                                  static_cast<unsigned char>(c));
                              }),
                            tbcoeffstring.end());
                          // Extract and assign third body coeff values
                          extractTBCoeff(tbcoeffstring, tbRN);
                        }
                    }
                }
            }
        }
      readFile.close();
    }
}

/*--------------------------------------------------------------------*/
//  Extract the species from the given string and assign the values for
//  nu' or nu''
/** \param[in]  a_str   String to extract values
 *  \param[out] a_str   String with current species removed
 *  \param[in]  a_ProdReact
 *                      Set to -1 if the string is of the reactants and
 *                      is set to 1 if the string is of the products
 *  \param[in]  a_cRN   The number of the current reaction in the .inp file
 *//*-----------------------------------------------------------------*/

void
FileThermParser::extractSpecies(std::string& a_str,
                                const int&   a_ProdReact,
                                const int&   a_cRN)
{
  const int initLength = a_str.length() + 1;
  // Check to see if species are duplicated
  int dupRef = -1;
  // Shows what reference location the current reaction holds
  int refRNum = a_cRN*m_thPhys.m_maxSpecPerReact;
  for (int spn = 0; spn != initLength; spn++)
    {
      // FIXME: This will not work with ionized species
      // Make sure the species is listed in the input file
      int refSpec = -1;
      // Check if there are multiple moles of a species, like 2H2
      Real numMol = 1.;
      char fc = a_str.front();
      if (std::isdigit(fc))
        {
          std::string numMolStr;
          while (std::isdigit(a_str.front()) || a_str.front() == '.')
            {
              numMolStr += a_str.front();
              // Now remove the number from the string
              a_str = a_str.substr(1,a_str.length());
            }
          numMol = std::atof(numMolStr.c_str());
        }
      std::string spec;
      size_t loc = a_str.find("+");
      // Check if this is the last species in reaction side
      if (loc == std::string::npos)
        {
          spec = a_str;
          // End the loop
          spn = initLength - 1;
        }
      else
        {
          spec = a_str.substr(0,loc);
          // Remove the part of the string we have already read
          a_str = a_str.substr(loc + 1, a_str.length());
        }
      for (int i = 0; i != CRDparam::g_numSpecies; ++i)
        {
          if (spec == CRDparam::g_speciesNames[i] ||
              spec == m_thPhys.m_altSpecNames[i])
            {
              // If it is a new species for this side of the reaction
              if (dupRef != i && dupRef != -1)
                {
                  refRNum += 1;
                  if (a_ProdReact < 0)
                    {
                      m_thPhys.m_numRects[a_cRN] += 1;
                    }
                  else
                    {
                      m_thPhys.m_numProds[a_cRN] += 1;
                    }
                }
              // If we are reading in reactants
              if (a_ProdReact < 0)
                {
                  m_thPhys.m_refNup[refRNum] = i;
                  m_thPhys.m_nup[refRNum] += numMol;
                }
              // If we are reading in products
              else
                {
                  m_thPhys.m_refNupp[refRNum] = i;
                  m_thPhys.m_nupp[refRNum] += numMol;
                }
              refSpec = i;
              dupRef = i;
            }
        }
      if (refSpec == -1)
        {
          CRD::msg << "FileThermParser::extractSpecies: " << spec
                   << " found in reaction " << a_cRN
                   << " but not specified in inputs!" << CRD::error;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Extract the species and coeff values for the third body reactions
/** \param[in]  a_str   String to extract values
 *  \param[out] a_str   String with current species removed
 *  \param[in]  a_tbRN  The number of the current TB reaction
 *//*-----------------------------------------------------------------*/

void
FileThermParser::extractTBCoeff(std::string&  a_str,
                                const int&    a_tbRN)
{
  const int initLength = a_str.length() + 1;
  // Loop should end with the break, not when i == initLength
  for (int len = 0; len != initLength; len++)
    {
      std::string spec;
      std::string coeffstring;
      size_t loc = a_str.find("/");
      int refSpec = -1;
  
      spec = a_str.substr(0,loc);
      // Remove the part from the string we have already read
      a_str = a_str.substr(loc + 1, a_str.length());
      loc = a_str.find("/");
      coeffstring = a_str.substr(0,loc);
      const Real coeffVal = std::atof(coeffstring.c_str());
      for (int i = 0; i != CRDparam::g_numSpecies; ++i)
        {
          if (spec == CRDparam::g_speciesNames[i] ||
              spec == m_thPhys.m_altSpecNames[i])
            {
              // Determine the component based on the reaction and species
              const int rcomp = i*m_thPhys.m_numTBReact + a_tbRN;
              m_thPhys.m_tbAlpha[rcomp] = coeffVal;
              refSpec = i;
            }
        }
      if (refSpec == -1)
        {
          CRD::msg << "ThermPhysics::extractTBCoeff: " << spec << " found in "
                   << "reaction but not specified in inputs!" << CRD::error;
        }
      // Cut the string again
      a_str = a_str.substr(loc + 1, a_str.length());
      if (a_str.length() <= 1)
        {
          break;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Extract the species and coeff values for the third body reactions
/** \param[in]  a_str   String to extract values
 *  \param[out] a_str   Empty string
 *  \param[in]  a_cRN   Current reaction number
 *  \param[in]  a_arbRN Current arbitrary reaction number
 *  \param[in]  a_FORD  T - looking at FORD rates, F - looking at RORD rates
 *//*-----------------------------------------------------------------*/

void
FileThermParser::extractARBCoeff(std::string&  a_str,
                                 const int&    a_cRN,
                                 const int&    a_arbRN,
                                 const bool    a_FORD,
                                 const bool    a_firstArb)
{
  int refRNum = a_cRN*m_thPhys.m_maxSpecPerReact;
  int refARBN = a_arbRN*m_thPhys.m_maxSpecPerReact;
  // If this is the first arbitrary order reaction line for this reaction
  // initialize the values as nu values
  if (a_firstArb)
    {
      for (int i = 0; i != m_thPhys.m_maxSpecPerReact; ++i)
        {
          m_thPhys.m_ford[refARBN + i] = m_thPhys.m_nup[refRNum + i];
          m_thPhys.m_rord[refARBN + i] = m_thPhys.m_nupp[refRNum + i];
        }
    }
  const int initLength = a_str.length() + 1;
  size_t startloc = a_str.find_first_not_of(" ");
  a_str = a_str.substr(startloc, initLength);
  size_t loc = a_str.find(" ");
  std::string spec = a_str.substr(0,loc);
  // Remove part of the string we have already read
  a_str = a_str.substr(loc+1, initLength);
  // Give delimiter that shouldn't ever be present
  Real coeff = extractValue(a_str, "*");
  // Number of species on this side of the reaction
  int numSpecs = m_thPhys.m_numRects[a_cRN];
  if (!a_FORD)
    {
      numSpecs = m_thPhys.m_numProds[a_cRN];
    }
  int refSpec = -1;
  for (int i = 0; i != numSpecs; ++i)
    {
      int curRefSpec = m_thPhys.m_refNup[refRNum + i];
      if (!a_FORD)
        {
          curRefSpec = m_thPhys.m_refNupp[refRNum + i];
        }
      std::string specName = CRDparam::g_speciesNames[curRefSpec];
      std::string altName = m_thPhys.m_altSpecNames[curRefSpec];
      if (spec == specName || spec == altName)
        {
          if (a_FORD)
            {
              m_thPhys.m_ford[refARBN + i] = coeff;
            }
          else
            {
              m_thPhys.m_rord[refARBN + i] = coeff;
            }
          refSpec = curRefSpec;
        }
    }
  if (refSpec == -1)
    {
      CRD::msg << "ThermPhysics::extractARBCoeff: " << spec
               << " found in reaction " << a_cRN
               << " but not specified in inputs!" << CRD::error;
    }
}

/*--------------------------------------------------------------------*/
//  We use the TB alpha values for the pressure-dependent reactions that
//  do not use TB
//  IE, if the PR reaction is not TB it would look like (+H2O).
//  In that case, the alpha values for this reaction will all be zero except
//  the one for H2O and the m_thirdBody[a_tbRN] = -1
/** \param[in]  a_str   String containing the species name
 *  \param[in]  a_tbRN  The number of the current TB reaction
 *//*-----------------------------------------------------------------*/

void
FileThermParser::extractPRSpeciesAlpha(const std::string a_str,
                                       const int&        a_tbRN)
{
  int refSpec = -1;
  for (int i = 0; i != CRDparam::g_numSpecies; ++i)
    {
      if (a_str == CRDparam::g_speciesNames[i] ||
          a_str == m_thPhys.m_altSpecNames[i])
        {
          // Determine the component based on the reaction and species
          const int rcomp = i*m_thPhys.m_numTBReact + a_tbRN;
          m_thPhys.m_tbAlpha[rcomp] = 1.;
        }
      else
        {
          const int rcomp = i*m_thPhys.m_numTBReact + a_tbRN;
          m_thPhys.m_tbAlpha[rcomp] = 0.;
        }
    }
  if (refSpec == -1)
    {
      CRD::msg << "ThermPhysics::extractPRSpeciesAlpha: " << a_str
               << " found in reaction but not specified in inputs!"
               << CRD::error;
    }
}

/*--------------------------------------------------------------------*/
//  Fill tables for forward and backward reaction rate lookup tables
/** 
 *//*-----------------------------------------------------------------*/

void
FileThermParser::fillReactTables(std::vector<Real>& a_preAF,
                                 std::vector<Real>& a_betai,
                                 std::vector<Real>& a_EAR,
                                 std::vector<Real>& a_REVpreAB,
                                 std::vector<Real>& a_REVbetai,
                                 std::vector<Real>& a_REVEAR)
{
  const int nRCT = m_thPhys.m_numReactions;
  const int numSpecies = CRDparam::g_numSpecies;
  Real factor1 = 1./20.;
  Real patmR = 101325.0*10./(m_thPhys.m_Rmol*1.E7);
  Real sixth = 1./6.;
  Real twelfth = 1./12.;
  std::vector<Real> gfval(numSpecies);
  m_thPhys.m_lookupKfwd.resize(nRCT);
  m_thPhys.m_lookupKbkwd.resize(nRCT);
  m_thPhys.m_lookupKeq.resize(nRCT);
  int curRev = 0;
  for (int rctcomp = 0; rctcomp != nRCT; ++rctcomp)
    {
      m_thPhys.m_lookupKfwd[rctcomp].assign(m_thPhys.m_lookupSize,0.);
      m_thPhys.m_lookupKbkwd[rctcomp].assign(m_thPhys.m_lookupSize,0.);
      m_thPhys.m_lookupKeq[rctcomp].assign(m_thPhys.m_lookupSize,0.);
      // Number of reactants for this reaction
      const int cnumr = m_thPhys.m_numRects[rctcomp];
      // Number of products for this reaction
      const int cnump = m_thPhys.m_numProds[rctcomp];
      // Current index for nup and nupp
      const int nucomp = m_thPhys.m_maxSpecPerReact*rctcomp;
      for (int i = 0; i != m_thPhys.m_lookupSize; ++i)
        {
          const Real T = m_thPhys.m_lookupLoT + (Real)(i)*m_thPhys.m_lookupDelT;
          Real Tinv = 1./T;
          Real bkconst = patmR*Tinv;
          Real Tlog = std::log(T);
          Real kfwd = a_preAF[rctcomp]*std::pow(T, a_betai[rctcomp])*
            std::exp(-a_EAR[rctcomp]*Tinv);
          m_thPhys.m_lookupKfwd[rctcomp][i] = kfwd;
          // Solve for normalized Gibbs free energy
          if (T >= m_thPhys.m_midLookupT)
            {
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  gfval[sp] =
                    -(m_thPhys.m_hnB1H[sp] + m_thPhys.m_hnA2H[sp] -
                      0.5*m_thPhys.m_hnA1H[sp]*Tinv
                      +Tlog*(-m_thPhys.m_hnA3H[sp]*T + m_thPhys.m_hnA2H[sp]) +
                      T*(m_thPhys.m_hnA3H[sp] - m_thPhys.m_hnB2H[sp] -
                         T*(m_thPhys.m_hnA4H[sp]*0.5 +
                            T*(m_thPhys.m_hnA5H[sp]*sixth +
                               T*(m_thPhys.m_hnA6H[sp]*twelfth +
                                  T*m_thPhys.m_hnA7H[sp]*factor1)))))*Tinv;
                }
            }
          else
            {
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  gfval[sp] =
                    -(m_thPhys.m_hnB1L[sp] + m_thPhys.m_hnA2L[sp] -
                      0.5*m_thPhys.m_hnA1L[sp]*Tinv
                      +Tlog*(-m_thPhys.m_hnA3L[sp]*T + m_thPhys.m_hnA2L[sp]) +
                      T*(m_thPhys.m_hnA3L[sp] - m_thPhys.m_hnB2L[sp] -
                         T*(m_thPhys.m_hnA4L[sp]*0.5 +
                            T*(m_thPhys.m_hnA5L[sp]*sixth +
                               T*(m_thPhys.m_hnA6L[sp]*twelfth +
                                  T*m_thPhys.m_hnA7L[sp]*factor1)))))*Tinv;
                }
            }
          if (m_thPhys.m_revReact[rctcomp] == 0)
            {
              m_thPhys.m_lookupKbkwd[rctcomp][i] = a_REVpreAB[curRev]*
                std::pow(T, a_REVbetai[curRev])*
                std::exp(-a_REVEAR[curRev]*Tinv);
            }
          else if (m_thPhys.m_revReact[rctcomp] > 0)
            {
              // dsmdh is the $\sum^{N_s}_{n=1} \nu_{n,r} -\frac{G_n}{R_u}$
              Real dsmdh = 0.;
              // sumnuval is the $\sum^{N_s}_{n=1} \nu_{n,r}$
              Real sumnuval = 0.;
              for (int offcomp = nucomp; offcomp != nucomp+cnumr; ++offcomp)
                {
                  const int spcomp1 = m_thPhys.m_refNup[offcomp];
                  sumnuval -= m_thPhys.m_nup[offcomp];
                  dsmdh -= m_thPhys.m_nup[offcomp]*gfval[spcomp1];
                }
              for (int offcomp = nucomp; offcomp != nucomp+cnump; ++offcomp)
                {
                  const int spcomp1 = m_thPhys.m_refNupp[offcomp];
                  sumnuval += m_thPhys.m_nupp[offcomp];
                  dsmdh += m_thPhys.m_nupp[offcomp]*gfval[spcomp1];
                }
              // Solve for equilibrium rate
              Real keq = std::exp(dsmdh)*std::pow(bkconst,sumnuval);
              m_thPhys.m_lookupKeq[rctcomp][i] = keq;
              // Replace kbkwd with $k_{b,r} \prod [X_n]^{\nu^{''}_{n,r}}
              m_thPhys.m_lookupKbkwd[rctcomp][i] = kfwd/keq;
            }
          else
            {
              m_thPhys.m_lookupKbkwd[rctcomp][i] = 0.;
            }
        }
       if (m_thPhys.m_revReact[rctcomp] == 0)
         {
           curRev++;
         }
    }
}

/*--------------------------------------------------------------------*/
//  Read in the transport data
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::readTransportFile(std::vector<std::array<Real, 4> >& a_muVecL,
                                   std::vector<std::array<Real, 4> >& a_muVecH,
                                   std::vector<std::array<Real, 4> >& a_kVecL,
                                   std::vector<std::array<Real, 4> >& a_kVecH)
{
  const int numSpecies = CRDparam::g_numSpecies;
  // Vector to check which species were not found
  std::vector<bool> propCheck(numSpecies, false);
  // Vector to check which species have improper ranges
  std::vector<bool> rangeCheck(numSpecies, false);
  std::string transFile = m_thPhys.m_transFile;
  if (!CH_System::fileExists(m_thPhys.m_transFile.c_str()))
    {
      CRD::msg << "Input: 'transport_file' " << m_thPhys.m_transFile
               << " does not exist!" << CRD::error;
    }
  int thLen = transFile.length();
  char *filename = new char[thLen + 2];
  std::strcpy(filename, transFile.c_str());
  ifstream readFile(filename);
  delete [] filename;
  // The number of chars in each value we are extracting is 15
  const int valLength = 15;
  if (readFile.is_open())
    {
      // The current number corresponding to the species
      int curSpec = -1;
      int curLine = -1;
      // Shows if we are reading the low values for thermal conductivity
      bool condLow = true;
      std::string readLine;
      while (std::getline(readFile, readLine))
        {
          std::string specstring =
            readLine.substr(0,readLine.find(" ",1));
          for (int i = 0; i != numSpecies; ++i)
            {
              if (specstring == CRDparam::g_speciesNames[i] ||
                  specstring == m_thPhys.m_altSpecNames[i])
                {
                  // Ensure it is not a binary interaction
                  std::string secspec = readLine.substr(16,17);
                  if (secspec[0] == ' ')
                    {
                      curSpec = i;
                      curLine = 0;
                      propCheck[i] = true;
                      rangeCheck[i] = true;
                      condLow = true;
                      break;
                    }
                }
            }
          if (curSpec == -1)
            {
              continue;
            }
          findReplace(readLine, "E ", "E+");
          if (curLine == 1)
            {
              std::string lowTS = readLine.substr(3,9);
              Real lowT = extractValue(lowTS, 6);
              // If the low temperature on the lower interval is too high,
              // Then we should set the lower interval as not found and
              // skip to the upper intervals
              if (lowT > 900.)
                {
                  curLine++;
                  condLow = false;
                  rangeCheck[curSpec] = false;
                }
            }
          if (curLine == 1)
            {
              std::string numLine = readLine.substr(20,readLine.size());
              a_muVecL[curSpec][0] = extractValue(numLine,valLength);
              a_muVecL[curSpec][1] = extractValue(numLine,valLength);
              a_muVecL[curSpec][2] = extractValue(numLine,valLength);
              a_muVecL[curSpec][3] = extractValue(numLine,valLength);
            }
          else if (curLine == 2)
            {
              std::string numLine = readLine.substr(20,readLine.size());
              a_muVecH[curSpec][0] = extractValue(numLine,valLength);
              a_muVecH[curSpec][1] = extractValue(numLine,valLength);
              a_muVecH[curSpec][2] = extractValue(numLine,valLength);
              a_muVecH[curSpec][3] = extractValue(numLine,valLength);
            }
          // If a third temperature interval is given, ignore it
          else if (readLine[1] == 'V')
            {
              continue;
            }
          else if (readLine[1] == 'C' && condLow == true)
            {
              std::string numLine = readLine.substr(20,readLine.size());
              a_kVecL[curSpec][0] = extractValue(numLine,valLength);
              a_kVecL[curSpec][1] = extractValue(numLine,valLength);
              a_kVecL[curSpec][2] = extractValue(numLine,valLength);
              a_kVecL[curSpec][3] = extractValue(numLine,valLength);
              condLow = false;
            }
          else if (readLine[1] == 'C')
            {
              std::string numLine = readLine.substr(20,readLine.size());
              a_kVecH[curSpec][0] = extractValue(numLine,valLength);
              a_kVecH[curSpec][1] = extractValue(numLine,valLength);
              a_kVecH[curSpec][2] = extractValue(numLine,valLength);
              a_kVecH[curSpec][3] = extractValue(numLine,valLength);
              curSpec = -1;
            }
          curLine++;
        }
    }
  // Check which species are not found and set them to be the same as N2
  for (int i = 0; i != numSpecies; ++i)
    {
      if (rangeCheck[i] == false)
        {
          a_muVecL[i][0] = 0.62526577;
          a_muVecL[i][1] = -0.31779652E2;
          a_muVecL[i][2] = -0.16407983E4;
          a_muVecL[i][3] = 0.17454992E1;
          a_kVecL[i][0] = 0.85439436;
          a_kVecL[i][1] = 0.10573224E3;
          a_kVecL[i][2] = -0.12347848E5;
          a_kVecL[i][3] = 0.47793128;
          if (propCheck[i] == false)
            {
              a_muVecH[i][0] = 0.87395209;
              a_muVecH[i][1] = 0.56152222E3;
              a_muVecH[i][2] = -0.17394809E6;
              a_muVecH[i][3] = -0.39335958;
              a_kVecH[i][0] = 0.88407146;
              a_kVecH[i][1] = 0.13357293E3;
              a_kVecH[i][2] = -0.1142964E5;
              a_kVecH[i][3] = 0.24417019;
              CRD::msg << "FileThermParser::readTransportFile: Transport "
                       << "data of " << CRDparam::g_speciesNames[i]
                       << " was not found in transport data file. "
                       << "Using transport data of N2!" << CRD::warn;
            }
          else
            {
              CRD::msg << "FileThermParser::readTransportFile: Transport "
                       << "data of " <<CRDparam::g_speciesNames[i]
                       << " was only partially transport data file. "
                       << "300 - 1000 K data is using transport data of N2!"
                       << CRD::warn;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set the values in the lookup tables for transport properties
/**
 *//*-----------------------------------------------------------------*/

void
FileThermParser::fillTransportTables(
  std::vector<std::array<Real, 4> >& a_muVecL,
  std::vector<std::array<Real, 4> >& a_muVecH,
  std::vector<std::array<Real, 4> >& a_kVecL,
  std::vector<std::array<Real, 4> >& a_kVecH)
{
  const int nS = CRDparam::g_numSpecies;
  std::vector<Real> initVec(nS,0.);
  m_thPhys.m_lookupMu.assign(nS,initVec);
  m_thPhys.m_lookupKappa.assign(nS,initVec);
  for (int sp = 0; sp != nS; ++sp)
    {
      m_thPhys.m_lookupMu[sp].assign(m_thPhys.m_lookupSize,0.);
      m_thPhys.m_lookupKappa[sp].assign(m_thPhys.m_lookupSize,0.);
      for (int i = 0; i != m_thPhys.m_lookupSize; ++i)
        {
          const Real T = m_thPhys.m_lookupLoT + (Real)(i)*m_thPhys.m_lookupDelT;
          if (T >= m_thPhys.m_midLookupT)
            {
              m_thPhys.m_lookupMu[sp][i] = exp(a_muVecH[sp][0]*log(T) +
                                               (a_muVecH[sp][2]/T +
                                                a_muVecH[sp][1])/T +
                                               a_muVecH[sp][3])*1.E-7;
              m_thPhys.m_lookupKappa[sp][i] = exp(a_kVecH[sp][0]*log(T) +
                                                  (a_kVecH[sp][2]/T +
                                                   a_kVecH[sp][1])/T +
                                                  a_kVecH[sp][3])*1.E-4;
            }
          else
            {
              m_thPhys.m_lookupMu[sp][i] = exp(a_muVecL[sp][0]*log(T) +
                                               (a_muVecL[sp][2]/T +
                                                a_muVecL[sp][1])/T +
                                               a_muVecL[sp][3])*1.E-7;
              m_thPhys.m_lookupKappa[sp][i] = exp(a_kVecL[sp][0]*log(T) +
                                                  (a_kVecL[sp][2]/T +
                                                   a_kVecL[sp][1])/T +
                                                  a_kVecL[sp][3])*1.E-4;
            }
        }
    }
}
