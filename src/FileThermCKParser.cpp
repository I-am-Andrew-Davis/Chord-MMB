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

#include "FileThermCKParser.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "DataTemp.H"

/*******************************************************************************
 *
 * Class FileThermCKParser: member definitions
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
FileThermCKParser::readThermFile()
{
  if (!CH_System::fileExists(m_thPhys.m_thermoFile.c_str()))
    {
      CRD::msg << "Input: 'thermo_file' " << m_thPhys.m_thermoFile
               << " does not exist!" << CRD::error;
    }
  const int numSpecies = CRDparam::g_numSpecies;
  // Note: Chemkin files do not specify the molar masses so they
  // must be hard-coded for each reference element
  // FIXME: The element and reference lists should be comprehensive
  // List of elements that could exist in the species
  // Note: two letter elements must be specified as upper then lower case,
  // like Be, not BE
  std::vector<std::string> elementList = {"He", "C", "H", "O", "N", "Ar"};
  // Molar masses for reference elements
  // FIXME: This is not comprehensive
  std::vector<Real> molMassEl =
    {4.002602, 12.0107, 1.00794, 15.999, 14.0067, 39.948};
  const int numRefElements = elementList.size();
  CH_assert(numRefElements == molMassEl.size());
  m_thPhys.m_molMass.assign(numSpecies, 0.);
  // Assign the molar masses
  for (int i = 0; i != numSpecies; ++i)
    {
      int specSize = CRDparam::g_speciesNames[i].length();
      for (int j = 0; j != numRefElements; ++j)
        {
          size_t locspec = CRDparam::g_speciesNames[i].find(elementList[j]);
          // Check if the species is a reference state
          if (CRDparam::g_speciesNames[i] == elementList[j])
            {
              m_thPhys.m_molMass[i] = molMassEl[j]/1000.;
              continue;
            }
          // Check if the species contains one of the elements
          else if (locspec != std::string::npos)
            {
              // Remove everything before th current species of interest
              std::string speciesName =
                CRDparam::g_speciesNames[i].substr(locspec, specSize);
              // Number of elements in species
              Real numEl = 0.;
              int sizeEl = elementList[j].length();
              for (int i = 0; i != speciesName.length(); ++i)
                {
                  std::string cutst = speciesName.substr(i,sizeEl);
                  if (cutst == elementList[j])
                    {
                      numEl++;
                      int diffLen = speciesName.length() - (i+sizeEl);
                      if (diffLen > 0)
                        {
                          char nc = speciesName.at(i+sizeEl);
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
                                  char nc2 = speciesName.at(i+sizeEl+1);
                                  if (std::isdigit(nc2))
                                    {
                                      numEl = numEl*10 + (nc2 - '0');
                                    }
                                }
                            }
                        }
                    }
                }
              m_thPhys.m_molMass[i] += numEl*molMassEl[j]/1000.;
            }
        }
      if (m_thPhys.m_molMass[i] == 0.)
        {
          CRD::msg << "Molar mass for species " << CRDparam::g_speciesNames[i]
                   << " is not known" << CRD::error;
        }
    }
  int thLen = m_thPhys.m_thermoFile.length();
  char *filename = new char[thLen + 2];
  std::strcpy(filename, m_thPhys.m_thermoFile.c_str());
  ifstream readFile(filename);
  delete [] filename;
  // The number of chars in each value we are extracting is 15
  const int valLength = 15;
  // Set up here to prevent compiler warning
  if (readFile.is_open())
    {
      // The current number corresponding to the species
      int curSpec = -1;
      int curLine = -1;
      std::string readLine;
      while (std::getline(readFile, readLine))
        {
          // Skip commented out lines
          if (readLine.find("!") <= 3)
            {
              continue;
            }
          if (readLine.find("REACTIONS") != std::string::npos)
            {
              break;
            }
          if (curSpec == -1)
            {
              if (readLine.find(" ") == 0)
                {
                  continue;
                }
              else
                {
                  std::string specstring =
                    readLine.substr(0,readLine.find(" ",1));
                  for (int i = 0; i != numSpecies; ++i)
                    {
                      if (specstring == CRDparam::g_speciesNames[i] ||
                          specstring == m_thPhys.m_altSpecNames[i])
                        {
                          curSpec = i;
                          curLine = 0;
                          break;
                        }
                    }
                }
            }
          if (curSpec == -1)
            {
              continue;
            }
          // Extract the first portion of the high enthalpy coefficients
          if (curLine == 1)
            {
              m_thPhys.m_hnA1H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA2H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA3H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA4H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA5H[curSpec] = extractValue(readLine,valLength);
            }
          else if (curLine == 2)
            {
              m_thPhys.m_hnA6H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA7H[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA1L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA2L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA3L[curSpec] = extractValue(readLine,valLength);
            }
          else if (curLine == 3)
            {
              m_thPhys.m_hnA4L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA5L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA6L[curSpec] = extractValue(readLine,valLength);
              m_thPhys.m_hnA7L[curSpec] = extractValue(readLine,valLength);
              curSpec = -1; // Reset values for next species
            }
          curLine++;
        }
    }
  readFile.close();
  // Check if any species were not found
  for (int i = 0; i != numSpecies; ++i)
    {
      // Heat of formation is assumed to be included in fits
      m_thPhys.m_H0[i] = 0.;
      m_thPhys.m_Rn[i] = m_thPhys.m_Rmol/m_thPhys.m_molMass[i];
      if (m_thPhys.m_hnA1L[i] == -1.)
        {
          CRD::msg << "FileThermCKParser::readThermFile: No inputs found for "
                   << CRDparam::g_speciesNames[i] << "!" << CRD::error;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Set the values in the lookup tables for thermodynamic properties
/**
 *//*-----------------------------------------------------------------*/

void
FileThermCKParser::fillThermoTables()
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
          if (T >= m_thPhys.m_midLookupT)
            {
              m_thPhys.m_lookupH[sp][i] = m_thPhys.m_Rn[sp]*
                (m_thPhys.m_hnA6H[sp] + T*(m_thPhys.m_hnA1H[sp] +
                                           T*(m_thPhys.m_hnA2H[sp]/2. +
                                     T*(m_thPhys.m_hnA3H[sp]/3. +
                                        T*(m_thPhys.m_hnA4H[sp]/4. +
                                           T*m_thPhys.m_hnA5H[sp]/5.)))));
            }
          else
            {
              m_thPhys.m_lookupH[sp][i] = m_thPhys.m_Rn[sp]*
                (m_thPhys.m_hnA6L[sp] + T*(m_thPhys.m_hnA1L[sp] +
                                  T*(m_thPhys.m_hnA2L[sp]/2. +
                                     T*(m_thPhys.m_hnA3L[sp]/3. +
                                        T*(m_thPhys.m_hnA4L[sp]/4. +
                                           T*m_thPhys.m_hnA5L[sp]/5.)))));
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Fill tables for forward and backward reaction rate lookup tables
/** 
 *//*-----------------------------------------------------------------*/

void
FileThermCKParser::fillReactTables(std::vector<Real>& a_preAF,
                                   std::vector<Real>& a_betai,
                                   std::vector<Real>& a_EAR,
                                   std::vector<Real>& a_REVpreAB,
                                   std::vector<Real>& a_REVbetai,
                                   std::vector<Real>& a_REVEAR)
{
  const int nRCT = m_thPhys.m_numReactions;
  const int numSpecies = CRDparam::g_numSpecies;
  Real patmR = 101325.0*10./(m_thPhys.m_Rmol*1.E7);
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
                  gfval[sp] = m_thPhys.m_hnA7H[sp] - m_thPhys.m_hnA6H[sp]*Tinv +
                    m_thPhys.m_hnA1H[sp]*(Tlog - 1.) +
                    T*(m_thPhys.m_hnA2H[sp]/2. + T*(m_thPhys.m_hnA3H[sp]/6. +
                                           T*(m_thPhys.m_hnA4H[sp]/12. +
                                              T*m_thPhys.m_hnA5H[sp]/20.)));
                }
            }
          else
            {
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  gfval[sp] = m_thPhys.m_hnA7L[sp] - m_thPhys.m_hnA6L[sp]*Tinv +
                    m_thPhys.m_hnA1L[sp]*(Tlog - 1.) +
                    T*(m_thPhys.m_hnA2L[sp]/2. + T*(m_thPhys.m_hnA3L[sp]/6. +
                                           T*(m_thPhys.m_hnA4L[sp]/12. +
                                              T*m_thPhys.m_hnA5L[sp]/20.)));
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
