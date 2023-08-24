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
 * \file CRDmsg.cpp
 *
 * \brief Member functions for CRDmsg
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

//----- Chombo Library -----//

#include "parstream.H"
#include "CH_assert.H"

//----- Internal -----//

#include "CRDmsg.H"


/*******************************************************************************
 *
 * Class Msg: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/** The constructor here forces definition of all statics in class but
 *  pout may not be defined yet so m_out must be explicitly set.
 *//*-----------------------------------------------------------------*/

CRD::Msg::Msg()
:
  m_lineSize(defaultLineSize),
  m_out(0),  // User will quickly notice if this is not explicitly set
  m_errorHeader("EE "),
  m_warnHeader ("WW "),
  m_header1    ("=> "),
  m_header2    ("-> "),
  m_bodyHeader ("   "),
  m_varHeaderFill('-'),
  m_varNameTerm  (" "),
  m_varHeaderTerm(" "),
  m_varNameSize(std::max((size_t)20, m_lineSize/2 - m_bodyHeader.size() -
                         m_varHeaderTerm.size())),
  m_defaultPrec(6),
  m_verbosity(0),
  m_terminateOnError(true)
{
  CH_assert(m_varNameSize > m_varNameTerm.size());
  m_defaultPrec = this->precision();
  char tmpVerbHeader[4];
  for (int i = 0; i != 5; ++i)
    {
      sprintf(tmpVerbHeader, "%1d> ", i);
      m_verbHeader[i] = tmpVerbHeader;
    }
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Modification of default headers and setup
/*--------------------------------------------------------------------*/

//
void
CRD::Msg::setLineSize(const size_t a_lineSize)
{
  CH_assert(a_lineSize > 0);
  m_lineSize = a_lineSize;
  m_varNameSize = std::max((size_t)20, m_lineSize/2 - m_bodyHeader.size() -
                           m_varHeaderTerm.size());
  CH_assert(m_varNameSize > m_varNameTerm.size());
}

//
void
CRD::Msg::setOutputStream(std::ostream& a_out)
{
  m_out = &a_out;
}

//
void
CRD::Msg::setErrorHeader(const char* a_header)
{
  m_errorHeader = a_header;
}

//
void
CRD::Msg::setWarnHeader(const char* a_header)
{
  m_warnHeader = a_header;
}

//
void
CRD::Msg::setHeader1(const char* a_header)
{
  m_header1 = a_header;
}

//
void
CRD::Msg::setHeader2(const char* a_header)
{
  m_header2 = a_header;
}

//
void
CRD::Msg::setBodyHeader(const char* a_header)
{
  m_bodyHeader = a_header;
}

//
void
CRD::Msg::setVerbosityHeader(const char* a_header, const int a_idxVerb)
{
  CH_assert(a_idxVerb >= 0 && a_idxVerb <= 4);
  m_verbHeader[a_idxVerb] = a_header;
}

/*--------------------------------------------------------------------*/
/// Set the default floating point precision
/** Note that this only changes storage of the default precision.  To
 *  change the precision used, you also need to call the
 *  'setDefaultFloat' member function
 *  \param[in]  a_prec  Precision
 *//*-----------------------------------------------------------------*/

void
CRD::Msg::setDefaultPrec(const int a_prec)
{
  m_defaultPrec = a_prec;
}

/*--------------------------------------------------------------------*/
/// Set the verbosity of the output
/** This should have the same value as CRDparam::g_verbosity.
 *  CRD::Msg does not use it directly to keep this low-level class
 *  separate.
 *  \param[in] a_verbosity
 *                      Level of verbosity (0 to 4)
 *//*-----------------------------------------------------------------*/

void
CRD::Msg::setVerbosity(const int a_verbosity)
{
  m_verbosity = a_verbosity;
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Write a message to describe a variable in form name:value
/*--------------------------------------------------------------------*/

void
CRD::Msg::writeVariable()
{
  const char* header = 0;
  std::string variable(str());
  std::string::size_type brk = variable.find_first_of('\n');
  std::string name(m_varNameSize, m_varHeaderFill);
  if (brk == std::string::npos)
    {
      header = m_bodyHeader.c_str();
    }
  else
    {
      const int nameFillLen = m_varNameSize - m_varNameTerm.size();
      const int nameFullLen =
        m_varNameSize + m_varHeaderTerm.size();
      if (brk > nameFillLen)
        {
          variable.replace(brk, 1, m_varNameTerm + m_varHeaderTerm);
          name = variable.substr(0, nameFullLen);
          variable.erase(0, nameFullLen);
        }
      else
        {
          name.replace(0,
                       brk + m_varNameTerm.size(),
                       variable.substr(0, brk) + m_varNameTerm);
          name.append(m_varHeaderTerm);
          variable.erase(0, brk+1);
        }
      name.insert(0, m_bodyHeader);
      header = name.c_str();
      str(variable);  // Reset the stringstream
    }
  linePrint(str().c_str(), *m_out, header, -1, m_lineSize);
}

/*******************************************************************************
 *
 * Class Msg: external related functions
 *
 ******************************************************************************/


/*============================================================================*/
/// Insert headers and line breaks into console output for nice printing
/** The primary purpose is to place line breaks in a string for nice printing.
 *  However, this routine also supports the placement of a header before each
 *  line with the option of indenting lines.  For even more flexible output,
 *  the number of lines printed can be specified and used in conjuction with
 *  the returned index in the string.
 *  \param[in]  a_message
 *                      The string to print (must be non-null)
 *  \param[in]  a_out   The output stream (default std::cout)
 *  \param[out] a_out   Modified stream
 *  \param[in]  a_header
 *                      Header to print on lines
 *  \param[in]  a_numLines
 *                      > 0 the number of lines to print
 *                      = 0 print all lines
 *                      < 0 line number to start indenting (instead of printing
 *                          header)
 *  \param[in]  a_lineSize
 *                      The size of a line.  Default is specified by
 *                      CRD::defaultLineSize
 *  \return             New index in the string
 *
 *  \note
 *  <ul>
 *    <li> The header takes up space in the line.  So if your header is too
 *         long, you will get some really short text lines.
 *    <li> Headers can wrap.  Indentations will not wrap but will match the size
 *         of the wrapped header.
 *    <li> If a space is not found for a line break, a word will be broken.
 *    <li> You can go all the way down to 1 character per line if you want too.
 *    <li> If a header is all spaces, +2 will be added to indentations
 *  </ul>
 *
 *  Example:
 *                                                                     \verbatim
#include<iostream>
#include "CRDmsg.H"

int main() {

   const char* text = "This text is simply used to test the output method.  "
      "You should not try to derive any kind of deep meaning from it.";

//--Typical printing

   std::cout << "Typical printing:\n\n";
   CRD::linePrint(text);

//--With header

   std::cout << "\nWith header:\n\n";
   CRD::linePrint(text, std::cout, "--->");

//--Change line size

   std::cout << "\nChange line size to 40 and indent after first line:\n\n";
   CRD::linePrint(text, std::cout, "->", -1, 40);

//--Get real fancy

   std::cout << "\nGet real fancy:\n\n";
   
   int j = CRD::linePrint(text, std::cout, "---\\ ", 1, 40);
   j += CRD::linePrint(text+j, std::cout, "----\\", 1, 40);
   j += CRD::linePrint(text+j, std::cout, "-----\\", 1, 40);
   CRD::linePrint(text+j, std::cout, "------\\", 0, 40);
}                                                
 *                                                                  \endverbatim
 *//*=========================================================================*/

size_t
CRD::linePrint(const char*   a_message,
               std::ostream& a_out,
               const char*   a_header,
               int           a_numLines,
               size_t        a_lLineSize)
{

//--Initialize input

  const size_t lenMessage(std::strlen(a_message));
  unsigned int headerSize(0);
  unsigned int lenHeader(0);
  if (a_lLineSize <= 0)
    {
      a_lLineSize = CRD::defaultLineSize;
    }
  // Get header size.  Remainder if longer than line size.
  if (a_header != 0)
    {
      lenHeader = std::strlen(a_header);
      headerSize = lenHeader % a_lLineSize;
      // std::cout << "HS: " << headerSize << std::endl;
    }
  size_t startIndentLines(0);  // Line to stop printing header and just indent
  if (a_numLines < 0)
    {
      startIndentLines = -a_numLines;
      a_numLines = 0;
    }

//--Local variables

  size_t lineNumber(0);
  size_t iss1(0);

//--Begin line output loop

  while (true)
    {
      ++lineNumber;
      // Find next non-blank character if not first line
      if (lineNumber != 1)
        {
          iss1 += std::strspn(a_message+iss1, " ");
          if (iss1 == lenMessage)
            {
              return iss1;
            }
        }
      // Find end of line
      int iss2 = std::min(lenMessage, iss1 + a_lLineSize - headerSize);
      if ((iss2 != lenMessage) && (a_message[iss2] != ' '))
        {
          // Find previous blank character (line break)
          size_t i = iss2 - 1;
          for (; ( (i != iss1) && (a_message[i] != ' ') ); --i);
          if (i != iss1)
            {
              iss2 = i + 1;  // Otherwise full line with no space
            }
        }
      // Eliminate blanks at end
      for (--iss2 ; a_message[iss2] == ' ' ; --iss2);
      ++iss2;

//--Now have line in substring between iss1 and iss2

      // Write the header or indentation
      if (a_header != 0)
        {
          if (startIndentLines != 0 && lineNumber > startIndentLines)
            {
              for ( size_t i = 0 ; i != headerSize ; ++i )
                {
                  a_out.put(' ');
                }
            }
          else
            {
              size_t j = 0;
              for (size_t i = 0 ; i != lenHeader ; ++i)
                {
                  a_out.put(a_header[i]);
                  ++j;
                  if (j == a_lLineSize)
                    {
                      a_out << std::endl;
                      j = 0;
                    }
                }
              // +2 to indentation if header is all blank
              if (startIndentLines != 0 && lineNumber == startIndentLines)
                {
                  if (std::strspn(a_header, " ") == lenHeader)
                    {
                      headerSize = (lenHeader + 2) % a_lLineSize;
                    }
                }
            }
        }
      // Write the message
      for (size_t i = iss1 ; i != iss2 ; ++i)
        {
          a_out.put(a_message[i]);
        }
      a_out << std::endl;
      iss1 = iss2;
      if (lineNumber == a_numLines)
        {
          return iss1;
        }
    }
}
