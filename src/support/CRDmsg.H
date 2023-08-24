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
 * \file CRDmsg.H
 *
 * \brief Advanced terminal output support
 *
 *//*+*************************************************************************/

#ifndef _CRDMSG_H_
#define _CRDMSG_H_

//----- Standard Library -----//

#include <ostream>
#include <sstream>
#include <limits>

//----- Chombo Library -----//

#include "REAL.H"


namespace CRD
{

/// Stream markers for terminating a message and specifying a type
enum IOTerm
{
  warn,                               ///< Writes a warning
  error,                              ///< Writes an error and exits
  abort,                              ///< Writes an error an aborts
  h1,                                 ///< Writes a header message with =>
  h2,                                 ///< Writes a header mesasge with ->
  body,                               ///< Writes a regular message
  var,                                ///< Writes a variable (name: value)
  v0,
  v1,
  v2,
  v3,
  v4,
  varend,                             ///< Same as var but no Chng. to verbosity
  end                                 ///< Generic end message
};

/// Stream manipulators that set the verbosity and type of message
enum IOFormat
{
  fwarn,
  ferror,
  fabort,
  fh1,
  fh2,
  fbody,
  fvar,
  fv0,
  fv1,
  fv2,
  fv3,
  fv4
};

/// Stream manipulators that inject special strings
enum IOPrint
{
  wasNumError,
  wasNumWarning
};
  
/// Verbosity manipulator
struct verb
{
  verb(const int a_msgVerb)
  :
  m_msgVerb(a_msgVerb)
    {
      CH_assert(a_msgVerb >= -1 && a_msgVerb <= 4);
    }
  int m_msgVerb;
};

const size_t defaultLineSize(80);     ///< Default maximum size of line when
                                      ///< printing message.

/// Insert headers and line breaks into console output for nice printing
size_t linePrint(const char*   a_message,
                 std::ostream& a_out = std::cout,
                 const char*   a_header = 0,
                 int           a_numLines = 0,
                 size_t        a_lLineSize = defaultLineSize);


/*******************************************************************************
 */
/// Output stream for writing a message
/**
 *  Similar functionality to 'parstream.H' and 'MayDay.H' from Chombo except
 *  manipulators or terminators can be used to define the type of message.
 *
 *  Example:
 *  For most headers and body text
 *                                                                         \code
 *  CRD::msg << "This is header 1" << CRD::h1;                          \endcode
 *
 *  For errors and warnings
 *                                                                         \code
 *  CRD::msg << "Invalid argument" << CRD::error;                       \endcode
 *
 *  For variables (variable and value), newline is used to separate entries
 *                                                                         \code
 *  CRD::msg << "CFL\n1.0" << CRD::var;                                 \endcode
 *
 *  Manipulators can also be used to set the verbosity and type of message.
 *  Type manipulators are the same as the terminators but begin with an 'f'.
 *  These let the class know in advance what the verbosity and type of the
 *  next message to print will be.    Since the class knows the type of
 *  message, the generic terminator 'CRD::end' can be used to print the message.
 *  Manipulators can also be used to force a low message verbosity so that it is
 *  always printed.  Note that type manipulators set a prescribed verbosity
 *  so the verbosity manipulator should be used after type is set.
 *                                                                         \code
 *  CRD::msg << CRD::fbody << CRD::verb(-1) << "Usage" << CRD::end;     \endcode
 *
 *  The above will always print.  Note that using 'CRD::body' as a terminator
 *  above would have reset the verbosity to the default for 'body' before
 *  printing the message.  All messages start out as CRD::body and default
 *  verbosity (i.e., that used to set level of printing).
 *
 *  If the message verbosity is set higher than allowed for printing, the stream
 *  will not even be buffered.  This may or may not have a performance impact.
 *  It is suggested to set the type and verbosity in advance for low priority
 *  messages
 *                                                                         \code
 *  CRD::msg << CRD::fv3 << "AMRLevel::define (level" << m_level << ')'
 *           << CRD::end;                                               \endcode
 *
 *  In the above, nothing would be buffered with the verbosity set less than 3.
 *  With verbosity set to 1, the following example would only print the header,
 *  "=> ", since nothing is buffered but the verbosity is changed to 0 (that
 *  of h1) at termination
 *                                                                         \code
 *  CRD::msg << CRD::fv3 << "AMRLevel::define (level" << m_level << ')'
 *           << CRD::h1;                                                \endcode
 *
 *  It is suggested to always use 'CRD::end' as a terminatory when message
 *  type has been specified using a manipulator.
 *
 *  IMPORTANT:  You have to set the stream 'setOutputStream' before using this.
 *              There is no checking to make sure it is set.
 *
 *//*+*************************************************************************/

class Msg : public std::ostringstream
{
public:

  // There should be only one instance of class Msg
  static Msg& getInstance()
    {
      static Msg instance;
      return instance;
    }

  // Use synthesized destructor

private:

  /// Default constructor
  Msg();

  // Copy and assignment prohibited
  Msg(const Msg&);
  Msg& operator=(const Msg&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  // Modification of default headers and setup
  void setLineSize(const size_t a_lineSize);
  void setOutputStream(std::ostream& a_out);
  void setErrorHeader(const char* a_header);
  void setWarnHeader(const char* a_header);
  void setHeader1(const char* a_header);
  void setHeader2(const char* a_header);
  void setBodyHeader(const char* a_header);
  void setVerbosityHeader(const char* a_header, const int a_idxVerb);
  void setDefaultPrec(const int a_int);
  void setVerbosity(const int a_verbosity);

  // Retrieval of values
  size_t lineSize() const { return m_lineSize; }
  std::ostream& pout() const { return *m_out; }

  /// Write a newline
  void newline(const int a_nlVerb = 1) const;

  /// Put bytes in units
  void putByteUnits(const std::size_t a_bytes, const int a_nlVerb = 1);

  // Formatting
  void setFloatDefault();
  void setFloatSN();
  void setMaxPrecFloatSN();
  void setPrecFloatSN(const int a_prec);
  void setFloatFP();
  void setMaxPrecFloatFP();
  void setPrecFloatFP(const int a_prec);

  /// Set termination flag
  void setTerminateOnError(const bool a_termFlag);

  /// Whether to terminate on error or abort message
  bool terminateOnError() const { return m_terminateOnError; }

  /// Get number of errors issued
  int numErrors() const { return m_numError; }

  /// Get number of warnings issued
  int numWarnings() const { return m_numWarning; }

  /// Reset count of number of errors
  void resetNumErrors() { m_numError = 0; }

  /// Reset count of number of warnings
  void resetNumWarnings() { m_numWarning = 0; }

  /// General buffering of stream
  template <typename T>
  Msg& operator<<(const T& a_val);

  ///  Set type of message (and verbosity) using manipulator
  Msg& operator<<(const CRD::IOFormat a_fmt);

  /// Modify the verbosity of the message from type default
  Msg& operator<<(const CRD::verb& a_fmt);

  /// Inject a special string into the stream
  Msg& operator<<(const CRD::IOPrint a_fmt);

  /// Termination (writes the line)
  std::ostringstream& operator<<(CRD::IOTerm a_term);


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the verbosity based on message type
  void classifyVerbosity(const CRD::IOTerm a_term);

  /// Write a message with error header
  void writeError();
  
  /// Write a message with warning header
  void writeWarning();

  /// Write a message with standard headers
  void writeMessage(const char* const a_header);

  /// Write a message to describe a variable in form name:value
  void writeVariable();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  size_t m_lineSize;                  ///< Size for line breaks
  std::ostream *m_out;                ///< Stream to use
  std::string m_errorHeader;          ///< Header for error messages
  std::string m_warnHeader;           ///< Header for warning messages
  std::string m_header1;              ///< Header 1
  std::string m_header2;              ///< Header 2
  std::string m_bodyHeader;           ///< Standard text header
  // Variables ({[m_bodyHeader]|name[m_varNameTerm][m_varHeaderFill]*|
  //             [m_varHeaderTerm]}variable) where {} denotes header and
  //             || denotes m_varNameSize characters.
  char        m_varHeaderFill;        ///< Fill character for variable header
  std::string m_varNameTerm;          ///< Termination of variable name
  std::string m_varHeaderTerm;        ///< Termination of variable header
  size_t m_varNameSize;               ///< Size of variable name in header
                                      ///< (includes size of varNameTerm but
                                      ///< excludes sizes of varHeaderTerm and
                                      ///< m_bodyHeader)
  std::string m_verbHeader[5];        ///< Verbosity headers (0-4)
  int m_defaultPrec;                  ///< Default precision for floating point
  int m_verbosity;                    ///< Amount of output requested
  int m_msgVerb;                      ///< Vebosity of a specific message.
                                      ///< Comparing this to m_verbosity
                                      ///< determines if the message should be
                                      ///< printed
  CRD::IOTerm m_msgType;              ///< Type specified by termination
  int m_numError;                     ///< Number of warnings printed
  int m_numWarning;                   ///< Number of errors printed

  bool m_terminateOnError;            ///< T - Terminate the program if error
                                      ///<     (default)
                                      ///< F - Continue running if error
};

//--Definition of Msg for writing messages

extern Msg& msg;

}  // namespace CRD


/*******************************************************************************
 *
 * Class Msg: inline member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Write a newline
/** Note that this doesn't put a newline in the buffered text, rather
 *  it places it directly on the stream.  It can be used between
 *  messages for nice whitespace.
 *  \param[in]  a_nlVerb
 *                      Verbosity at which to print (default is 1)
 *//*-----------------------------------------------------------------*/

inline void
CRD::Msg::newline(const int a_nlVerb) const
{
  if (a_nlVerb <= m_verbosity)
    {
      (*m_out) << std::endl;
    }
}

/*--------------------------------------------------------------------*/
//  Put bytes in units
/** \param[in]  a_byets Number of bytes
 *  \param[in]  a_nlVerb
 *                      Verbosity at which to print (default is 1)
 *//*-----------------------------------------------------------------*/

inline void
CRD::Msg::putByteUnits(const std::size_t a_bytes, const int a_nlVerb)
{
  if (a_nlVerb <= m_verbosity)
    {
      static constexpr const char* units[] =
      {
        " B",    // 0
        " KiB",  // 1
        " MiB",  // 2
        " GiB",  // 3
        " TiB"   // 4
      };
      
      double bytes;
      int lblUnits;
      if (a_bytes < 1024ull)
        {
          bytes = a_bytes;
          lblUnits = 0;
        }
      else if (a_bytes < 1048576ull)
        {
          bytes = static_cast<double>(a_bytes)/1024.;
          lblUnits = 1;
        }
      else if (a_bytes < 1073741824ull)
        {
          bytes = static_cast<double>(a_bytes)/1048576.;
          lblUnits = 2;
        }
      else if (a_bytes < 1099511627776ull)
        {
          bytes = static_cast<double>(a_bytes)/1073741824.;
          lblUnits = 3;
        }
      else
        {
          bytes = static_cast<double>(a_bytes)/1099511627776.;
          lblUnits = 4;
        }
      (*static_cast<std::ostringstream*>(this)) << bytes << units[lblUnits];
    }
}

/*--------------------------------------------------------------------*/
/// Set default formatting for floating point output
/** The default formatting is general
 *//*--------------------------------------------------------------------*/

inline void
CRD::Msg::setFloatDefault()
{
  this->setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  this->precision(m_defaultPrec);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to scientific notation
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setFloatSN()
{
  this->setf(std::ios_base::scientific, std::ios_base::floatfield);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to scientific notation at maximum precision
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setMaxPrecFloatSN()
{
  this->setf(std::ios_base::scientific, std::ios_base::floatfield);
  this->precision(std::numeric_limits<Real>::digits10);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to scientific notation at specified precision
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setPrecFloatSN(const int a_prec)
{
  this->setf(std::ios_base::scientific, std::ios_base::floatfield);
  this->precision(a_prec);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to fixed point notation
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setFloatFP()
{
  this->setf(std::ios_base::fixed, std::ios_base::floatfield);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to fixed point notation at maximum precision
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setMaxPrecFloatFP()
{
  this->setf(std::ios_base::fixed, std::ios_base::floatfield);
  this->precision(std::numeric_limits<Real>::digits10);
}

/*--------------------------------------------------------------------*/
/// Set FP formatting to fixed point notation at specified precision
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::setPrecFloatFP(const int a_prec)
{
  this->setf(std::ios_base::fixed, std::ios_base::floatfield);
  this->precision(a_prec);
}

/*--------------------------------------------------------------------*/
//  Set termination flag
/** \param[in]  a_termFlag
 *                      T - Terminate on error or abort message
 *                      F - Continue running on error or abort message
 *//*-----------------------------------------------------------------*/

inline void
CRD::Msg::setTerminateOnError(const bool a_termFlag)
{
  m_terminateOnError = a_termFlag;
}

/*--------------------------------------------------------------------*/
//  General buffering of stream
/** \tparam     T       Type to buffer
 *  \param[in]  a_val   Value to buffer
 *  \note
 *  <ul>
 *    <li> I don't know why, but if we don't put it in we get
 *         ambiguous operator<< errors.
 *  </ul>
 *//*-----------------------------------------------------------------*/

template <typename T>
inline CRD::Msg&
CRD::Msg::operator<<(const T& a_val)
{
  // Don't even bother accumulating in buffer if the verbosity is too high
  // for printing
  if (m_msgVerb <= m_verbosity)
    {
      (*static_cast<std::ostringstream*>(this)) << a_val;
    }
  return *this;
}

/*--------------------------------------------------------------------*/
//  Set type of message (and verbosity) using manipulator
/** \param[in]  a_fmt   Manipulator
 *//*-----------------------------------------------------------------*/

inline CRD::Msg&
CRD::Msg::operator<<(const CRD::IOFormat a_fmt)
{
  m_msgType = static_cast<CRD::IOTerm>(a_fmt);
  classifyVerbosity(m_msgType);
  return *this;
}

/*--------------------------------------------------------------------*/
//  Modify the verbosity of the message from type default
/** \param[in]  a_fmt   Manipulator setting verbosity
 *//*-----------------------------------------------------------------*/

inline CRD::Msg&
CRD::Msg::operator<<(const CRD::verb& a_fmt)
{
  m_msgVerb = a_fmt.m_msgVerb;
  return *this;
}


/*--------------------------------------------------------------------*/
//  Inject a special string into the stream
/** \param[in]  a_fmt   Type of special string
 *//*-----------------------------------------------------------------*/

inline CRD::Msg&
CRD::Msg::operator<<(CRD::IOPrint a_fmt)
{
  static const char *const plural[][2] =
    {
      { "were "    , "was "     },
      { " errors"  , " error"   },
      { " warnings", " warning" }
    };
  // Don't even bother accumulating in buffer if the verbosity is too high for
  // printing
  if (m_msgVerb <= m_verbosity)
    {
      switch (a_fmt)
        {
        case CRD::wasNumError:
        {
          const int oneError = (int)(m_numError == 1);
          (*static_cast<std::ostringstream*>(this)) << plural[0][oneError];
          (*static_cast<std::ostringstream*>(this)) << m_numError;
          (*static_cast<std::ostringstream*>(this)) << plural[1][oneError];
          break;
        }
        case CRD::wasNumWarning:
        {
          const int oneWarn = (int)(m_numWarning == 1);
          (*static_cast<std::ostringstream*>(this)) << plural[0][oneWarn];
          (*static_cast<std::ostringstream*>(this)) << m_numWarning;
          (*static_cast<std::ostringstream*>(this)) << plural[2][oneWarn];
          break;
        }
        }
    }
  return *this;
}

/*--------------------------------------------------------------------*/
//  Termination (writes the line)
/** This is where the magic happens
 *  \param[in]  a_term  Type of terminated or 'CRD::end' for generic
 *                      termination
 *//*-----------------------------------------------------------------*/

inline std::ostringstream&
CRD::Msg::operator<<(CRD::IOTerm a_term)
{
  if (a_term == CRD::end)  // Used stored message type (possibly set by
    {                      // manipulator
      a_term = m_msgType;
    }
  else                     // Use verbosity for type
    {
      classifyVerbosity(a_term);
    }
  if (m_msgVerb <= m_verbosity)
    {
      const char* header = 0;
      switch (a_term)
        {
        case CRD::warn:
          writeWarning();
          goto reset;
          break;
        case CRD::error:
          writeError();
          if (terminateOnError())
            {
              // ::exit(1);
              ::abort();
            }
          goto reset;
          break;
        case CRD::abort:
          writeError();
          if (terminateOnError())
            {
              ::abort();
            }
          goto reset;
          break;
        case CRD::h1:
          header = m_header1.c_str();
          break;
        case CRD::h2:
          header = m_header2.c_str();
          break;
        case CRD::body:
          header = m_bodyHeader.c_str();
          break;
        case CRD::var:
        case CRD::varend:
          writeVariable();
          goto reset;
          break;
        case CRD::v0:
        case CRD::v1:
        case CRD::v2:
        case CRD::v3:
        case CRD::v4:
          header = m_verbHeader[(int)a_term - (int)CRD::v0].c_str();
          break;
        case CRD::end:  // To avoid warnings
          break;
        }
      writeMessage(header);
    }
  reset: ;
  str("");
  m_msgVerb = m_verbosity;  // Force buffering by default
  m_msgType = CRD::body;    // Will be a body message by default
  return *this;
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the verbosity based on message type
/** \param[in]  a_type  Message type
 *//*-----------------------------------------------------------------*/

inline void
CRD::Msg::classifyVerbosity(const CRD::IOTerm a_type)
{
  switch (a_type)
    {
    case CRD::warn:
    case CRD::error:
    case CRD::abort:
      m_msgVerb = -1;
      break;
    case CRD::h1:
    case CRD::v0:
      m_msgVerb = 0;
      break;
    case CRD::h2:
    case CRD::body:
    case CRD::var:
    case CRD::v1:
      m_msgVerb = 1;
      break;
    case CRD::v2:
      m_msgVerb = 2;
      break;
    case CRD::v3:
      m_msgVerb = 3;
      break;
    case CRD::v4:
      m_msgVerb = 4;
      break;
    case CRD::varend:
    case CRD::end:  // To avoid warnings
      break;
    }
}

/*--------------------------------------------------------------------*/
//  Write a message with error header
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::writeError()
{
  ++m_numError;
  linePrint(str().c_str(), *m_out, m_errorHeader.c_str(), -1, m_lineSize);
}

/*--------------------------------------------------------------------*/
//  Write a message with warning header
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::writeWarning()
{
  ++m_numWarning;
  linePrint(str().c_str(), *m_out, m_warnHeader.c_str(), -1, m_lineSize);
}

/*--------------------------------------------------------------------*/
//  Write a message with standard headers
/*--------------------------------------------------------------------*/

inline void
CRD::Msg::writeMessage(const char *const a_header)
{
  linePrint(str().c_str(), *m_out, a_header, -1, m_lineSize);
}

#endif /* ! defined _CRDMSG_H_ */
