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
 * \file TagLevelFactory.cpp
 *
 * \brief Member functions for TagLevelFactory
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <map>

//----- Chombo Library -----//

#include "ParmParse.H"

#include "UsingNamespace.H"

//----- Internal -----//

#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "CRDparam.H"
#include "CRDmsg.H"
#include "CRDPhysics.H"
#include "TagMethodBuffer.H"
#include "TagMethodPhysBox.H"
#include "TagMethodBaseBox.H"
#include "TagMethodBoundary.H"
#include "TagMethodValue.H"
#include "TagMethodGradient.H"
#include "TagMethodTBL.H"


/*******************************************************************************
 *
 * Class TagLevelFactory: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Construct by reading the input file
/**
 *//*-----------------------------------------------------------------*/

TagLevelFactory::TagLevelFactory(const int a_tagBufferSize)
{
  CH_assert(CRDparam::maxAMRLevel() < 4096);  // For safe sprintf

  if (CRDparam::maxAMRLevel() == 0)
    // No tagging is needed.  Create a default method of only adding a buffer.
    {
      TagLevel* tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      m_tagLevelVec.resize(1);
      m_levelMapVec.resize(1);
      m_tagLevelVec[0] = RefCountedPtr<TagLevel>(tagLevel);
      m_levelMapVec[0] = 0;
      return;
    }

  CRD::msg.setTerminateOnError(false);
  CRD::msg << "Parameters for 'tagging'" << CRD::h2;
  ParmParse ppTAG("tag");

  std::map<int, TagLevel*> tagLevels;
  if (ppTAG.contains("all_prepend") || ppTAG.contains("all_append"))
    {
      for (int level = 0; level < CRDparam::maxAMRLevel(); ++level)
        {
          tagLevels.insert({level, new TagLevel});
        }
    }
  // Prepend method to all levels
  if (ppTAG.contains("all_prepend"))
    {
      const int numMethod = ppTAG.countname("all_prepend");
      for (int idxMethod = 0; idxMethod != numMethod; ++idxMethod)
        {
          PP_String pp_string;
          int status = ppTAG.queryval(
            "all_prepend", ParmParse::ppString, &pp_string, 0, idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Error reading tag method for prepending to "
                "all levels" << CRD::error;
            }
          const char* methodName = pp_string.c_str();
          readTagMethod(ppTAG, "all_prepend", -1, idxMethod, methodName,
                        tagLevels);
        }
    }
  // Specific methods for specific levels
  char tagLevelName[32] = { '\0' };
  for (int level = 0; level < CRDparam::maxAMRLevel(); ++level)
    {
      sprintf(tagLevelName, "level_%d", level);
      const int numMethod = ppTAG.countname(tagLevelName);
      if (numMethod > 0)
        {
          if (tagLevels.find(level) == tagLevels.end())
            {
              tagLevels.insert({level, new TagLevel});
            }
        }
      for (int idxMethod = 0; idxMethod != numMethod; ++idxMethod)
        {
          PP_String pp_string;
          int status = ppTAG.queryval(
            tagLevelName, ParmParse::ppString, &pp_string, 0, idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Error reading tag method for level " << level
                       << '!' << CRD::error;
            }
          const char* methodName = pp_string.c_str();
          readTagMethod(ppTAG, tagLevelName, level, idxMethod, methodName,
                        tagLevels);
        }
    }
  // Append method to all levels
  if (ppTAG.contains("all_append"))
    {
      const int numMethod = ppTAG.countname("all_append");
      for (int idxMethod = 0; idxMethod != numMethod; ++idxMethod)
        {
          PP_String pp_string;
          int status = ppTAG.queryval(
            "all_append", ParmParse::ppString, &pp_string, 0, idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Error reading tag method for appending to "
                "all levels" << CRD::error;
            }
          const char* methodName = pp_string.c_str();
          readTagMethod(ppTAG, "all_append", -1, idxMethod, methodName,
                        tagLevels);
        }
    }
  // Append a tag-buffer method to all levels
  if (a_tagBufferSize > 0)
    {
      for (auto& iter : tagLevels)
        {
          iter.second->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
          CRD::msg << "Tag (level"  << iter.first << ") method buffer\n"
                   << a_tagBufferSize << CRD::var;
        }
    }

  // Preliminary analysis
  // There should be at least 1 tag level
  const int numTagLevel = tagLevels.size();
  if (numTagLevel < 1)
    {
      CRD::msg << "Input: Tagging was not defined on any level" << CRD::error;
    }
  // There should be a tagLevel for the base grid otherwise there is no way to
  // create higher levels
  if (tagLevels.begin()->first != 0)
    {
      CRD::msg << "Input: Tagging must be defined for level 0" << CRD::error;
    }

  // Set to terminate on error
  CRD::msg.setTerminateOnError(true);
  CRD::msg.newline();
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the input file for "
        "cell tag information.  Aborting now." << CRD::error;
    }

  m_tagLevelVec.resize(numTagLevel);
  m_levelMapVec.resize(numTagLevel);
  int c = 0;
  for (auto& iter : tagLevels)
    {
      m_levelMapVec[c] = iter.first;
      m_tagLevelVec[c] = RefCountedPtr<TagLevel>(iter.second);
      m_tagLevelVec[c]->FIXMEsetTagBufferSize(a_tagBufferSize);
      ++c;
    }
}

/*--------------------------------------------------------------------*/
//  Constructor for single TagLevel object
/** The TagLevel object is created externally and copied here.  This
 *  constructor creates only one tagLevel object.
 *  \param[in]  a_tagLevel
 *                      Object that tags cells on an AMR level.
 *                      Allocate with new and do not delete.
 *//*-----------------------------------------------------------------*/

TagLevelFactory::TagLevelFactory(TagLevel *const a_tagLevel)
  :
  m_tagLevelVec(1),
  m_levelMapVec(1)
{
  m_tagLevelVec[0] = RefCountedPtr<TagLevel>(a_tagLevel);
  m_levelMapVec[0] = 0;
}

/*--------------------------------------------------------------------*/
//  Constructor for multiple TagLevel objects
/** The TagLevel object is created externally and copied here.
 *  \param[in]  a_tagLevelVec
 *                      Vector of TagLevel objects.  These will be
 *                      assigned to do tagging on levels according to
 *                      the map
 *  \param[in]  a_levelMapVec
 *                      A vector of AMR level indices.  If an AMR
 *                      level index (idx) is in the range,
 *                      (a_levelMapVec[i] <= idx < a_levelMapVec[i+1])
 *                      then a_tagLevelVec[i] is assigned to that
 *                      level.
 *  \note
 *  <ul>
 *    <li> The TagLevel objects should be allocated with new and not
 *         deleted
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory::TagLevelFactory(
  const std::vector<TagLevel*> a_tagLevelVec,
  const std::vector<int>       a_levelMapVec)
  :
  m_tagLevelVec(a_tagLevelVec.size()),
  m_levelMapVec(a_levelMapVec.size())
{
  CH_assert(a_tagLevelVec.size() > 0);
  CH_assert(a_levelMapVec.size() > 0);
  CH_assert(a_tagLevelVec.size() == a_levelMapVec.size());
  CH_assert(a_levelMapVec[0] == 0);
  for (int i = 0; i != a_tagLevelVec.size(); ++i)
    {
      m_tagLevelVec[i] = RefCountedPtr<TagLevel>(a_tagLevelVec[i]);
      m_levelMapVec[i] = a_levelMapVec[i];
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagLevelFactory::~TagLevelFactory()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Create a new TagLevel object
/** In the default implementation, all levels will receive a pointer
 *  to the same object.
 *  \param[in]  a_level Index of this level with coarsest defined
 *                      as 0.
 *  \param[in]  a_fnRefRatio
 *                      Refinement ratio between this level and the
 *                      next finer level.
 *  The default implementation uses the copy constructor to create
 *  a copy of the member TagLevel object.
 *//*-----------------------------------------------------------------*/

const RefCountedPtr<TagLevel>&
TagLevelFactory::new_tagLevel(const int a_level,
                              const int a_fnRefRatio) const
{
  CH_assert(m_levelMapVec[0] == 0);  // Checked in constructor but check again
                                     // since following loop may be infinite
                                     // otherwise
  /// Find 'i' such that a_levelMapVec[i] <= a_level < a_levelMapVec[i+1]
  int i = m_levelMapVec.size() - 1;
  while (a_level < m_levelMapVec[i]) --i;
  return m_tagLevelVec[i];
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Read details of a tag method and add to appropriate level(s)
/** \param[in]  a_ppTAG ParmParse object
 *  \param[in]  a_level Level for tag method (-1 means all levels)
 *  \param[in]  a_idxMethod
 *                      Index of the method (as determined by
 *                      ParmParse)
 *  \param[in]  a_methodName
 *                      Name of tagging method from input file
 *  \param[in]  a_tagLevels
 *                      Level a_level must have been inserted
 *  \param[out] a_tagLevels
 *                      New tag method inserted
 *//*-----------------------------------------------------------------*/

void
TagLevelFactory::readTagMethod(ParmParse&                a_ppTAG,
                               const char*               a_tagLevelName,
                               const int                 a_level,
                               const int                 a_idxMethod,
                               const char*               a_methodName,
                               std::map<int, TagLevel*>& a_tagLevels)
{

/*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 *
 * If adding new methods, spend time to make sure the processing of input is
 * correct and feedback is clear.  Tagging can be confusing.
 *
 *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*/

//--Method: Physical box

  if (strcmp(a_methodName, "physical_box") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      int physVecSize = SpaceDim;
      bool haveLastArg = (numArgs % 2 == 0);
      // If numArgs = 7,8 assume 3D vectors
      // If numArgs = 5,6 assume 2D vectors (require SpaceDim = 2)
      // If numArgs = 3,4 assume 1D vectors (require SpaceDIm = 1)
      switch (numArgs)
        {
        case 8:
        case 7:
          physVecSize = 3;
          break;
        case 6:
        case 5:
          physVecSize = 2;
          if (SpaceDim > 2)
            {
              CRD::msg << "Input: Tag method 'physical_box', insufficient "
                "vector sizes in arguments" << CRD::error;
            }
          break;
        case 4:
        case 3:
          physVecSize = 1;
          if (SpaceDim > 1)
            {
              CRD::msg << "Input: Tag method 'physical_box', insufficient "
                "vector sizes in arguments" << CRD::error;
            }
          break;
        case 2:
        case 1:
          CRD::msg << "Input: Tag method 'physical_box', insufficient "
            "arguments" << CRD::error;
          break;
        default:
          // High-dim users are more reliable?
          break;
        }
      // Low corner of box in physical space
      stc::Vector<double, 3> lo;
      int status = a_ppTAG.queryarr(
        a_tagLevelName, ParmParse::ppDouble, lo.dataPtr(),
        1,
        physVecSize,
        a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'physical_box', error reading lower "
            "corner" << CRD::error;
        }
      // High corner of box in physical space
      stc::Vector<double, 3> hi;
      status = a_ppTAG.queryarr(
        a_tagLevelName, ParmParse::ppDouble, hi.dataPtr(),
        1 + physVecSize,
        physVecSize,
        a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'physical_box', error reading upper "
            "corner" << CRD::error;
        }
      // Check all cells or fast exclusion of boxes
      bool fastExclusion = true;
      if (haveLastArg)
        {
          PP_String checkAllCells;
          status = a_ppTAG.queryval(
            a_tagLevelName, ParmParse::ppString, &checkAllCells,
            1 + 2*physVecSize,
            a_idxMethod);
          if (status != 0 &&
              strcmp(checkAllCells.c_str(), "check_all_cells") == 0)
            {
              fastExclusion = false;
            }
        }
      const RealVect rlo(lo);
      const RealVect rhi(hi);
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodPhysBox(rlo, rhi, fastExclusion));
            }
          CRD::msg << "Tag (level" << level << ") physical box\n";
          CRD::msg << rlo << ':' << rhi;
          if (!fastExclusion) CRD::msg << " check all cells";
          CRD::msg << CRD::var;
        }
    }

//--Method: Base box

  else if (strcmp(a_methodName, "base_box") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      int vecSize = SpaceDim;
      bool haveLastArg = (numArgs % 2 == 0);
      // If numArgs = 7,8 assume 3D vectors
      // If numArgs = 5,6 assume 2D vectors (require SpaceDim = 2)
      // If numArgs = 3,4 assume 1D vectors (require SpaceDIm = 1)
      switch (numArgs)
        {
        case 8:
        case 7:
          vecSize = 3;
          break;
        case 6:
        case 5:
          vecSize = 2;
          if (SpaceDim > 2)
            {
              CRD::msg << "Input: Tag method 'base_box', insufficient "
                "vector sizes in arguments" << CRD::error;
            }
          break;
        case 4:
        case 3:
          vecSize = 1;
          if (SpaceDim > 1)
            {
              CRD::msg << "Input: Tag method 'base_box', insufficient "
                "vector sizes in arguments" << CRD::error;
            }
          break;
        case 2:
        case 1:
          CRD::msg << "Input: Tag method 'base_box', insufficient "
            "arguments" << CRD::error;
          break;
        default:
          // High-dim users are more reliable?
          break;
        }
      // Low corner of base box space
      stc::Vector<int, 3> lo;
      int status = a_ppTAG.queryarr(
        a_tagLevelName, ParmParse::ppInt, lo.dataPtr(),
        1,
        vecSize,
        a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'base_box', error reading lower "
            "corner" << CRD::error;
        }
      // High corner of box in physical space
      stc::Vector<int, 3> hi;
      status = a_ppTAG.queryarr(
        a_tagLevelName, ParmParse::ppInt, hi.dataPtr(),
        1 + vecSize,
        vecSize,
        a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'base_box', error reading upper "
            "corner" << CRD::error;
        }
      TagMethodBaseBox::Mode mode = TagMethodBaseBox::ModeAddTagsInsideBox;
      if (haveLastArg)
        {
          PP_String modeStr;
          status = a_ppTAG.queryval(
            a_tagLevelName, ParmParse::ppString, &modeStr,
            1 + 2*vecSize,
            a_idxMethod);
          if (status != 0)
            {
              if (strcmp(modeStr.c_str(), "add_inside") == 0)
                {
                  mode = TagMethodBaseBox::ModeAddTagsInsideBox;
                }
              else if (strcmp(modeStr.c_str(), "exclude_outside") == 0)
                {
                  mode = TagMethodBaseBox::ModeExcludeTagsOutsideBox;
                }
              else
                {
                  CRD::msg << "Input: Tag method 'base_box', unknown mode: "
                           << modeStr.c_str() << CRD::error;
                }
            }
          else
            {
              CRD::msg << "Input: Tag method 'base_box', error reading mode"
                       << CRD::error;
            }
        }
      const IntVect ilo(lo);
      const IntVect ihi(hi);
      const Box tagBox(ilo, ihi);
      if (tagBox.isEmpty())
        {
          CRD::msg << "Input: Tag method 'base_box', empty tag box: " << tagBox
                   << CRD::error;
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodBaseBox(tagBox, mode));
            }
          CRD::msg << "Tag (level" << level << ") base box\n";
          switch (mode)
            {
            case TagMethodBaseBox::ModeAddTagsInsideBox:
              CRD::msg << "add cells inside ";
              break;
            case TagMethodBaseBox::ModeExcludeTagsOutsideBox:
              CRD::msg << "exclude cells outside ";
              break;
            }
          CRD::msg << ilo << ':' << ihi;
          CRD::msg << CRD::var;
        }
    }

//--Method: Boundary

  else if (strcmp(a_methodName, "boundary") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      if (numArgs < 4 || numArgs > 5)
        {
          CRD::msg << "Input: Tag method 'boundary', insufficient arguments "
                   << "(block_name dir side [grow_normal])"
                   << CRD::error;
        }
      // Block index from block name
      PP_String stringPP;
      int status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 1, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading block name"
                   << CRD::error;
        }
      std::string blockName(stringPP.c_str());
      const int idxBlk = CRDparam::g_coordSys->whichBlock(blockName);
      if (idxBlk == -1)
        {
          CRD::msg << "Input: Tag method 'boundary', unable to find block with "
            "name " << blockName << CRD::error;
          CRD::msg << "Known names:" << CRD::end;
          for (int i = 0, i_end = CRDparam::g_coordSys->numBlocks(); i != i_end;
               ++i)
            {
              CRD::msg << "  " << std::setw(4) << i << ": "
                       << CRDparam::g_coordSys->blockInfo(i).m_name << CRD::end;
            }
        }
      // Direction
      int dir = 0;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppInt, &dir, 2, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading direction"
                   << CRD::error;
        }
      if (dir < 0 || dir >= SpaceDim)
        {
          CRD::msg << "Input: Tag method 'boundary', requires 0 <= dir < "
                   << SpaceDim << " (read " << dir << ')' << CRD::error;
        }
      // Side
      Side::LoHiSide side = Side::Lo;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 3, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading side"
                   << CRD::error;
        }
      if (std::strncmp(stringPP.c_str(), "lo", 2) == 0 ||
          std::strncmp(stringPP.c_str(), "min", 3) == 0 ||
          std::strncmp(stringPP.c_str(), "0", 1) == 0 ||
          std::strncmp(stringPP.c_str(), "-1", 2) == 0)
        {
          side = Side::Lo;
        }
      else if (std::strncmp(stringPP.c_str(), "hi", 2) == 0 ||
               std::strncmp(stringPP.c_str(), "max", 3) == 0 ||
               std::strncmp(stringPP.c_str(), "1", 1) == 0)
        {
          side = Side::Hi;
        }
      else
        {
          CRD::msg << "Input: Tag method 'boundary', " << stringPP.c_str()
                   << " unrecognized as a term for side" << CRD::error;
        }
      // Growth in normal direction
      int growNormal = 1;
      if (numArgs == 5)
        {
          status = a_ppTAG.queryval(
            a_tagLevelName, ParmParse::ppInt, &growNormal, 4, a_idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Tag method 'boundary', error reading growth "
                "in normal direction" << CRD::error;
            }
          if (growNormal < 1)
            {
              CRD::msg << "Input: Tag method 'boundary', requires "
                "'grow_normal' < 1 (read " << growNormal << ')' << CRD::error;
              growNormal = 1;
            }
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      static constexpr const char* lblSide[] = { " low", " high" };
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodBoundary(idxBlk, dir, side, growNormal));
            }
          CRD::msg << "Tag (level" << level << ") boundary\n"
                   << "block " << blockName << ", face " << dir << lblSide[side]
                   << ", grow " << growNormal << CRD::var;
        }
    }

//--Method: Buffer

  else if (strcmp(a_methodName, "buffer") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      if (numArgs != 2)
        {
          CRD::msg << "Input: Tag method 'buffer', insufficient arguments"
                   << "(buffer_size)"
                   << CRD::error;
        }
      int tagBufferSize = 0;
      int status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppInt, &tagBufferSize, 1, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'buffer', error reading buffer size"
                   << CRD::error;
        }
      if (tagBufferSize < 0)
        {
          CRD::msg << "Input: Tag method 'buffer', requires buffer size >= 0"
                   << CRD::error;
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodBuffer(tagBufferSize));
            }
          CRD::msg << "Tag (level" << level << ") buffer\n"
                   << tagBufferSize << CRD::var;
        }
    }

//--Method: Value

  else if (strcmp(a_methodName, "value") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      if (numArgs < 4 || numArgs > 5)
        {
          CRD::msg << "Input: Tag method 'value', insufficient arguments"
                   << "(comp lower_threshold upper_threshold [outside])"
                   << CRD::error;
        }
      // Component (with some support for names)
      PP_String stringPP;
      int status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 1, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'value', error reading component"
                   << CRD::error;
        }
      int idxComp = -1;
      TagMethodValue::ValueType valueType =
        TagMethodValue::ValueTypeConservative;
      std::string compName(stringPP.c_str());
      if (compName == "vorticity")
        {
          valueType = TagMethodValue::ValueTypeVorticity;
        }
      else if (compName == "hydrocarbon_flame")
        {
          valueType = TagMethodValue::ValueTypeHydrocarbonFlame;
        }
      else
        {
          for (int cU = 0, cU_end = CRDparam::g_CRDPhysics->numConservative();
               cU != cU_end; ++cU)
            {
              if (compName == CRDparam::g_CRDPhysics->consvStateName(cU))
                {
                  idxComp = cU;
                  valueType = TagMethodValue::ValueTypeConservative;
                  break;
                }
            }
          if (idxComp == -1)
            {
              for (int cW = 0, cW_end = CRDparam::g_CRDPhysics->numPrimitive();
                   cW != cW_end; ++cW)
                {
                  if (compName == CRDparam::g_CRDPhysics->primStateName(cW))
                    {
                      idxComp = cW;
                      valueType = TagMethodValue::ValueTypePrimitive;
                      break;
                    }
                }
            }
          if (idxComp == -1)
            // Try to convert to a numeric index (interpreted as primitive
            // state)
            {
              valueType = TagMethodValue::ValueTypePrimitive;
              try { idxComp = std::stoi(compName); }
              catch (...)
                {
                  CRD::msg << "Input: Tag method 'value', unrecognized "
                    "component" << CRD::error;
                }
              compName.clear();
              if (idxComp < 0 ||
                  idxComp >= CRDparam::g_CRDPhysics->numPrimitive())
                {
                  CRD::msg << "Input: Tag method 'value', component not in "
                    "range of primitive variables 0:"
                           << CRDparam::g_CRDPhysics->numPrimitive() - 1
                           << CRD::error;
                }
              else
                {
                  compName = CRDparam::g_CRDPhysics->primStateName(idxComp);
                }
            }
        }
      // Lower threshold
      double lowerThresh = 0.;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppDouble, &lowerThresh, 2, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'value', error reading lower threshold"
                   << CRD::error;
        }
      // Upper threshold
      double upperThresh = 0.;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppDouble, &upperThresh, 3, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'value', error reading upper threshold"
                   << CRD::error;
        }
      if (lowerThresh > upperThresh)
        {
          CRD::msg << "Input: Tag method 'value', requires lower threshold <= "
            "upper threshold" << CRD::error;
        }
      // Inside or outside
      bool inside = true;
      if (numArgs == 5)
        {
          status = a_ppTAG.queryval(
            a_tagLevelName, ParmParse::ppString, &stringPP, 4, a_idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Tag method 'value', error reading threshold "
                "side" << CRD::error;
            }
          if (std::strncmp(stringPP.c_str(), "out", 3) == 0)
            {
              inside = false;
            }
          else if (std::strncmp(stringPP.c_str(), "in", 2) != 0)
            {
          CRD::msg << "Input: Tag method 'value', unrecognized side '"
                   << stringPP.c_str() << "'" << CRD::error;
            }
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodValue(idxComp,
                                   lowerThresh,
                                   upperThresh,
                                   inside,
                                   valueType));
            }
          static constexpr const char* lblSide[] = { " outside", " inside" };
          CRD::msg << "Tag (level" << level << ") value\n"
                   << compName << lblSide[inside] << ' ' << lowerThresh
                   << ':' << upperThresh << CRD::var;
        }
    }

//--Method: Gradient

  else if (strcmp(a_methodName, "gradient") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      if (numArgs < 3 || numArgs > 4)
        {
          CRD::msg << "Input: Tag method 'gradient', insufficient arguments"
                   << "(comp threshold [low_value])"
                   << CRD::error;
        }
      // Component (with some support for names)
      PP_String stringPP;
      int status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 1, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'gradient', error reading component"
                   << CRD::error;
        }
      int idxComp = -1;
      bool usePrimitive = true;
      std::string compName(stringPP.c_str());
      for (int cU = 0, cU_end = CRDparam::g_CRDPhysics->numConservative();
           cU != cU_end; ++cU)
        {
          if (compName == CRDparam::g_CRDPhysics->consvStateName(cU))
            {
              idxComp = cU;
              usePrimitive = false;
              break;
            }
        }
      if (idxComp == -1)
        {
          for (int cW = 0, cW_end = CRDparam::g_CRDPhysics->numPrimitive();
               cW != cW_end; ++cW)
            {
              if (compName == CRDparam::g_CRDPhysics->primStateName(cW))
                {
                  idxComp = cW;
                  break;
                }
            }
        }
      if (idxComp == -1)
        // Try to convert to a numeric index (interpreted as primitive
        // state)
        {
          try { idxComp = std::stoi(compName); }
          catch (...)
            {
              CRD::msg << "Input: Tag method 'gradient', unrecognized "
                "component" << CRD::error;
            }
          compName.clear();
          if (idxComp < 0 ||
              idxComp >= CRDparam::g_CRDPhysics->numPrimitive())
            {
              CRD::msg << "Input: Tag method 'gradient', component not in "
                "range of primitive variables 0:"
                       << CRDparam::g_CRDPhysics->numPrimitive() - 1
                       << CRD::error;
            }
        }
      // Threshold
      double threshold = 0.;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppDouble, &threshold, 2, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'gradient', error reading threshold"
                   << CRD::error;
        }
      // Lower value
      double lowVal = 0.;
      if (numArgs == 4)
        {
          double lowVal = 0.;
          status = a_ppTAG.queryval(
            a_tagLevelName, ParmParse::ppDouble, &lowVal, 3, a_idxMethod);
          if (status == 0)
            {
              CRD::msg << "Input: Tag method 'gradient', error reading lower "
                "value" << CRD::error;
            }
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodGradient(idxComp,
                                      threshold,
                                      usePrimitive,
                                      lowVal));
            }
          if (compName.size() == 0)
            {
              compName = CRDparam::g_CRDPhysics->primStateName(idxComp);
            }
          // static constexpr const char* lblSide[] = { " outside", " inside" };
          CRD::msg << "Tag (level" << level << ") gradient\n"
                   << compName << ' ' << threshold;
          if (lowVal != 0.)
            {
              CRD::msg << ' ' << lowVal;
            }
          CRD::msg << CRD::var;
        }
    }

//--Method: Boundary-Layer

  else if (strcmp(a_methodName, "boundary_layer") == 0)
    {
      const int numArgs = a_ppTAG.countval(a_tagLevelName, a_idxMethod);
      if (numArgs != 7)
        {
          CRD::msg << "Input: Tag method 'boundary_layer', "
                   << "insufficient arguments "
                   << "(block_name dir side streamwise_dir start_delta "
                   << "Reynolds_number)"
                   << CRD::error;
        }
      // Block index from block name
      PP_String stringPP;
      int status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 1, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading block name"
                   << CRD::error;
        }
      std::string blockName(stringPP.c_str());
      const int idxBlk = CRDparam::g_coordSys->whichBlock(blockName);
      if (idxBlk == -1)
        {
          CRD::msg << "Input: Tag method 'boundary', unable to find block with "
            "name " << blockName << CRD::error;
        }
      // Direction
      int dir = 0;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppInt, &dir, 2, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading direction"
                   << CRD::error;
        }
      if (dir < 0 || dir >= SpaceDim)
        {
          CRD::msg << "Input: Tag method 'boundary', requires 0 <= dir < "
                   << SpaceDim << " (read " << dir << ')' << CRD::error;
        }
      // Side
      Side::LoHiSide side = Side::Lo;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppString, &stringPP, 3, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading side"
                   << CRD::error;
        }
      if (std::strncmp(stringPP.c_str(), "lo", 2) == 0 ||
          std::strncmp(stringPP.c_str(), "min", 3) == 0 ||
          std::strncmp(stringPP.c_str(), "0", 1) == 0 ||
          std::strncmp(stringPP.c_str(), "-1", 2) == 0)
        {
          side = Side::Lo;
        }
      else if (std::strncmp(stringPP.c_str(), "hi", 2) == 0 ||
               std::strncmp(stringPP.c_str(), "max", 3) == 0 ||
               std::strncmp(stringPP.c_str(), "1", 1) == 0)
        {
          side = Side::Hi;
        }
      else
        {
          CRD::msg << "Input: Tag method 'boundary', " << stringPP.c_str()
                   << " unrecognized as a term for side" << CRD::error;
        }
      // Streamwise direction
      int streamwiseDir = 0;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppInt, &streamwiseDir, 4, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading "
                   << "streamwiseDir" << CRD::error;
        }
      if (streamwiseDir < 0 || streamwiseDir >= SpaceDim)
        {
          CRD::msg << "Input: Tag method 'boundary', requires "
                   << "0 <= streamwiseDir < " << SpaceDim << " (read " << dir
                   << ')' << CRD::error;
        }
      // Starting boundary layer thickness
      double delta = 0;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppDouble, &delta, 5, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading initial "
            "boundary-layer thickness" << CRD::error;
        }
      if (delta < 0)
        {
          CRD::msg << "Input: Tag method 'boundary', requires "
            "'delta' >= 0 (read " << delta << ')' << CRD::error;
          delta = 0.;
        }
      // Starting Reynolds number (x-position Reynolds number)
      double Re = 0;
      status = a_ppTAG.queryval(
        a_tagLevelName, ParmParse::ppDouble, &Re, 6, a_idxMethod);
      if (status == 0)
        {
          CRD::msg << "Input: Tag method 'boundary', error reading initial "
            "Reynolds-number" << CRD::error;
        }
      if (Re < 0)
        {
          CRD::msg << "Input: Tag method 'boundary', requires "
            "'Re' >= 0 (read " << delta << ')' << CRD::error;
          Re = 0.;
        }
      int level = a_level;
      int level_end = a_level + 1;
      if (a_level == -1)
        {
          level = 0;
          level_end = CRDparam::maxAMRLevel();
        }
      static constexpr const char* lblSide[] = { " low", " high" };
      for (; level != level_end; ++level)
        {
          if (CRD::msg.numErrors() == 0)
            {
              a_tagLevels[level]->appendTagMethod(
                new TagMethodTBL(idxBlk, dir, side, streamwiseDir, delta, Re));
            }
          CRD::msg << "Tag (level" << level << ") boundary\n"
                   << "block " << blockName << ", face " << dir << lblSide[side]
                   << ", streamwiseDir " << streamwiseDir << ", delta "
                   << delta << ", ReynoldsNumber " << Re << CRD::var;
        }
    }

//--Method: Unrecognized

  else
    {
      CRD::msg << "Input: Unrecognized tag method '" << a_methodName << "' on ";
      if (a_level == -1)
        {
          CRD::msg << "all levels";
        }
      else
        {
          CRD::msg << "level " << a_level << CRD::error;
        }
      CRD::msg << CRD::error;
    }
}
