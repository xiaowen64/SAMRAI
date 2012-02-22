/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Parallel I/O class buffer to manage parallel ostreams output
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Set the prefix that begins every new output line.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void ParallelBuffer::setPrefixString(
   const std::string& text)
{
   d_prefix = text;
}

/*
 *************************************************************************
 *
 * Set the primary output stream.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void ParallelBuffer::setOutputStream1(
   std::ostream* stream)
{
   d_ostream1 = stream;
}

/*
 *************************************************************************
 *
 * Set the secondary output stream.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void ParallelBuffer::setOutputStream2(
   std::ostream* stream)
{
   d_ostream2 = stream;
}

/*
 *************************************************************************
 *
 * Output a string to the output stream by invoking the
 * outputString(string,length) method.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void ParallelBuffer::outputString(
   const std::string& text)
{
   outputString(text, static_cast<int>(text.length()));
}

}
}
