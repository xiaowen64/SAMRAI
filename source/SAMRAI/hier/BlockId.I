/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BlockId::BlockId():
   d_value(invalidId().d_value) {
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BlockId::BlockId(
   const BlockId& other):
   d_value(other.d_value) {
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BlockId::BlockId(
   const int& value):
   d_value(value) {
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BlockId::~BlockId()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   d_value = s_invalid_id.d_value;
#endif
}

/*
 ********************************************************************************
 * Assignment operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
BlockId& BlockId::operator = (const BlockId& rhs)
{
   d_value = rhs.d_value;
   return *this;
}

/*
 ********************************************************************************
 * Assignment operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BlockId::setId(const int& rhs)
{
   d_value = rhs;
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::isValid() const
{
   return d_value != s_invalid_id.d_value;
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const int& BlockId::getBlockValue() const
{
   return d_value;
}

/*
 ********************************************************************************
 * Equality operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator == (const BlockId& rhs) const
{
   return d_value == rhs.d_value;
}

/*
 ********************************************************************************
 * Equality operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator == (const int& rhs) const
{
   return d_value == rhs;
}

/*
 ********************************************************************************
 * Inequality operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator != (const BlockId& rhs) const
{
   return !(*this == rhs);
}

/*
 ********************************************************************************
 * Less-than operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator < (const BlockId& rhs) const
{
   return d_value < rhs.d_value;
}

/*
 ********************************************************************************
 * Greater-than operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator > (const BlockId& rhs) const
{
   return d_value > rhs.d_value;
}

/*
 ********************************************************************************
 * Less-than-or-equal-to operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator <= (const BlockId& rhs) const
{
   return d_value <= rhs.d_value;
}

/*
 ********************************************************************************
 * Greater-thanor-equal-to operator.
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BlockId::operator >= (const BlockId& rhs) const
{
   return d_value >= rhs.d_value;
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const BlockId& BlockId::zero()
{
   return s_zero_id;
}

/*
 ********************************************************************************
 ********************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const BlockId& BlockId::invalidId()
{
   return s_invalid_id;
}

}
}
