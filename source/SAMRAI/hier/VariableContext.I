/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple integer id and namestring variable context 
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
int VariableContext::getCurrentMaximumInstanceNumber()
{
   return s_instance_counter;
}

SAMRAI_INLINE_KEYWORD
int VariableContext::getIndex() const
{
   return d_index;
}

SAMRAI_INLINE_KEYWORD
const std::string& VariableContext::getName() const
{
   return d_name;
}

SAMRAI_INLINE_KEYWORD
bool VariableContext::operator == (
   const VariableContext& other) const
{
   return d_index == other.d_index;
}

}
}
