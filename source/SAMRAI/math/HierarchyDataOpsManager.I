/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton manager for hierarchy data operation objects.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

SAMRAI_INLINE_KEYWORD
HierarchyDataOpsManager*
HierarchyDataOpsManager::getManager()
{
   if (!s_pdat_op_manager_instance) {
      s_pdat_op_manager_instance = new HierarchyDataOpsManager();
   }
   return s_pdat_op_manager_instance;
}

SAMRAI_INLINE_KEYWORD
void
HierarchyDataOpsManager::shutdownCallback()
{
   if (s_pdat_op_manager_instance) {
      delete s_pdat_op_manager_instance;
   }
   s_pdat_op_manager_instance = ((HierarchyDataOpsManager *)NULL);
}

SAMRAI_INLINE_KEYWORD
void
HierarchyDataOpsManager::registerSingletonSubclassInstance(
   HierarchyDataOpsManager* subclass_instance)
{
   if (!s_pdat_op_manager_instance) {
      s_pdat_op_manager_instance = subclass_instance;
   } else {
      TBOX_ERROR("HierarchyDataOpsManager internal error...\n"
         << "Attempting to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

}
}
