/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An restart manager singleton class
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Returns d_is_from_restart.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
bool RestartManager::isFromRestart()
{
   return d_is_from_restart;
}

/*
 *************************************************************************
 *
 * Returns the root of the restart database which was opened at the
 * start of the simulation.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
boost::shared_ptr<Database> RestartManager::getRootDatabase()
{
   return d_database_root;
}

/*
 *************************************************************************
 *
 * Sets the root database
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void RestartManager::setRootDatabase(
   const boost::shared_ptr<Database>& database)
{
   if (!database) {
      d_database_root.reset();
      d_is_from_restart = false;
   } else {
      d_database_root = database;
      d_is_from_restart = true;
   }
}

/*
 *************************************************************************
 *
 * Sets the database factory
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void RestartManager::setDatabaseFactory(
   const boost::shared_ptr<DatabaseFactory>& database_factory)
{
   d_database_factory = database_factory;
}

/*
 *************************************************************************
 *
 * writeRestartFile
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void RestartManager::writeRestartFile(
   const std::string& root_dir_name)
{
   writeRestartFile(root_dir_name, 0);
}

SAMRAI_INLINE_KEYWORD
void RestartManager::registerSingletonSubclassInstance(
   RestartManager* subclass_instance)
{
   if (!s_manager_instance) {
      s_manager_instance = subclass_instance;
   } else {
      TBOX_ERROR("RestartManager internal error...\n"
         << "Attemptng to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Remove all items from the restart item list.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void RestartManager::clearRestartItems()
{
   d_restart_items_list.clearItems();
}

/*
 *************************************************************************
 *
 * Write simulation state to root database
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void RestartManager::writeRestartToDatabase()
{
   if (d_database_root) {
      writeRestartFile(d_database_root);
   } else {
      TBOX_ERROR("writeRestartToDatabase has no database to write to"
         << std::endl);
   }
}

}
}
