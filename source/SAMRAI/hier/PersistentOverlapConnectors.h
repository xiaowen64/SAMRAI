/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Manager of Connectors incident from a common BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_PersistentOverlapConnectors
#define included_hier_PersistentOverlapConnectors

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>

namespace SAMRAI {
namespace hier {

class Connector;
class BoxLevel;

/*!
 * @brief A managager of overlap Connectors incident from a
 * BoxLevel, used to store and, if needed, generate overlap
 * Connectors in BoxLevel.
 *
 * PersistantOverlapConnectors provides a mechanism for objects using
 * the BoxLevel to look up overlap Connectors for the
 * BoxLevel.  Connectors created or returned by this class are
 * complete overlap Connectors that always contain the correct overlap
 * neighbor data.
 *
 * For improved scalability, Connectors can be constructed externally
 * and copied into the collection.  Connectors can also be
 * automatically computed using a non-scalable global search.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 * <b>bool check_accessed_connectors:</b> When true, check Connectors when
 * they are accessed.  The check is an non-scalable operation and is
 * meant for debugging.
 *
 * <b>string implicit_connector_creation_rule:</b> How to proceed when
 * findConnector() cannot find any suitable overlap Connector.  Values
 * can be "ERROR", "WARN" (default) or "SILENT".  If "SILENT",
 * silently get a globalized version of the head BoxLevel and look for
 * overlaps.  If "WARN", do the same thing but write a warning to the
 * log.  If "ERROR", exit with an error.
 *
 * @note Creating overlap Connectors by global search is not scalable.
 * Nevertheless, the default is "WARN", so that application
 * development need not worry about missing overlap Connectors.  To
 * selectively enable automatic Connector generation, set this to
 * "ERROR" and use findOrCreateConnector() instead of findConnector()
 * where you are unsure if the Connector has been created.
 *
 * @see findConnector()
 * @see findOrCreateConnector()
 * @see hier::Connector
 */
class PersistentOverlapConnectors
{

public:
   /*!
    * @brief Deletes all Connectors to and from this object
    */
   ~PersistentOverlapConnectors();

   /*!
    * @brief Create an overlap Connector, computing relationships by
    * globalizing data.
    *
    * The base will be the BoxLevel that owns this object.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see hier::Connector
    * @see hier::Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width);

   /*!
    * @brief Create an overlap Connector with its transpose, computing
    * relationships by globalizing data.
    *
    * The base will be the BoxLevel that owns this object.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see hier::Connector
    * @see hier::Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    * @param[in] transpose_connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& connector_width,
      const IntVector& transpose_connector_width);

   /*!
    * @brief Create an overlap Connector using externally
    * computed relationships.
    *
    * Create the Connector initialized with the arguments.
    * The base will be the BoxLevel that owns this object.
    *
    * @see hier::Connector
    * @see hier::Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    * @param[in] relationships
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width,
      const Connector& relationships);

   /*!
    * @brief Cache the supplied overlap Connector.
    *
    * The head will be the supplied head and the base will be the
    * BoxLevel that owns this object.
    *
    * @param[in] head
    * @param[in] connector
    *
    * @pre myBoxLevel().isInitialized()
    * @pre connector
    */
   void
   cacheConnector(
      const BoxLevel& head,
      boost::shared_ptr<Connector>& connector);

   /*!
    * @brief Find an overlap Connector with the given head and minimum
    * Connector width.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a single Connector is
    * arbitrary and should be re-examined.
    *
    * @par Assertions
    * If no Connector fits the criteria and @c
    * implicit_connector_creation_rule is "ERROR", an assertion is
    * thrown.  To automatically create the Connector instead, use
    * findOrCreateConnector() or set @c
    * implicit_connector_creation_rule to "WARN" or "SILENT".
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Find an overlap Connector with its transpose with the given head
    * and minimum Connector widths.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a single Connector is
    * arbitrary and should be re-examined.
    *
    * @par Assertions
    * If no Connector fits the criteria and @c
    * implicit_connector_creation_rule is "ERROR", an assertion is
    * thrown.  To automatically create the Connector instead, use
    * findOrCreateConnector() or set @c
    * implicit_connector_creation_rule to "WARN" or "SILENT".
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] transpose_min_connector_width Find the transpose overlap
    *      Connector satisfying this minimum Connector width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      const IntVector& transpose_min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Find or create an overlap Connectors with the
    * given head and minimum Connector width.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a
    * single Connector is arbitrary and should be re-examined.
    *
    * If no Connector fits the criteria, a new one is created using
    * global search for edges.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findOrCreateConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Find or create an overlap Connectors with its transpose with the
    * given head and minimum Connector widths.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a
    * single Connector is arbitrary and should be re-examined.
    *
    * If no Connector fits the criteria, a new one is created using
    * global search for edges.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    * @param[in] transpose_min_connector_width Find the transpose overlap
    *      Connector satisfying this minimum ghost cell width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findOrCreateConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      const IntVector& transpose_min_connector_width,
      bool exact_width_only = false);

   /*!
    * @brief Returns whether the object has overlap
    * Connectors with the given head and minimum Connector
    * width.
    *
    * TODO:  does the following comment mean that this must be called
    * before the call to findConnector?
    *
    * If this returns true, the Connector fitting the specification
    * exists and findConnector() will not throw an assertion.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    * @param[in] exact_width_only If true, reject Connectors that do not
    *      match the requested width exactly.
    *
    * @return True if a Connector is found, otherwise false.
    */
   bool
   hasConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only = false) const;

   /*!
    * @brief Delete stored Connectors.
    */
   void
   clear();

   const BoxLevel&
   myBoxLevel()
   {
      return d_my_box_level;
   }

private:
   // Unimplemented default constructor.
   PersistentOverlapConnectors();

   //@{ @name Methods meant only for BoxLevel to use.

   /*!
    * @brief Constructor, to be called from the BoxLevel
    * allocating the object.
    *
    * @param my_box_level The BoxLevel served by this
    * object.
    */
   explicit PersistentOverlapConnectors(
      const BoxLevel& my_box_level);

   /*
    * Read from the input database.
    */
   void
   getFromInput();

   /*
    * Private method used by different public versions of findConnector.
    */
   boost::shared_ptr<Connector>
   privateFindConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only);

   /*
    * Private method used by different public versions of
    * findOrCreateConnector.
    */
   boost::shared_ptr<Connector>
   privateFindOrCreateConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      bool exact_width_only);

   /*
    * Private method used by cacheConnector.
    */
   void
   privateCacheConnector(
      const BoxLevel& head,
      boost::shared_ptr<Connector>& connector);

   //@}

   //@{
   /*!
    * @brief Only BoxLevels are allowed to construct a
    * PersistentOverlapConnectors.
    */
   friend class BoxLevel;
   //@}

   typedef tbox::Array<boost::shared_ptr<Connector> > ConVect;

   /*!
    * @brief Persistent overlap Connectors incident from me.
    */
   ConVect d_cons_from_me;

   /*!
    * @brief Persistent overlap Connectors incident to me.
    */
   ConVect d_cons_to_me;

   /*!
    * @brief Reference to the BoxLevel served by this object.
    */
   const BoxLevel& d_my_box_level;

   /*!
    * @brief Whether to check overlap Connectors when they are created.
    */
   static char s_check_created_connectors;

   /*!
    * @brief Whether to check overlap Connectors when they are accessed.
    */
   static char s_check_accessed_connectors;

   /*!
    * @brief Whether to force Connector finding functions to create
    * connectors that are missing.
    *
    * See input parameter implicit_connector_creation_rule.
    */
   static char s_implicit_connector_creation_rule;

};

}
}

#endif // included_hier_PersistentOverlapConnectors
