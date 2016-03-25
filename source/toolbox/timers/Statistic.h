//
// File:        Statistic.h
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    \f$       \f$
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Class to record statistics during program execution.
//

#ifndef included_tbox_Statistic
#define included_tbox_Statistic

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_AbstractStream
#include "tbox/AbstractStream.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

namespace SAMRAI {
 namespace tbox {

class Statistician;

/**
 * Class Statistic defines a simple object that can be used to record
 * information generated during the course of a simulation for post-
 * processing later.  Each statistic object is created by the singleton
 * Statistian object and is defined by a name string
 * identifier and is characterized by the sort of information it may record.
 * Depending on how the object is created, it may record processor
 * information (i.e., a separate value for each processor), or patch
 * information (i.e., a separate value for each patch on each processor).
 * An example of the former may be the total number of cells on each
 * processor.  An example of the second may be the number of cells on each
 * patch.  Each recorded data item may be any numerical value, but it 
 * will always be stored as a double for simplicity.  The string identifier 
 * for a processor stat is "PROC_STAT" and the string identifier for a 
 * patch stat is "PATCH_STAT".
 * 
 * An example use of a Statistic to record the number of gridcells on each
 * processor is as follows:
 *
 *    Pointer<Statistic> stat_num_gridcells =
 *        Statistician::getStatistician()->
 *        getStatistic("NumberGridcells", "PROC_STAT");
 *    ...
 *    stat_num_gridcells->recordProcStat(num_cells_on_proc);
 *    ...
 *
 * The type of the statistic restricts the way in which the statistic 
 * object may be used to record information.  For a "processor" stat
 * only the recordProcStat() functions can be used.  For a "patch"
 * stat only the recordPatchStat() functions can be used.
 * 
 * Typically, the information is recorded to generate a time sequence of
 * values.  But this need not be the case always.  An optional time
 * stamp may be provided for each value as it is recorded.  In any case,
 * the sequence order of the values is determined by the recording order.
 *
 * Also, the Statistician class is used to manage Statistic 
 * objects.  It provided a global point of access for creating and accessing
 * statistic objects and supports post-processing statistic information 
 * in parallel.
 *
 * In some cases, it may be desirable to record information for each
 * level in a calculation; e.g., the number of cells on each processor
 * on level zero, level 1, etc.  In this case, one can cimply create a
 * separate statistic object for each level. 
 *
 * @see tbox::Statistician
 */

class Statistic : public DescribedClass
{
   friend class Statistician;
public:
   /**
    * Virtual destructor destroys recorded object data.
    */
   virtual ~Statistic();

   /**
    * Return string name identifier for statistic object.
    */
   string getName() const;

   /**
    * Return string statistic type identifier for statistic object.
    */
   string getType() const;

   /**
    * Return integer instance identifier for statistic object.
    */
   int getInstanceId() const;

   /**
    * Return integer length of list of statistic sequence records.
    * This value is either the length of the processor statistic list
    * or the patch statistic list, whichever corresponds to the statistic 
    * type.
    */
   int getStatSequenceLength() const;

   /**
    * Reset the state of the statistic information.
    */
   void reset();

   /**
    * Record double processor statistic value. The optional sequence number 
    * argument identifies where in timestep sequence the value should be.
    * If the sequence number is not specified, an internal counter will 
    * determine the appropriate sequence number. When assertion checking 
    * is active, an unrecoverable exception will result if this function 
    * is called and "PATCH_STAT" was specified in the constructor.
    */
   void recordProcStat(double value,
                       int seq_num = -1);

   /**
    * Record double patch statistic value.  The patch number refers to
    * the global patch number on a level.  The sequence number 
    * argument identifies where in timestep sequence the value should be.
    * The sequence number MUST be explicitly specified because the number 
    * of patches on each processor will generally be different at
    * each sequence step.  When assertion checking is active, an 
    * unrecoverable exception will result if this function is called and 
    * "PROC_STAT" was specified in the constructor.
    */
   void recordPatchStat(int patch_num,
                        double value,
                        int seq_num);

   /**
    * Return true if size of stream required to pack all statistic
    * data can be determined for all processors without exchanging
    * any details of structure of statistic data.  Otherwise, return false.
    */
   bool canEstimateDataStreamSize();

   /**
    * Return integer number of bytes needed to stream the statistic data.
    * This is the amount needed by the stat transaction class.
    */
   int getDataStreamSize();

   /**
    * Pack contents of statistic data structure into message stream.
    */
   void packStream(AbstractStream& stream);

   /**
    * Unpack contents of statistic data structure from message stream.
    */
   void unpackStream(AbstractStream& stream);

   /**
    * Print statistic data to given output stream.  Floating point precision
    * can be specified (default is 12).
    */
   void printClassData(ostream& stream, int precision = 12) const;

   /**
    * Write statistic data members to database. When assertion checking 
    * is on, the database pointer must be non-null.
    */
   virtual void putToDatabase( Pointer<Database> db );

   /**
    * Read restarted times from restart database.  When assertion checking 
    * is on, the database pointer must be non-null.
    */
   virtual void getFromRestart(Pointer<Database> db);

   /*
    * These structures are used to store statistic data entries.
    * They need to be declared public for the Sun CC compiler.
    */
   struct ProcStat {
      double value;        // stat record value
   };

   struct PatchStatRecord {
      int    patch_id;      // global patch number
      double value;         // stat record value
   };

   struct PatchStat {
#ifdef LACKS_NAMESPACE_IN_DECLARE
      List<PatchStatRecord> patch_records; // stat record
#else
      List<Statistic::PatchStatRecord> patch_records; // stat record
#endif
   };

protected:
   /**
    * The constructor for the Statistic class sets the name string
    * and the statistic type for a statistic object.
    */
   Statistic(const string& name,
                  const string& stat_type,
                  int instance_id);

   /**
    * Return const reference to list of processor records.
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   const Array<ProcStat>& getProcStatSeqArray() const;
#else
   const Array<Statistic::ProcStat>& getProcStatSeqArray() const;
#endif

   /**
    * Return const reference to list of patch records.
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   const Array<PatchStat>& getPatchStatSeqArray() const;
#else
   const Array<Statistic::PatchStat>& getPatchStatSeqArray() const;
#endif

private:

   /*
    * Static double value used to indicate when a particular sequence entry
    * is skipped.
    */
   static double s_empty_seq_tag_entry;

   /*
    * Assess whether the processor or patch stat arrays need to be resized.
    */
   void checkArraySizes(int seq_num);

   
   // The following two members are not implemented
   Statistic(const Statistic&);
   void operator=(const Statistic&);

   /*
    * The enumerated type maps statistic types to integers.
    */
   enum STATISTIC_RECORD_TYPE{ PROC_STAT = 0, PATCH_STAT = 1 };

   /*
    * Name, instance id, and type identifier for this statistic object.
    */
   string d_object_name;
   int    d_instance_id;
   int    d_stat_type;         // see STATISTIC_RECORD_TYPE above.

   /*
    * Arrays of records.  Note that one of these will always be empty.
    * Integer sequence length refers to length of list corresponding
    * to stat type.
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   Array<ProcStat> d_proc_array;
   Array<PatchStat> d_patch_array;
#else
   Array<Statistic::ProcStat> d_proc_array;
   Array<Statistic::PatchStat> d_patch_array;
#endif

   /*
    * Sequence and patch counters (NOTE: patch counter use for patch stats
    * only) and high-water-mark array sizes for proc and patch stats.
    */
   int d_seq_counter;
   int d_total_patch_entries;
   int d_proc_stat_array_size;
   int d_patch_stat_array_size;
};

}  
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Statistic.I"
#endif
#endif
