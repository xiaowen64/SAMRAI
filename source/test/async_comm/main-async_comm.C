/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/async_comm/main-async_comm.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Test program for asynchromous communication classes
 */

#include "SAMRAI_config.h"

#include <iomanip>

#include "tbox/AsyncCommGroup.h"
#include "tbox/AsyncCommStage.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/TimerManager.h"


#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
using namespace tbox;
#endif

using namespace std;


/*
************************************************************************
*                                                                      *
* An implementation of relaunchable job that prints out its index      *
* when it completes.                                                   *
*                                                                      *
************************************************************************
*/
struct JobImplementation : public RelaunchableJob {
public:
   void continueJob() { TBOX_ERROR("Unimplenented method."); }
   JobState getJobState()
      { TBOX_ERROR("Unimplenented method."); return JOB_IS_COMPLETED; }
   AsyncCommGroup *getCommunicationGroup()
      { TBOX_ERROR("Unimplenented method."); return NULL; }
   int local_index; // Index in array of active groups.
};

/*
************************************************************************
*                                                                      *
* This program tests the asynchronous communication classes:           *
*   AsyncCommGroup                                                     *
*   AsyncCommStage                                                     *
*                                                                      *
* 1. Group the processors.  See code for heuristic rule for            *
*    defining groups.                                                  *
*                                                                      *
* 2. Perform asynchronous communication within each group.             *
*                                                                      *
* 3. Check results.                                                    *
*                                                                      *
*************************************************************************
*/

int main( int argc, char *argv[] )
{

   int pass_count = 0;
   int fail_count = 0;

   /*
    * Initialize MPI, SAMRAI.
    */
   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::startup();

   /* 
    * Create block so that smart pointer references are removed before 
    * SAMRAIManger is shutdown.  If this is not done memory leaks are
    * introduced.
    */
   {

      /*
       * Process command line arguments.  For each run, the input 
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
       */
      string input_filename;

      if (argc != 2)  {
	 TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n" 
		    << "  options:\n"
		    << "  none at this time" << endl);
      } else {
	 input_filename = argv[1];
      }


      /*
       * Make sure all processes are alive and well before running test.
       */
      const int rank = SAMRAI_MPI::getRank();
      SAMRAI_MPI::barrier();
      cout << "Process " << setw(5) << rank << " is ready." << endl;


      /*
       * Created a separate communicator for testing,
       * to avoid possible interference with other communications
       * by SAMRAI library.
       */
      tbox::SAMRAI_MPI::comm isolated_communicator;
#ifdef HAVE_MPI
      MPI_Comm_dup( MPI_COMM_WORLD, &isolated_communicator );
      cout << "Process " << setw(5) << rank
	   << " duplicated MPI_COMM_WORLD." << endl;
#endif


      /*
       * Create input database and parse all data in input file.
       */

      Pointer<Database> input_db = new InputDatabase("input_db");
      InputManager::getManager()->parseInputFile(input_filename, input_db);


      /*
       * Set up the timer manager.
       */
      if ( input_db->isDatabase("TimerManager") ) {
	 TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      Pointer<Database> main_db = input_db->getDatabase("Main");
      string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);


      /*
       * Start logging.
       */
      const string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
      if (log_all_nodes) {
	 PIO::logAllNodes(log_file_name);
      } else {
	 PIO::logOnlyNodeZero(log_file_name);
      }

      plog << "********************* Note! *********************\n"
	   << "* The asychronous communication classes are meant for\n"
	   << "* large processor counts.\n"
	   << "*\n"
	   << "* For this test to be significant, you should run it on\n"
	   << "* lots of processors.  I recommend running on a\n"
	   << "* 'massively parallel' processor count (however you\n"
	   << "* would like to define massively parallel).\n"
	   << "*\n"
	   << "* This program should not be used for performance\n"
	   << "* testing because performance may be context-sensitive\n"
	   << "* and these groups are rather contrived.";

      plog << "\n\n\n";



      const int sync_bcast_cycles =
	 main_db->getIntegerWithDefault("sync_bcast_cycles", 1);
      const int sync_sumreduce_cycles =
	 main_db->getIntegerWithDefault("sync_sumreduce_cycles", 1);
      const int asyncany_bcast_cycles =
	 main_db->getIntegerWithDefault("asyncany_bcast_cycles", 1);
      const int asyncany_sumreduce_cycles =
	 main_db->getIntegerWithDefault("asyncany_sumreduce_cycles", 1);
      const int asyncsome_bcast_cycles =
	 main_db->getIntegerWithDefault("asyncsome_bcast_cycles", 1);
      const int asyncsome_sumreduce_cycles =
	 main_db->getIntegerWithDefault("asyncsome_sumreduce_cycles", 1);

      int sync_bcast_count = 0;
      int sync_sumreduce_count = 0;
      int asyncany_bcast_count = 0;
      int asyncany_sumreduce_count = 0;
      int asyncsome_bcast_count = 0;
      int asyncsome_sumreduce_count = 0;



      const int def_num_groups = (SAMRAI_MPI::getNodes()+1)/2;
      plog << "Default num groups: " << def_num_groups << endl;
      const int num_groups =
	 main_db->getIntegerWithDefault("num_groups", def_num_groups);
      plog << "Num groups: " << num_groups << endl;
      const int num_children =
	 main_db->getIntegerWithDefault("num_children", 2);


      if ( SAMRAI_MPI::getRank() == 0 ) {
	 cout << "Running num_groups = " << num_groups << endl;
	 cout << "Running num_children = " << num_children << endl;
      }



      int gi; // Group index.
      int ai; // Active group index.



      /*
       * Output data for stage methods.
       */


      int count = 0;
      while ( ( sync_bcast_count          < sync_bcast_cycles          ) ||
	      ( sync_sumreduce_count      < sync_sumreduce_cycles      ) ||
	      ( asyncany_bcast_count      < asyncany_bcast_cycles      ) ||
	      ( asyncany_sumreduce_count  < asyncany_sumreduce_cycles  ) ||
	      ( asyncsome_bcast_count     < asyncsome_bcast_cycles     ) ||
	      ( asyncsome_sumreduce_count < asyncsome_sumreduce_cycles ) ) {


	 if ( SAMRAI_MPI::getRank() == 0 ) {
	    cout << " Starting cycle number " << count << endl;
	 }


	 plog << "\n\n\n***************** Beginning Cycle Number "
	      << count << " *******************\n\n";


	 Array<Array<int> > group_ids(num_groups);
	 Array<int> owners(num_groups);
	 Array<int> active_flags(num_groups);

	 Array<int> active_groups(num_groups);
	 int num_active_groups = 0;


	 /*
	  * Define groups.
	  * Group n includes all process whose rank divisible by n+1 -- with
	  * a slight variation.  With each testing cycle, the "rank" is
	  * increased by one.
	  * Set owner of each group to roughly the processor in the middle
	  * of the group.
	  */
	 for ( int n=0; n<num_groups; ++n ) {

	    int gsize = (SAMRAI_MPI::getNodes()+n)/(n+1);
	    group_ids[n].resizeArray( gsize );
	    active_flags[n] = false;
	    for ( int i=0; i<gsize; ++i ) {
	       group_ids[n][i] = i*(n+1);
	       group_ids[n][i] = ( group_ids[n][i] + count ) % SAMRAI_MPI::getNodes();
	       if ( group_ids[n][i] == rank ) {
		  active_groups[num_active_groups++] = n;
		  active_flags[n] = true;
	       }
	    }

	    owners[n] = group_ids[n][ gsize/2 ];

	 }
	 active_groups.resizeArray(num_active_groups);


	 /*
	  * Write out group data.
	  */
	 plog << "Group definitions (" << num_groups << " groups):\n";
	 plog << "(Groups with '*' contains the local process, "
	      << rank << ".)\n";
	 plog << "(Groups with '**' is owned by the local process\n\n";
	 plog << " ID  size owner members...\n";
	 for ( int n=0; n<num_groups; ++n ) {
	    plog << setw(3) << n
		 << setw(5) << group_ids[n].size()
		 << setw(4) << owners[n]
		 << (active_flags[n]?'*':' ')
		 << (owners[n]==rank?'*':' ') << ':';
	    for ( int i=0; i<group_ids[n].size(); ++i ) {
	       plog << "  " << group_ids[n][i];
	    }
	    plog << '\n';
	 }
	 plog << '\n';
      
	 plog << "Active groups (" << num_active_groups << " groups):";
	 for ( ai=0; ai<num_active_groups; ++ai ) {
	    plog << "  " << active_groups[ai];
	 }
	 plog << "\n\n";



	 /*
	  * Initialize data for sum-reduce tests.
	  * Compute the correct sum for comparison.
	  */

	 Array<int> sum(num_active_groups);
	 Array<int> correct_sum(num_active_groups);

	 for ( ai=0; ai<num_active_groups; ++ai ) {
	    sum[ai] = 1 + rank;
	    correct_sum[ai] = 0;
	    Array<int> &g = group_ids[active_groups[ai]];
	    for ( int j=0; j<g.size(); ++j ) {
	       correct_sum[ai] += 1 + g[j];
	    }
	 }


	 /*
	  * Initialize data for broadcast test.
	  * Broadcast data is 1001 + the group index.
	  */
	 Array<int> bcdata(num_active_groups);
	 Array<int> correct_bcdata(num_active_groups);
	 for ( ai=0; ai<num_active_groups; ++ai ) {
	    gi = active_groups[ai];
	    bcdata[ai] = rank == owners[gi] ? 1001 + gi : -1;
	    correct_bcdata[ai] = 1001 + gi;
	 }



	 /*
	  * Create the communication stage and groups.
	  * Each group uses its group index as the MPI tag.
	  */
	 JobImplementation *jobs = new JobImplementation[num_active_groups];
	 AsyncCommStage comm_stage;
	 Array<AsyncCommGroup*> comm_groups(num_active_groups);
	 for ( ai=0; ai<num_active_groups; ++ai ) {
	    gi = active_groups[ai];
	    jobs[ai].local_index = ai;
	    plog << "Allocating group " << gi << "\n";
	    comm_groups[ai] = comm_stage.allocateCommGroup(num_children,
							   &jobs[ai]);
	    comm_groups[ai]->setGroupAndRootRank(group_ids[gi], owners[gi]);
	    comm_groups[ai]->setMPICommunicator( isolated_communicator );
	    comm_groups[ai]->setMPITag(1000000*count+gi);
	    comm_groups[ai]->setUseBlockingSendToParent(false);
	    comm_groups[ai]->setUseBlockingSendToChildren(false);
	 }

	 if ( sync_bcast_count < sync_bcast_cycles ) {
	    /*
	     * For the synchronous (groupwise) broadcast test,
	     * each group broadcasts the its group id.
	     */
	    plog << "\n\n\n*********** Synchronous Broadcast "
		 << sync_bcast_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai )
	       if ( rank != owners[active_groups[ai]] ) bcdata[ai] = -1;
	    plog << "Job Group Result Correct  Note\n";
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       AsyncCommGroup &comm_group = *comm_groups[ai];
	       comm_group.beginBcast( &bcdata[ai], 1 );
	       comm_group.waitOperation();
	       TBOX_ASSERT( comm_group.isDone() );
	       gi = active_groups[ai];
	       plog << setw(3) << ai
		    << setw(5) << gi
		    << setw(8) << bcdata[ai]
		    << setw(8) << correct_bcdata[ai]
		  ;
	       plog << "  Bcast difference = "
		    << bcdata[ai]-correct_bcdata[ai];
	       if ( bcdata[ai] != correct_bcdata[ai] ) {
		  plog << "  Error!";
		  cout << "Error in bcast result for group "
		       << gi << endl;
		  ++fail_count;
	       }
	       else ++pass_count;
	       plog << endl;
	    }
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT( comm_stage.isDone() );
	    ++sync_bcast_count;
	 }




	 if ( sync_sumreduce_count < sync_sumreduce_cycles ) {
	    /*
	     * For the sum advanceSome reduce test,
	     * each group sums up the ranks of its members, plus 1.
	     */
	    plog << "\n\n\n*********** Synchronous Sum Reduce "
		 << sync_sumreduce_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai ) sum[ai] = 1 + rank;
	    plog << "Job Group Result Correct  Note\n";
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       AsyncCommGroup &comm_group = *comm_groups[ai];
	       comm_group.beginSumReduce( &sum[ai], 1 );
	       comm_group.waitOperation();
	       TBOX_ASSERT( comm_group.isDone() );
	       gi = active_groups[ai];
	       plog << setw(3) << ai
		    << setw(5) << gi
		    << setw(8) << sum[ai]
		    << setw(8) << correct_sum[ai]
		  ;
	       if ( rank == owners[gi] ) {
		  plog << "  Sum reduce difference = "
		       << sum[ai]-correct_sum[ai];
		  if ( sum[ai] != correct_sum[ai] ) {
		     plog << "  Error!";
		     cout << "Error in sum reduce result for group "
			  << gi << endl;
		     ++fail_count;
		  }
		  else ++pass_count;
	       }
	       else {
		  plog << "  Not owner (not checking)";
	       }
	       plog << endl;
	    }
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT ( comm_stage.isDone() );
	    ++sync_sumreduce_count;
	 }




	 if ( asyncany_bcast_count < asyncany_bcast_cycles ) {
	    /*
	     * For the advanceSome broadcast test,
	     * each group broadcasts the its group id.
	     */
	    plog << "\n\n\n*********** advanceAny Broadcast "
		 << asyncany_bcast_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai )
	       if ( rank != owners[active_groups[ai]] ) bcdata[ai] = -1;
	    plog << "Job Group Result Correct  Note\n";
	    ai = 0;
	    int counter = 0;
	    while ( counter < num_active_groups || ! comm_stage.isDone() ) {
	       JobImplementation *job = NULL;
	       if ( counter < num_active_groups ) {
		  if ( comm_groups[counter]->beginBcast( &bcdata[counter], 1 ) ) {
		     TBOX_ASSERT( comm_groups[counter]->isDone() );
		     job = &jobs[counter];
		  }
		  ++counter;
	       }
	       else {
		  RelaunchableJob *base_job = comm_stage.advanceAny();
		  job = dynamic_cast<JobImplementation*>(base_job);
		  TBOX_ASSERT( job != NULL );
	       }
	       if ( job != NULL ) {
		  ai = job->local_index;
		  gi = active_groups[ai];
		  plog << setw(3) << ai
		       << setw(5) << gi
		       << setw(8) << bcdata[ai]
		       << setw(8) << correct_bcdata[ai]
		     ;
		  plog << "  Bcast difference = "
		       << bcdata[ai]-correct_bcdata[ai];
		  if ( bcdata[ai] != correct_bcdata[ai] ) {
		     plog << "  Error!";
		     cout << "Error in bcast result for group "
			  << gi << endl;
		     ++fail_count;
		  }
		  else ++pass_count;
		  plog << endl;
		  TBOX_ASSERT( comm_groups[ai]->isDone() );
	       }
	       ++ai;
	    }
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT( comm_stage.isDone() );
	    ++asyncany_bcast_count;
	 }




	 if ( asyncany_sumreduce_count < asyncany_sumreduce_cycles ) {
	    /*
	     * For the advanceSome broadcast test,
	     * each group broadcasts the its group id.
	     */
	    plog << "\n\n\n*********** advanceAny Sum Reduce "
		 << asyncany_sumreduce_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai ) sum[ai] = 1 + rank;
	    plog << "Job Group Result Correct  Note\n";
	    ai = 0;
	    int counter = 0;
	    while ( counter < num_active_groups || ! comm_stage.isDone() ) {
	       JobImplementation *job = NULL;
	       if ( counter < num_active_groups ) {
		  if ( comm_groups[counter]->beginSumReduce( &sum[counter], 1 ) ) {
		     TBOX_ASSERT( comm_groups[counter]->isDone() );
		     job = &jobs[counter];
		  }
		  ++counter;
	       }
	       else {
		  RelaunchableJob *base_job = comm_stage.advanceAny();
		  job = dynamic_cast<JobImplementation*>(base_job);
		  TBOX_ASSERT( job != NULL );
	       }
	       if ( job != NULL ) {
		  ai = job->local_index;
		  gi = active_groups[ai];
		  plog << setw(3) << ai
		       << setw(5) << gi
		       << setw(8) << sum[ai]
		       << setw(8) << correct_sum[ai]
		     ;
		  if ( rank == owners[gi] ) {
		     plog << "  Sum reduce difference = "
			  << sum[ai]-correct_sum[ai];
		     if ( sum[ai] != correct_sum[ai] ) {
			plog << "  Error!";
			cout << "Error in sum reduce result for group "
			     << gi << endl;
			++fail_count;
		     }
		     else ++pass_count;
		  }
		  else {
		     plog << "  Not owner (not checking)";
		  }
		  plog << endl;
		  TBOX_ASSERT( comm_groups[ai]->isDone() );
	       }
	       ++ai;
	    }
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT( comm_stage.isDone() );
	    ++asyncany_sumreduce_count;
	 }




	 if ( asyncsome_bcast_count < asyncsome_bcast_cycles ) {
	    /*
	     * For the advanceSome broadcast test,
	     * each group broadcasts the its group id.
	     */
	    plog << "\n\n\n*********** advanceSome Broadcast "
		 << asyncsome_bcast_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai )
	       if ( rank != owners[active_groups[ai]] ) bcdata[ai] = -1;
	    Array<RelaunchableJob*> completed(num_active_groups);
	    int num_complete = 0;
	    num_complete = 0;
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       AsyncCommGroup &comm_group = *comm_groups[ai];
	       bool fin = comm_group.beginBcast( &bcdata[ai], 1 );
	       if ( fin ) {
		  completed[num_complete++] = &jobs[ai];
	       }
	    }
	    plog << "Job Group Result Correct  Note\n";
	    do {
	       for ( int n=0; n<num_complete; ++n ) {
		  JobImplementation *job =
		     dynamic_cast<JobImplementation*>(completed[n]);
		  TBOX_ASSERT( job != NULL );
		  ai = job->local_index;
		  gi = active_groups[ai];
		  plog << setw(3) << ai
		       << setw(5) << gi
		       << setw(8) << bcdata[ai]
		       << setw(8) << correct_bcdata[ai]
		     ;
		  plog << "  Bcast difference = "
		       << bcdata[ai]-correct_bcdata[ai];
		  if ( bcdata[ai] != correct_bcdata[ai] ) {
		     plog << "  Error!";
		     cout << "Error in bcast result for group "
			  << gi << endl;
		     ++fail_count;
		  }
		  else ++pass_count;
		  plog << endl;
		  TBOX_ASSERT( comm_groups[ai]->isDone() );
	       }
	       num_complete = comm_stage.advanceSome( completed );
	    } while ( num_complete != 0 );
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT( comm_stage.isDone() );
	    ++asyncsome_bcast_count;
	 }




	 if ( asyncsome_sumreduce_count < asyncsome_sumreduce_cycles ) {
	    /*
	     * For the sum advanceSome reduce test,
	     * each group sums up the ranks of its members, plus 1.
	     */
	    plog << "\n\n\n*********** advanceSome Sum Reduce "
		 << asyncsome_sumreduce_count << " ************\n";
	    for ( ai=0; ai<num_active_groups; ++ai ) sum[ai] = 1 + rank;
	    Array<RelaunchableJob*> completed(num_active_groups);
	    int num_complete = 0;
	    num_complete = 0;
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       AsyncCommGroup &comm_group = *comm_groups[ai];
	       bool fin = comm_group.beginSumReduce( &sum[ai], 1 );
	       if ( fin ) {
		  completed[num_complete++] = &jobs[ai];
	       }
	    }
	    plog << "Job Group Result Correct  Note\n";
	    do {
	       for ( int n=0; n<num_complete; ++n ) {
		  JobImplementation *job =
		     dynamic_cast<JobImplementation*>(completed[n]);
		  TBOX_ASSERT( job != NULL );
		  ai = job->local_index;
		  gi = active_groups[ai];
		  plog << setw(3) << ai
		       << setw(5) << gi
		       << setw(8) << sum[ai]
		       << setw(8) << correct_sum[ai]
		     ;
		  if ( rank == owners[gi] ) {
		     plog << "  Sum reduce difference = "
			  << sum[ai]-correct_sum[ai];
		     if ( sum[ai] != correct_sum[ai] ) {
			plog << "  Error!";
			cout << "Error in sum reduce result for group "
			     << gi << endl;
			++fail_count;
		     }
		     else ++pass_count;
		  }
		  else {
		     plog << "  Not owner (not checking)";
		  }
		  plog << endl;
		  TBOX_ASSERT( comm_groups[ai]->isDone() );
	       }
	       num_complete = comm_stage.advanceSome( completed );
	    } while ( num_complete != 0 );
	    for ( ai=0; ai<num_active_groups; ++ai ) {
	       TBOX_ASSERT( comm_groups[ai]->isDone() );
	    }
	    TBOX_ASSERT ( comm_stage.isDone() );
	    ++asyncsome_sumreduce_count;
	 }

	 delete[] jobs;

	 ++count;
      }



      plog << '\n';
      plog << "pass_count = " << pass_count << endl;
      plog << "fail_count = " << fail_count << endl;
      plog << "\n************** Test completed **************\n" << endl;
      input_db->printClassData(tbox::plog);


      /*
       * Clean up and exit.
       */

      TimerManager::getManager()->print(plog);

#if defined(HAVE_MPI)
      MPI_Comm_free( &isolated_communicator );
#endif

      if ( fail_count == 0 ) {
	 tbox::pout << "\nPASSED:  async_comm" << endl;
      }

      cout << "Process " << setw(5) << rank << " exiting." << endl;

   }

   tbox::SAMRAIManager::shutdown();
   SAMRAI_MPI::finalize();


   return(fail_count); 
}
