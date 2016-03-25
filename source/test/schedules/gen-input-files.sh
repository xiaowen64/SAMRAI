# This script generates input files that are used to 
# test the different options available in 
# xfer_RefineSchedule.  Here's what it does:
#
#   1. Appends "append-[n].input" to test-2d.input
#      to create the inputs test-2d-append-[n].input.
#      Here, test-2d.input is pulled directly from
#      the test/applications/Euler directory.
#   2. Runs the 2D Euler executable on the new
#      set of input files.
#
# The different "n" in the append files signifies the
# different test cases.  For example, append-1.input might
# specify "refine_schedule_generation_method = ORIG_NSQUARED",     
# append-2.input would specify a different option, etc.
#
# This script should be invoked from within 
# source/test/schedules/Makefile as:
#
# new-algorithms:
#	sh gen-input-files.sh 
#

for n in 1 2 3
do
  #
  # construct "original" input file - cat input file (test-2d.input) 
  # with append-n.input
  #
  echo "test-2d.input append-$n.input > test-2d-append-$n.input" 
  cat test-2d.input append-$n.input > test-2d-append-$n.input 

done
