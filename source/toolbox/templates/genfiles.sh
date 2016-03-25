##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: shell script to create SAMRAI template files in the repository
##

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# set up PERL program name and the make-templates.pl command
#

PERL=${PERL:-perl}
MT="$PERL ../../scripts/make-template.pl"

#
# create a new scratch temporary directory
#

rm -rf ./tmp
mkdir ./tmp 

#
# Create empty filename files for each type 
#
for t in default Complex Double Float Integer bool dcomplex char float double int; do
    touch ./tmp/$t.filenames
done


for t in bool dcomplex char float double int; do
   ${MT} default.filenames ./tmp tbox tbox::MathUtilities $t
done
#
# basic arrays for built-in types
#

for t in bool dcomplex char "char*" float double int string; do
   ${MT} default.filenames ./tmp tbox tbox::Array $t
done

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array bool
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array double
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array int

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array bool
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array double
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array int

${MT} default.filenames ./tmp tbox tbox::Array tbox::List int

#
# other array templates
#

${MT} default.filenames ./tmp tbox tbox::Array tbox::DatabaseBox
${MT} default.filenames ./tmp tbox tbox::Array tbox::Schedule::ScheduleMessageStream
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer tbox::Statistic
${MT} default.filenames ./tmp tbox tbox::Array tbox::Statistic::PatchStat
${MT} default.filenames ./tmp tbox tbox::Array tbox::Statistic::ProcStat

${MT} default.filenames ./tmp tbox tbox::Array tbox::List tbox::Pointer tbox::Transaction

${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer tbox::HDFDatabase
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer tbox::Timer
${MT} default.filenames ./tmp tbox tbox::Array tbox::AsyncCommGroup*
${MT} default.filenames ./tmp tbox tbox::Array tbox::AsyncCommStage::StagedGroup*
${MT} default.filenames ./tmp tbox tbox::Array tbox::RelaunchableJob*
${MT} default.filenames ./tmp tbox tbox::List tbox::RelaunchableJob*

${MT} default.filenames ./tmp tbox tbox::Array tbox::SAMRAI_MPI::request
${MT} default.filenames ./tmp tbox tbox::Array tbox::SAMRAI_MPI::status

#
# basic lists for built-in types
#

for t in double int string; do
   ${MT} default.filenames ./tmp tbox tbox::List $t
done

#
# other list templates
#

${MT} default.filenames ./tmp tbox tbox::List tbox::HDFDatabase::KeyData
${MT} default.filenames ./tmp tbox tbox::List tbox::InputDatabase::KeyData
${MT} default.filenames ./tmp tbox tbox::List tbox::Parser::ParseData
${MT} default.filenames ./tmp tbox tbox::List tbox::RestartManager::RestartItem
${MT} default.filenames ./tmp tbox tbox::List tbox::Statistic::PatchStatRecord

${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer tbox::Database
${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer tbox::Transaction
${MT} default.filenames ./tmp tbox tbox::List tbox::Timer*

#
# pointer templates
#

${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Arena
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Database
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::FixedArena
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::HDFDatabase
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::InputDatabase
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::MessageStream
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::NullDatabase
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Schedule
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Statistic
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Timer
${MT} default.filenames ./tmp tbox tbox::Pointer tbox::Transaction

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticNond ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticNond
sh ../../scripts/depend

rm -rf ./tmp
