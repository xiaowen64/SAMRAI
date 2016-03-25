/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dlbg/DLBGTest.h $
  Copyright:	(c) 1997-2000 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	DLBGTest class declaration
*/

#ifndef included_DLBGTest
#define included_DLBGTest


#include <string>
using namespace std;

#include "tbox/Pointer.h"
#include "tbox/Database.h"


/*
  SAMRAI classes
*/
#include "CartesianVizamraiDataWriter.h"
#include "VisItDataWriter.h"
#include "VisDerivedDataStrategy.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "CartesianRobinBcHelper.h"
#include "RobinBcCoefStrategy.h"
#include "SinusoidalFrontTagger.h"


using namespace SAMRAI;


/*!
  @brief Class to test new PIND algorithm.
*/
template<int DIM>
class DLBGTest
  : public appu::VisDerivedDataStrategy<DIM>
{

public:

  /*!
    @brief Constructor.
  */
  DLBGTest(
    /*! Ojbect name */
    const string &object_name ,
    /*! Input database */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database );

  ~DLBGTest();



  mesh::StandardTagAndInitStrategy<DIM> *getStandardTagAndInitObject();




  //@{ @name SAMRAI::appu::VisDerivedDataStrategy<DIM> virtuals

  virtual bool packDerivedDataIntoDoubleBuffer(
    double* buffer ,
    const hier::Patch<DIM> &patch ,
    const hier::Box<DIM> &region ,
    const string &variable_name ,
    int depth_id) const;

  //@}



public:

  /*
    Deallocate patch data allocated by this class.
  */
   void computeHierarchyData( hier::PatchHierarchy<DIM> &hierarchy,
                              double time );

  /*!
    @brief Deallocate internally managed patch data on level.
  */
  void deallocatePatchData( hier::PatchLevel<DIM> &level );

  /*!
    @brief Deallocate internally managed patch data on hierarchy.
  */
  void deallocatePatchData( hier::PatchHierarchy<DIM> &hierarchy );

  /*!
    @brief Tell a Vizamrai plotter which data to write for this class.
  */
  int registerVariablesWithPlotter(
    tbox::Pointer<appu::CartesianVizamraiDataWriter<DIM> > writer
    );
  /*!
    @brief Tell a VisIt plotter which data to write for this class.
  */
  int registerVariablesWithPlotter(
    tbox::Pointer<appu::VisItDataWriter<DIM> > writer
    );



private:

  string d_name;
  tbox::Pointer<hier::PatchHierarchy<DIM> > d_hierarchy;

  SinusoidalFrontTagger<DIM> d_tagger;

  /*!
    @brief Front time.
  */
  double d_time;

};


#endif	// included_DLBGTest
