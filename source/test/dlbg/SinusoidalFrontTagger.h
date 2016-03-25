/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dlbg/SinusoidalFrontTagger.h $
  Copyright:	(c) 1997-2000 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	SinusoidalFrontTagger class declaration
*/

#ifndef included_SinusoidalFrontTagger
#define included_SinusoidalFrontTagger


#include <string>
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
#include "StandardTagAndInitStrategy.h"
#include "CellData.h"
#include "NodeData.h"


using namespace std;
using namespace SAMRAI;


/*!
  @brief Class to tag a sinusoidal "front" in given domain.
*/
template<int DIM>
class SinusoidalFrontTagger
  : public mesh::StandardTagAndInitStrategy<DIM>,
    public appu::VisDerivedDataStrategy<DIM>
{

public:

  /*!
    @brief Constructor.
  */
  SinusoidalFrontTagger(
    /*! Ojbect name */
    const string &object_name ,
    /*! Input database */
    SAMRAI::tbox::Database *database=NULL );

  ~SinusoidalFrontTagger();


  //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy<DIM> virtuals

public:

  /*!
    @brief Allocate and initialize data for a new level
    in the patch hierarchy.

    This is where you implement the code for initialize data on the
    grid.  Nevermind when it is called or where in the program that
    happens.  All the information you need to initialize the grid
    are in the arguments.

    @see SAMRAI::mesh::StandardTagAndInitStrategy<DIM>::initializeLevelData()
  */
  virtual void initializeLevelData (
    /*! Hierarchy to initialize */
    const tbox::Pointer< SAMRAI::hier::BasePatchHierarchy<DIM> > hierarchy ,
    /*! Level to initialize */
    const int level_number ,
    const double init_data_time ,
    const bool can_be_refined ,
    /*! Whether level is being introduced for the first time */
    const bool initial_time ,
    /*! Level to copy data from */
    const tbox::Pointer< SAMRAI::hier::BasePatchLevel<DIM> > old_level
      = tbox::Pointer< SAMRAI::hier::PatchLevel<DIM> >((0)) ,
    /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true );

  virtual void resetHierarchyConfiguration (
    /*! New hierarchy */
    tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<DIM> > new_hierarchy ,
    /*! Coarsest level */ int coarsest_level ,
    /*! Finest level */ int finest_level );

  virtual void applyGradientDetector(
    const tbox::Pointer< hier::BasePatchHierarchy<DIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation );

  //@}



   void initializePatchData (
      hier::Patch<DIM> &patch ,
      const double init_data_time ,
      const bool initial_time ,
      const bool allocate_data );

  bool packDerivedDataIntoDoubleBuffer(double *buffer,
				       const hier::Patch<DIM> &patch,
				       const hier::Box<DIM> &region,
				       const string &variable_name,
				       int depth_index) const;


public:

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
    appu::CartesianVizamraiDataWriter<DIM> &writer );
  /*!
    @brief Tell a VisIt plotter which data to write for this class.
  */
  int registerVariablesWithPlotter(
    appu::VisItDataWriter<DIM> &writer );


   /*
     Compute patch data allocated by this class, on a hierarchy.
   */
   void computeHierarchyData( hier::PatchHierarchy<DIM> &hierarchy,
                              double time );

  /*!
    @brief Compute distance and tag data for a level.
  */
  void computeLevelData(
    const hier::PatchHierarchy<DIM> &hierarchy,
    const int ln,
    const double time,
    const int dist_id,
    const int tag_id,
    const tbox::Pointer<hier::PatchLevel<DIM> > &old_level = tbox::Pointer<hier::PatchLevel<DIM> >() ) const;

  /*!
    @brief Compute distance and tag data for a patch.
  */
  void computePatchData(
    const hier::Patch<DIM> &patch,
    const double time,
    pdat::NodeData<DIM,double> *dist_data,
    pdat::CellData<DIM,int> *tag_data) const;

  /*!
    @brief Copy tag src to tag dst only where tag src is not 0.
  */
  void overlayTagData( const pdat::CellData<DIM,int> &src,
                       pdat::CellData<DIM,int> &dst ) const;



private:
  string d_name;
  tbox::Pointer<hier::PatchHierarchy<DIM> > d_hierarchy;

  /*!
    @brief Period of sinusoid.
  */
  double d_period;

   /*!
     @brief Initial displacement.
   */
  double d_init_disp[DIM];

   /*!
     @brief Front velocity.
   */
  double d_velocity[DIM];

  /*!
    @brief Amplitude of sinusoid.
  */
  double d_amplitude;

  /*!
    @brief Number of cells to tag around cells intersecting the front.
  */
  int d_adaption_buffer;

  tbox::Pointer<hier::VariableContext> d_context;

  /*!
    @brief Distance from the front in the x direction.
  */
  int d_dist_id;
  /*!
    @brief Value of tag based on distance from front.
   */
  int d_tag_id;

  /*!
    @brief Whether to allocate data on the mesh.
  */
  bool d_allocate_data;

  /*!
    @brief Front time.
  */
  double d_time;

};


#endif	// included_ssup_SinusoidalFrontTagger
