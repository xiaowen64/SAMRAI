##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisVolInterface.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Isosurface Interface
##

class cvisVolInterface {
    inherit cvisInterface

    variable Cropping 0
    variable CroppingRegionPlanes

    variable Scaling 1.0
    variable RealBounds

    common Number 0
    
    variable Input ""

    variable ColorFileName ""
    variable ColorFileTypes {{ColorMaps {.vcmap}}}

    variable AlphaFileName ""
    variable AlphaFileTypes {{ColorMaps {.vamap}}}

    constructor {name} {
	cvisInterface::constructor $name
    } {
	foreach plane {x y z} {
	    set CroppingRegionPlanes($plane,min) 0
	    set CroppingRegionPlanes($plane,max) 1
	    set RealBounds($plane,min) 0.0
	    set RealBounds($plane,max) 1.0
	}
    }

    method Save {file} {

	puts $file "$this SetCroppingRegionPlanes \
		$CroppingRegionPlanes(x,min) \
		$CroppingRegionPlanes(x,max) \
		$CroppingRegionPlanes(y,min) \
		$CroppingRegionPlanes(y,max) \
		$CroppingRegionPlanes(z,min) \
		$CroppingRegionPlanes(z,max)"
	puts $file "$this SetCropping $Cropping"
	if {[string length $ColorFileName]} {
	    puts $file "$this SetColorFileName \"$ColorFileName\""
	}
	if {[string length $AlphaFileName]} {
	    puts $file "$this SetAlphaFileName \"$AlphaFileName\""
	}
    }

    variable ColorFileBrowser ".null"
    method SetColorFileName {colorfilename} {
	set ColorFileName $colorfilename
	foreach plane $Input {
	    $plane SetColorFileName $ColorFileName
	}

    }

    method LoadColorMap {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $ColorFileBrowser]]} {
	    $ColorFileBrowser SetOpenPrompt
	} {
	    set ColorFileBrowser [cvisFileBrowser #auto]
	    $ColorFileBrowser SetTitle "Select ColorMap File"
	    $ColorFileBrowser SetFileTypes $ColorFileTypes
	    $ColorFileBrowser SetOpenPrompt
	}

	set filename [$ColorFileBrowser Prompt]
	SetColorFileName $filename

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    variable AlphaFileBrowser ".null"
    method SetAlphaFileName {colorfilename} {
	set AlphaFileName $colorfilename
	foreach plane $Input {
	    $plane SetAlphaFileName $AlphaFileName
	}

    }

    method LoadAlphaMap {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $AlphaFileBrowser]]} {
	    $AlphaFileBrowser SetOpenPrompt
	} {
	    set AlphaFileBrowser [cvisFileBrowser #auto]
	    $AlphaFileBrowser SetTitle "Select AlphaMap File"
	    $AlphaFileBrowser SetFileTypes $AlphaFileTypes
	    $AlphaFileBrowser SetOpenPrompt
	}

	set filename [$AlphaFileBrowser Prompt]
	SetAlphaFileName $filename

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method AddInput {iso} {
	lappend Input $iso

	$iso SetInterface $this
    }
    
    method ExecuteSetCropping {} {
	set frame $NotebookPage.vol
	set cs [$frame childsite]
	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    
	    if {$Cropping} {
		$plane_frame.max configure -state normal
		$plane_frame.min configure -state normal
	    } {
		$plane_frame.max configure -state disabled
		$plane_frame.min configure -state disabled
	    }
	}
	
	foreach item $Input {
	    $item SetCropping $Cropping
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetCropping {cropping} {
	if {$Cropping != $cropping } {
	    set Cropping $cropping


	    ExecuteSetCropping
	}
    }

    method GetCropping {} {
	return $Cropping
    }

    method ExecuteSetCroppingRegionPlanes {} {

	foreach item $Input {

	    $item SetCroppingRegionPlanes \
		    [expr $CroppingRegionPlanes(x,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(x,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,max) * $Scaling]
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetCroppingRegionPlanes {min_x max_x min_y max_y min_z max_z} {

	set CroppingRegionPlanes(x,min) $min_x
	set CroppingRegionPlanes(x,max) $max_x

	set CroppingRegionPlanes(y,min) $min_y
	set CroppingRegionPlanes(y,max) $max_y

	set CroppingRegionPlanes(z,min) $min_z
	set CroppingRegionPlanes(z,max) $max_z

	ExecuteSetCroppingRegionPlanes
    }

    method Update {} {

	set item [lindex $Input 0]

	set Scaling [$item GetScaling]

	scan [$item GetBounds] "%f %f %f %f %f %f" \
		RealBounds(x,min) RealBounds(x,max) \
		RealBounds(y,min) RealBounds(y,max) \
		RealBounds(z,min) RealBounds(z,max)

	# convert from screen to real coords
	foreach plane {x y z} {
	    set RealBounds($plane,min) [expr $RealBounds($plane,min) / $Scaling]
	    set RealBounds($plane,max) [expr $RealBounds($plane,max) / $Scaling]
	}

	set frame $NotebookPage.vol
	set cs [$frame childsite]
	foreach plane {x y z} {
	    set plane_frame $cs.$plane

	    $plane_frame.mintitle configure -text "Min $plane ($RealBounds($plane,min))"
	    $plane_frame.maxtitle configure -text "Max $plane ($RealBounds($plane,max))"
	}
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# IsoValue controls
	#-----------------------------------------------------------------
	set frame $NotebookPage.vol
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Cropping Region"
	set cs [$frame childsite]

	checkbutton $cs.cropping \
		-text "Crop Data" \
		-variable [scope Cropping] \
		-command "$this ExecuteSetCropping"

	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    frame $plane_frame
	    
	    label $plane_frame.mintitle -text "Min $plane ($RealBounds($plane,min))"
	    entry $plane_frame.min \
		-textvariable [scope CroppingRegionPlanes($plane,min)] 
	    bind $plane_frame.min <Return> "$this ExecuteSetCroppingRegionPlanes"

	    label $plane_frame.maxtitle -text "Max $plane ($RealBounds($plane,max))"
	    entry $plane_frame.max \
		    -textvariable [scope CroppingRegionPlanes($plane,max)] 
	    bind $plane_frame.max <Return> "$this ExecuteSetCroppingRegionPlanes"
	    
	    if {$Cropping} {
		$plane_frame.max configure -state normal
		$plane_frame.min configure -state normal
	    } {
		$plane_frame.max configure -state disabled
		$plane_frame.min configure -state disabled
	    }
	}

	#------------------------------------------------------------------
	# Color Scaling
	#------------------------------------------------------------------
	set frame $NotebookPage.scaleframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Color Mapping"
	set cs [$frame childsite]

	set frame $cs.cmap
	frame $frame
	button $frame.colorfilename_button -text "Load ColorMap" \
		-command  "$this LoadColorMap" 

	entry $frame.colorfilename \
		-textvariable [scope ColorFileName]

	set frame $cs.amap
	frame $frame
	button $frame.alphafilename_button -text "Load AlphaMap" \
		-command  "$this LoadAlphaMap" 

	entry $frame.alphafilename \
		-textvariable [scope AlphaFileName]

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.vol
	set cs [$frame childsite]

	pack $cs.cropping -side left

	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    pack $plane_frame.mintitle -side left
	    pack $plane_frame.min -side left -padx 4 -pady 4
	    pack $plane_frame.maxtitle -side left
	    pack $plane_frame.max -side left -padx 4 -pady 4
	    pack $plane_frame -fill x
	}

 	pack $cs -fill x
  	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.scaleframe
	pack $frame -fill x
	set cs [$frame childsite]

	set frame $cs.cmap
	pack $frame.colorfilename_button -side left -anchor w -padx 4 -pady 8
	pack $frame.colorfilename -side left \
		-fill x -expand yes -pady 4
	pack $frame -fill x

	set frame $cs.amap
	pack $frame.alphafilename_button -side left -anchor w -padx 4 -pady 4
	pack $frame.alphafilename -side left \
		-fill x -expand yes -pady 4
	pack $frame -fill x

	pack $cs -fill x
    }   
}

