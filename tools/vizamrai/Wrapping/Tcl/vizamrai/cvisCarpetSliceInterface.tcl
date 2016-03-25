##
## File:        cvisCarpetSliceInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Carpet Slice Interface
##

class cvisCarpetSliceInterface {
    inherit cvisInterface

    common Number 0
    
    variable Slice 0
    variable SlicePlane "z"
    variable Input ""

    variable IndexScaling
    variable IndexBounds

    variable RealBounds 

    variable Scaling 1.0

    variable ScaleFactor 0.5
    variable Offset 0.0

    variable MaxScaleFactor 2.0
    variable MaxBounds

    variable PatchBoundaryIsVisible 1
    variable PatchBoundaryWidth 1
    variable BoundaryEdgesInterior 1

    variable CellToNode 0

    variable NumPatchColors 8
    variable PatchColors
    variable PatchVisible

    variable DataIsVisible 1

    method Save {file} {
	set tracervar [::itcl::local cvisTrace #auto]
	puts $file "$this SetSlicePlane $SlicePlane"
	puts $file "$this SetSlice $Slice"
	puts $file "$this SetScaleFactor $ScaleFactor"
	puts $file "$this SetOffset $Offset"
	puts $file "$this SetCellToNode $CellToNode"
	puts $file "$this SetPatchBoundaryWidth $PatchBoundaryWidth"
	puts $file "$this SetPatchBoundaryVisibility $PatchBoundaryIsVisible"
	puts $file "$this SetBoundaryEdgesInterior $BoundaryEdgesInterior"
	set level 0
	while {$level < $NumPatchColors} {
	    puts $file "$this SetPatchColor $level $PatchColors($level)"
	    puts $file "$this SetPatchVisibility $level $PatchVisible($level)"
	    incr level
	}
    }

    method SetPatchColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	set PatchColors($level) "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]

	set frame $NotebookPage.boundaries
	set cs_main [$frame childsite]
	set frame $cs_main.patchboundary
	set cs [$frame childsite]
	set frame $cs.level_$level
	$frame.level_color configure -bg $color
	
	foreach plane $Input {
	    eval $plane SetBorderColor $level $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method ExecuteSetPatchColor {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $PatchColors($level)
	
	set value [color Prompt]
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetPatchColor $level $r $g $b
	    
	}
	
	delete object color 
    }

    constructor {name} {
	cvisInterface::constructor $name
    } {
	foreach plane {x y z} {
	    set RealBounds(min,$plane) 0.0
	    set RealBounds(max,$plane) 1.1

	    set IndexBounds(min,$plane) 0
	    set IndexBounds(max,$plane) 1
	    set IndexScaling($plane) 1.0
	}

	set NumPatchColors 8
	# green
	set PatchColors(0) "0 1 0"
	# yellow
	set PatchColors(1) "1 1 0"
	# light blue
	set PatchColors(2) "0 1 1"
	# purple
	set PatchColors(3) "1 0 1"
	# pink 
	set PatchColors(4) "1 0.75 0.80"
	# orange
	set PatchColors(5) "1 0.65 0"
	# red
	set PatchColors(6) "1 0 0"
	# blue
	set PatchColors(7) "0 0 1"
	
	set PatchVisible(0) "1"
	set PatchVisible(1) "1"
	set PatchVisible(2) "1"
	set PatchVisible(3) "1"
	set PatchVisible(4) "1"
	set PatchVisible(5) "1"
	set PatchVisible(6) "1"
	set PatchVisible(7) "1"
    }

    method SetBoundaryEdgesInterior {boundaryedgesinterior} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BoundaryEdgesInterior != $boundaryedgesinterior} {
	    set BoundaryEdgesInterior $boundaryedgesinterior
	    ExecuteSetBoundaryEdgesInterior
	}
    }
    
    method GetBoundaryEdgesInterior {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $BoundaryEdgesInterior
    }

    method ExecuteSetBoundaryEdgesInterior { } {
	set tracervar [::itcl::local cvisTrace #auto]

  	foreach plane $Input {
	    $plane SetBoundaryEdgesInterior $BoundaryEdgesInterior
  	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetPatchVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$PatchVisible($level) != $visible} {
	    set PatchVisible($level) $visible
	    ExecuteSetPatchVisibility $level
	}
    }

    method GetPatchVisibility {level} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $PatchVisible($level)
    }

    method ExecuteSetPatchVisibility {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetBorderVisibility $level $PatchVisible($level)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetPatchBoundaryVisibility {patchboundaryvisible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$PatchBoundaryIsVisible != $patchboundaryvisible} {
	    set PatchBoundaryIsVisible $patchboundaryvisible
	    ExecuteSetPatchBoundaryVisibility
	}
    }
    
    method GetPatchBoundaryVisibility {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $PatchBoundaryIsVisible
    }

    method ExecuteSetPatchBoundaryVisibility { } {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetBoundaryEdgesVisibility $PatchBoundaryIsVisible
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }


    method SetPatchBoundaryWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$PatchBoundaryWidth != $width} {
	    set PatchBoundaryWidth $width
	    ExecuteSetPatchBounaryWidth
	}
    }
    
    method GetPatchBoundaryWidth {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $PatchBoundaryWidth
    }

    method ExecuteSetPatchBoundaryWidth { } {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetBoundaryEdgeWidth $PatchBoundaryWidth
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method UpdateScaleFactor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetScaleFactor $ScaleFactor
	}
    }

    method ExecuteSetScaleFactor { } {
	set tracervar [::itcl::local cvisTrace #auto]

	UpdateScaleFactor

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetScaleFactor {scalefactor} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$ScaleFactor != $scalefactor} {
	    set ScaleFactor $scalefactor
	    ExecuteSetScaleFactor
	}
    }

    method GetScaleFactor {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $ScaleFactor
    }

    method SetOffset {offset} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Offset != $offset} {
	    set Offset $offset
	    ExecuteSetOffset
	}
    }

    method GetOffset {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Offset
    }

    method ExecuteSetOffset {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetOffset $Offset
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetCellToNode {celltonode} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$CellToNode != $celltonode} {
	    set CellToNode $celltonode
	    ExecuteSetCellToNode
	}
    }

    method GetCellToNode {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $CellToNode
    }

    method ExecuteSetCellToNode {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetCellToNode $CellToNode
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method AddInput {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]

	# If this is the first plane then set defaults to a visible
	# plane

  	if {[string length $Input]} {
	    set setsliceflag 0
  	} {
	    set setsliceflag 1
	}

	lappend Input $sliceplane

	# if no input then set the values to be what this input 
	# is

	Update

	$sliceplane SetInterface $this
	
	if {$setsliceflag} {
  	    set Slice 0
	    ExecuteSetSlice
	}

	ExecuteSetSlicePlane
    }

    method ExecuteSetSlice {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    # Put slice in the middle of the cell

	    set real [expr (($Slice + 0.5) * $IndexScaling($SlicePlane)) \
		    + $RealBounds(min,$SlicePlane) ]
	    set screen [expr $real * $Scaling]

	    $plane SetSlice $screen

	    # Update real value in user display
	    set frame $NotebookPage.scale
	    set cs [$frame childsite]
	    $cs.realvalue configure -text $real
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetSlicePlane {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$SlicePlane != $sliceplane} {
	    set SlicePlane $sliceplane
	    ExecuteSetSlicePlane
	}
    }

    method GetSlicePlane {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $SlicePlane
    }

    method ExecuteSetSlicePlane {} {
	set tracervar [::itcl::local cvisTrace #auto]
	Update

	foreach plane $Input {
	    $plane SetSlicePlane $SlicePlane
	}

	ExecuteSetSlice
    }

    method GetSliceInScreenCoord {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set real [expr (($Slice + 0.5) * $IndexScaling($SlicePlane)) \
		+ $RealBounds(min,$SlicePlane) ]
	set screen [expr $real * $Scaling]
	return $screen
    }

    method SetDataVisible {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	set DataIsVisible $visible
	DataVisibleChange
    }

    method DataVisibleChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    # SGS this does not work yet
	    $plane SetDataVisibility $DataIsVisible
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set MaxScaleFactor 2.0
	set ScaleFactorDelta 0.01

	# Initialize values

	set item [lindex $Input 0]

	set Scaling [$item GetScaling]    

	scan [$item GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z)
	
	scan [$item GetBounds] "%f %f %f %f %f %f" \
		RealBounds(min,x) RealBounds(max,x) \
		RealBounds(min,y) RealBounds(max,y) \
		RealBounds(min,z) RealBounds(max,z)
	
	scan [$item GetIndexBounds] "%f %f %f %f %f %f" \
		IndexBounds(min,x) IndexBounds(max,x) \
		IndexBounds(min,y) IndexBounds(max,y) \
		IndexBounds(min,z) IndexBounds(max,z)

       # convert from screen to real coords
	foreach plane {x y z} {
	    set RealBounds(min,$plane) [expr $RealBounds(min,$plane) / $Scaling]      
	    set RealBounds(max,$plane) [expr $RealBounds(max,$plane) / $Scaling]
	}
	    
	set frame $NotebookPage.scaleframe
	set cs [$frame childsite]

	set frame $NotebookPage.scale
	set cs [$frame childsite]
	$cs.scale configure \
		-from $IndexBounds(min,$SlicePlane) \
		-to $IndexBounds(max,$SlicePlane)

	# ScaleFactor can change do due autoscaling if 
	# switching variables (hence MaxScaleFactor changes)
	
	# SGS Sept 5, 2002.
	# This looks like an old hack when we were doing
	# scaling factor calc here rather than in module itself.
#	UpdateScaleFactor 

	# SGS Sept 5, 2002.
	# Why do we need this here, patch width is not changed?
#	ExecuteSetPatchBoundaryWidth
    }

    method SetSlice {slice} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Slice $slice

	# Force value to be in the domain
	if {$Slice < $IndexBounds(min,$SlicePlane)} {
	    set Slice $IndexBounds(min,$SlicePlane)
	}
	    
	if {$Slice > $IndexBounds(max,$SlicePlane)} {
	    set Slice $IndexBounds(max,$SlicePlane)
	}

	ExecuteSetSlice
    }

    # SGS need to have methods to set the scaler precision etc

    # SGS method to set the min and max scalar or always do the
    # add input method?

    method SetNotebookPage {notebook} {
	set tracervar [::itcl::local cvisTrace #auto]
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# Slice Plane controls
	#-----------------------------------------------------------------
	#-----------------------------------------------------------------
	# Visibility
	#-----------------------------------------------------------------

	set frame $NotebookPage.visibility
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Visible" -variable [scope DataIsVisible] \
		-command "$this DataVisibleChange"

	#-----------------------------------------------------------------
	# X/Y/Z Plane Selector
	#-----------------------------------------------------------------
	set frame $NotebookPage.plane

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Slice Axis"
	set cs [$frame childsite]


	radiobutton $cs.x -text "X" \
		-variable [scope SlicePlane] \
		-value "x" \
		-command "$this ExecuteSetSlicePlane"

	radiobutton $cs.y -text "Y" \
		-variable [scope SlicePlane] \
		-value "y" \
		-command "$this ExecuteSetSlicePlane"

	radiobutton $cs.z -text "Z" \
		-variable [scope SlicePlane] \
		-value "z" \
		-command "$this ExecuteSetSlicePlane"


	#-----------------------------------------------------------------
	# Slice value
	#-----------------------------------------------------------------
	set frame $NotebookPage.scale
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Slice Location"
	set cs [$frame childsite]

	label $cs.realvalue \
		-text "[expr ($Slice + 0.5) * $IndexScaling($SlicePlane) + $RealBounds(min,$SlicePlane)]"

	scale $cs.scale \
		-length 100 \
		-from $IndexBounds(min,$SlicePlane) \
		-to $IndexBounds(max,$SlicePlane) \
		-variable [scope Slice] \
		-orient horizontal -highlightthickness 0

	# Need to add button binding here so when slider changes
	# we update the plane
	bind $cs.scale <ButtonRelease> "$this ExecuteSetSlicePlane"

	#-----------------------------------------------------------------
	# Height map Scaling 
	#-----------------------------------------------------------------
	set frame $NotebookPage.scaleframe

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Height Displacement"
	set cs [$frame childsite]

	label $cs.scalelabel -text "Scaling Factor"

	set MaxScaleFactor 2.0
	set ScaleFactorDelta 0.01

	scale $cs.scalefactor \
		-variable [scope ScaleFactor] \
		-length 100 \
		-from 0 -to $MaxScaleFactor \
		-resolution $ScaleFactorDelta -digits 3 \
		-orient horizontal -highlightthickness 0

	bind $cs.scalefactor <ButtonRelease> \
		"$this ExecuteSetScaleFactor"

	label $cs.offsettitle -text "Offset"
	entry $cs.offset \
		-textvariable [scope Offset] 
	bind $cs.offset <Return> "$this ExecuteSetOffset"

	#-----------------------------------------------------------------
	# CellToNode
	#-----------------------------------------------------------------

	set frame $NotebookPage.celltonode
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Cell To Node Conversion"
	set cs [$frame childsite]

	radiobutton $cs.cell \
		-variable [scope CellToNode] \
		-text "leave data at cell center" \
	        -value 0 \
		-command "$this ExecuteSetCellToNode"

	radiobutton $cs.node \
		-variable [scope CellToNode] \
		-text "average data to nodes" \
	        -value 1 \
		-command "$this ExecuteSetCellToNode"

	#-----------------------------------------------------------------
	# Patch Boundary 
	#-----------------------------------------------------------------

	set frame $NotebookPage.boundaries
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Boundary Drawing"
	set cs_main [$frame childsite]

	set frame $cs_main.patchboundary
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Patch"
	set cs [$frame childsite]

	label $cs.edgelabel -text "Edge Width"

	scale $cs.edgewidth \
		-length 20 \
		-from 1 -to 20 \
		-variable [scope PatchBoundaryWidth] \
		-orient horizontal -highlightthickness 0

	bind $cs.edgewidth <ButtonRelease> \
		"$this ExecuteSetPatchBoundaryWidth"

	checkbutton $cs.visible \
		-text "Visible" -variable [scope PatchBoundaryIsVisible] \
		-command "$this ExecuteSetPatchBoundaryVisibility"

	set level 0
	while {$level < $NumPatchColors} {
	    set frame $cs.level_$level
	    frame $frame
	    checkbutton $frame.level_visible \
		    -text "Level $level" -variable [scope PatchVisible($level)] \
		    -command "$this ExecuteSetPatchVisibility $level"
	    scan $PatchColors($level) "%f %f %f" r g b
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]
	    button $frame.level_color -bg $color \
		    -command "$this ExecuteSetPatchColor $level"
	    incr level
	}

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------

	set frame $NotebookPage.visibility
	set cs [$frame childsite]
	pack $cs.visible -anchor w 
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.plane
	set cs [$frame childsite]
	pack $cs.x -side left -pady 4
	pack $cs.y -side left
	pack $cs.z -side left
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.scale
	set cs [$frame childsite]
	pack $cs.realvalue -side right 
  	pack $cs.scale -pady 2 -expand yes -fill x 
	pack $cs -fill x
	pack $frame -fill x


	#-----------------------------------------------------------------
	set frame $NotebookPage.scaleframe
	set cs [$frame childsite]

	pack $cs.scalelabel
	pack $cs.scalefactor -expand yes -fill x

	pack $cs.offsettitle -side left
	pack $cs.offset -side left -padx 4 -pady 4

	pack $cs -fill x
	pack $frame -fill x
	
	#-----------------------------------------------------------------
	set frame $NotebookPage.celltonode
	set cs [$frame childsite]

	pack $cs.cell -anchor w -pady 4
	pack $cs.node -anchor w

  	pack $cs -fill x
  	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.boundaries
	pack $frame -fill x
	set cs_main [$frame childsite]
	pack $cs_main -fill x
	#-----------------------------------------------------------------
	set frame $cs_main.patchboundary
	set cs [$frame childsite]

	pack $frame -pady 4 -side left

	pack $cs.edgelabel 
	pack $cs.edgewidth -expand yes -fill x

	pack $cs.visible -anchor w

	set level 0
	while {$level < $NumPatchColors} {
	    set frame $cs.level_$level
	    pack $frame.level_visible -side left
	    pack $frame.level_color -side right
	    incr level
	    pack $frame -anchor w
	}
    }   
}

