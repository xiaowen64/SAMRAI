class cvisDataBrowserInterface {
    inherit cvisInterface

    variable RenderWindowIsVisible 1

    variable DataIsVisible 1
    variable LabelIsVisible 1

    variable CellBoundaryWidth 3
    variable CellBoundaryIsVisible 1

    variable Input ""

    variable MaxBounds

    variable LabelFormat "%1.2e"

    constructor {name} {
	cvisInterface::constructor $name
    } {
    }

    method SetFormat {} {
	foreach plane $Input {
	    $plane SetLabelFormat $LabelFormat
	}
    }

    method AddInput {browser} {
	lappend Input $browser
	$browser SetInterface $this

	foreach plane "x y z" { 
	    set Bounds(min,$plane) $::cvis::math::FLT_MAX
	    set Bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	}

	foreach item $Input {
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    Bounds(min,x) Bounds(max,x) \
		    Bounds(min,y) Bounds(max,y) \
		    Bounds(min,z) Bounds(max,z)

	    set max_bounds [expr -$::cvis::math::FLT_MAX]
	    foreach i "x y z" {
		set t [expr $Bounds(max,$i) - $Bounds(min,$i)]
		if {$max_bounds < $t} { 
		    set max_bounds $t
		}
	    }
	}
	

	set MaxBounds $max_bounds
	SetCellBoundaryWidth
    }

    method SetCellBoundaryWidth {} {
	foreach plane $Input {
	    $plane SetCellBoundaryEdgeWidth \
		    [expr $CellBoundaryWidth*$MaxBounds*0.001]
	}
    }

    method SetCellBoundaryVisible {} {
	foreach plane $Input {
	    $plane SetCellBoundaryVisibility $CellBoundaryIsVisible
	}
    }

    method SetRenderWindowVisibility {} {
	foreach plane $Input {
	    $plane SetRenderWindowVisibility $RenderWindowIsVisible
	}
    }

    method SetDataVisibility {} {
	foreach plane $Input {
	    $plane SetDataVisibility $DataIsVisible
	}
    }

    method SetLabelVisibility {} {
	foreach plane $Input {
	    $plane SetLabelVisibility $LabelIsVisible
	}
    }

    method Update {} {
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook
	#-----------------------------------------------------------------
	# Cell Boundary 
	#-----------------------------------------------------------------
	set frame $NotebookPage.window
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Data Browser Window"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	checkbutton $frame.visible \
		-text "Window Visible" -variable [scope RenderWindowIsVisible] \
		-command "$this SetRenderWindowVisibility"

	checkbutton $frame.datavisible \
		-text "Data Visible" -variable [scope DataIsVisible] \
		-command "$this SetDataVisibility"



	#-----------------------------------------------------------------
	# Text Label
	#-----------------------------------------------------------------
	set frame $NotebookPage.label
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Data Text Labels"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken


	label $frame.labelformattext \
		-text "Label Format"
	
	entry $frame.labelformat \
		-textvariable [scope LabelFormat]

	bind $frame.labelformat <Return> "$this SetFormat"

	checkbutton $frame.labelvisible \
		-text "Text Label Visible" -variable [scope LabelIsVisible] \
		-command "$this SetLabelVisibility"

	#-----------------------------------------------------------------
	# Cell Boundary 
	#-----------------------------------------------------------------
	set frame $NotebookPage.cellboundary
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Cell Boundary "
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	label $frame.edgelabel -text "Edge Width"

	scale $frame.edgewidth \
		-length 20 \
		-from 1 -to 20 \
		-variable [scope CellBoundaryWidth] \
		-orient horizontal -highlightthickness 0

	bind $frame.edgewidth <ButtonRelease> \
		"$this SetCellBoundaryWidth"

	checkbutton $frame.visible \
		-text "Visible" -variable [scope CellBoundaryIsVisible] \
		-command "$this SetCellBoundaryVisible"


	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.window
	pack $frame.title
	pack $frame.titleseperator -fill x

	pack $frame.visible -anchor w
	pack $frame.datavisible -anchor w 


	pack  $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.label

	pack $frame.title
	pack $frame.titleseperator -fill x

	pack $frame.labelvisible -anchor w -side bottom

	pack $frame.labelformattext -anchor w -side left
	pack $frame.labelformat -anchor w 

	pack  $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.cellboundary
	pack $frame.title
	pack $frame.titleseperator -fill x

	pack $frame.edgelabel
	pack $frame.edgewidth -expand yes -fill x

	pack $frame.visible -anchor w

	pack  $frame -fill x
    }   
}

