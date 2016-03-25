##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisTabbedNotebook.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Tabbed Notebook
##


class cvisTabbedNotebook {
    variable ParentWindow ""
    variable TabbedNotebook ""

    variable MaxHeight 0
    variable MaxWidth 0

    variable Frames ""

    constructor { parent } {
	set ParentWindow $parent
	set TabbedNotebook [iwidgets::tabnotebook $parent.tnb]
	$TabbedNotebook configure -tabpos w
	pack $TabbedNotebook  -expand yes -fill both -padx 4 -pady 4
    }

    method Pack { } {

	$TabbedNotebook select 0

	update idletasks

	# Try to compute the size for the window
	foreach frame $Frames {
	    set page $frame.canvas.frame
	    set w [winfo reqwidth $page]
	    if {$w > $MaxWidth} {
		set MaxWidth $w
	    }
	    
	    set h [winfo reqheight $page]
	    if {$h > $MaxHeight} {
		set MaxHeight $h
	    }
	}

	# SGS This is a hack, when switching to using iwidget 
	# labeled frames the above code no longer worked
	# What's going on?
	set MaxHeight 750

	foreach frame $Frames {
	    $frame.canvas config -scrollregion "0 0 $MaxWidth $MaxHeight" \
		    -width $MaxWidth -height $MaxHeight

	    set page $frame.canvas.frame
	    $frame.canvas create window 0 0 -anchor nw -window $page \
		    -width $MaxWidth -height $MaxHeight
	}

	set MaxWidth [expr 20 + $MaxWidth + \
		[winfo reqwidth $TabbedNotebook.canvas.tabset]]

	set h [winfo reqheight $TabbedNotebook.canvas.tabset]
	if {$h > $MaxHeight} {
	    set MaxHeight $h
	}
	
	$TabbedNotebook configure -width $MaxWidth -height $MaxHeight
    }
    
    public method AddInterface {interface} {
	
	set frame [$TabbedNotebook add -label [$interface GetName]]

	frame $frame.scrollframe
	set frame $frame.scrollframe

	# 340x520
	canvas $frame.canvas \
		-width 10 \
		-height 10 \
		-yscrollcommand [list $frame.yscrollbar set] \
		-xscrollcommand [list $frame.xscrollbar set] 

	scrollbar $frame.xscrollbar -orient horizontal \
		-command [list $frame.canvas xview]
	scrollbar $frame.yscrollbar -orient vertical \
		-command [list $frame.canvas yview]

	pack $frame.xscrollbar -side bottom -fill x 
	pack $frame.yscrollbar -side right -fill y
	pack $frame.canvas -expand y -fill both
	pack $frame -expand y -fill both

	set page $frame.canvas.frame
	frame $page -relief flat
	# SGS can we calc the width/height 

	$interface SetNotebookPage $page

	lappend Frames $frame
    }
}


