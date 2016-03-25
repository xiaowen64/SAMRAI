##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisBaseClasses.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Base classes and global variables
##

namespace eval cvis {
    namespace eval options {
	variable DEFAULT_CAMERA "default.cam"
	variable LOWMEM 0
    }

    namespace eval math {
	variable DEGREESTORADIANS 0.0174533
	variable FLT_MAX 3.402823466E+38
	variable FLT_MIN 1.175494351E-38
	variable FLT_EPSILON 1.192092896E-07
	variable INT_MAX 2147483647 
	variable INT_MIN -2147483648

	set BOUNDARY_JITTER 0.0001
	set BOUNDARY_TOLERANCE 0.00001

	set JITTER 0.02
    }

    namespace eval RestartFile {
	variable ObjectList ""

	proc Open {filename} {
	    # SGS add checking for file etc
	    # SGS turn off rendering so we don't redraw 
	    #     overly much
	    ::cvis::queue::SetPipeline 0
	    source $filename
	    ::cvis::queue::SetPipeline 1
	}
	
	proc Save {filename} {
	    set file [open $filename "w"]
	    foreach item $::cvis::RestartFile::ObjectList {
		$item Save $file
	    }
	    close $file
	}
	
	proc Add {object} {
	    lappend ::cvis::RestartFile::ObjectList $object
	}
    }

    namespace eval queue {
	variable queue {}
	variable state 1

	proc GetPipeline {} {
	    return $::cvis::queue::state
	}

	# Turn execution on/off
	proc SetPipeline { state } {
	    set tracervar [::itcl::local cvisTrace #auto]
	    set ::cvis::queue::state $state
	    $tracervar print COMM "State = $state $::cvis::queue::state"
	    ::cvis::queue::Execute
	}

	# Add an event to execute
	proc Add { event } {
	    set tracervar [::itcl::local cvisTrace #auto]
	    $tracervar print COMM "Adding Event= $event"
	    lappend ::cvis::queue::queue $event
	    $tracervar print COMM "PostAdd Queue "
	    foreach event $::cvis::queue::queue {
		$tracervar print COMM $event
	    }
	}
	
	# Execute all of the queued events
	proc Execute { } {
	    set tracervar [::itcl::local cvisTrace #auto]

	    $tracervar print COMM "State = $::cvis::queue::state"

	    if {$::cvis::queue::state} {

		::cvis::BusyExecute {
		    set trimed_queue ""
		    
		    $tracervar print COMM "PreTrim Queue "
		    
		    foreach event $::cvis::queue::queue {
			$tracervar print COMM $event
		    }
		    
		    set queue $::cvis::queue::queue
		    set ::cvis::queue::queue {}
		    
		    while {[llength $queue]} {
			set event [lindex $queue 0]
			set queue [lrange $queue 1 end]
			
			set object [lindex $event 0]
			set method [lindex $event 1]
			set args [lrange $event 2 end]
			
			set found_equiv 0
			foreach event $queue {
			    if {![string compare $object [lindex $event 0]]} {
				if {![string compare $method [lindex $event 1]]} {
				    set found_equiv 1
				}
			    }
			}
			
			if {!$found_equiv} {
			    $tracervar print COMM "Adding $object $method $args"
			    lappend trimed_queue "$object $method $args"
			}
		    }
		    
		    
		    $tracervar print COMM "PostTrimQueue "
		    foreach event $trimed_queue {
			$tracervar print COMM $event
		    }
		    
		    # 
		    foreach event $trimed_queue {
			$tracervar print COMM "Executing task $event"
			eval $event
		    }
		}
	    }
	}
    }

    namespace eval gui {} {

	variable PanelHeaderPosition nw
	
	proc About {} {
	    set tracervar [::itcl::local cvisTrace #auto]
	    if {[winfo exists .about]} {
		.about activate
	    } {
		global cvisName
		global cvisRelease
		global cvisDate
		global cvis_library
		
		iwidgets::dialogshell .about -title "About"
		set frame [.about childsite]
		.about add dismiss -text "Dismiss" -command {.about deactivate}
		
		button $frame.vtk -text "VTK" \
			-image [image create photo -file \
			[file join $cvis_library images vtk.gif]]
		
		button $frame.tcltk -text "TCLTK" \
			-image [image create photo -file \
			[file join $cvis_library images tcltk.gif]]

		label $frame.vizamraiTitle \
			-anchor w \
			-justify left \
			-font {-family times -size 16 -weight bold} \
			-text "Vizamrai"
		pack $frame.vizamraiTitle -expand yes -fill both
		
		label $frame.vizamraiInfo \
			-anchor w \
			-justify left \
			-font {-family times -size 12 -weight bold} \
			-text "Copyright 1997-2001\nLawrence Livermore National Security, LLC\n$cvisRelease  $cvisDate"
		pack $frame.vizamraiInfo -expand yes -fill both
		
		frame $frame.copyright
		text $frame.copyright.text \
			-yscrollcommand "$frame.copyright.scrollbar set"
		
		scrollbar $frame.copyright.scrollbar -orient vert -command "$frame.copyright.text yview"
		
		$frame.copyright.text insert end "SAMRAI code release number: UCRL-CODE-2002-004
		
Copyright Notice
----------------
		
(c) 1997-2001
Lawrence Livermore National Security, LLC.
All rights reserved.
	
This work was produced at the University of California, Lawrence Livermore
National Laboratory (UC LLNL) under contract no. W-7405-ENG-48 (Contract 48)
between the U.S. Department of Energy (DOE) and The Regents of the University
of California (University) for the operation of UC LLNL.  The rights of the 
Federal Government are reserved under Contract 48 subject to the restrictions 
agreed upon by the DOE and University as allowed under DOE Acquisition Letter 
97-1.

Disclaimer
----------

This work was prepared as an account of work sponsored by an agency of
the United States Government.  Neither the United States Government nor
the University of California nor any of their employees makes any warranty,
express or implied, or assumes any legal liability or responsibility for
the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not
infringe privately-owned rights.  Reference  herein to any specific
commercial products, process, or service by trade name, trademark,
manufacturer, or otherwise, does not necessarily constitute or imply
its  endorsement, recommendation, or favoring by the United States
Government or the University of California.  The views and opinions
of authors expressed herein do not necessarily state or reflect those
of the United States Government or the University of California, and
shall not be used for advertising or product endorsement purposes.

Notification of Commercial Use
------------------------------

Commercialization of this product is prohibited without notifying the 
Department of Energy (DOE) or Lawrence Livermore National Laboratory (LLNL)."

		pack $frame.copyright.scrollbar -side right -fill y
		pack $frame.copyright.text -side left -fill both -expand true
		pack $frame.copyright

		pack $frame.vtk -side left
		pack $frame.tcltk -side left

		vtkVersion version
		set vtk_version [version GetVTKVersion]
		version Delete
		
		label $frame.vtkInfo \
			-anchor w \
			-text "VTK  Version $vtk_version"
		pack $frame.vtkInfo -expand yes -fill both -padx 4 -pady 4
		
		global tk_patchLevel
		label $frame.tclInfo \
			-anchor w \
			-text "TCK/TK  Version: $tk_patchLevel"
		pack $frame.tclInfo -expand yes -fill both -padx 4 -pady 4
		
		
		.about activate
	    }
	}
    }

    variable exit_flag 0
    proc Exit {} {
	set ::cvis::exit_flag 1
    }

    proc WaitForExit {} {
	vwait ::cvis::exit_flag
    }

    # Usage: withBusyCursor { script ... }
    #
    # This is from ASPN (aspn.activestate.com) 
    # Written by Joe English
    # License:
    # The source code in this recipe is in the public domain.
    # 
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    # ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
    # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
    # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY 
    # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
    # USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    proc BusyExecute {body} {
	global errorInfo errorCode
	set busy {}
	set list {.}
	# Traverse the widget hierarchy to locate widgets with 
	# a nondefault -cursor setting.
	#
	while {$list != ""} {
	    set next {}
	    foreach w $list {
		catch {set cursor [$w cget -cursor]}
		if {[winfo toplevel $w] == $w || $cursor != ""} {
		    lappend busy $w $cursor
		    set cursor {}
		}
		set next [concat $next [winfo children $w]]
	    }
	    set list $next
	}
	
	# Change the cursor:
	#
	foreach {w _} $busy {
	    catch {$w configure -cursor watch}
	}
	update 

	# Execute the script body.
	#
	set rc [catch {uplevel 1 $body} result]
	set ei $errorInfo
	set ec $errorCode

	# Restore the original cursor settings.
	#
	foreach {w cursor} $busy {
	    catch {$w configure -cursor $cursor}
	}
	
	# Return script result to caller.
	#
	return -code $rc -errorinfo $ei -errorcode $ec $result
    }
}

class cvisDisplayObject {
    variable RenderWindow ""
    variable Name
    variable Create 1
    variable Interface

    constructor { name} {
	set Name $name
    }

    method GetName {} {
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
    }

    method GetRenderWindow { } {
	return $RenderWindow
    }
    
    method SetInterface {interface} {
	set Interface interface
    }

    method GetInterface {} {
	return $Interface
    }

    method SetInitialize { } {
	set Create 1
    }
}

class cvisInterface {
    variable NotebookPage ""
    variable Name ""
    variable RenderWindow ""

    constructor {name} {
	set Name $name

	# Add All interfaces to the restart file
	::cvis::RestartFile::Add $this
    }

    # Default restart save function does nothing
    # should implement commands to reset state for 
    # each interface
    method Save {file} {
    }

    method GetName {} {
    }

    method SetNotebookPage {notebookpage} {
	set NotebookPage $notebookpage
    }

    method GetNotebookPage { } {
	return $NotebookPage
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
    }

    method GetRenderWindow { } {
	return $RenderWindow
    }
}

class cvisText {
    variable Color "1 1 1"
    variable FontFamily 0
    variable Bold 1
    variable Italic 1
    variable Shadow 1

    method SetColor {r g b} {
	set $Color "$r $g $b"
    }

    method GetColor {color} {
	return $Color
    }

    method SetFontFamily {font} {
	set FontFamily $font
    }
    
    method GetFontFamily {} {
	return $FontFamily
    }

    method SetBold {bold} {
	set Bold $bold
    }

    method GetBold {} {
	return $Bold
    }

    method SetItalic {italic} {
	set Italic $italic
    }

    method GetItalic {} {
	return $Italic
    }

    method SetShadow {shadow} {
	set Shadow $shadow
    }
    
    method GetShadow {} {
	return $Shadow
    }
}

