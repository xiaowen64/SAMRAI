##
## File:        cvisSlicePlaneInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Slice plane interface
##

class cvisOperatorInterface {
    inherit cvisInterface

    variable Input ""

    variable FileName "temp.ascii"
    variable FileNameProp
    variable FileTypes {{ascii {.ascii}}}
    variable AutoSave 0

    variable ImageFileBrowser ".null"

    constructor {name} {
	cvisInterface::constructor $name
    } {
	set FileNameProp [cvisProperty $name.FileName]
	$name.FileName AddCallback "$this UpdateFileName"
    }

    method Save {file} {
	if {[string length $FileName]} {
	    puts $file "$this SetFileName \"$FileName\""
	}
	puts $file "$this SetAutoSave $AutoSave"
    }

    method AddInput {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]

	set Slice 0

	lappend Input $sliceplane

	Update

	$sliceplane SetInterface $this
    }

    method GetFileName { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$FileNameProp GetValue]
    }
    
    method GetFileNameProp {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$FileNameProp GetReference]
    }

    method SetFileName { filename } {
	set tracervar [::itcl::local cvisTrace #auto]
	set FileName $filename
	$FileNameProp SetValue $filename
    }

    method SetEquivWindowName { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$FileNameProp SetEquiv $property
    }

    method UpdateFileName {filename} {
	set tracervar [::itcl::local cvisTrace #auto]
	set temp "[file rootname $filename].ascii"
	set FileName $temp

	foreach item $Input {
	    $item SetFileName $FileName
	}
    }

    method SetAutoSave {autosave} {
	if {$AutoSave != $autosave} {
	    set AutoSave $autosave
	    ExecuteSetAutoSave
	}
    }

    method ExecuteSetAutoSave {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach item $Input {
	    $item SetAutoSave $AutoSave
	}
    }

    method SaveFileAsASCII { } {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach item $Input {
	    $item SaveFileAsASCII 
	}
    }

    method FileNameBrowser { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $ImageFileBrowser]]} {
	    $ImageFileBrowser SetSavePrompt
	} {
	    set ImageFileBrowser [cvisFileBrowser #auto]
	    $ImageFileBrowser SetTitle "Select Image File"
	    $ImageFileBrowser SetFileTypes $FileTypes
	    $ImageFileBrowser SetSavePrompt
	}

	set FileName [$ImageFileBrowser Prompt]

	SetFileName $FileName
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method SetNotebookPage {notebook} {
	set tracervar [::itcl::local cvisTrace #auto]
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# Construct controls for operator output
	#-----------------------------------------------------------------
	set frame $NotebookPage.filenameframe

	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Operator Output"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	set frame $NotebookPage.file
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Output File"
	set cs [$frame childsite]

	button $cs.save \
		-text "Save Current Data" \
		-command "$this SaveFileAsASCII"
	
	button $cs.selectfile \
		-text "FileName" \
		-command  "$this FileNameBrowser" 

	entry $cs.filename -textvariable [scope FileName]

  	checkbutton $cs.autosave \
  		-text "Auto Save on New Load" -variable [scope AutoSave] \
  		-command "$this ExecuteSetAutoSave"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.filenameframe
	pack $frame.title
	pack $frame.titleseperator -fill x

	set frame $NotebookPage.file
	set cs [$frame childsite]

	pack $cs.save -side top -anchor w \
		-padx 4 -pady 8
	pack $cs.autosave  -side bottom -anchor w -pady 4
	pack $cs.selectfile -side left \
		-padx 4 -pady 4
	pack $cs.filename -side left \
		-fill x -expand yes -pady 4
	pack $cs -fill x 

	pack $frame -fill x
   }   
}

