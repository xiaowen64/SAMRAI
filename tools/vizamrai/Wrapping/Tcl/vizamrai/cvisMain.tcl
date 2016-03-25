##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisMain.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Main Window
##

class cvisMain {
    variable ParentWindow ""
    variable TabbedNotebook ""

    variable WindowNameProp

    constructor { parent } {
	set ParentWindow $parent
	iwidgets::menubar $ParentWindow.menuBar -menubuttons {
	    menubutton file -text "File" -menu {
		options -tearoff false

		command Open -label "Open" -command {$this Open} \
			-helpstr "Open saved state file"

		command SaveAs -label "SaveAs" -command {$this SaveAs} \
			-helpstr "Save state to a file"

		command exit -label "Exit" -command {$this ExitApplication} \
			-helpstr "Exit Vizamrai"

	    }
	    menubutton help -text "Help" -menu {
		options -tearoff false

		command help -label "Help Contents" \
			-command {$this HelpContents} \
			-helpstr "Help Contents"
		
		command about -label "About" \
			-command {$this About} \
			-helpstr "About Vizamrai"
	    }
	}

	pack $ParentWindow.menuBar -side top -fill x
	set TabbedNotebook [cvisTabbedNotebook notebook $ParentWindow]

	set WindowNameProp [cvisProperty $this.WindowName]
	$WindowNameProp AddCallback "$this UpdateWindowName"
    }

    variable FileName ""
    variable FileTypes {{{Vizamrai State File} {.vstate}}}

    destructor {
  	delete object $WindowNameProp
    }

    method Open {} {
	set ImageFileBrowser [cvisFileBrowser #auto]
	$ImageFileBrowser SetTitle "Select Vizamrai Save State File"
	$ImageFileBrowser SetFileTypes $FileTypes
	$ImageFileBrowser SetOpenPrompt
	set FileName [$ImageFileBrowser Prompt]
	::cvis::RestartFile::Open $FileName
    }

    method SaveAs {} {
	set ImageFileBrowser [cvisFileBrowser #auto]
	$ImageFileBrowser SetTitle "Select Vizamrai Save State File"
	$ImageFileBrowser SetFileTypes $FileTypes
	$ImageFileBrowser SetSavePrompt
	set FileName [$ImageFileBrowser Prompt]
	::cvis::RestartFile::Save $FileName
    }

    method HelpContents { } {
	set tracervar [::itcl::local cvisTrace #auto]
	global cvis_library 
	if { [winfo exists .help] } {
	    .help activate 
	} {
	    iwidgets::hyperhelp .help \
		    -title "Vizamrai Help" \
		    -modality none \
		    -topics {Contents} \
		    -helpdir [file join $cvis_library help]
	    .help showtopic Contents
	    .help activate
	}
    }
    
    method About { } {
	set tracervar [::itcl::local cvisTrace #auto]
	::cvis::gui::About
    }

    method Pack { } {
	set tracervar [::itcl::local cvisTrace #auto]

	$TabbedNotebook Pack
    }
    
    public method AddInterface {interface} {
	set tracervar [::itcl::local cvisTrace #auto]
	$TabbedNotebook AddInterface $interface
    }

    method ExitApplication {} {
	::cvis::Exit
    }

    method GetWindowName {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$WindowNameProp GetValue]
    }

    method GetWindowNameProp {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$WindowNameProp GetReference]
    }

    method SetEquivWindowName { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$WindowNameProp SetEquiv $property
    }

    method UpdateWindowName { windowname } {
	set tracervar [::itcl::local cvisTrace #auto]
	wm title . "Vizamrai CP: $windowname"
    }

}


