##
## File:        cvisFileBrowser.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: File Browser 
##


class cvisFileBrowser {
    variable FileName ""
    variable Directory "."
    variable Title ""
    variable FileTypes ""
    variable DefaultExtension ""

    variable IsOpen 1

    method SetOpenPrompt {} {
	set IsOpen 1
    }

    method SetSavePrompt {} {
	set IsOpen 0
    }

    method SetFileName {name} {
	set FileName $name
    }
    
    method GetFileName {} {
	return $FileName
    }

    method SetDirectory {directory} {
	set Directory $directory
    }

    method GetDirectory { } {
	return $Directory
    }

    method SetTitle {title} {
	set Title $title
    }
    
    method GetTitle { } {
	return $Title
    }
    
    method SetFileTypes {types} {
	set FileTypes $types
    }

    method GetFileTypes {} {
	return $FileTypes
    }

    method SetDefaultExtension {ext} {
	set DefaultExtension $ext
    }

    method GetDefaultExtension {} {
	return $DefaultExtension
    }

    method Prompt {} {
	if $IsOpen {
	    set filename [tk_getOpenFile -title $Title \
		    -filetypes $FileTypes \
		    -defaultextension $DefaultExtension \
		    -initialdir $Directory]
	} {
	    set filename [tk_getSaveFile -title $Title \
		    -filetypes $FileTypes \
		    -defaultextension $DefaultExtension \
		    -initialdir $Directory]
	}

	# If user did not select new name then keep old name
	if {[string length $filename]} {
	    set FileName $filename
	    set Directory [file dirname $filename]
	}

	return $FileName
    }
}
