##
## File:        cvisProperty.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Property for "joining" variables in seperate components
##

class cvisProperty { 
    variable Value ""
    variable CallbackList 

    constructor { } {
	set CallbackList [[cvisList #auto] GetReference]
    }
    
    method SetValue {value} {
	set Value $value
	DoCallbacks
    }

    method GetValue {} {
	return $Value
    }
    
    method GetCallbackList {} {
	return $CallbackList
    }
    
    method AddCallback {cmd} {
	$CallbackList AddElement $cmd
    }

    method SetEquiv { property } {
	$CallbackList Merge [$property GetCallbackList]
    }

    method Print { } {
	$CallbackList Print
    }

    method GetReference { } {
	return $this
    }

    method DoCallbacks { } {
	$CallbackList InitTraversal
	set item [$CallbackList GetNextItem] 
	while {[string length $item]} {
	    eval [$item GetData] {$Value}
	    set item [$CallbackList GetNextItem]
	}
    }
}






