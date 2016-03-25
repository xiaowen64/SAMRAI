##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisProperty.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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






