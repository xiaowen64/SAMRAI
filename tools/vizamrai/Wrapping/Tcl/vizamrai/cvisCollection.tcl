##
## File:        cvisCollection.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Collection class
##
class cvisCollection {
    variable Collection ""
    variable Position 0

    method AddItem {item} {
	lappend Collection $item
    }

    method RemoveItem {item} {
	set remove_location [lsearch $Collection $item]
	set Collection [lreplace $Collection $remove_location $remove_location]
    }

    method GetNumberOfItems {} {
	return [llength $Collection]
    }

    method InitTraversal {} {
	set Position 0
    }

    method GetNextItem {} {
	if { $Position > [llength $Collection] } {
	    set value ""
	} {
	    set value [lindex $Collection $Position]
	    incr Position
	}
	return $value
    }
}
