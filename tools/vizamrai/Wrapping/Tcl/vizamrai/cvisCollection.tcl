##
## File:        cvisCollection.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
