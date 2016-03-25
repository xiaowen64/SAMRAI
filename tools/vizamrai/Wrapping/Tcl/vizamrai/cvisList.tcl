##
## File:        cvisList.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Linked List
##


class cvisListElement {
    variable Next ""
    variable Previous ""
    
    variable Data

    method GetNext { } {
	return $Next
    }

    method SetNext {next} {
	set Next $next
    }

    method GetPrevious { } {
	return $Previous
    }

    method SetPrevious {previous} {
	set Previous $previous
    }

    method GetData { } {
	return $Data
    }

    method SetData {data} {
	set Data $data
    }
    
    method GetReference { } {
	return $this
    }
}

class cvisListHead {
    variable Head ""
    variable Tail ""
    
    method GetTail { } {
	return $Tail
    }
    
    method SetTail {tail } {
	set Tail $tail
    }
    
    method GetHead { } {
	return $Head
    }
    
    method AddElement { data } {
	set new_element [[cvisListElement #auto] GetReference]

	$new_element SetPrevious $Tail
	$new_element SetNext ""
	$new_element SetData $data
	
	if {[string length $Tail]} {
	    $Tail SetNext $new_element
	    set Tail $new_element
	} {
	    set Head $new_element
	    set Tail $new_element
	}
	
    }

    method Merge {mergelist } {
	set tail [$mergelist GetTail]
	if {[string length $tail]} {
	    $tail SetNext $Head
	}
	$mergelist SetTail $Tail
    }
    
    method GetReference { } {
	return $this
    }
}

class cvisList { 
    
    variable Current ""
    variable ListHead ""
    
    constructor { } {
	set ListHead [[cvisListHead #auto] GetReference]
	return $this
    }
    
    method GetReference { } {
	return $this
    }

    method GetListHead { } {
	return $ListHead
    }
    
    # SGS How do we delete?
    
    method Print { } {
	InitTraversal
	set item [GetNextItem]
	while {[string length $item]} {
	    puts "$item --[$item GetData]"
	    set item [GetNextItem]
	}
    }
    
    method AddElement { data } {    
	$ListHead AddElement $data
    }
    
    method Merge { mergelist } {
	# Make the two lists one and the same by appending this list
	# to the end of the argument list

	$ListHead Merge [$mergelist GetListHead]

	set ListHead [$mergelist GetListHead]
    }
    
    method InitTraversal { } {
	# Initialize the list traversal
	set Current [$ListHead GetHead]
    }
    
    method GetNextItem { } {
	if {[string length $Current]} {
	    # Return the current and point to next for next iteration
	    set ret $Current
	    set Current [$Current GetNext]
	} {
	    # If we are at the end of the list return null
	    set Current ""
	    set ret ""
	}
	
	return $ret
    }

}




