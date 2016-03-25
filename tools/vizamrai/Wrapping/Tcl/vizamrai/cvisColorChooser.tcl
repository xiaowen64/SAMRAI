##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisColorChooser.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Release:     $$
## Revision:    $$
## Modified:    $$
## Description: Color chooser widget class
##

class cvisColorChooser {
    variable rcolor 0
    variable gcolor 0
    variable bcolor 0
    
    method SetColor {r g b} {
	set rcolor $r
	set gcolor $g
	set bcolor $b
    }

    method GetColor {} {
	return "$rcolor $gcolor $bcolor"
    }

    method Prompt {} {
	set color [format "#%02x%02x%02x" \
		[expr int($rcolor * 255.0)] \
		[expr int($gcolor * 255.0)] \
		[expr int($bcolor * 255.0)]]

	set value [tk_chooseColor -initialcolor $color]
	if {[string length $value]} {
	    scan $value "#%2x%2x%2x" r g b
	    set rcolor [expr $r / 255.0]
	    set gcolor [expr $g / 255.0]
	    set bcolor [expr $b / 255.0]
	} {
	    return ""
	}
	return "$rcolor $gcolor $bcolor"
    }
}
