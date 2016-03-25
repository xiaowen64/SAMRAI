##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisColorBar.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: ColorBar to show data ranges
##

class cvisColorBar {
    inherit cvisDisplayObject cvisText

    variable LookupTable ""

    variable Title " "

    variable PositionX 0.01
    variable PositionY 0.05
    variable Width 0.15
    variable Height 1.0
    variable Orientation 1

    variable TitleSize 0.5
    variable LabelSize 0.5

    variable NumberOfLabels 5

    variable Visible 1

    variable Exists 0

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
    }

    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]

	::cvis::queue::Add "$this ExecuteSetTableRange"
    }

    method ExecuteSetTableRange {} {
	}
    }

    method SetVisibility {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Visible != $visible} {
	    set Visible $visible
	    }
	}
    }

    method GetColorBarActor { } {
	set tracervar [::itcl::local cvisTrace #auto]

	} {
	    return ""
	}
    }

    method SetLookupTable {lut} {
	set tracervar [::itcl::local cvisTrace #auto]
	$lut AddDepends $this
	set LookupTable $lut
	UpdateLookupTable
    }
    
    method UpdateLookupTable {} {
	set tracervar [::itcl::local cvisTrace #auto]
	}
    }
    
    method GetLookupTable { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return $LookupTable 
    }

    method SetTitle {title} { 
	# If title is "" then we get a core dump under Mesa
	# so force at least a space
	if {[string length $title]} {
	    set Title $title
	} {
	    set Title " "
	}

	}
    }

    method GetTitle {} {
	return $Title
    }

    variable LabelFormat "%1.2e"
    method SetLabelFormat {format} {
	set LabelFormat $format
	}
    }

    method GetLabelFormat {} {
	return $LabelFormat
    }

    # Methods needed to support cvisText API
    method SetColor {r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Color "$r $g $b"
	}
    }

    method SetFontFamily {font} {
	set tracervar [::itcl::local cvisTrace #auto]
	set FontFamily $font
	}
    }

    method SetBold {bold} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Bold $bold
	}
    }

    method SetItalic {italic} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Italic $italic
	}
    }

    method SetShadow {shadow} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Shadow $shadow
	}
    }

    method ExecuteCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {![string length $LookupTable]} {
	    return
	}

	if {$Exists} {
	    return
	}

	set Exists 1





	
    }

    method ExecuteDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Exits 0
    }

    method Update {} { 
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {![string length $RenderWindow] } { 
	    return
	}

	if {!$Exists} {
	    ::cvis::queue::Add "$this ExecuteCreate"
	}
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
	Update
    }
    
    method SetWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]

	set Width $width
	
	}
    }

    method GetWidth { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return $Width 
    }

    method SetHeight {height} {
	set tracervar [::itcl::local cvisTrace #auto]

	set Height $height
	
	}
    }

    method GetHeight { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return $Height
    }

    method SetPosition {x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	set PositionX $x
	set PositionY $y
	
	}
    }

    method GetPosition { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return "$PositionX $PositionY"
    }

    method SetOrientation {orientation} {
	set tracervar [::itcl::local cvisTrace #auto]

	set Orientation $orientation
	
	}
    }

    method GetOrientation { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return "$Orientation"
    }

    method SetLabelSize {labelsize} {
	set tracervar [::itcl::local cvisTrace #auto]

	set LabelSize $labelsize
	
	}
    }

    method GetLabelSize { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return "$LabelSize"
    }

    method SetTitleSize {titlesize} {
	set tracervar [::itcl::local cvisTrace #auto]

	set TitleSize $titlesize
	
	}
    }

    method GetTitleSize { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return "$TitleSize"
    }

    method SetNumberOfLabels {number} {
	set tracervar [::itcl::local cvisTrace #auto]

	set NumberOfLabels $number
	
	}
    }

    method GetNumberOfLabels { } {
	set tracervar [::itcl::local cvisTrace #auto]

	return "$NumberOfLabels"
    }

}










