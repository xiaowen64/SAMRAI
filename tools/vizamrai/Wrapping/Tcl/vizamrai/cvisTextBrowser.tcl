class cvisTextBrowser {
    variable Color "1 1 1"
    variable FontFamily 0
    variable Bold 0
    variable Italic 0
    variable Shadow 0

    variable Title "Select Font"

    variable RenderWindow ""    
    variable TextObject ""

    variable Topwindow ""

    constructor {} {
	# Save this in the restart file
	::cvis::RestartFile::Add $this
    } {
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
    }

    method GetRenderWindow {} {
	return $RenderWindow
    }

    method SetInput {input} {
	set TextObject $input
    }

    method GetInput {} {
	return $TextObject
    }

    method Save {file} {
	puts $file "$this SetColor $Color"
	puts $file "$this SetFontFamily $FontFamily"
	puts $file "$this SetBold $Bold"
	puts $file "$this SetItalic $Italic"
	puts $file "$this SetShadow $Shadow"
	puts $file "$this SetTitle {$Title} "
    }

    method OK {} {
	wm withdraw $Topwindow
    }

    method SetTitle {title} {
	set Title $title
	if {[winfo exists $Topwindow]} {
	    wm title $Topwindow "$Title"
	}
    }

    method GetTitle {} {
	return $Title
    }

    method SetColor {r g b} {
	set $Color "$r $g $b"
	$TextObject SetColor $r $g $b
	$RenderWindow Render

	set frame $Topwindow.color
	if {[winfo exists $frame.color]} {
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]

	    $frame.color configure -bg $color
	}
    }

    method GetColor {color} {
	return $Color
    }

    method SelectColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $Color
	
	set value [color Prompt]
	
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b
	    SetColor $r $g $b
	}	

	delete object color 
    }

    method SetFontFamily {font} {
	set FontFamily $font
	$TextObject SetFontFamily $FontFamily
	$RenderWindow Render
    }
    
    method GetFontFamily {} {
	return $FontFamily
    }

    method SetBold {bold} {
	set Bold $bold
	ExecuteSetBold
    }

    method GetBold {} {
	return $Bold
    }

    method ExecuteSetBold {} {
	$TextObject SetBold $Bold
	$RenderWindow Render
    }

    method SetItalic {italic} {
	set Italic $italic
	ExecuteSetItalic
    }

    method GetItalic {} {
	return $Italic
    }

    method ExecuteSetItalic {} {
	$TextObject SetItalic $Italic
	$RenderWindow Render
    }

    method SetShadow {shadow} {
	set Shadow $shadow
	ExecuteSetShadow
    }
    
    method GetShadow {} {
	return $Shadow
    }

    method ExecuteSetShadow {} {
	$TextObject SetShadow $Shadow
	$RenderWindow Render
    }

    method Prompt {} {
	if {[winfo exists $Topwindow]} {
	    wm deiconify $Topwindow
	} {
	    # Make a valid window name out of $this
	    regsub -all "::" $this "" Topwindow
	    regsub -all "\\\." $Topwindow "" Topwindow
	    set Topwindow .textbrowser$Topwindow

	    toplevel $Topwindow
	    
	    wm title $Topwindow "$Title"
	    
	    set frame $Topwindow
	    label $frame.title -text "Font Options"
	    # Force this to be a little wider than default width
	    $frame.title configure -width 25
	    frame $frame.titleseperator \
		    -height 2 -borderwidth 1 -relief sunken
	    
	    # Font
	    set frame $Topwindow.font
	    frame $frame -relief sunken -borderwidth 2
	    
	    label $frame.title -text "Font Family"
	    frame $frame.titleseperator \
		    -height 2 -borderwidth 1 -relief sunken
	    
	    frame $frame.type
	    radiobutton $frame.type.arial -text "Arial" \
		    -variable [scope FontFamily] \
		    -value "0" \
		    -command "$this SetFontFamily 0"
	    
	    radiobutton $frame.type.courier -text "Courier" \
		    -variable [scope FontFamily] \
		    -value "1" \
		    -command "$this SetFontFamily 1"
	    
	    radiobutton $frame.type.times -text "Times" \
		    -variable [scope FontFamily] \
		    -value "2" \
		    -command "$this SetFontFamily 2"
	    
	    # Bold
	    checkbutton $frame.bold \
		    -text "Bold" -variable [scope Bold] \
		    -command "$this ExecuteSetBold"
	    
	    # Italic
	    checkbutton $frame.italic \
		    -text "Italic" -variable [scope Italic] \
		    -command "$this ExecuteSetItalic"
	    
	    
	    # Shadow
	    checkbutton $frame.shadow \
		    -text "Shadow" -variable [scope Shadow] \
		    -command "$this ExecuteSetShadow"
	    
	    # Color
	    set frame $Topwindow.color
	    frame $frame -relief sunken -borderwidth 2
	    
	    label $frame.label -text "Color"
	    scan $Color "%f %f %f" r g b
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]
	    button $frame.color -bg $color \
		    -command "$this SelectColor"
	    
	    
	    # ok
	    set frame $Topwindow.control
	    frame $frame -relief sunken -borderwidth 2
	    button $frame.ok \
		    -text "OK" \
		    -command "$this OK"
	    
	    #-----------------------------------------------------------------
	    # Pack widgets
	    #-----------------------------------------------------------------
	    
	    set frame $Topwindow
	    pack $frame.title
	    pack $frame.titleseperator -fill x
	    
	    set frame $Topwindow.font
	    pack $frame.title
	    pack $frame.titleseperator -fill x
	    
	    pack $frame.type -anchor w
	    pack $frame.type.arial -anchor w
	    pack $frame.type.courier -anchor w
	    pack $frame.type.times -anchor w
	    
	    pack $frame.bold -anchor w
	    pack $frame.italic -anchor w
	    pack $frame.shadow -anchor w
	    
	    pack $frame -fill x
	    
	    set frame $Topwindow.color
	    pack $frame.label -anchor w -side left
	    pack $frame.color -anchor w
	    pack $frame -fill x

	    set frame $Topwindow.control
	    pack $frame.ok
	    pack $frame -fill x

	} 
    }
}






