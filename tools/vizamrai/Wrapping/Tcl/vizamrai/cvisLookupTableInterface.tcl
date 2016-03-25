##
## File:        cvisLookupTableInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: LookupTable interface
##


class cvisLookupTableInterface {
    inherit cvisInterface

    variable Input ""

    variable Title " "

    variable Visible 1

    variable Position
    variable Size
    variable Orientation 1

    variable TitleSize 0.5
    variable LabelSize 0.5

    variable LabelFormat "%1.2e"

    variable NumberOfLabels 5

    variable Type 0

    variable AutoRangeScaling 1
    
    variable ScalarRangeMin 0
    variable ScalarRangeMax 1

    variable ColorFileName ""
    variable ColorFileTypes {{ColorMaps {.vcmap}}}

    variable AlphaFileName ""
    variable AlphaFileTypes {{ColorMaps {.vamap}}}

    variable K 5.0
    variable K0 -0.001
    variable K1 0.5

    constructor {name} {
	cvisInterface::constructor $name
	cvisTextBrowser $this.textbrowser 

	set Position(x) 0.01
	set Position(y) 0.05

	set Size(x) 0.15
	set Size(y) 1.0
    } {
    }
    
    method Save {file} {
	puts $file "$this SetScalarRange $ScalarRangeMin $ScalarRangeMax"
	puts $file "$this SetType $Type"
	puts $file "$this SetK $K"
	puts $file "$this SetK0 $K0"
	puts $file "$this SetK1 $K1"
	puts $file "$this SetTitle {$Title}"
	puts $file "$this SetTitleSize $TitleSize"
	puts $file "$this SetLabelSize $LabelSize"
	puts $file "$this SetLabelFormat \"[GetLabelFormat]\""
	puts $file "$this SetNumberOfLabels $NumberOfLabels"
	puts $file "$this SetVisibility $Visible"
	puts $file "$this SetAutoRangeScaling $AutoRangeScaling"
	puts $file "$this SetScalarRange $ScalarRangeMin $ScalarRangeMax"
	if {[string length $ColorFileName]} {
	    puts $file "$this SetColorFileName \"$ColorFileName\""
	}
	if {[string length $AlphaFileName]} {
	    puts $file "$this SetAlphaFileName \"$AlphaFileName\""
	}
	puts $file "$this SetPosition $Position(x) $Position(y)"
	puts $file "$this SetSize     $Size(x) $Size(y)"
	puts $file "$this SetOrientation $Orientation"
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
	$this.textbrowser SetRenderWindow $renderwindow
    }


    variable ColorFileBrowser ".null"
    method SetColorFileName {colorfilename} {
	set ColorFileName $colorfilename
	foreach item $Input {
	    $item SetColorFileName $ColorFileName
	}

    }

    method LoadColorMap {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $ColorFileBrowser]]} {
	    $ColorFileBrowser SetOpenPrompt
	} {
	    set ColorFileBrowser [cvisFileBrowser #auto]
	    $ColorFileBrowser SetTitle "Select ColorMap File"
	    $ColorFileBrowser SetFileTypes $ColorFileTypes
	    $ColorFileBrowser SetOpenPrompt
	}

	set filename [$ColorFileBrowser Prompt]
	SetColorFileName $filename

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    variable AlphaFileBrowser ".null"
    method SetAlphaFileName {colorfilename} {
	set AlphaFileName $colorfilename
	foreach item $Input {
	    $item SetAlphaFileName $AlphaFileName
	}

    }

    method LoadAlphaMap {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $AlphaFileBrowser]]} {
	    $AlphaFileBrowser SetOpenPrompt
	} {
	    set AlphaFileBrowser [cvisFileBrowser #auto]
	    $AlphaFileBrowser SetTitle "Select AlphaMap File"
	    $AlphaFileBrowser SetFileTypes $AlphaFileTypes
	    $AlphaFileBrowser SetOpenPrompt
	}

	set filename [$AlphaFileBrowser Prompt]
	SetAlphaFileName $filename

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetScalarRange {min max} {
	set doit 0
	if { $ScalarRangeMin != $min } {
	    set ScalarRangeMin $min
	    set doit 1
	}

	if { $ScalarRangeMax != $max } {
	    set ScalarRangeMax $max
	    set doit 1
	}

	if {$doit} {
	    ExecuteSetScalarRange
	}
    }

    method GetScalarRange {} {
	return "$ScalarRangeMin $ScalarRangeMax"
    }




    method ExecuteSetScalarRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach item $Input {
	    $item SetTableRange $ScalarRangeMin $ScalarRangeMax
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetAutoRangeScaling {autorangescaling} {
	if {$AutoRangeScaling != $autorangescaling} {
	    set AutoRangeScaling $autorangescaling
	    ExecuteSetAutoRangeScaling
	}
    }

    method ExecuteSetAutoRangeScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach item $Input {
	    $item SetAutoTableRangeScaling $AutoRangeScaling
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetTitle {title} {
	set Title $title
	ExecuteSetTitle
    }

    method GetTitle {} {
	return $Title
    }
    

    method ExecuteSetTitle {} {
	foreach item $Input {
	    [$item GetColorBar] SetTitle $Title
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetVisibility {visible} {
	if {$Visible != $visible} {
	    set Visible $visible
	    ExecuteSetVisibility
	}
    }
    
    method GetVisibibility {} {
	return $Visible
    }

    method ExecuteSetVisibility {} {
	foreach item $Input {
	    $item SetVisibility $Visible
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetType {type} {
	if {$type != $Type} {
	    set Type $type
	    SetTypeChange
	}
    }
    
    method GetType {} {
	return $Type
    }
    
    method SetTypeChange {} {
	foreach item $Input {
	    $item SetType $Type
	}
	
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetK {k} {
	if {$k != $K} {
	    set K $k
	    ExecuteSetScaleOptions
	}
    }
    method GetK {} {
	return $K
    }

    method GetK0 {} {
	return $K0
    }
    method SetK0 {k} {
	if {$k != $K0} {
	    set K0 $k
	    ExecuteSetScaleOptions
	}
    }
    
    method GetK1 {} {
	return $K1
    }
    method SetK1 {k} {
	if {$k != $K1} {
	    set K1 $k
	    ExecuteSetScaleOptions
	}
    }



    method ExecuteSetScaleOptions {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach item $Input {
	    $item SetK $K
	    $item SetK0 $K0
	    $item SetK1 $K1
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method AddInput {input} {
	lappend Input $input

	$input SetInterface $this

	$this.textbrowser SetInput [$input GetColorBar]
	
	Update
    }

    method Update {} {
    }

    method LegendOptions {} {
	# What about multiple items?
	if {![string length [info commands $this.textbrowser]]} {
	    set item [lindex $Input 0]
	    $this.textbrowser SetTitle "Colorbar Legend Options"
	}
	$this.textbrowser Prompt 
    }

    method SetSize {x y} {
	set Size(x) $x
	set Size(y) $y
	ExecuteSetSize
    }
    
    method GetSize {} {
	return $Size(x) $Size(y)
    }

    method ExecuteSetSize {} {
	foreach item $Input {
	    [$item GetColorBar] SetWidth $Size(x)
	    [$item GetColorBar] SetHeight $Size(y)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}

    }

    method SetPosition {x y} {
	set Position(x) $x
	set Position(y) $y
	ExecuteSetPosition
    }
    
    method GetPosition {} {
	return $Position(x) $Position(y)
    }

    
    method ExecuteSetPosition {} {
	foreach item $Input {
	    [$item GetColorBar] SetPosition $Position(x) $Position(y)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetOrientation {orientation} {
	set Orientation $orientation
	ExecuteSetOrientation
    }
    
    method GetOrientation {} {
	return $Orientation
    }

    
    method ExecuteSetOrientation {} {
	foreach item $Input {
	    [$item GetColorBar] SetOrientation $Orientation
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetTitleSize {size} {
	set TitleSize $size
	ExecuteSetTitleSize
    }

    method GetTitleSize {} {
	return $TitleSize
    }
    

    method ExecuteSetTitleSize {} {
	foreach item $Input {
	    [$item GetColorBar] SetTitleSize $TitleSize
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetLabelSize {size} {
	set LabelSize $size
	ExecuteSetLabelSize
    }

    method GetLabelSize {} {
	return $LabelSize
    }
    

    method ExecuteSetLabelSize {} {
	foreach item $Input {
	    [$item GetColorBar] SetLabelSize $LabelSize
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method ExecuteSetLabelFormat {} {
	foreach item $Input {
	    [$item GetColorBar] SetLabelFormat $LabelFormat
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }
    
    method SetLabelFormat {format} {
	set LabelFormat $format
	ExecuteSetLabelFormat
    }

    method GetLabelFormat {} {
	return $LabelFormat
    }

    method SetNumberOfLabels {number} {
	set NumberOfLabels $number
	ExecuteSetNumberOfLabels
    }

    method GetNumberOfLabels {} {
	return $NumberOfLabels
    }
    

    method ExecuteSetNumberOfLabels {} {
	foreach item $Input {
	    [$item GetColorBar] SetNumberOfLabels $NumberOfLabels
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }


    method SetNotebookPage {notebook} {
	set NotebookPage $notebook
	#-----------------------------------------------------------------
	# Visibility
	#-----------------------------------------------------------------
	set frame $NotebookPage.visible
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Color Bar Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Visible" -variable [scope Visible] \
		-command "$this ExecuteSetVisibility"

	#-----------------------------------------------------------------
	# Text
	#-----------------------------------------------------------------
	set frame $NotebookPage.legend
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Legend"
	set cs [$frame childsite]

	set frame $cs.label
	frame $frame

	label $frame.title -text "Title"
	entry $frame.label \
		-textvariable [scope Title] 
	bind $frame.label <Return> "$this ExecuteSetTitle"

	set frame $cs.format
	frame $frame 

	label $frame.labelformattext \
		-text "Label Format"
	
	entry $frame.labelformat \
		-textvariable [scope LabelFormat]

	bind $frame.labelformat <Return> "$this ExecuteSetLabelFormat"

	button $cs.textbutton -text "Legend Font Selector" \
		-command  "$this LegendOptions"

	set frame $cs.numberoflabels
	frame $frame
	label $frame.label \
		-text "Number of Labels"

	scale $frame.size \
		-from 2 -to 20 \
		-variable [scope NumberOfLabels] \
		-orient horizontal \
		-highlightthickness 0
	    
	# Need to add button binding here so when slider changes
	# we update the plane
	bind $frame.size <ButtonRelease> \
		"$this ExecuteSetNumberOfLabels"

	set frame $cs.titlesize
	frame $frame
	label $frame.label \
		-text "Title Size"

	scale $frame.size \
		-from 0 -to 1.0 -resolution 0.01 \
		-variable [scope TitleSize] \
		-orient horizontal \
		-highlightthickness 0
	    
	# Need to add button binding here so when slider changes
	# we update the plane
	bind $frame.size <ButtonRelease> \
		"$this ExecuteSetTitleSize"

	set frame $cs.labelsize
	frame $frame
	label $frame.label \
		-text "Label Size"

	scale $frame.size \
		-from 0 -to 1.0 -resolution 0.01 \
		-variable [scope LabelSize] \
		-orient horizontal \
		-highlightthickness 0
	    
	# Need to add button binding here so when slider changes
	# we update the plane
	bind $frame.size <ButtonRelease> \
		"$this ExecuteSetLabelSize"


	#-----------------------------------------------------------------
	# Orientation
	#-----------------------------------------------------------------
	set frame $NotebookPage.orientation
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Orientation"
	set cs [$frame childsite]

	frame $cs.type

	radiobutton $cs.type.horizontal -text "Horizontal" \
		-variable [scope Orientation] \
		-value "0" \
		-command "$this ExecuteSetOrientation"

	radiobutton $cs.type.verticle -text "Verticle" \
		-variable [scope Orientation] \
		-value "1" \
		-command "$this ExecuteSetOrientation"

	#-----------------------------------------------------------------
	# Position/Size 
	#-----------------------------------------------------------------
	set frame $NotebookPage.colorbar
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Color Bar Position and Size"
	set cs [$frame childsite]

	foreach plane {x y} {
	    set frame $cs.$plane
	    frame $frame
	    label $frame.label \
		    -text "$plane"

	    scale $frame.position \
		    -from 0 -to 1.0 -resolution 0.01 \
		    -variable [scope Position($plane)] \
		    -orient horizontal \
		    -highlightthickness 0
	    
	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $frame.position <ButtonRelease> \
		    "$this ExecuteSetPosition"

	    scale $frame.size \
		    -from 0 -to 1.0 -resolution 0.01 \
		    -variable [scope Size($plane)] \
		    -orient horizontal \
		    -highlightthickness 0
	    
	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $frame.size <ButtonRelease> \
		    "$this ExecuteSetSize"

	}

	#-----------------------------------------------------------------
	# Type
	#-----------------------------------------------------------------
	set frame $NotebookPage.type
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Color Mapping Type"
	set cs [$frame childsite]

	frame $cs.type

	radiobutton $cs.type.linear -text "Linear" \
		-variable [scope Type] \
		-value "0" \
		-command "$this SetTypeChange"

	radiobutton $cs.type.log -text "Log" \
		-variable [scope Type] \
		-value "1" \
		-command "$this SetTypeChange"

	radiobutton $cs.type.schlieren -text "Schlieren" \
		-variable [scope Type] \
		-value "2" \
		-command "$this SetTypeChange"

	frame $cs.options

	label $cs.options.ktitle -text "k"
	entry $cs.options.k \
		-textvariable [scope K] 
	bind $cs.options.k <Return> "$this ExecuteSetScaleOptions"

	label $cs.options.k0title -text "k0"
	entry $cs.options.k0 \
		-textvariable [scope K0] 
	bind $cs.options.k0 <Return> "$this ExecuteSetScaleOptions"

	label $cs.options.k1title -text "k1"
	entry $cs.options.k1 \
		-textvariable [scope K1] 
	bind $cs.options.k1 <Return> "$this ExecuteSetScaleOptions"

	#------------------------------------------------------------------
	# Color Scaling
	#------------------------------------------------------------------
	set frame $NotebookPage.scaleframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Color Map Range"
	set cs [$frame childsite]

	checkbutton $cs.autoscale \
		-text "Use Min/Max from dataset" \
		-variable [scope AutoRangeScaling] \
		-command "$this ExecuteSetAutoRangeScaling"
	
	set frame $cs.minmax
	frame $frame

	label $frame.mintitle -text "Min"
	entry $frame.min \
		-textvariable [scope ScalarRangeMin] 
	bind $frame.min <Return> "$this ExecuteSetScalarRange"

	label $frame.maxtitle -text "Max"
	entry $frame.max \
		-textvariable [scope ScalarRangeMax] 
	bind $frame.max <Return> "$this ExecuteSetScalarRange"

	set frame $cs.cmap
	frame $frame
	button $frame.colorfilename_button -text "Load ColorMap" \
		-command  "$this LoadColorMap" 

	entry $frame.colorfilename \
		-textvariable [scope ColorFileName]

	set frame $cs.amap
	frame $frame
	button $frame.alphafilename_button -text "Load AlphaMap" \
		-command  "$this LoadAlphaMap" 

	entry $frame.alphafilename \
		-textvariable [scope AlphaFileName]

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------

	#-----------------------------------------------------------------
	set frame $NotebookPage.visible
	set cs [$frame childsite]
	pack $cs.visible -anchor w -pady 4
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.legend
	pack $frame -fill x
	set cs [$frame childsite]
	set frame $cs.label

	pack $frame.title -anchor w -side left
	pack $frame.label -fill x -pady 8
	pack $frame -fill x 

	set frame $cs.format
	pack $frame.labelformattext -anchor w -side left
	pack $frame.labelformat -anchor w 
	pack $frame -fill x

	pack $cs.textbutton -fill x

	set frame $cs.numberoflabels
	pack $frame.label -side left
	pack $frame.size -side left -expand yes -fill x
	pack $frame -fill x

	set frame $cs.titlesize
	pack $frame.label -side left
	pack $frame.size -side left -expand yes -fill x
	pack $frame -fill x

	set frame $cs.labelsize
	pack $frame.label -side left
	pack $frame.size -side left -expand yes -fill x
	pack $frame -fill x

	pack $cs -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.orientation
	set cs [$frame childsite]
	pack $cs.type -expand yes -fill x 
	pack $cs.type.horizontal -side left -pady 4
	pack $cs.type.verticle -side left
  	pack $cs -fill x
  	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.colorbar
	pack $frame -fill x
	set cs [$frame childsite]

	foreach plane {x y} {
	    set frame $cs.$plane
	    pack $frame.label -side left
	    pack $frame.position -pady 2 -side left -expand yes -fill x 
	    pack $frame.size -side left -expand yes -fill x 

	    pack $frame -fill x
	}

	pack $cs -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.type
	set cs [$frame childsite]

	pack $cs.type -expand yes -fill x 
	pack $cs.type.linear -side left -pady 4
	pack $cs.type.log -side left
	pack $cs.type.schlieren -side left

	pack $cs.options.ktitle -side left
	pack $cs.options.k -side left -padx 4 -pady 4

	pack $cs.options.k0title -side left
	pack $cs.options.k0 -side left -padx 4 -pady 4

	pack $cs.options.k1title -side left
	pack $cs.options.k1 -side left -padx 4 -pady 4
	
	pack $cs.options -fill x

	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.scaleframe
	set cs [$frame childsite]
	pack $frame -fill x

	pack $cs.autoscale  -side top -padx 4 -pady 4 -anchor w

	set frame $cs.minmax
	pack $frame.mintitle -side left 
	pack $frame.min -side left -padx 4 -pady 4
	pack $frame.maxtitle -side left
	pack $frame.max -side left -padx 4 -pady 4
	pack $frame -fill x

	set frame $cs.cmap
	pack $frame.colorfilename_button -side left -anchor w -padx 4 -pady 4
	pack $frame.colorfilename -side left \
		-fill x -expand yes -pady 4
	pack $frame -fill x

	set frame $cs.amap
	pack $frame.alphafilename_button -side left -anchor w -padx 4 -pady 4
	pack $frame.alphafilename -side left \
		-fill x -expand yes -pady 4
	pack $frame -fill x

	pack $cs -fill x

    }   
}


