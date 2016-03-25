##
## File:        cvisReadSAMRAIVectorInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: User interface for file reader for Vizamrai data
##

class cvisReadSAMRAIVectorInterface {
    inherit cvisInterface

    variable AnimationFileNames ""

    variable FileName ""
    variable FileTypes {{SAMRAI {.vis}}}

    # Keep track of the modification time for the current file
    # and only load if newer file is found.
    variable FileTime 0

    variable VariableNames ""
    variable VariableName ""
    variable Variable 0

    variable VariableMenu ""

    variable Reader ""

    variable StopAnimation 0
    variable CurrentFileIndex 0

    variable AutoAnimation 0
    variable AutoAnimationTimeInterval 5

    variable Stride

    variable Cropping 0
    variable CroppingRegionPlanes

    variable CropBoxIsVisible 0
    variable CropBoxColor "1.0 1.0 1.0"

    variable AnimationFileListBox ""

    constructor {name} {
	cvisInterface::constructor $name
    } {
	foreach plane {x y z} {
	    set Stride($plane) 1
	    set CroppingRegionPlanes($plane,min) 0
	    set CroppingRegionPlanes($plane,max) 1
	}
    }

    method FileInfo {} {

	set reader [$Reader GetReader]
	puts "FileName [$reader GetFileName]"
	puts "Number of Patches [$reader GetNumberOfPatchBoundaries]"
	puts "Number of Levels [$reader GetNumberOfPatchBoundaries]"
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    set patchlevel [$Reader GetPatchLevel $patch]
	    
	    puts "$patchlevel $item"

	    puts "    dims [$item GetNumberOfCells]"
	    
	    set item [$collection GetNextItem]
	    incr patch
	}
    }

    method Save {file} {
	if {[string length $FileName]} {
	    puts $file "$this SetFileName \"$FileName\""
	}
	if {[string length $VariableName]} {
	    puts $file "$this SetCurrentVariableName \"$VariableName\""
	}
	puts $file "$this UpdateAnimationFileList"
	puts $file "$this SetStride [GetStride]"
	puts $file "$this SetCroppingRegionPlanes \
		$CroppingRegionPlanes(x,min) \
		$CroppingRegionPlanes(x,max) \
		$CroppingRegionPlanes(y,min) \
		$CroppingRegionPlanes(y,max) \
		$CroppingRegionPlanes(z,min) \
		$CroppingRegionPlanes(z,max)"
	puts $file "$this SetCropping $Cropping"
	puts $file "$this SetCropBoxVisibility $CropBoxIsVisible"
	puts $file "$this SetCropBoxColor $CropBoxColor"
	puts $file "$this SetAutoAnimationTimeInterval $AutoAnimationTimeInterval"
    }

    method SetAutoAnimation {autoanimation} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$AutoAnimation != $autoanimation} {
	    set AutoAnimation $autoanimation
	    ExecuteSetAutoAnimation
	}
    }

    method GetAutoAnimation {} {
	return $AutoAnimation
    }

    method ExecuteSetAutoAnimation {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$AutoAnimation} {
	    $this DoAutoAnimation
	}
    }

    method DoAutoAnimation {} {
	if {$AutoAnimation} {
	    $this UpdateAnimationFileList
	    $this GotoEndAnimation
	    after [expr $AutoAnimationTimeInterval * 1000] \
		    "$this DoAutoAnimation"
	}
    }

    method SetAutoAnimationTimeInterval {time} {
	set AutoAnimationTimeInterval $time
    }
    
    method GetAutoAnimationTimeInterval {} {
	return $AutoAnimationTimeInterval
    }

    method ExecuteSetStride {} {
	set tracervar [::itcl::local cvisTrace #auto]	
	$Reader SetStride $Stride(x) $Stride(y) $Stride(z)
    }

    method SetStride {x y z} {
	set Stride(x) $x
	set Stride(y) $y
	set Stride(z) $z
    }

    method GetStride {} {
	return "$Stride(x) $Stride(y) $Stride(z)"
    }

    method SetCurrentVariableName { variablename } {
	set tracervar [::itcl::local cvisTrace #auto]
	set VariableName $variablename
	UpdateCurrentVariableName
    }

    method UpdateCurrentVariableName {} {
	set tracervar [::itcl::local cvisTrace #auto]
	$Reader SetCurrentVariableName $VariableName
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}

	set NumVariables [$Reader GetNumberOfVariables]
	for {set i 0} { $i < $NumVariables } {incr i} {
	    set varname [$Reader GetVariableName $i]
	    if { ![string compare $VariableName $varname] } {
		set Variable $i
	    }
	}
    }

    method GetCollection { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetCollection]
    }

    method GetReader { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Reader
    }

    method SetInput { reader } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Reader $reader
    }

    method AddOutput {list} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Outputs $list
    }

    method GetFileName { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $FileName
    }

    # Force a new file load no matter what the timestamp
    # claims.
    method SetFileNameAlwaysLoad {filename} {
	set tracervar [::itcl::local cvisTrace #auto]
	set FileTime 0
	$this SetFileName $filename
    }

    method SetFileName { filename } {
	set tracervar [::itcl::local cvisTrace #auto]

	# Get the modification time of the new file
	set filetime [file mtime $filename]

	# If this is a new filename or if the file is newer 
	# then load
	if {[string compare $filename $FileName] || \
		($filetime > $FileTime)} {

	    set FileName $filename
	    set FileTime $filetime

	    set oldvariablename $VariableName

	    $Reader SetFileName $filename

	    # Remove the old entries from the variable chooser
	    $VariableMenu delete 0 last
	    set VariableNames ""
	    set VariableName  ""

	    # Need to set the variable names from the file
	    # Need to set default variable name
	    set NumVariables 0
	    for {set i 0} { $i < [$Reader GetNumberOfVariables] } {incr i} {
		set depth [$Reader GetDepth $i]

		if {$depth == 3 || $depth == 2} {
		    set varname [$Reader GetVariableName $i]
		
		    # If we find one with the previous name then default to
		    # what the user had set
		    if { ![string compare $oldvariablename $varname] } {
			set VariableName $varname
			set Variable $i
		    }
		    
		    lappend VariableNames $varname
		    $VariableMenu add radiobutton -label $varname \
			    -variable [scope VariableName]\
			    -command "$this UpdateCurrentVariableName"
		    
		    incr NumVariables 
		}
	    }
	    
	    # If we did not find a match with previous variable name set to
	    # default, first one in the list
	    if {![string length $VariableName]} {
		set VariableName [lindex $VariableNames 0]
		set Variable 0
	    }
	    
	    SetCurrentVariableName $VariableName

	    set AnimationFileNames \
		    [file rootname [file rootname $filename]].*.vis

	    UpdateAnimationFileList
	    
	    Update
	} 
    }

    variable FileBrowser ".null"
    method FileNameBrowser { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $FileBrowser]]} {
	} {
	    set FileBrowser [cvisFileBrowser #auto]
	    $FileBrowser SetTitle "Select FileName"
	    $FileBrowser SetFileTypes $FileTypes
	}
	
	set filename [$FileBrowser Prompt]

	# If not filename keep the existing name
	if {[string length $filename]} { 
	    SetFileName $filename
	    
	} 
    }

    method UpdateAnimationFileList { } {
	set tracervar [::itcl::local cvisTrace #auto]

	set filelist [glob $AnimationFileNames]
	set filelist [lsort $filelist]

	$AnimationFileListBox delete 0 end

	set CurrentFileIndex 0
	set index  0
	
	foreach file $filelist { 
	    $AnimationFileListBox insert end $file

	    if { ![string compare $file $FileName]} {
		set CurrentFileIndex $index
	    }
	    incr index
	}

	UpdateAnimation
    }

    method AnimationSelection { } {
	set tracervar [::itcl::local cvisTrace #auto]
	
	set CurrentFileIndex [$AnimationFileListBox curselection]

	UpdateAnimation
	$this SetFileName [$AnimationFileListBox get $CurrentFileIndex]
    }

    method UpdateAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]

	# Make the file the current selection in file listbox
	# and move into view
	$AnimationFileListBox selection clear 0 end
	$AnimationFileListBox selection set $CurrentFileIndex $CurrentFileIndex

	set visible_y [$AnimationFileListBox cget -visibleitems]

	set visible_y [string range $visible_y \
		[expr [string last x $visible_y] + 1] end]

	# Scroll so selection stays in window
	set top [expr $CurrentFileIndex - ($visible_y / 2)]
	if { $top < 0 } {
	    set top 0
	}
	
	$AnimationFileListBox yview $top
    }

    method StepForwardAnimation { } { 
	set tracervar [::itcl::local cvisTrace #auto]

	if { $CurrentFileIndex < [expr [$AnimationFileListBox size] - 1] } { 
	    incr CurrentFileIndex
	    UpdateAnimation
	    $this SetFileName [$AnimationFileListBox get $CurrentFileIndex]
	}
    }

    method StepBackwardAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]

	if { $CurrentFileIndex > 0 } { 
	    incr CurrentFileIndex -1
	    UpdateAnimation
	    $this SetFileName [$AnimationFileListBox get $CurrentFileIndex]
	}
    }

    method PlayAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]

	set t [time {
	    set oldpipelinestate [::cvis::queue::GetPipeline]
	    ::cvis::queue::SetPipeline 0
	    while { $CurrentFileIndex < [$AnimationFileListBox size] } {
		# Test if user wants to stop
		update
		if {$StopAnimation} {
		    set StopAnimation 0
		    break
		}
		
		# View next file
		UpdateAnimation
		$this SetFileNameAlwaysLoad [$AnimationFileListBox get $CurrentFileIndex]
		
		$RenderWindow RenderNow
		incr CurrentFileIndex
	    }
	    ::cvis::queue::SetPipeline $oldpipelinestate
	    
	    if { $CurrentFileIndex == [$AnimationFileListBox size] } {
		incr CurrentFileIndex -1
	    }
	}
	]
	    
        # puts "Rendering time $t"
    }

    method StopAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]
	set StopAnimation 1
    }

    method GotoStartAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]

	set CurrentFileIndex 0
	if {[$AnimationFileListBox size]} {
	    UpdateAnimation
	    $this SetFileName [$AnimationFileListBox get $CurrentFileIndex]
	}
    }

    method GotoEndAnimation { } {
	set tracervar [::itcl::local cvisTrace #auto]


	if {[$AnimationFileListBox size]} {
	    set CurrentFileIndex [expr [$AnimationFileListBox size] - 1]
	    UpdateAnimation
	    $this SetFileName [$AnimationFileListBox get $CurrentFileIndex]
	} {
	    set CurrentFileIndex 0
	}
    }

    method ExecuteSetCropping {} {
	$Reader SetCropping $Cropping
    }

    method SetCropping {cropping} {
	if {$Cropping != $cropping } {
	    set Cropping $cropping
	    ExecuteSetCropping
	}
    }

    method GetCropping {} {
	return $Cropping
    }

    method SetRenderWindow {renderwindow} {

	if {[string length [info commands $this.actor]]} {
	    $RenderWindow DeleteActor $this.actor
	} {
	    vtkOutlineSource $this.outlineSource

	    vtkPolyDataMapper $this.outlineMapper
	    $this.outlineMapper SetImmediateModeRendering \
		    $::cvis::options::LOWMEM
	    $this.outlineMapper SetInput [$this.outlineSource GetOutput]
	    
	    cvisActor $this.actor
	    $this.actor SetMapper $this.outlineMapper
	    $this.actor SetVisibility $CropBoxIsVisible
	    eval [$this.actor GetProperty] SetColor $CropBoxColor

	    $this.actor SetPickable 0
		
	}

	::cvisInterface::SetRenderWindow $renderwindow

	$RenderWindow AddActor $this.actor 

	if {[string length $Reader]} {
	    set Scaling [$Reader GetScaling]
	} {
	    set Scaling 1.0
	}

	eval $this.outlineSource SetBounds \
		[expr $CroppingRegionPlanes(x,min) * $Scaling] \
		[expr $CroppingRegionPlanes(x,max) * $Scaling] \
		[expr $CroppingRegionPlanes(y,min) * $Scaling] \
		[expr $CroppingRegionPlanes(y,max) * $Scaling] \
		[expr $CroppingRegionPlanes(z,min) * $Scaling] \
		[expr $CroppingRegionPlanes(z,max) * $Scaling]
	
	$RenderWindow Render
    }

    method ExecuteSetCroppingRegionPlanes {} {

	set Scaling [$Reader GetScaling]

  	$Reader SetCroppingRegionPlanes  \
		    [expr $CroppingRegionPlanes(x,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(x,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,max) * $Scaling]

	if {[string length [info commands $this.outlineSource]]} {
	    eval $this.outlineSource SetBounds \
		    [expr $CroppingRegionPlanes(x,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(x,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(y,max) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,min) * $Scaling] \
		    [expr $CroppingRegionPlanes(z,max) * $Scaling]
	    
	    $RenderWindow Render
	}
    }

    method SetCroppingRegionPlanes {min_x max_x min_y max_y min_z max_z} {

	set CroppingRegionPlanes(x,min) $min_x
	set CroppingRegionPlanes(x,max) $max_x

	set CroppingRegionPlanes(y,min) $min_y
	set CroppingRegionPlanes(y,max) $max_y

	set CroppingRegionPlanes(z,min) $min_z
	set CroppingRegionPlanes(z,max) $max_z

	ExecuteSetCroppingRegionPlanes
    }

    method SetCropBoxColor {r g b} {
	set CropBoxColor "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	
	set frame $NotebookPage.cropping
	set cs [$frame childsite]

	$cs.control.color.color configure -bg $color

	if {[string length $RenderWindow]} {
	    if {[string length [info commands $this.actor]]} {
		eval [$this.actor GetProperty] SetColor $CropBoxColor
	    }
	    $RenderWindow Render
	}
    }

    method SelectCropBoxColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $CropBoxColor
	
	set value [color Prompt]
	
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b
	    SetCropBoxColor $r $g $b
	}	

	delete object color 
    }

    method ExecuteSetCropBoxVisibility {} {
	if {[string length $RenderWindow]} {
	    $this.actor SetVisibility $CropBoxIsVisible
	    $RenderWindow Render
	}
    }

    method SetCropBoxVisibility {visible} {
	if {$visible != $CropBoxIsVisible} {
	    set CropBoxIsVisible $visible
	    ExecuteSetCropBoxVisibility
	}
    }
    
    method GetCropBoxVisibility {} {
	return $CropBoxIsVisible
    }

    method Update {} {
	set frame $NotebookPage.cropping
	set cs [$frame childsite]

	set Scaling [$Reader GetScaling]

	scan [[$Reader GetReader] GetBounds] "%f %f %f %f %f %f" \
		Bounds(x,min) Bounds(x,max) \
		Bounds(y,min) Bounds(y,max) \
		Bounds(z,min) Bounds(z,max)

	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    set Bounds($plane,min) [expr $Bounds($plane,min) / \
		    $Scaling]

	    set value [format "%1.3e" $Bounds($plane,min)]
	    $plane_frame.mintitle configure -text "Min $plane ($value)"
	    set Bounds($plane,max) [expr $Bounds($plane,max) / \
		    $Scaling]
	    set value [format "%1.3e" $Bounds($plane,max)]
	    $plane_frame.maxtitle configure -text "Max $plane ($value)"
	}
	
    }

    method SetNotebookPage {notebookpage} {
	set tracervar [::itcl::local cvisTrace #auto]
	set NotebookPage $notebookpage

	global cvis_library

	#------------------------------------------------------------------
	# Title
	#------------------------------------------------------------------
	set frame $NotebookPage.titleframe
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Input File"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	#------------------------------------------------------------------
	# Filename
	#------------------------------------------------------------------

	set frame $NotebookPage.animationframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "File"
	set cs [$frame childsite]

	frame $cs.filenameframe

	button $cs.filenameframe.selectfile -text "FileName" \
		-command  "$this FileNameBrowser" 

	entry $cs.filenameframe.filename \
		-textvariable [scope FileName]

	#------------------------------------------------------------------
	# Animaion
	#------------------------------------------------------------------

	iwidgets::scrolledlistbox $cs.animationfilelist -selectmode single \
		-visibleitems 20x5 \
		-selectioncommand "$this AnimationSelection"

	set AnimationFileListBox $cs.animationfilelist

	frame $cs.animationcontrols -relief sunken

	button $cs.animationcontrols.update -text "Update Files" \
		-image [image create photo -file \
		[file join $cvis_library images makelist.gif]] \
		-command  "$this UpdateAnimationFileList" 

	button $cs.animationcontrols.stop \
		-image [image create photo -file \
		[file join $cvis_library images stop.gif]] \
		-command  "$this StopAnimation"

	button $cs.animationcontrols.play  \
		-image [image create photo -file \
		[file join $cvis_library images play.gif]] \
		-command  "$this PlayAnimation"

	button $cs.animationcontrols.first \
		-image [image create photo -file \
		[file join $cvis_library images start.gif]] \
		-command  "$this GotoStartAnimation"

	button $cs.animationcontrols.last \
		-image [image create photo -file \
		[file join $cvis_library images end.gif]] \
		-command  "$this GotoEndAnimation"

	button $cs.animationcontrols.forward \
		-image [image create photo -file \
		[file join $cvis_library images next.gif]] \
		-command  "$this StepForwardAnimation"

	button $cs.animationcontrols.backward  \
		-image [image create photo -file \
		[file join $cvis_library images previous.gif]] \
		-command  "$this StepBackwardAnimation"

	frame $cs.autoanimation

	checkbutton $cs.autoanimation.onoff \
		-text "Automatically check for new files" -variable [scope AutoAnimation] \
		-command "$this ExecuteSetAutoAnimation"

	label $cs.autoanimation.title \
		-text "Check interval (seconds)" 
	entry $cs.autoanimation.value \
		-textvariable [scope AutoAnimationTimeInterval]

	label $cs.filetitle -text "Animation Filenames (use wildcard '*')" 
	entry $cs.filenames -textvariable [scope AnimationFileNames]
	
	#------------------------------------------------------------------
	# Variable 
	#------------------------------------------------------------------

	set frame $NotebookPage.variableframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Variable"
	set cs [$frame childsite]

	set VariableMenu [tk_optionMenu $cs.variablename \
		[scope VariableName] ""]
	$cs.variablename configure -width 10 

	#-----------------------------------------------------------------
	# Stride
	#-----------------------------------------------------------------

	set frame $NotebookPage.strideframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Downsize"
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    frame $frame
	    label $frame.slicename \
		    -text "$plane"

	    scale $frame.scale \
		    -from 1 \
		    -to 20 \
		    -variable [scope Stride($plane)] \
		    -orient horizontal \
		    -highlightthickness 0
	    
	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $frame.scale <ButtonRelease> \
		    "$this ExecuteSetStride"
	}

	#-----------------------------------------------------------------
	# Cropping
	#-----------------------------------------------------------------

	set frame $NotebookPage.cropping
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Cropping Region"
	set cs [$frame childsite]

	frame $cs.control

	checkbutton $cs.control.cropping \
		-text "Crop Data" \
		-variable [scope Cropping] \
		-command "$this ExecuteSetCropping"

	if {[string length $Reader]} {
	    scan [$Reader GetIndexBounds] "%f %f %f %f %f %f" \
		    Bounds(x,min) Bounds(x,max) \
		    Bounds(y,min) Bounds(y,max) \
		    Bounds(z,min) Bounds(z,max)
	} {
	    foreach plane {x y z} {
		set Bounds($plane,min) 0
		set Bounds($plane,max) 1
	    }
	}

	checkbutton $cs.control.visible \
		-text "CropBox Visible" -variable [scope CropBoxIsVisible] \
		-command "$this ExecuteSetCropBoxVisibility"

	frame $cs.control.color

	label $cs.control.color.label -text CropBoxColor
	
	scan $CropBoxColor "%f %f %f" r g b
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	button $cs.control.color.color -bg $color \
		-command "$this SelectCropBoxColor"

	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    frame $plane_frame
	    
	    set value [format "%1.3e" $Bounds($plane,min)]
	    label $plane_frame.mintitle -text "Min $plane ($value)"
	    entry $plane_frame.min \
		-textvariable [scope CroppingRegionPlanes($plane,min)] 
	    bind $plane_frame.min <Return> "$this ExecuteSetCroppingRegionPlanes"

	    set value [format "%1.3e" $Bounds($plane,max)]
	    label $plane_frame.maxtitle -text "Max $plane ($value)"
	    entry $plane_frame.max \
		    -textvariable [scope CroppingRegionPlanes($plane,max)] 
	    bind $plane_frame.max <Return> "$this ExecuteSetCroppingRegionPlanes"
	    
	}

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.titleframe
	pack $frame.title
	pack $frame.titleseperator -fill x

	set frame $NotebookPage.animationframe
	set cs [$frame childsite]

	pack $cs.filenameframe.selectfile -side left -padx 4 -pady 4
	pack $cs.filenameframe.filename -side left \
		-fill x -expand yes
	pack $cs.filenameframe -fill x -pady 4

	#------------------------------------------------------------------
	pack $cs.filetitle -anchor w
	pack $cs.filenames -padx 4 -pady 4 -fill x -expand yes
	pack $cs.animationfilelist -padx 4 -pady 4 -fill x -expand yes 
	pack $cs.animationcontrols
	pack $cs.animationcontrols.update \
		-side left -anchor w 
	pack $cs.animationcontrols.first \
		-side left -anchor w 
	pack $cs.animationcontrols.backward \
		-side left -anchor w
	pack $cs.animationcontrols.play \
		-side left -anchor w 
	pack $cs.animationcontrols.forward \
		-side left -anchor w
	pack $cs.animationcontrols.last \
		-side left -anchor w
	pack $cs.animationcontrols.stop \
		-side left -anchor w
	pack $cs.autoanimation.onoff -anchor w -pady 4
	pack $cs.autoanimation.title -anchor w -side left
	pack $cs.autoanimation.value -padx 4 -pady 4 \
		-side left
	pack $cs.autoanimation -anchor w 
	pack $cs -fill x
	pack $frame -fill x

	#------------------------------------------------------------------
	set frame $NotebookPage.variableframe

	set cs [$frame childsite]
	pack $cs.variablename -side left -anchor w -fill x -expand yes -pady 6
	pack $cs -fill x 
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.strideframe
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    pack $frame.slicename -side left
	    pack $frame.scale -expand yes -fill x 
	    pack $frame -fill x
	}
	pack $NotebookPage.strideframe -fill x 
	
	#-----------------------------------------------------------------
	set frame $NotebookPage.cropping
	set cs [$frame childsite]

	pack $cs.control.cropping -side left
	pack $cs.control.visible -side left
	pack $cs.control.color.label -side left 
	pack $cs.control.color.color -side left
	pack $cs.control.color -side left 
	pack $cs.control -fill x

	foreach plane {x y z} {
	    set plane_frame $cs.$plane
	    pack $plane_frame.mintitle -side left
	    pack $plane_frame.min -side left  -pady 4
	    pack $plane_frame.maxtitle -side left
	    pack $plane_frame.max -side left  -pady 4
	    pack $plane_frame -fill x
	}
	
	pack $cs -fill x
  	pack $frame -fill x

    }	
}



