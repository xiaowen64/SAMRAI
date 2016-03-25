# KeyFrame.tcl - Keyframes for vtk
# 
# These procs provide a simple (yet powerful) keyframing capability
# for vtk tcl applications
# A keyframe contains a time sorted ordering of methods
# and arguments for an object. Key locations, colors, and
# other parameters are stored and can be interpolated
# at intermediate points.
# The protocol for keyframing is illustrated in the following
# example to keyframe the position of the camera camera1:
# the renderer ren1:
# 1) KeyNew CameraPosition camera1 SetPosition
#    create a keyframe for camera1 that will use the
#    SetPosition method of the camera
# 2) KeyAdd CameraPosition camera1 [camera1 GetPosition]
#    adds the current camera position as a key frame
# 3) Repeat 2, changing the camera position for each step
# 4) KeyRun CameraPosition 30
#    runs key frame for 30 steps. Steps that lie
#    between keyframes are interpolated with Cubic Splines.
# After each step is interpolated, the proc KeyRender
# is invoked. This proc can be redefined by the use to
# do something more sophisticated. The default proc does
# a: renWin Render
#
############################################
# Create a new keyframe for object and method
proc KeyNew {keyFrame object method} {
    upvar $keyFrame a
    set a(object) $object
    set a(method) $method
    set a(counter) 0
    set a(debug) 0
}

#
# Reset the keyframe count to 0
proc KeyReset {keyFrame} {
  upvar $keyFrame a
  set a(counter) 0
  if { $a(debug) == 1} {puts "Resetting $keyFrame"}
}

#
# Add a new keyframe with contents value
proc KeyAdd {keyFrame value} {
  upvar $keyFrame a
  set a($a(counter)) $value
  incr a(counter)
}

#
# Run a keyframe for "frames" frames
proc KeyRun {keyFrame frames} {
  upvar $keyFrame a
  set object $a(object)
  set method $a(method)
  set depends [llength $a(0)]
  #
  # create splines if they do not exist
  #
  for {set j 0} {$j < $depends} {incr j} {
    set spline $object${method}Spline$j
    if {[info commands $spline] == ""} {
      vtkKochanekSpline $spline
    }
    $spline RemoveAllPoints
  }
  #
  # add points to the splines
  #
  for {set i 0} {$i < $a(counter)} {incr i} {
      for {set j 0} {$j < $depends} {incr j} {
	$object${method}Spline$j AddPoint $i [lindex $a($i) $j]
      }
  }
  #
  # evaluate splines at key frames
  #
  global KeyCurrentFrameNumber 
  global KeyNumberOfFrames 
  set KeyNumberOfFrames $frames
  for {set i 0} {$i < $frames} {incr i} {
      set KeyCurrentFrameNumber $i
      set t [expr ( $a(counter) - 1.0 ) / ( $frames - 1) * $i]
      KeyGoto a $t
  }
}

# Goto keyframe #
proc KeyGoto {keyFrame t} {
  upvar $keyFrame a
  set object $a(object)
  set method $a(method)
  set depends [llength $a(0)]
  set keyCommand "$object $method "
  for {set j 0} {$j < $depends} {incr j} {
    lappend keyCommand "[$object${method}Spline$j Evaluate $t]"
  } 
    if {$a(debug) == 1} {puts "$keyCommand"}
  eval $keyCommand
  KeyRender
}

# Save a keyframe in "file"
proc KeySaveInFile {keyFrame fileName} {
  upvar $keyFrame a
  set filePtr [open $fileName w]
  set frames $a(counter)
  puts $filePtr "KeyNew $keyFrame $a(object) \"$a(method)\""
    for {set i 0} {$i < $frames} {incr i} {
    puts $filePtr "KeyAdd $keyFrame \"$a($i)\""
  }
  close $filePtr
}

proc GetInterpolate { camera } {
    set values ""
    eval lappend values [$camera GetPosition]
    eval lappend values [$camera GetFocalPoint]
    eval lappend values [$camera GetViewUp]
}

proc SetInterpolate { Camera px py pz fx fy fz vx vy vz} {
    $Camera SetPosition $px $py $pz
    $Camera SetFocalPoint $fx $fy $fz
    $Camera SetViewUp $vx $vy $vz
}

# Called after keyframe is executed
set KeyCurrentFrameNumber 0
set KeyNumberOfFrames 4
set KeyNumberOfInputFiles 51
set KeyInputFileName "eng3d.%04d.vis"


#  proc KeyRender {} {
#      global KeyCurrentFrameNumber
#      global KeyNumberOfFrames
#      set number_of_files 14
#      set filename "data/room2damr.%04d.vis"

#      set input_file_number [expr round((double($KeyCurrentFrameNumber) / \
#  	    ($KeyNumberOfFrames-1)) * ($number_of_files-1))] 

#      puts $input_file_number
#  # To Read in a file at each time step
#  #    readSamrai SetFileName [format $filename $input_file_number]
#  #    readSamrai Update

#      [renderWindow GetRenderer] ResetCameraClippingRange
#      puts "Rendering $KeyCurrentFrameNumber"
#      renderWindow Render
#      set filename  [format "image.%05d.ppm" $KeyCurrentFrameNumber]
#      renderWindow SetFileName $filename
#      renderWindow SaveImageAsPPM
#      #    exec gzip $filename &
#  }

#  # Called after keyframe is executed
#  proc KeyRender {} {
#      [ren1 GetLights] InitTraversal
#      set light [[ren1 GetLights] GetNextItem]
#      set camera [ren1 GetActiveCamera]
#      $camera SetViewUp 0 -1 0
#      $camera SetClippingRange 10 5000
#      eval $light SetPosition [$camera GetPosition]
#      eval $light SetFocalPoint [$camera GetFocalPoint]
#      renWin Render
#  }

#  # Called after keyframe is executed
#  proc KeyRender {} {
#    renWin Render
#  }


