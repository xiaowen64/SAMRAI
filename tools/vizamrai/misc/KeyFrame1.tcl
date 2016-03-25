source ../misc/KeyFrame.tcl

proc KeyRender {} {
    global KeyCurrentFrameNumber
    global KeyNumberOfFrames
    global KeyNumberOfInputFiles
    global KeyInputFileName

    set input_file_number [expr round((double($KeyCurrentFrameNumber) / \
	    ($KeyNumberOfFrames-1)) * ($KeyNumberOfInputFiles-1))] 

    puts $input_file_number
# To Read in a file at each time step
#    readSamrai SetFileName [format $KeyInputFileName $input_file_number]
#    readSamrai Update

    [renderWindow GetRenderer] ResetCameraClippingRange
    puts "Rendering $KeyCurrentFrameNumber"
    renderWindow Render
    set filename  [format "seq1.%05d.ppm" $KeyCurrentFrameNumber]
    renderWindow SetFileName $filename
    renderWindow SaveImageAsPPM
}


