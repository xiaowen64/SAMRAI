set VTK_INSTALL_PATH $env(VTK_INSTALL_PATH)

proc Execute {command} {
    global VTK_INSTALL_PATH
    global argv 
    exec $VTK_INSTALL_PATH/bin/$command &
}

button .threeslice \
	-text "Three Slice" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images threeslice.gif]] \
	-command {Execute vizamrai.3slice}

button .singleslice \
	-text "Single Slice" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images singleslice.gif]] \
	-command {Execute vizamrai.1slice}

button .carpet \
	-text "Surface" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images surface.gif]] \
	-command {Execute vizamrai.carpet}

button .iso \
	-text "IsoSurface" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images isosurface.gif]] \
	-command {Execute vizamrai.iso}

button .vector \
	-text "Volume" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images vector.gif]] \
	-command {Execute vizamrai.vector}


button .volume \
	-text "Volume" \
	-image [image create photo -file \
	        [file join $VTK_INSTALL_PATH lib/vtk/tcl/vizamrai/images volume.gif]] \
	-command {Execute vizamrai.vol}


button .exit \
	-text "Exit" \
	-command {exit}

pack .threeslice -fill x
pack .singleslice -fill x
pack .carpet -fill x
pack .iso -fill x
pack .vector -fill x
pack .volume -fill x
pack .exit -fill x
