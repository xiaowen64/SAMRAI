proc MakeBuilding {name bounds texture texture_scaling} {
    set scaling [readSamrai GetScaling]



    for {set i 0} { $i < 6} {incr i} {
	set bounds [lreplace $bounds $i $i [expr [lindex $bounds $i] * $scaling]]
    }
    vtkCubeSource building.$name.source
    eval building.$name.source SetBounds $bounds

    vtkPNMReader building.$name.tiff
    building.$name.tiff SetFileName $texture

    vtkTexture building.$name.atext
    building.$name.atext SetInput [building.$name.tiff GetOutput]
    building.$name.atext InterpolateOn

    vtkTransformTextureCoords building.$name.xform
    eval building.$name.xform SetScale $texture_scaling
    building.$name.xform SetInput [building.$name.source GetOutput]

    vtkPolyDataMapper building.$name.mapper
    building.$name.mapper SetInput [building.$name.xform GetOutput]

    vtkActor building.$name.actor
    building.$name.actor SetMapper building.$name.mapper
    building.$name.actor SetTexture building.$name.atext

    renderWindow AddActor building.$name.actor
}

proc DeleteBuilding {name} {
    renderWindow RemoveActor building.$name.actor
    building.$name.actor Delete
    building.$name.mapper Delete
    building.$name.tiff Delete
    building.$name.atext Delete
    building.$name.xform Delete
    building.$name.source Delete
}

MakeBuilding p1  { 12.0 24.0 12.0 24.0 0.0 48.0 } \
	../misc/Bumpplat.ppm {0.01 0.01 0.01}

MakeBuilding p2  { 60.0 72.0 12.0 24.0 0.0 48.0 } \
	../misc/Bumpplat.ppm {0.01 0.01 0.01}

MakeBuilding p3  { 36.0 48.0 28.0 40.0 24.0 72.0} \
	../misc/Bumpplat.ppm {0.01 0.01 0.01}

MakeBuilding p4 { 12.0 24.0 44.0 56.0 0.0 48.0 } \
	../misc/Bumpplat.ppm {0.01 0.01 0.01}

MakeBuilding p5 { 60.0 72.0 44.0 56.0 0.0 48.0 } \
	../misc/Bumpplat.ppm {0.01 0.01 0.01}

