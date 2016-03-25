proc MakeBuilding {name bounds} {
    set scaling [readSamrai GetScaling]

    for {set i 0} { $i < 6} {incr i} {
	set bounds [lreplace $bounds $i $i [expr [lindex $bounds $i] * $scaling]]
    }
    vtkCubeSource building.$name.source
    eval building.$name.source SetBounds $bounds

    vtkPolyDataMapper building.$name.mapper
    building.$name.mapper SetInput [building.$name.source GetOutput]

    vtkActor building.$name.actor
    building.$name.actor SetMapper building.$name.mapper

    renderWindow AddActor building.$name.actor
}

proc DeleteBuilding {name} {
    renderWindow RemoveActor building.$name.actor
    building.$name.actor Delete
    building.$name.mapper Delete
    building.$name.source Delete
}

MakeBuilding 1  { 8.0 32.0 32.0 48.0 0.0 40.0 }
MakeBuilding 2  {88.0 112.0 16.0 40.0 0.0 56.0}
MakeBuilding 3a {48.0 72.0  8.0 16.0 0.0 32.0}
MakeBuilding 3b {56.0 64.0 16.0 24.0 0.0 32.0}
MakeBuilding 4a {56.0 64.0 32.0 40.0 0.0 32.0}
MakeBuilding 4b {48.0 72.0 40.0 48.0 0.0 32.0}


