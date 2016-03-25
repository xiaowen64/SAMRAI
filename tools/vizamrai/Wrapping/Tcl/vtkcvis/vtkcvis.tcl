if {[info commands vtkSamraiStructuredPointsReader] != "" ||
    [::vtk::load_component vtkcvisTCL] == ""} {

    package provide vtkcvis 1.0
}

