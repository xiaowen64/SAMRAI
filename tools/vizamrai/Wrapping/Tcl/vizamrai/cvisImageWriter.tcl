##
## File:        cvisRenderWindow.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2001 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Render Window for 3D data
##

class cvisImageWriter {
    variable FileTypes { \
	    {PNG {.png}} \
	    {BMP {.bmp}} \
	    {JPEG {.jpg}} \
	    {PPM {.ppm}} \
	    {PostScript {.ps}} \
	    {TIFF {.tif}} }

    variable WriterTable

    variable Type "png"

    variable Input ""

    variable FileName "" 
    
    variable Writer ""

    method constructor { } {
	set WriterTable(bmp) vtkBMPWriter
	set WriterTable(jpg) vtkJPEGWriter
	set WriterTable(png) vtkPNGWriter
	set WriterTable(ppm) vtkPNMWriter
	set WriterTable(ps)  vtkPostScriptWriter
	set WriterTable(tif) vtkPostTIFFWriter
    }

    method destructor { } {
	if {[string length $Writer]} {
	    $Writer Delete
	    set Writer ""
	}
    }
    
    method GetFileTypes {} {
	return $FileTypes
    }

    method SetType {type} {
	set Type $type
    }

    method GetType {} {
	return $Type
    }

    method SetInput {input} {
	set Input $input
	if {[string length $Writer]} {
	    $Writer SetInput $input
	}
    }

    method GetInput {} {
	return $Input
    }

    method SetFileName {filename} {
	set FileName $filename
	set Type [string trimleft [file extension $FileName] "."]

	if {[string length $Writer]} {
	    if {[string compare $WriterTable($Type) [$Writer GetClassName]]} {
		$Writer Delete
		set Writer ""
	    
	    }
	}

	if {![string length $Writer]} {
	    set Writer [$WriterTable($Type) $this.writer]
	}

	$Writer SetFileName $FileName
	if {[string length $Input]} {
	    $Writer SetInput $Input
	}
    }

    method GetFileName {} {
	return $FileName
    }

    method Write {} {
	if {[string length $Writer]} {
	    $Writer Write
	}
    }
}