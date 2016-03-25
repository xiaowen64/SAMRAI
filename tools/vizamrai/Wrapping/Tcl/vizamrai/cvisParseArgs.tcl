##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisParseArgs.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Parse common command line arguments
##


# Parse the command line and default files for settings
# Returns list of unused arguments
proc cvisParseArgs {arguments } {
    upvar $arguments parsed_args
    global argv 

    # Force loading of base stuff
    cvisDisplayObject changethis "changethis"
    cvisTrace foo

    set optlist {
	boolean lowmem
	string camera
	string restart
    }

    while { ![typedopts $argv $optlist found values argv ]} {
	puts stderr "Error: parsing command line $values(_ERROR_)"
    }

    if { $found(lowmem) } {
	set ::cvis::options::LOWMEM 1
    }

    if { $found(camera) } {
	if { [file exists $values(camera)] } {
	    set ::cvis::options::DEFAULT_CAMERA $values(camera)
	    set parsed_args(camera) $values(camera)
	} {
	    puts stderr "Error: camera file <$values(camera)> does not exist"
	}
    }

    if { $found(restart) } {
	if { [file exists $values(restart)] } {
	    set parsed_args(restart) $values(restart)
	} {
	    puts stderr "Error: restart file <$values(restart)> does not exist"
	}
    }

    if [string length "$argv"] {
	set parsed_args(filename) "$argv"
    }
    
    return
}
