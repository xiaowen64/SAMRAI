#set start 533333 
#set end  -533333
set start 64
set end  [expr -$start]

set cropping(x,min) 0
set cropping(x,max) 320

set cropping(y,min) -128
set cropping(y,max) 128

set cropping(z,min) -128

set slices 200

set spacing [expr ($end - $start) / $slices]

for {set plane 0} {$plane < $slices} {incr plane 1} {
    puts "Rendering start $plane"

    set cropping(z,max) [expr $start + $plane * $spacing]

    vol SetCroppingRegionPlanes \
	    $cropping(x,min) $cropping(x,max) \
	    $cropping(y,min) $cropping(y,max) \
	    $cropping(z,min) $cropping(z,max) 

    puts [vol GetCroppingRegionPlanes]
	    
    renderWindow SetFileName frames/frame.[format "%05d" $plane].ppm
    renderWindow Render
    renderWindow SaveImageAsPPM
    puts "Rendering done $plane"
}
