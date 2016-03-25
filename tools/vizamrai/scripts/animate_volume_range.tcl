set start 1
set end   1e14



set slices 10

set spacing [expr log($end - $start) / $slices]

for {set plane 0} {$plane < $slices} {incr plane 1} {
    set min_range [expr $start + exp($plane * $spacing)]
    set max_range [expr $start + exp(($plane+4) * $spacing)]

    puts "$min_range $max_range"

    colorTable SetTableRange $min_range $max_range

    renderWindow SetFileName frames/frame.[format "%05d" $plane].ppm
    renderWindow Render
    renderWindow SaveImageAsPPM
}
