for {set plane 0} {$plane < 160} {incr plane 1} {
    slicePlane.x SetSlice $plane
    renderWindow SetFileName frames/frame.[format "%05d" $plane].ppm
    renderWindow Render
    renderWindow SaveImageAsPPM
}
