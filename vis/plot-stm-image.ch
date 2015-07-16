set terminal postscript enhanced color eps size 3,3
set output 'surface_plot.eps'
set view 0,0

set size ratio 1

unset surface
set dgrid3d 110,110
set cbrange [0.0:2]
set xrange [0:111]
set yrange [0:111]
unset xtics
unset ytics
set colorbox horizontal
set colorbox user origin 0.3,0.15 size 0.4,0.1
set cbtics ("Pit" 0, "No defect" 2)
set palette defined ( 0 "#2d0000",0.1 "#931c00", 1 "#d68800", 2 "#d68800" )
set pm3d at b corners2color c2
splot 'forplot.txt' using 1:2:3 notitle
exit
