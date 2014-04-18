set xrange [0:300]
set yrange [0:360]
set xlabel "k [MeV/c]"
set ylabel "angk [deg]"
set dgrid3d 300/5,360/5,1

#set contour base #set contour at the base
#set cntrparam levels incremental -1, 0.1, 1
#set pm3d at b   #set color at base, no contour line
set pm3d map  # set a 2D map
#set palette defined (-0.3 "black",-0.15 "blue",0 "white", 0.15 "green", 0.3 "red")
set palette defined ( 0.0  "#000000",\
                      0.1  "#FF00FF",\
                      0.2  "#0000FF",\
                      0.3  "#00FFFF",\
                      0.4  "#00FF33",\
                      0.5  "#00FF00",\
                      0.6  "#33FF00",\
                      0.7  "#FFFF00",\
                      0.8  "#FF6600",\
                      0.9  "#FF0000",\
                      1.0  "#FFFFFF" )

#========= multiplot
set multiplot layout 3, 6 title filename[9:strlen(filename)-4]
set tmargin 0
set bmargin 0
set lmargin 5
set rmargin 0

set cbrange [0:0.15]
set cbtics 0.05
set format cb "%5.2f"
set title "DWIA 1s1/2 x 1  " ;splot filename using 1:2:($11*1) title "" pal
set title "DWIA 1p3/2 x 4  " ;splot filename using 1:2:($15*4) title "" pal
set title "DWIA 1p1/2 x 4  " ;splot filename using 1:2:($19*4) title "" pal
set title "DWIA 1d5/2 x 9  " ;splot filename using 1:2:($23*9) title "" pal
set title "DWIA 2s1/2 x 0.5" ;splot filename using 1:2:($27*0.5) title "" pal
set title "DWIA 1d3/2 x 9  " ;splot filename using 1:2:($31*9) title "" pal

set cbrange [-0.5:0.5]
set cbtics 0.1        
set title "A00n0 1s1/2" ;splot filename using 1:2:12 title "" pal
set title "A00n0 1p3/2" ;splot filename using 1:2:16 title "" pal
set title "A00n0 1p1/2" ;splot filename using 1:2:20 title "" pal
set title "A00n0 1d5/2" ;splot filename using 1:2:24 title "" pal
set title "A00n0 2s1/2" ;splot filename using 1:2:28 title "" pal
set title "A00n0 1d3/2" ;splot filename using 1:2:32 title "" pal

set cbrange [-0.5:0.5]
set title "A00n0 1s1/2-2s1/2";splot filename using 1:2:($12-$28) title "" pal
set title "A00n0 1p3/2-1p1/2";splot filename using 1:2:($16-$20) title "" pal
set title ""; splot 0.1
set title "A00n0 1d5/2-1d3/2";splot filename using 1:2:($24-$32) title "" pal

unset multiplot
