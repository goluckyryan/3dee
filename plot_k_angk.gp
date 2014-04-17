## the bash command is
# $ gnuplot -e "variable=value" -p plot_kp.gp

##### command line -e input ##################

#if (!exists("pngIO")) pngIO = 0 ; 
#print 'pngIO = '.pngIO

##### SETTING ################################
filename = "paraOut_16O_JA1.0_JB1.5_angNN90_Sp10.0.dat"
sizeX = 600
sizeY = 600
##### Plotting script ########################
reset #reset all graphic setting
print '============================================'
print 'file : ',filename
pngName = filename[9:strlen(filename)-4].'.png'

set term x11 size sizeX,sizeY
set style data lines
#set key autotitle columnhead
set xrange [0:300]
set yrange [0:180]
set xlabel "k [MeV/c]"
set ylabel "angk [deg]"
set dgrid3d 300/5,180/5,1

#set contour base #set contour at the base
#set cntrparam levels incremental -1, 0.1, 1
#set pm3d at b   #set color at base, no contour line
set pm3d map  # set a 2D map
#set palette maxcolors 100

#========= multiplot
set multiplot layout 3, 2
set tmargin 0
set bmargin 0
set lmargin 2
set rmargin 0

#set cbrange [0:0.04]
set title "DWIA 1p3/2"  ;splot filename using 1:2:15 title ""  pal
set title "DWIA  1p1/2" ;splot filename using 1:2:19 title "" pal

set cbrange [-0.5:0.5]
set title "A00n0 1p3/2" ;splot filename using 1:2:16 title ""  pal
set title "A00n0 1p1/2" ;splot filename using 1:2:20 title "" pal

set cbrange [-0.6:0.6]
set title "A00n0 1p3/2-1p1/2";splot filename using 1:2:($16-$20) title "" pal

unset multiplot
#---------- pause and wait for eneter to stop

pause mouse keypress (" ==== press s to save, any for exit \n")
#print 'keypress ASCII = '.MOUSE_KEY

if (MOUSE_KEY == 115){
   set term png small size sizeX, sizeY
   set output pngName
   print 'png output -> '.pngName
   
}
if (MOUSE_KEY != 115){
   print 'no output'
}