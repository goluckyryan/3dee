## the bash command is
# $ gnuplot -e "variable=value" -p plot_kp.gp

##### command line -e input ##################

#if (!exists("pngIO")) pngIO = 0 ; 
#print 'pngIO = '.pngIO

##### SETTING ################################
filename = "paraOut_16O_JA1.0_JB1.5_angk60_angNN70_Sp10.0.dat"
sizeX = 600
sizeY = 800
##### Plotting script ########################
reset #reset all graphic setting
print '============================================'
print 'file : ',filename
pngName = filename[9:strlen(filename)-4].'.png'

set term x11 size sizeX,sizeY
set grid   #set graphic grid
set style data lines
set key autotitle columnhead
#set key box linestyle -1
set xrange [0:300]
unset xtics

#-------------- set mulit plot <row>, <col>
set multiplot layout 2, 1 title filename
set tmargin 0
#set bmargin 0
set lmargin 10
set rmargin 3
set grid

#--- 1st plot
set origin 0,0.75
set size 1, 0.25
unset xlabel
set ylabel 'DWIA'
set yrange[0:0.25]
set xtics 0, 20, 300
set ytics 0,0.05,0.25
#set xtics ("" tics(0), ""tics(1),""tics(2),""tics(3),""tics(4))

# the sytex (A?B:C) is read as (if A then B else C)
plot filename using 1:($10<10?$10:0)     title "1 x 1s1/2" lc rgb "#FF0000" lw 2, \
           "" using 1:($14<10?4*$14:0)   title "4 x 1p3/2" lc rgb "#FF6600" lw 2, \
           "" using 1:($18<10?4*$18:0)   title "4 x 1p1/2" lc rgb "#999900" lw 2, \
           "" using 1:($22<10?9*$22:0)   title "9 x 1d5/2" lc rgb "#00FF00" lw 2, \
           "" using 1:($26<10?0.5*$26:0) title "0.5 x 2s1/2" lc rgb "#00FFFF" lw 2, \
           "" using 1:($30<10?9*$30:0)   title "9 x 1d3/2" lc rgb "#0000FF" lw 2  

#--- 2nd plot
set origin 0, 0.5
set size 1, 0.25
set ylabel 'A00n0'
set xtics 0,20,300.
set yrange [-1:1]
set ytics -1,0.2,1

plot  0 lc rgb "#000000", \
     filename using 1:($11<10?$11:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($15<10?$15:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($19<10?$19:0) lc rgb "#999900" lw 2, \
           "" using 1:($23<10?$23:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($27<10?$27:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($31<10?$31:0) lc rgb "#0000FF" lw 2

#--- 3rd plot
set origin 0,0.25
set size 1, 0.25
set ylabel 'Pn000'
set xtics 0,20,300.
set yrange [-1:1]
set ytics -1,0.2,1

plot  0 lc rgb "#000000", \
     filename using 1:($12<10?$12:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($16<10?$16:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($20<10?$20:0) lc rgb "#999900" lw 2, \
           "" using 1:($24<10?$24:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($28<10?$28:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($32<10?$32:0) lc rgb "#0000FF" lw 2

#--- 4th plot
set origin 0,0
set size 1, 0.25
set xlabel 'k [MeV/c]'
set ylabel 'P0n00'
set xtics 0,20,300.
set yrange [-1:1]
set ytics -1,0.2,1

plot  0 lc rgb "#000000", \
     filename using 1:($13<10?$13:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($17<10?$17:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($21<10?$21:0) lc rgb "#999900" lw 2, \
           "" using 1:($25<10?$25:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($29<10?$29:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($33<10?$33:0) lc rgb "#0000FF" lw 2

unset multiplot


#---------- pause and wait for eneter to stop

pause mouse keypress (" ==== press s to save, any for exit \n")
#print 'keypress ASCII = '.MOUSE_KEY

if (MOUSE_KEY == 115){
   set term png small size sizeX, sizeY
   set output pngName
   print 'png output -> '.pngName

   
   #-------------- set mulit plot <row>, <col>
   set multiplot layout 2, 1 title filename
   set tmargin 0
   #set bmargin 0
   set lmargin 10
   set rmargin 3
   set grid

   #--- 1st plot
   set origin 0,0.75
   set size 1, 0.25
   unset xlabel
   set ylabel 'DWIA'
   set yrange[0:0.25]
   set xtics 0, 20, 300
   set ytics 0,0.05,0.25
   plot filename using 1:($10<10?$10:0)     title "1 x 1s1/2" lc rgb "#FF0000" lw 2, \
           "" using 1:($14<10?4*$14:0)   title "4 x 1p3/2" lc rgb "#FF6600" lw 2, \
           "" using 1:($18<10?4*$18:0)   title "4 x 1p1/2" lc rgb "#999900" lw 2, \
           "" using 1:($22<10?9*$22:0)   title "9 x 1d5/2" lc rgb "#00FF00" lw 2, \
           "" using 1:($26<10?0.5*$26:0) title "0.5 x 2s1/2" lc rgb "#00FFFF" lw 2, \
           "" using 1:($30<10?9*$30:0)   title "9 x 1d3/2" lc rgb "#0000FF" lw 2  

   #--- 2nd plot
   set origin 0, 0.5
   set size 1, 0.25
   set ylabel 'A00n0'
   set xtics 0,20,300.
   set yrange [-1:1]
   set ytics -1,0.2,1
   
   plot  0 lc rgb "#000000", \
         filename using 1:($11<10?$11:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($15<10?$15:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($19<10?$19:0) lc rgb "#999900" lw 2, \
           "" using 1:($23<10?$23:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($27<10?$27:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($31<10?$31:0) lc rgb "#0000FF" lw 2

   #--- 3rd plot
   set origin 0,0.25
   set size 1, 0.25
   set ylabel 'Pn000'
   set xtics 0,20,300.
   set yrange [-1:1]
   set ytics -1,0.2,1

   plot  0 lc rgb "#000000", \
         filename using 1:($12<10?$12:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($16<10?$16:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($20<10?$20:0) lc rgb "#999900" lw 2, \
           "" using 1:($24<10?$24:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($28<10?$28:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($32<10?$32:0) lc rgb "#0000FF" lw 2

   #--- 4th plot
   set origin 0,0
   set size 1, 0.25
   set xlabel 'k [MeV/c]'
   set ylabel 'P0n00'
   set xtics 0,20,300.
   set yrange [-1:1]
   set ytics -1,0.2,1

   plot  0 lc rgb "#000000", \
        filename using 1:($13<10?$13:0) lc rgb "#FF0000" lw 2, \
           "" using 1:($17<10?$17:0) lc rgb "#FF6600" lw 2, \
           "" using 1:($21<10?$21:0) lc rgb "#999900" lw 2, \
           "" using 1:($25<10?$25:0) lc rgb "#00FF00" lw 2, \
           "" using 1:($29<10?$29:0) lc rgb "#00FFFF" lw 2, \
           "" using 1:($33<10?$33:0) lc rgb "#0000FF" lw 2

   unset multiplot
   
   
}
if (MOUSE_KEY != 115){
   print 'no output'
}