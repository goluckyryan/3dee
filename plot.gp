file1 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testBS.dat"
file2 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testBS_nL.dat"


#file1 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testDirac_nD_nL.dat"
#file2 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testDirac_nD_L.dat"
#file3 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testDirac_D_nL.dat"
#file4 = "~/3dee_program/sample/result/3d_23F_Ti0289_Sp13.3_Tc005_ang100_phi100_testDirac_D_L.dat"

set xlabel "T1 [MeV]"
set ylabel "triple cross section [ub/sr/sr]"

set grid

plot file1 u 1:19 every ::1 with lp title "nD_nL",\
     file2 u 1:19 every ::1 with lp title "nD_L"
#     file3 u 1:19 every ::1 with lp title "D_nL",\
#     file4 u 1:19 every ::1 with lp title "D_L"
     
pause -1  "Hit return to continue"  # Wait until a carriage return is hit
