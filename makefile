CC=g++

all: 3DeeGen_k.o 3DeeGen_angk_angNN.o 3DeeGen_k_angk.o 3DeeGen_k_angk_angNN.o 3DeeGen_k_angk_angNN_Lab.o 3DeeGen_T1_ang1_T2_ang2.o 3DeeGen_Tc_angc_angd.o 3DeeGen_Tc_angc_angd_betad.o KECal_Tc_angc_angd_betad.o

3DeeGen_k.o: 3DeeGen_k.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k.cpp -o 3DeeGen_k.o

3DeeGen_angk_angNN.o: 3DeeGen_angk_angNN.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_angk_angNN.cpp -o 3DeeGen_angk_angNN.o

3DeeGen_k_angk.o: 3DeeGen_k_angk.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k_angk.cpp -o 3DeeGen_k_angk.o

3DeeGen_k_angk_angNN.o: 3DeeGen_k_angk_angNN.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k_angk_angNN.cpp -o 3DeeGen_k_angk_angNN.o

3DeeGen_k_angk_angNN_Lab.o: 3DeeGen_k_angk_angNN_Lab.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k_angk_angNN_Lab.cpp -o 3DeeGen_k_angk_angNN_Lab.o

3DeeGen_T1_ang1_T2_ang2.o: 3DeeGen_T1_ang1_T2_ang2.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_T1_ang1_T2_ang2.cpp -o 3DeeGen_T1_ang1_T2_ang2.o

3DeeGen_Tc_angc_angd.o: 3DeeGen_Tc_angc_angd.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd.cpp -o 3DeeGen_Tc_angc_angd.o

3DeeGen_Tc_angc_angd_betad.o: 3DeeGen_Tc_angc_angd_betad.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd_betad.cpp -o 3DeeGen_Tc_angc_angd_betad.o

KECal_Tc_angc_angd_betad.o: KECal_Tc_angc_angd_betad.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) KECal_Tc_angc_angd_betad.cpp -o KECal_Tc_angc_angd_betad.o

clean:
	rm -rfv *.o 
