CC=g++

all: make_infile.o read_outfile.o 3DeeGen_k.o 3DeeGen_angk_angNN.o 3DeeGen_k_angk.o 3DeeGen_k_angk_angNN.o 3DeeGen_k_angk_angNN_Lab.o 3DeeGen_T1_ang1_T2_ang2.o 3DeeGen_Tc_angc_angd.o

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

3DeeGen_Tc_angc_angd.o: 3DeeGen_Tc_angc_angd.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd.cpp -o 3DeeGen_Tc_angc_angd.o

make_infile.o: make_infile.cpp
	$(CC) make_infile.cpp -o make_infile.o

read_outfile.o: read_outfile.cpp
	$(CC) read_outfile.cpp -o read_outfile.o

clean:
	rm -rfv *.o 
