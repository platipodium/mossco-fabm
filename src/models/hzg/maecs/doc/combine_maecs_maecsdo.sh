#!/bin/bash

#move the original maecs.F90 as maecs.F90_ORIGINAL
mv ../maecs.F90 ../maecs_F90

#create a temporary file with some documentation
echo "!> @file maecs_maecsdo_combined.F90" > ../maecs.F90_1
#echo -en "\n" >> ../maecs.F90_1
echo "!> @brief A temporary file built by combining maecs.F90 and maecs_do.F90 for documentation purposes">> ../maecs.F90_1
echo -en "\n" >> ../maecs.F90_1

#append the original file and remove the line containing 'end module fabm_hzg_maecs' 
cat ../maecs_F90 >> ../maecs.F90_1
sed -i '/end module fabm_hzg_maecs/d' ../maecs.F90_1

#append the contents of maecs_do inside the intermediate file
cp ../maecs.F90_1 ../maecs.F90_2
cat ../maecs_do.F90 >> ../maecs.F90_2 

#append the line 'end module fabm_hzg_maecs' to the end of the file
cp ../maecs.F90_2 ../maecs.F90_3
echo -en "\n" >> ../maecs.F90_3
echo "end module fabm_hzg_maecs" >> ../maecs.F90_3
#echo -en "\n" >> ../maecs.F90_3

#copy the resulting maecs.F90_3 file as an .F90 code so that it's parsed by doxygen
cp ../maecs.F90_3 maecs_maecsdo_combined.F90

#rename the maecs_do.F90 as a non-fortran file so that it's NOT parsed by doxygen
mv ../maecs_do.F90 ../maecs_do_F90

#remove the intermediate files
rm ../maecs.F90_*
