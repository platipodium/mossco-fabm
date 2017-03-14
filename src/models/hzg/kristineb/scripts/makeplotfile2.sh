#! /bin/bash
#
#author: Michael Bengfort (michael.bengfort@hzg.de)
# run with "./makeplotfile.sh 'testscript.pdf' 'output.dat' 13 'low' 'x11'"
plotfile='plotfile.plt'
Obsfolder='../Observations'
outputfile=${1}
inputfile1=${2}
inputfile2=${3}
phyto_num=${4}
terminal=${5}
#########################
rm ${plotfile}
echo 'reset' >> ${plotfile}
echo 'unset multiplot' >> ${plotfile}
if [ ${terminal} = 'x11' ]; then
  echo "set terminal x11 enhanced" >> ${plotfile}
else
  echo 'set terminal pdf enhanced color fontscale 0.2' >> ${plotfile}
  echo "set output '$outputfile'" >> ${plotfile}
fi
echo 'set multiplot layout 2,3' >> ${plotfile}
echo 'set xdata time' >> ${plotfile}
echo 'set timefmt "%Y-%m-%d %H:%M:%S"' >> ${plotfile}
echo 'set format x "%m-%d"' >> ${plotfile}
echo 'set ytics nomirror' >> ${plotfile}
#Define Linecolors
echo "set style line 2  lc rgb '#0025ad' lt 1 lw 1.5 # --- blue" >>${plotfile}
echo "set style line 3  lc rgb '#0042ad' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 4  lc rgb '#0060ad' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 5  lc rgb '#007cad' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 6  lc rgb '#0099ad' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 7  lc rgb '#00ada4' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 8  lc rgb '#00ad88' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 9  lc rgb '#00ad6b' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 10 lc rgb '#00ad4e' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 11 lc rgb '#00ad31' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 12 lc rgb '#00ad14' lt 1 lw 1.5 #      ." >>${plotfile}
echo "set style line 13 lc rgb '#09ad00' lt 1 lw 1.5 # --- green" >>${plotfile}
echo "set style line 14 lc rgb '#08ad00' lt 1 lw 1.5 # ---?" >>${plotfile}
######

# ##################################
let row=3*${phyto_num}+3
echo "set key top right" >> ${plotfile}
echo "set ylabel 'N ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
printf "plot '${inputfile1}'  u 1:${row}  w l t 'low CO2' lc rgb 'blue' lw 1.5,'${Obsfolder}/lowCO2.dat' u 1:3 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ", '${inputfile2}'  u 1:${row}  w l t 'high CO2' lc rgb 'red' lw 1.5,'${Obsfolder}/highCO2.dat' u 1:3 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}

let row=3*${phyto_num}+4
echo "set key top right" >> ${plotfile}
echo "set ylabel 'P ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
printf "plot '${inputfile1}' u 1:${row} w l t 'low CO2' lc rgb 'blue' lw 2, '${Obsfolder}/lowCO2.dat' u 1:4 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ",'${inputfile2}' u 1:${row} w l t 'high CO2' lc rgb 'red' lw 2,'${Obsfolder}/highCO2.dat' u 1:4 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}

let row=${row}+3
echo "unset y2tics" >> ${plotfile}
echo "unset y2label" >> ${plotfile}
echo "set ylabel 'Chl_a {/Symbol m}gL^{-1}'" >> ${plotfile}
printf "plot '${inputfile1}' u 1:${row} w l lw 2 lc rgb 'blue' t 'low CO2', '${Obsfolder}/lowCO2.dat' u 1:7 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >>${plotfile}
printf ", '${inputfile2}' u 1:${row} w l lw 2 lc rgb 'red' t 'high CO2', '${Obsfolder}/highCO2.dat' u 1:7 w p pt 6 ps 0.5 lc rgb 'red' notitle" >>${plotfile}
printf "\n" >> ${plotfile}
###
let row=${row}+3
##
echo "set ylabel '<L> ({/Symbol m}m)'" >>${plotfile}
printf "plot '${inputfile1}' u 1:${row} w l lc rgb 'blue' lw 2 t 'low CO2','${Obsfolder}/lowCO2.dat' u 1:8 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ", '${inputfile2}' u 1:${row} w l lc rgb 'red' lw 2 t 'high CO2','${Obsfolder}/highCO2.dat' u 1:8 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}
let tmp=$row
###
echo "set title 'low CO2'" >>${plotfile}
let row=${tmp}+1+6*${phyto_num}+2+3
echo "unset key" >> ${plotfile}
echo "set ylabel 'day of experiment'" >> ${plotfile}
echo "set xlabel 'size class'" >> ${plotfile}
echo "set tics out"  >> ${plotfile}
echo "unset xtics"  >> ${plotfile}
echo "set pm3d map" >> ${plotfile}
echo "set xdata" >> ${plotfile}
echo "set format x '%1.0f'" >> ${plotfile}
echo "set ticslevel 0" >> ${plotfile}
let endblock=${row}-1
echo "splot '${inputfile1}' every ::${endblock} matrix w pm3d" >> ${plotfile}
##
echo "set title 'high CO2'" >>${plotfile}
let row=${tmp}+1+6*${phyto_num}+2+3
echo "unset key" >> ${plotfile}
echo "set ylabel 'day of experiment'" >> ${plotfile}
echo "set xlabel 'size class'" >> ${plotfile}
echo "set pm3d map" >> ${plotfile}
echo "set xdata" >> ${plotfile}
echo "set format x '%1.0f'" >> ${plotfile}
echo "set ticslevel 0" >> ${plotfile}
let endblock=${row}-1
echo "splot '${inputfile2}' every ::${endblock} matrix w pm3d" >> ${plotfile}

###
echo "unset multiplot" >> ${plotfile}
echo "set terminal x11" >> ${plotfile}