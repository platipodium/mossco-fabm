#! /bin/bash
#
#author: Michael Bengfort (michael.bengfort@hzg.de)
# run with "./makeplotfile.sh 'testscript.pdf' 'output.dat' 17 'low' 'x11'"
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
echo 'set multiplot layout 3,3' >> ${plotfile}
#echo 'set xdata time' >> ${plotfile}
#echo 'set timefmt "%Y-%m-%d %H:%M:%S"' >> ${plotfile}
#echo 'set format x "%m-%d"' >> ${plotfile}
#echo 'set ytics nomirror' >> ${plotfile}
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
echo "set xlabel 'day of experiment'" >> ${plotfile}
# ##################################
let row=3*${phyto_num}+3
echo "set key top right" >> ${plotfile}
echo "set ylabel 'N ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
#printf "plot '${inputfile1}'  u 1:${row}  w l t 'low CO2' lc rgb 'blue' lw 1.5,'${Obsfolder}/lowCO2.dat' u 1:3 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
#printf ", '${inputfile2}'  u 1:${row}  w l t 'high CO2' lc rgb 'red' lw 1.5,'${Obsfolder}/highCO2.dat' u 1:3 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf "plot '${inputfile1}'  u (\$0*2):${row}  w l t 'low CO2' lc rgb 'blue' lw 1.5,'${Obsfolder}/lowCO2.dat' u (\$0*2):3 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ", '${inputfile2}'  u (\$0*2):${row}  w l t 'high CO2' lc rgb 'red' lw 1.5,'${Obsfolder}/highCO2.dat' u (\$0*2):3 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}

let row=3*${phyto_num}+4
echo "set key top right" >> ${plotfile}
echo "set ylabel 'P ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
#printf "plot '${inputfile1}' u 1:${row} w l t 'low CO2' lc rgb 'blue' lw 2, '${Obsfolder}/lowCO2.dat' u 1:4 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
#printf ",'${inputfile2}' u 1:${row} w l t 'high CO2' lc rgb 'red' lw 2,'${Obsfolder}/highCO2.dat' u 1:4 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf "plot '${inputfile1}' u (\$0*2):${row} w l t 'low CO2' lc rgb 'blue' lw 2, '${Obsfolder}/lowCO2.dat' u (\$0*2):4 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ",'${inputfile2}' u (\$0*2):${row} w l t 'high CO2' lc rgb 'red' lw 2,'${Obsfolder}/highCO2.dat' u (\$0*2):4 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}

let row=${row}+3
echo "unset y2tics" >> ${plotfile}
echo "unset y2label" >> ${plotfile}
echo "set ylabel 'Chl_a {/Symbol m}gL^{-1}'" >> ${plotfile}
#printf "plot '${inputfile1}' u 1:${row} w l lw 2 lc rgb 'blue' t 'low CO2', '${Obsfolder}/lowCO2.dat' u 1:7 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >>${plotfile}
#printf ", '${inputfile2}' u 1:${row} w l lw 2 lc rgb 'red' t 'high CO2', '${Obsfolder}/highCO2.dat' u 1:7 w p pt 6 ps 0.5 lc rgb 'red' notitle" >>${plotfile}
printf "plot '${inputfile1}' u (\$0*2):${row} w l lw 2 lc rgb 'blue' t 'low CO2', '${Obsfolder}/lowCO2.dat' u (\$0*2):7 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >>${plotfile}
printf ", '${inputfile2}' u (\$0*2):${row} w l lw 2 lc rgb 'red' t 'high CO2', '${Obsfolder}/highCO2.dat' u (\$0*2):7 w p pt 6 ps 0.5 lc rgb 'red' notitle" >>${plotfile}
printf "\n" >> ${plotfile}
###
let row=${row}+3
##
echo "set ylabel '<L> ({/Symbol m}m)'" >>${plotfile}
#printf "plot '${inputfile1}' u 1:${row} w l lc rgb 'blue' lw 2 t 'low CO2','${Obsfolder}/lowCO2.dat' u 1:9 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
#printf ", '${inputfile2}' u 1:${row} w l lc rgb 'red' lw 2 t 'high CO2','${Obsfolder}/highCO2.dat' u 1:9 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf "plot '${inputfile1}' u (\$0*2):${row} w l lc rgb 'blue' lw 2 t 'low CO2','${Obsfolder}/lowCO2.dat' u (\$0*2):9 w p pt 6 ps 0.5 lc rgb 'blue' notitle" >> ${plotfile}
printf ", '${inputfile2}' u (\$0*2):${row} w l lc rgb 'red' lw 2 t 'high CO2','${Obsfolder}/highCO2.dat' u (\$0*2):9 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
printf '\n' >> ${plotfile}
##
let tmp=$row
echo "set title 'low CO2'">>${plotfile}
echo "set ylabel ''">>${plotfile}
let row=${row}+2+6*${phyto_num}+2+3
let row2=${row}+1
let row3=${row2}+1
let row4=${row3}+1
let row5=${row4}+1
echo "set autoscale y" >> ${plotfile}
printf "plot '${inputfile1}' u (\$0*2):(0*\$0):${row} w filledcu title 'growth', '${inputfile1}' u (\$0*2):(0*\$0) lt -1 notitle, '${inputfile1}' u (\$0*2):${row} lt -1 notitle" >> ${plotfile}
printf ", '${inputfile1}' u (\$0*2):(0*\$0):${row2} w filledcu title 'respiration', '${inputfile1}' u (\$0*2):(0*\$0) lt -1 notitle, '${inputfile1}' u (\$0*2):${row2} lt -1 notitle" >> ${plotfile}
printf ",'${inputfile1}' u (\$0*2):${row2}:(\$${row2}+\$${row3}) w filledcu title 'sinking', '${inputfile1}' u (\$0*2):${row2} lt -1 notitle, '${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}) lt -1 notitle" >> ${plotfile}
printf ",'${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}):(\$${row2}+\$${row3}+\$${row4}) w filledcu title 'aggregation', '${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}) lt -1 notitle, '${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}) lt -1 notitle" >> ${plotfile}
printf ",'${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}):(\$${row2}+\$${row3}+\$${row4}+\$${row5}) w filledcu title 'grazing', '${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}) lt -1 notitle, '${inputfile1}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}+\$${row5}) lt -1 notitle" >> ${plotfile}
printf "\n" >> ${plotfile}
##
echo "set title 'high CO2'">>${plotfile}
let row=${tmp}+2+6*${phyto_num}+2+3
let row2=${row}+1
let row3=${row2}+1
let row4=${row3}+1
let row5=${row4}+1
echo "set autoscale y" >> ${plotfile}
printf "plot '${inputfile2}' u (\$0*2):(0*\$0):${row} w filledcu title 'growth', '${inputfile2}' u (\$0*2):(0*\$0) lt -1 notitle, '${inputfile2}' u (\$0*2):${row} lt -1 notitle" >> ${plotfile}
printf ", '${inputfile2}' u (\$0*2):(0*\$0):${row2} w filledcu title 'respiration', '${inputfile2}' u (\$0*2):(0*\$0) lt -1 notitle, '${inputfile2}' u (\$0*2):${row2} lt -1 notitle" >> ${plotfile}
printf ",'${inputfile2}' u (\$0*2):${row2}:(\$${row2}+\$${row3}) w filledcu title 'sinking', '${inputfile2}' u (\$0*2):${row2} lt -1 notitle, '${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}) lt -1 notitle" >> ${plotfile}
printf ",'${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}):(\$${row2}+\$${row3}+\$${row4}) w filledcu title 'aggregation', '${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}) lt -1 notitle, '${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}) lt -1 notitle" >> ${plotfile}
printf ",'${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}):(\$${row2}+\$${row3}+\$${row4}+\$${row5}) w filledcu title 'grazing', '${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}) lt -1 notitle, '${inputfile2}' u (\$0*2):(\$${row2}+\$${row3}+\$${row4}+\$${row5}) lt -1 notitle" >> ${plotfile}
printf "\n" >> ${plotfile}
let row=${row5}+1
let tmp=$row
###
echo "set title 'biomass (low CO2)'" >>${plotfile}
#let row=${tmp} +1+6*${phyto_num}+2+3
echo "unset key" >> ${plotfile}
echo "set ylabel 'size (log_{ESD})'" >> ${plotfile}
echo "set ylabel offset -1.5,0" >> ${plotfile}
#echo "set tics out"  >> ${plotfile}
#echo "unset xtics"  >> ${plotfile}
echo "set pm3d map" >> ${plotfile}
#echo "set xdata" >> ${plotfile}
#echo "set format x '%1.0f'" >> ${plotfile}
#echo "set ticslevel 0" >> ${plotfile}
let endblock=${row}-1
echo "set xtics out">> ${plotfile}
echo "set yrange [0:${phyto_num}]" >> ${plotfile}
echo "set ytics ('-0.5' 0, '0.0' 1, '0.5' 2, '1.0' 3, '1.3' 4, '1.6' 5, '1.9' 6, '2.2' 7, '2.5' 8, '3.0' 9, '3.5' 10, '4.0' 11, '4.5' 12, '5.0' 13, '5.5' 14, '6.0' 15, '6.2' 16)" >> ${plotfile}
echo "set ytics offset 0,0.5" >> ${plotfile}
echo "splot '${inputfile1}' matrix u (\$2*2):(\$1-${endblock}):3 every ::${endblock} w pm3d" >> ${plotfile} 
##
echo "set title 'biomass (high CO2)'" >>${plotfile}
#let row=${tmp} #+1+6*${phyto_num}+2+3
echo "unset key" >> ${plotfile}
echo "set ylabel 'size (log_{ESD})'" >> ${plotfile}
echo "set ylabel offset -1.5,0" >> ${plotfile}
echo "set pm3d map" >> ${plotfile}
#echo "set xdata" >> ${plotfile}
#echo "set format x '%1.0f'" >> ${plotfile}
#echo "set ticslevel 0" >> ${plotfile}
let endblock=${row}-1
echo "set yrange [0:${phyto_num}]" >> ${plotfile}
echo "splot '${inputfile2}' matrix u (\$2*2):(\$1-${endblock}):3 every ::${endblock} w pm3d" >> ${plotfile}
####
###
echo "unset multiplot" >> ${plotfile}
echo "set terminal x11" >> ${plotfile}
