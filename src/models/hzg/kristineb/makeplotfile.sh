#! /bin/bash
#
# run with "./makeplotfile.sh 'testscript.pdf' 'output.dat' 13 'low' 'x11'"
plotfile='plotfile.plt'
Obsfolder='Observations'
outputfile=${1}
inputfile=${2}
phyto_num=${3}
loworhigh=${4}
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
echo 'set multiplot layout 4,3' >> ${plotfile}
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
echo "unset key" >> ${plotfile}
echo "set ylabel 'biomass ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
echo "set xlabel 'time (d)'" >> ${plotfile}
let lines=2
printf  "plot '${inputfile}' u 1:3 w l ls ${lines} t 'Phy_1'" >> ${plotfile}
for (( k = 1 ; k < ${3} ; k++))
	do
	  let row=3+$k*3
	  let species=1+$k
	  let lines=${lines}+1
	  printf ",'${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'" >> ${plotfile}
done
printf "\n" >> ${plotfile}
echo "set key" >> ${plotfile}
# ##################################
let row=${row}+3
echo 'set y2tics' >> ${plotfile}
echo "set ylabel 'N ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
echo "set y2label 'P ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
printf "plot '${inputfile}' u 1:${row} w l t 'N' lc rgb 'red'" >> ${plotfile}
let row=${row}+1
printf ",'${inputfile}'  u 1:${row} axes x1y2 w l  t 'P' lc rgb 'green','${Obsfolder}/${loworhigh}CO2.dat' u 1:3 axes x1y1 w p pt 6 ps 0.5 lc rgb 'red' notitle,'${Obsfolder}/${loworhigh}CO2.dat' u 1:4 axes x1y2 w p pt 6 ps 0.5 lc rgb 'green' notitle\n" >> ${plotfile}
printf '\n' >> ${plotfile}
echo "set ylabel 'D_N ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
echo "set y2label 'D_P ({/Symbol m}mol-C L^{-1})'" >> ${plotfile}
let row=${row}+1
printf "plot '${inputfile}' u 1:${row} w l t 'D_N'" >> ${plotfile}
let row=${row}+1
printf ", '${inputfile}' u 1:${row} axes x1y2 w l  t 'D_P'\n" >> ${plotfile}
####Focing
###
let row=${row}+1
echo "unset y2tics" >> ${plotfile}
echo "unset y2label" >> ${plotfile}
echo "set ylabel 'Chl_a {/Symbol m}gCL^{-1}'" >> ${plotfile}
echo "plot '${inputfile}' u 1:${row} w l t 'Chl_a', '${Obsfolder}/${loworhigh}CO2.dat' u 1:7 w p pt 6 ps 0.5 lc rgb 'black' notitle" >>${plotfile}
##
echo "set y2tics" >> ${plotfile}
echo "set ylabel '{/Symbol m}m'" >>${plotfile}
let row=${row}+1
printf "plot '${inputfile}' u 1:${row} axes x1y1 w l lc rgb 'red' t 'Mean cell size','${Obsfolder}/${loworhigh}CO2.dat' u 1:8 axes x1y1 w p pt 6 ps 0.5 lc rgb 'red' notitle" >> ${plotfile}
let row=${row}+1
printf ",'${inputfile}' u 1:${row} axes x1y2 w l t 'size diverstiy'\n" >> ${plotfile}
#
echo "unset key" >> ${plotfile}
echo "unset y2tics" >> ${plotfile}
echo '#set xrange ["2013-03-10":"2013-03-17"]' >> ${plotfile}
echo "set ylabel 'PAR-Forcing'" >> ${plotfile}
let row=${row}+1
let species=1
let lines=2
printf "plot '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'" >> ${plotfile}
for (( k=1 ; k < ${3} ; k++))
      do
         let row=${row}+1
         let species=${species}+1
         let lines=${lines}+1
         printf ", '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'" >> ${plotfile}
done
printf '\n'>> ${plotfile}
echo "set autoscale x" >> ${plotfile}
echo "set ylabel 'CO2-Forcing'" >> ${plotfile}
let row=${row}+1
let species=1
let lines=2
printf "plot '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'">> ${plotfile}
for (( k=1 ; k < ${3} ; k++))
      do
         let row=${row}+1
         let species=${species}+1
	 let lines=${lines}+1
	 printf ", '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'" >> ${plotfile}
done
printf '\n' >> ${plotfile}
echo "set ylabel 'grazing-Forcing'" >> ${plotfile}
echo "set logscale y" >> ${plotfile}
let row=${row}+1
let species=1
let lines=2
printf "plot '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'">> ${plotfile}
for (( k=1 ; k < ${3} ; k++))
      do
         let row=${row}+1
         let species=${species}+1
         let lines=${lines}+1
printf ", '${inputfile}' u 1:${row} w l ls ${lines} t 'Phy_${species}'" >> ${plotfile}
done
printf '\n' >> ${plotfile}
echo "set key" >> ${plotfile}
echo "unset logscale y" >> ${plotfile}
printf "\n" >>${plotfile}
##
echo "set key left" >> ${plotfile}
echo "set ylabel 'Temperature forcing'" >> ${plotfile}
let row=${row}+1
printf "plot '${inputfile}' u 1:${row} w l t 'TF_{Phy}'" >> ${plotfile}
let row=${row}+1
printf ", '${inputfile}' u 1:${row} w l t 'TF_{Zoo}'\n" >> ${plotfile}
#Uptake####
echo "set logscale y" >> ${plotfile}
echo "unset key" >>${plotfile}
echo "set ylabel 'uptake N'" >> ${plotfile}
let row=${row}+1
let species=1
let lines=2
printf "plot '${inputfile}' u 1:${row} w l ls ${lines} t 'uptake_N_${species}'" >> ${plotfile}
for (( k=1 ; k < ${3} ; k++))
      do
       let row=${row}+1
       let species=${species}+1
       let lines=${lines}+1
       printf ", '${inputfile}' u 1:${row} w l ls ${lines} t 'uptake_N_${species}'" >> ${plotfile}
done
printf '\n' >> ${plotfile}
echo "set ylabel 'uptake P'" >> ${plotfile}
let row=${row}+1
let species=1
let lines=2
printf "plot '${inputfile}' u 1:${row} w l ls ${lines} t 'uptake_P_${species}'" >> ${plotfile}
for (( k=1 ; k < ${3} ; k++))
      do
       let row=${row}+1
       let species=${species}+1
       let lines=${lines}+1
       printf ", '${inputfile}' u 1:${row} w l ls ${lines} t 'uptake_P_${species}'" >> ${plotfile}
done
printf '\n' >> ${plotfile}
###
echo "unset multiplot" >> ${plotfile}
echo "set terminal x11" >> ${plotfile}