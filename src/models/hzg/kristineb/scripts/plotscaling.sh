#! /bin/bash
#
#author: Michael Bengfort (michael.bengfort@hzg.de)
#This script reads scaling parameters from kristineb_scal.nml and creat a 
#gnuplot-scipt for plotting the scaling relations
# run with "./plotscaling.sh 'x11'"
terminal=${1}
plotfile='plotscaling.plt'
input='../kristineb_scal.nml'
tmpfile='tmpfile1'
plotfile='scal.plt'
rm ${plotfile}
echo 'reset' >> ${plotfile}
echo 'unset multiplot' >> ${plotfile}
tail -n 26 ${input} >> ${tmpfile}
head -n 25 ${tmpfile} >> ${plotfile}
rm ${tmpfile}
sed -i "s/,/ /g " ${plotfile} 
sed -i "s/!/#/g " ${plotfile} 
if [ ${terminal} = 'x11' ]; then
  echo "set terminal x11 enhanced" >> ${plotfile}
else
  echo 'set terminal pdf enhanced color fontscale 0.2' >> ${plotfile}
  echo "set output 'scaling.pdf'" >> ${plotfile}
fi
echo 'set multiplot layout 2,3' >> ${plotfile}
#convert to ESD
echo " b_qmin_N2=10.0**b_qmin_N" >>${plotfile}
echo " b_qmax_N2=10.0**b_qmax_N" >>${plotfile}
echo " b_vmax_N2=10.0**b_vmax_N" >>${plotfile}
echo " b_qmax_P2=10.0**b_qmax_P" >>${plotfile}
echo " b_qmin_P2=10.0**b_qmin_P" >>${plotfile}
echo " b_vmax_P2=10.0**b_vmax_P" >>${plotfile}
echo " b_carbon=10.0**b_carbon" >> ${plotfile}
echo " b_mumax2=10.0**b_mumax" >> ${plotfile}
echo " b_Mumax2=b_mumax2*(pi/6.)**a_mumax" >> ${plotfile}
echo " a_Mumax2=3*a_mumax" >> ${plotfile}
echo " b_mumax_large2=10.0**b_mumax_large" >> ${plotfile}
echo " b_Mumax_large2=b_mumax_large2*(pi/6.)**a_mumax_large" >> ${plotfile}
echo " a_Mumax_large2=3*a_mumax_large" >> ${plotfile} 
echo " b_Mumax_small2=10.0**b_mumax_small" >> ${plotfile} 
echo " b_Mumax_small2=b_Mumax_small2*(pi/6.)**a_mumax_small" >> ${plotfile}
echo " a_Mumax_small2=3*a_mumax_small" >> ${plotfile} 
echo " b_kn_N=10**b_kn_N" >> ${plotfile}
echo " b_kn_N=b_kn_N*(pi/6.)**a_kn_N" >> ${plotfile}
echo " a_kn_N=3*a_kn_N" >> ${plotfile}
echo " b_kn_P=10**b_kn_P" >> ${plotfile}
echo " b_kn_P=b_kn_P*(pi/6.)**a_kn_P" >> ${plotfile}
echo " a_kn_P=3*a_kn_P" >> ${plotfile}
echo " mumax_kr=2.26" >> ${plotfile}
#-
echo " b_qmin_N_tmp=(b_qmin_N2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_qmin_N=b_qmin_N_tmp*((pi/6.0)**(a_qmin_N-a_carbon))" >> ${plotfile}
echo " a_qmin_N=3.0*(a_qmin_N-a_carbon)" >> ${plotfile}
#
echo " b_qmax_N_tmp=(b_qmax_N2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_qmax_N=b_qmax_N_tmp*((pi/6.0)**(a_qmax_N-a_carbon))" >> ${plotfile}
echo " a_qmax_N=3.0*(a_qmax_N-a_carbon)" >> ${plotfile}
#
echo " b_qmax_P_tmp=(b_qmax_P2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_qmax_P=b_qmax_P_tmp*((pi/6.0)**(a_qmax_P-a_carbon))" >> ${plotfile}
echo " a_qmax_P=3.0*(a_qmax_P-a_carbon)" >> ${plotfile}
#
echo " b_qmin_P_tmp=(b_qmin_P2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_qmin_P=b_qmin_P_tmp*((pi/6.0)**(a_qmin_P-a_carbon))" >> ${plotfile}
echo " a_qmin_P=3.0*(a_qmin_P-a_carbon)" >> ${plotfile}
#
echo " b_vmax_P_tmp=(b_vmax_P2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_vmax_P=b_vmax_P_tmp*((pi/6.0)**(a_vmax_P-a_carbon))" >> ${plotfile}
echo " a_vmax_P=3.0*(a_vmax_P-a_carbon)" >> ${plotfile}
#
echo " b_vmax_N_tmp=(b_vmax_N2/b_carbon)*12.0*10.0**6  #moleX/moleC" >> ${plotfile}
echo " b_vmax_N=b_vmax_N_tmp*((pi/6.0)**(a_vmax_N-a_carbon))" >> ${plotfile}
echo " a_vmax_N=3.0*(a_vmax_N-a_carbon)" >> ${plotfile}
#----



echo "set logscale y" >> ${plotfile}
echo "set xrange [0:6]" >> ${plotfile}

echo "set title '{/Symbol m}_{max}'" >> ${plotfile}
echo "f(x)=b_Mumax_small2*exp(a_Mumax_small2*x)" >> ${plotfile}
echo "g(x)=b_Mumax_small2*exp(a_Mumax_small2*mumax_kr)*exp(a_Mumax_large2*(x-mumax_kr))" >> ${plotfile}
echo "set xlabel 'log_{ESD}'" >> ${plotfile}
echo "plot x < mumax_kr ? f(x) : g(x) t sprintf(\"{/Symbol a}_{s}= %1.2f , {/Symbol b}_{s}= %1.2f  ) \n {/Symbol a}_{l}= %1.2f , {/Symbol b}_{l}=-  )\", a_Mumax_small2, b_Mumax_small2, a_Mumax_large2)" >> ${plotfile}


echo "set title 'k_N'" >> ${plotfile}
echo "set key" >> ${plotfile}
echo "f(x)=b_kn_N*exp(a_kn_N*x)" >> ${plotfile}
echo "g(x)=b_kn_P*exp(a_kn_P*x)" >> ${plotfile}
echo "plot f(x) t sprintf(\"N ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_kn_N, b_kn_N), g(x) t sprintf(\"P ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_kn_P, b_kn_P)" >> ${plotfile}

echo "set title 'Carbon'" >> ${plotfile}
echo "f(x)=b_carbon*exp(a_carbon*x)" >> ${plotfile}
echo "plot f(x) t sprintf(\"{/Symbol a}= %1.2f , {/Symbol b}= %1.2f  \", a_carbon, b_carbon)" >> ${plotfile}

echo "set title 'Q_{min}'" >> ${plotfile}
echo "f(x)=b_qmin_N*exp(a_qmin_N*x)" >> ${plotfile}
echo "g(x)=b_qmin_P*exp(a_qmin_P*x)" >> ${plotfile}
echo "plot f(x) t sprintf(\"N ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_qmin_N, b_qmin_N), g(x) t sprintf(\"P ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_qmin_P, b_qmin_P)" >> ${plotfile}

echo "set title 'Q_{max}'" >> ${plotfile}
echo "f(x)=b_qmax_N*exp(a_qmax_N*x)" >> ${plotfile}
echo "g(x)=b_qmax_P*exp(a_qmax_P*x)" >> ${plotfile}
echo "plot f(x) t sprintf(\"N ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_qmax_N, b_qmax_N), g(x) t sprintf(\"P ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_qmax_P, b_qmax_P)" >> ${plotfile}

echo "set title 'v_{max}'" >> ${plotfile}
echo "f(x)=b_vmax_N*exp(a_vmax_N*x)" >> ${plotfile}
echo "g(x)=b_vmax_P*exp(a_vmax_P*x)" >> ${plotfile}
echo "plot f(x) t sprintf(\"N ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_vmax_N, b_vmax_N), g(x) t sprintf(\"P ({/Symbol a}= %1.2f , {/Symbol b}= %1.2f  )\", a_vmax_P, b_vmax_P)" >> ${plotfile}
#
#
gnuplot ${plotfile}