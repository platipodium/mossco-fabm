#! /bin/bash
#
#author: Michael Bengfort (michael.bengfort@hzg.de)
#This script runs the model for low and high CO2-pressure.
sed -i s/low/high/g run.nml
sed -i s/'co2_low = .true.,'/'co2_low = .false.,'/g kristineb_switch.nml
./fabm0d
sed -i s/high/low/g run.nml
sed -i s/'co2_low = .false.,'/'co2_low = .true.,'/g kristineb_switch.nml
./fabm0d
