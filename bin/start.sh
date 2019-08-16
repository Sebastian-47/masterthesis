#!/bin/bash

#1Pi/16
sed -i '11s/.*/angle_phi 0.0625/' options.dat
./SplitStep 6 >> rational1.dat

#2Pi/16
sed -i '11s/.*/angle_phi 0.125/' options.dat
./SplitStep 6 >> rational2.dat

#3Pi/16
sed -i '11s/.*/angle_phi 0.1875/' options.dat
./SplitStep 6 >> rational3.dat

#4Pi/16
sed -i '11s/.*/angle_phi 0.25/' options.dat
./SplitStep 6 >> rational4.dat

#5Pi/16
sed -i '11s/.*/angle_phi 0.3125/' options.dat
./SplitStep 6 >> rational5.dat

#6Pi/16
sed -i '11s/.*/angle_phi 0.375/' options.dat
./SplitStep 6 >> rational6.dat

#7Pi/16
sed -i '11s/.*/angle_phi 0.4375/' options.dat
./SplitStep 6 >> rational7.dat

#rand1
sed -i '11s/.*/angle_phi 0.07343529162169780/' options.dat
./SplitStep 6 >> rand1.dat

#rand2
sed -i '11s/.*/angle_phi 0.1355479824884209/' options.dat
./SplitStep 6 >> rand2.dat

#rand3
sed -i '11s/.*/angle_phi 0.1798402032968870/' options.dat
./SplitStep 6 >> rand3.dat

#rand4
sed -i '11s/.*/angle_phi 0.2165189820519779/' options.dat
./SplitStep 6 >> rand4.dat

#rand5
sed -i '11s/.*/angle_phi 0.2764071674820491/' options.dat
./SplitStep 6 >> rand5.dat

#rand6
sed -i '11s/.*/angle_phi 0.3864458249608038/' options.dat
./SplitStep 6 >> rand6.dat

#rand7
sed -i '11s/.*/angle_phi 0.4355641972973841/' options.dat
./SplitStep 6 >> rand7.dat

#golden ratio
sed -i '11s/.*/angle_phi 0.3819660112501052/' options.dat
./SplitStep 6 >> golden.dat


