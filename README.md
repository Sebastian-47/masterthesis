# master thesis - Sebastian Kinnewig

---------------------------------------------------
## 1.) Overview

This project belongs to belongs to my master thesis and is dedicated to compute the dynamics of the QWZ model.

---------------------------------------------------
## 2.) Dependencies

1. cmake

2. fftw3
   (for parallelizing it has to be installed with --enable-openmp)

3. openMP (optional only needed for parallelizing)

---------------------------------------------------
## 3.) Installation

1. mkdir bin && cd bin

2. cmake ../splitstep

3. make

---------------------------------------------------
## 4.) Installation if fftw was installed with --prefix

1. Change the header in the following files corresponding to example.cc 
  /splitstep/main.cc 
  /splitstep/arithmetics.cc 
  /splitstep/arithmetics.h
  /splitstep/wavefunction/wavefunction.cc
  /splitstep/wavefunction/wavefunction.h
  
2. Complete the installation pathes in CMakeLists_splitstep.txt and CMakeLists_wavefunction.txt

3. Replace the /splitstep/CMakeLists.txt with CMakeLists_splitstep.txt (just rename CMakeLists_splitstep.txt to CMakeLists.txt)

4. Replace the /splitstep/wavefunction/CMakeLists.txt with CMakeLists_wavefunction.txt (just rename CMakeLists_wavefunction.txt to CMakeLists.txt)

5. mkdir bin && cd bin

6. cmake ../splitstep

7. make

