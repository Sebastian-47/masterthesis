#include <iostream>

#include <complex>
#include <fstream>
#include <functional> 
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath> 

//Libary: FFTW (Fast Fourie Transformation)
#include <fftw3.h>

#include "wavefunction.h"

#define index(i,j) ((i)*size_m+(j))

//Contructor:
wavefunction::wavefunction():
  size_m(0),
  size_n(0),
  potential(NULL),
  psi_0(NULL),
  psi_1(NULL),
  //options for printing the wavefunction
  frame_rate(1),
  frame_number(0),
  saving_path(""),
  saving_name("")
{}

wavefunction::wavefunction(const wavefunction& wf){
  size_n = wf.size_n;
  size_m = wf.size_m;
  potential = new double[size_n * size_m];
  frame_rate = wf.frame_rate;
  frame_number = wf.frame_number;
  saving_path = wf.saving_path;
  saving_name = wf.saving_name;
  
  psi_0 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);
  psi_1 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);
  for(int i = 0; i < size_n * size_m; i++){
    potential[i] = wf.potential[i];
    psi_0[i] = wf.psi_0[i];
    psi_1[i] = wf.psi_1[i];
  }
}

wavefunction::wavefunction(int size_n_, int size_m_, double* potential_){
  //matrix image_matrix = image_to_matrix(name);
  size_n = size_n_;
  size_m = size_m_;
  potential = new double[size_n * size_m];
 
  psi_0 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);
  psi_1 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);

  frame_rate = 1;
  frame_number = 0;
  saving_path = "";
  saving_name = "";

  for(int i = 0; i < size_n * size_m; i++){
    potential[i] = potential_[i];
    psi_0[i] = 0;
    psi_1[i] = 0;
  }
  int half = ((size_n / 2) * size_m) + (size_m / 2);
  psi_0[half] = complex_d(1 / sqrt(2), 0);
  psi_1[half] = complex_d(1 / sqrt(2), 0);
}

//destructor:
wavefunction::~wavefunction(){
  delete[] potential;
  fftw_free(psi_0);
  fftw_free(psi_1);
}

//saving options:
void wavefunction::saving_options(int frame_rate_, std::string saving_path_, std::string saving_name_){
  frame_rate = frame_rate_;
  saving_path = saving_path_;
  saving_name = saving_name_;
}

void wavefunction::set_frame_rate(int frame_rate_){
  frame_rate = frame_rate_;
}

void wavefunction::set_saving_path(std::string saving_path_){
  saving_path = saving_path_;
}

void wavefunction::reset_frames(){
  frame_number = 0;
}

std::string wavefunction::generate_path(){
  std::stringstream strs;
  std::string out = saving_path;
  out += saving_name;
  out += ".csv.";
  strs.str("");
  strs << frame_number;
  out += strs.str();
  return out;
}

//functions:
void wavefunction::reset(){
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = 0;
    psi_1[i] = 0;
  }
}

void wavefunction::reset(int i, int j){
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = 0;
    psi_1[i] = 0;
  }
  psi_0[index(i,j)] = complex_d(1 / sqrt(2), 0);
  psi_1[index(i,j)] = complex_d(1 / sqrt(2), 0);
}

void wavefunction::lie_trotter_timestep_schroedinger(double dt, int steps, double direction, bool create_csv, int fftwthreads){
  //check if the files can be written 
  if(create_csv){
    if(saving_path == ""){
      std::cout << "No saving path defined, to save the ouput define a saving path!" << std::endl;
      create_csv = false;
    } 
    if(saving_name == ""){
      std::cout << "No saving name defined, to save the ouput define a saving name!" << std::endl;
      create_csv = false;
    }
  }

  //build the plan:
  //before we build the plan we need to save psi (creating a plan can delete the context of psi)
  complex_d* psi_0_backup = new complex_d[size_n * size_m];
  complex_d* psi_1_backup = new complex_d[size_n * size_m];
  for(int i = 0; i < size_n * size_m; i++){
    psi_0_backup[i] = psi_0[i]; 
    psi_1_backup[i] = psi_1[i]; 
  }
  //create the plans
  fftw_plan_with_nthreads(fftwthreads);
  std::cout << "Creating the fastest plan to do the fourie transformation ...";
  fftw_plan plan_forwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), +1, FFTW_PATIENT);
  fftw_plan plan_forwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), +1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  std::cout << "Creating the fastest plan to do the backward fourie transformation ...";
  fftw_plan plan_backwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), -1, FFTW_PATIENT);
  fftw_plan plan_backwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), -1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  //restore psi
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0_backup[i]; 
    psi_1[i] = psi_1_backup[i]; 
  }
  //clean
  delete[] psi_0_backup;  
  delete[] psi_1_backup;  

  //some constants needed during the computation
  const complex_d imag_i(0.0, 1.0); 
  double sqrt_mn = std::sqrt(size_m * size_n);
  double dt_n = dt / (size_n * size_n);
  double dt_m = dt / (size_m * size_m);
  double cos_dt = std::cos(dt);
  double sin_dt = std::sin(dt);
  int half_n = size_n / 2;
  int half_m = size_m / 2;
  complex_d psi_0_temp;
  complex_d psi_1_temp;
    
  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
          std::cout << "=";
        } else if(prog == pos){
          std::cout << ">";
        } else {
          std::cout << " ";
        }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }

    std::cout.flush();   
    std::cout << step_number << "\r";
    //Apply the potential (in position Space)
    for(int i = 0; i < size_n * size_m; i++){
      psi_0_temp = psi_0[i];
      psi_1_temp = psi_1[i];
      if(potential[i] == 255){
        psi_0[i] = (cos_dt * psi_0_temp) - (direction * imag_i * psi_1_temp * sin_dt);
        psi_1[i] = (- direction * imag_i * psi_0_temp * sin_dt) + (cos_dt * psi_1_temp);
      }else{
        psi_0[i] = (cos_dt * psi_0_temp) + (direction * imag_i * psi_1_temp * sin_dt);
        psi_1[i] = (direction * imag_i * psi_0_temp * sin_dt) + (cos_dt * psi_1_temp);
      }
    }

    //Fourie Transform: Position space -> Momentum space
    fftw_execute(plan_forwards_0);
    fftw_execute(plan_forwards_1);
      
    //Apply the gradient && and normalize the output (in momentum space)
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index(i,j)] = (psi_0[index(i,j)] * std::exp(direction * imag_i * ((dt_n * std::pow((i - half_n), 2)) + (dt_m * std::pow(j - half_m, 2))))) / sqrt_mn;
        psi_1[index(i,j)] = (psi_1[index(i,j)] * std::exp(direction * imag_i * ((dt_n * std::pow((i - half_n), 2)) + (dt_m * std::pow(j - half_m, 2))))) / sqrt_mn;
      }
    }
    
    //Fourie Transform: Momentum space -> Position Space
    fftw_execute(plan_backwards_0);
    fftw_execute(plan_backwards_1);

    //normalize the output
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }
  }
  std::cout << std::endl;
    
  //destroyes the plans
  fftw_destroy_plan(plan_forwards_0);
  fftw_destroy_plan(plan_forwards_1);
  fftw_destroy_plan(plan_backwards_0);
  fftw_destroy_plan(plan_backwards_1);
}

void wavefunction::lie_trotter_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv, int fftwthreads){
  //check if the files can be written 
  if(create_csv){
    if(saving_path == ""){
      std::cout << "No saving path defined, to save the ouput define a saving path!" << std::endl;
      create_csv = false;
    } 
    if(saving_name == ""){
      std::cout << "No saving name defined, to save the ouput define a saving name!" << std::endl;
      create_csv = false;
    }
  }

  //build the plan:
  //before we build the plan we need to save psi (creating a plan can delete the context of psi)
  complex_d* psi_0_backup = new complex_d[size_n * size_m];
  complex_d* psi_1_backup = new complex_d[size_n * size_m];
  for(int i = 0; i < size_n * size_m; i++){
    psi_0_backup[i] = psi_0[i]; 
    psi_1_backup[i] = psi_1[i]; 
  }
  //create the plans
  fftw_plan_with_nthreads(fftwthreads);
  std::cout << "Creating the fastest plan to do the fourie transformation ...";
  fftw_plan plan_forwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), +1, FFTW_PATIENT);
  fftw_plan plan_forwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), +1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  std::cout << "Creating the fastest plan to do the backward fourie transformation ...";
  fftw_plan plan_backwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), -1, FFTW_PATIENT);
  fftw_plan plan_backwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), -1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  //restore psi
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0_backup[i]; 
    psi_1[i] = psi_1_backup[i]; 
  }
  //clean
  delete[] psi_0_backup;  
  delete[] psi_1_backup;  

  //some constants needed during the computation
  const complex_d imag_i(0.0, 1.0); 
  const double E = 2.718281828459045; 
  const double PI = 3.141592653589793;
  double p1_pos;
  double p2_pos;
  double dt = direction * dt_in;
  double sqrt_mn = std::sqrt(size_m * size_n);
  int half_n = size_n / 2;
  int half_m = size_m / 2;
  complex_d psi_0_temp;
  complex_d psi_1_temp;
  

  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
        std::cout << "=";
      } else if(prog == pos){
        std::cout << ">";
      } else {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    //this->renormalize();
    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }

    std::cout.flush();   
    std::cout << step_number << "\r";

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index(i,j)] = exp( imag_i * (dt / 1.0) * uu) * psi_0[index(i,j)];
        psi_1[index(i,j)] = exp(- imag_i * (dt / 1.0) * uu) * psi_1[index(i,j)];
      }
    }

    //Fourie Transform: Position space -> Momentum space
    fftw_execute(plan_forwards_0);
    fftw_execute(plan_forwards_1);

    //normalize the output
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p1_pos = (2.0 * PI* (j - half_m)) / (size_m * 1.0);
        p1_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];

        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt) * ((1.0 + exp(2.0 * imag_i * dt)) * psi_0_temp - (exp(2.0 * imag_i * dt) - 1.0) * psi_0_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt) - 1.0) * psi_1_temp * sin(p1_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt) * ((1.0 + exp(2.0 * imag_i * dt)) * psi_1_temp + (exp(2.0 * imag_i * dt) - 1.0) * psi_1_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt) - 1.0) * psi_0_temp * sin(p1_pos));
      }
    }

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p2_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        p2_pos = (2.0 * PI * (j - half_m)) / (size_m * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];
        
        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt) * ((1.0 + exp(2.0 * imag_i * dt)) * psi_0_temp - (exp(2.0 * imag_i * dt) - 1.0) * psi_0_temp * cos(p2_pos) + imag_i * (exp(2.0 * imag_i * dt) - 1.0) * psi_1_temp * sin(p2_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt) * ((1.0 + exp(2.0 * imag_i * dt)) * psi_1_temp + (exp(2.0 * imag_i * dt) - 1.0) * psi_1_temp * cos(p2_pos) - imag_i * (exp(2.0 * imag_i * dt) - 1.0) * psi_0_temp * sin(p2_pos));
      }
    }

    //Fourie Transform: Momentum space -> Position Space
    fftw_execute(plan_backwards_0);
    fftw_execute(plan_backwards_1);

    //normalize the output
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }

    //apply the potential
    for(int i = 0; i < size_n * size_m; i++){
      if(potential[i] == 0){
        psi_0[i] = complex_d(0.0, 0.0); 
        psi_1[i] = complex_d(0.0, 0.0); 
      }
    }

  }
  std::cout << std::endl;

  //destroyes the plans
  fftw_destroy_plan(plan_forwards_0);
  fftw_destroy_plan(plan_forwards_1);
  fftw_destroy_plan(plan_backwards_0);
  fftw_destroy_plan(plan_backwards_1);
}

void wavefunction::strang_timestep_schroedinger(double dt, int steps, double direction, bool create_csv, int fftwthreads){
  //check if the files can be written 
  if(create_csv){
    if(saving_path == ""){
      std::cout << "No saving path defined, to save the ouput define a saving path!" << std::endl;
      create_csv = false;
    } 
    if(saving_name == ""){
      std::cout << "No saving name defined, to save the ouput define a saving name!" << std::endl;
      create_csv = false;
    }
  }

  // build the plan:
  //before we build the plan we need to save psi (creating a plan can delete the context of psi)
  complex_d* psi_0_backup = new complex_d[size_n * size_m];
  complex_d* psi_1_backup = new complex_d[size_n * size_m];
  for(int i = 0; i < size_n * size_m; i++){
    psi_0_backup[i] = psi_0[i]; 
    psi_1_backup[i] = psi_1[i]; 
  }
  //creating the plans  
  fftw_plan_with_nthreads(fftwthreads);
  std::cout << "Creating the fastest plan to do the fourie transformation ...";
  fftw_plan plan_forwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), +1, FFTW_PATIENT);
  fftw_plan plan_forwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), +1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  std::cout << "Creating the fastest plan to do the backward fourie transformation ...";
  fftw_plan plan_backwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), -1, FFTW_PATIENT);
  fftw_plan plan_backwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), -1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  //restoring psi
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0_backup[i]; 
    psi_1[i] = psi_1_backup[i]; 
  }
  //cleaning
  delete[] psi_0_backup;  
  delete[] psi_1_backup;  


  //some constants needed during the computation
  double PI = 3.14159265358979324;

  const complex_d imag_i(0.0, 1.0); 
  double sqrt_mn = std::sqrt(size_m * size_n);
  double dt_m =  dt / (size_m * size_m);
  double dt_n =  dt / (size_n * size_n);
  double cos_dt = std::cos(dt);
  double sin_dt = std::sin(dt);
  int half_n = size_n / 2;
  int half_m = size_m / 2;
  complex_d psi_0_temp;
  complex_d psi_1_temp;

  //Fourie Transform: Position space -> Momentum space
  fftw_execute(plan_forwards_0);
  fftw_execute(plan_forwards_1);
  
  //normalize the output
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0[i] / sqrt_mn;
    psi_1[i] = psi_1[i] / sqrt_mn;
  }

  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
          std::cout << "=";
        } else if(prog == pos){
          std::cout << ">";
        } else {
          std::cout << " ";
        }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    //Apply the first half of the gradient && and normalize the output (in momentum space)
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index(i,j)] = (psi_0[index(i,j)] * std::exp((direction * imag_i * ((dt_n * std::pow(i - half_n, 2)) + (dt_m * std::pow(j - half_m, 2)))) ));
        psi_1[index(i,j)] = (psi_1[index(i,j)] * std::exp((direction * imag_i * ((dt_n * std::pow(i - half_n, 2)) + (dt_m * std::pow(j - half_m, 2)))) ));
      }
    }

    //Fourie Transform: Momentum space -> Position Space
    fftw_execute(plan_backwards_0);
    fftw_execute(plan_backwards_1);
    //normalize the output
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }

    //Apply the potential (in position Space)
    for(int i = 0; i < size_n * size_m; i++){
      psi_0_temp = psi_0[i];
      psi_1_temp = psi_1[i];
      if(potential[i] == 255){
        psi_0[i] = (cos_dt * psi_0_temp) - (direction * imag_i * psi_1_temp * sin_dt);
        psi_1[i] = (- direction * imag_i * psi_0_temp * sin_dt) + (cos_dt * psi_1_temp);
      }else{
        psi_0[i] = (cos_dt * psi_0_temp) + (direction * imag_i * psi_1_temp * sin_dt);
        psi_1[i] = (direction * imag_i * psi_0_temp * sin_dt) + (cos_dt * psi_1_temp);
      }
    }
 
    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }

    //Fourie Transform: Position space -> Momentum space
    fftw_execute(plan_forwards_0);
    fftw_execute(plan_forwards_1);
    //normalize the output
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }
      
    //Apply the secound half of the gradient && and normalize the output (in momentum space)
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index(i,j)] = (psi_0[index(i,j)] * std::exp((direction * imag_i * ((dt_n * std::pow(i - half_n, 2)) + (dt_m * std::pow(j - half_m, 2)))) * 2.0 * PI ));
        psi_1[index(i,j)] = (psi_1[index(i,j)] * std::exp((direction * imag_i * ((dt_n * std::pow(i - half_n, 2)) + (dt_m * std::pow(j - half_m, 2)))) * 2.0 * PI ));
      }
    }
  }
  std::cout << std::endl;

  //Fourie Transform: Momentum space -> Position Space
  fftw_execute(plan_backwards_0);
  fftw_execute(plan_backwards_1);
  //normalize the output
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0[i] / sqrt_mn;
    psi_1[i] = psi_1[i] / sqrt_mn;
  }

  //destroyes the plans
  fftw_destroy_plan(plan_forwards_0);
  fftw_destroy_plan(plan_forwards_1);
  fftw_destroy_plan(plan_backwards_0);
  fftw_destroy_plan(plan_backwards_1);
}

void wavefunction::strang_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv, int fftwthreads){
  //check if the files can be written 
  if(create_csv){
    if(saving_path == ""){
      std::cout << "No saving path defined, to save the ouput define a saving path!" << std::endl;
      create_csv = false;
    } 
    if(saving_name == ""){
      std::cout << "No saving name defined, to save the ouput define a saving name!" << std::endl;
      create_csv = false;
    }
  }

  //build the plan:
  //before we build the plan we need to save psi (creating a plan can delete the context of psi)
  complex_d* psi_0_backup = new complex_d[size_n * size_m];
  complex_d* psi_1_backup = new complex_d[size_n * size_m];
  for(int i = 0; i < size_n * size_m; i++){
    psi_0_backup[i] = psi_0[i]; 
    psi_1_backup[i] = psi_1[i]; 
  }
  //create the plans
  fftw_plan_with_nthreads(fftwthreads);
  std::cout << "Creating the fastest plan to do the fourie transformation ...";
  fftw_plan plan_forwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), +1, FFTW_PATIENT);
  fftw_plan plan_forwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), +1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  std::cout << "Creating the fastest plan to do the backward fourie transformation ...";
  fftw_plan plan_backwards_0 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_0), reinterpret_cast<fftw_complex*>(psi_0), -1, FFTW_PATIENT);
  fftw_plan plan_backwards_1 = fftw_plan_dft_2d(size_n, size_m, reinterpret_cast<fftw_complex*>(psi_1), reinterpret_cast<fftw_complex*>(psi_1), -1, FFTW_PATIENT);
  std::cout << " done!" << std::endl;
  //restore psi
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0_backup[i]; 
    psi_1[i] = psi_1_backup[i]; 
  }
  //clean
  delete[] psi_0_backup;  
  delete[] psi_1_backup;  

  //some constants needed during the computation
  const complex_d imag_i(0.0, 1.0); 
  const double E = 2.718281828459045; 
  const double PI = 3.141592653589793;
  double p1_pos;
  double p2_pos;
  double dt = direction * dt_in;
  double dt_2 = dt / 2.0;
  double sqrt_mn = std::sqrt(size_m * size_n);
  int half_n = size_n / 2;
  int half_m = size_m / 2;
  complex_d psi_0_temp;
  complex_d psi_1_temp;
  
  //Fourie Transform: Position space -> Momentum space
  fftw_execute(plan_forwards_0);
  fftw_execute(plan_forwards_1);
  //normalize the output
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0[i] / sqrt_mn;
    psi_1[i] = psi_1[i] / sqrt_mn;
  }

  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
        std::cout << "=";
      } else if(prog == pos){
        std::cout << ">";
      } else {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    std::cout.flush();   
    std::cout << step_number << "\r";

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p1_pos = (2.0 * PI* (j - half_m)) / (size_m * 1.0);
        p1_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];

        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_0_temp - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * sin(p1_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_1_temp + (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * sin(p1_pos));
      }
    }

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p2_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        p2_pos = (2.0 * PI * (j - half_m)) / (size_m * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];
        
        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_0_temp - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * cos(p2_pos) + imag_i * (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * sin(p2_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_1_temp + (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * cos(p2_pos) - imag_i * (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * sin(p2_pos));
      }
    }

    //Fourie Transform: Momentum space -> Position Space
    fftw_execute(plan_backwards_0);
    fftw_execute(plan_backwards_1);

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index(i,j)] = exp( imag_i * (dt / 1.0) * uu) * psi_0[index(i,j)] / sqrt_mn;
        psi_1[index(i,j)] = exp(- imag_i * (dt / 1.0) * uu) * psi_1[index(i,j)] / sqrt_mn;
      }
    }

    //this->renormalize();
    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }

    //apply the potential
    for(int i = 0; i < size_n * size_m; i++){
    if(potential[i] == 0){
        psi_0[i] = complex_d(0.0, 0.0); 
        psi_1[i] = complex_d(0.0, 0.0); 
      }
    }

    //Fourie Transform: Position space -> Momentum space
    fftw_execute(plan_forwards_0);
    fftw_execute(plan_forwards_1);
    for(int i = 0; i < size_n * size_m; i++){
      psi_0[i] = psi_0[i] / sqrt_mn;
      psi_1[i] = psi_1[i] / sqrt_mn;
    }

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p2_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        p2_pos = (2.0 * PI * (j - half_m)) / (size_m * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];
        
        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_0_temp - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * cos(p2_pos) + imag_i * (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * sin(p2_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_1_temp + (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * cos(p2_pos) - imag_i * (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * sin(p2_pos));
      }
    }

    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //p1_pos = (2.0 * PI* (j - half_m)) / (size_m * 1.0);
        p1_pos = (2.0 * PI * (i - half_n)) / (size_n * 1.0);
        psi_0_temp = psi_0[index(i,j)];
        psi_1_temp = psi_1[index(i,j)];

        psi_0[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_0_temp - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * sin(p1_pos));
        psi_1[index(i,j)] = (1.0/2.0) * exp(- imag_i * dt_2) * ((1.0 + exp(2.0 * imag_i * dt_2)) * psi_1_temp + (exp(2.0 * imag_i * dt_2) - 1.0) * psi_1_temp * cos(p1_pos) - (exp(2.0 * imag_i * dt_2) - 1.0) * psi_0_temp * sin(p1_pos));
      }
    }

  }
  std::cout << std::endl;

  fftw_execute(plan_backwards_0);
  fftw_execute(plan_backwards_1);
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0[i] / sqrt_mn;
    psi_1[i] = psi_1[i] / sqrt_mn;
  }

  //destroyes the plans
  fftw_destroy_plan(plan_forwards_0);
  fftw_destroy_plan(plan_forwards_1);
  fftw_destroy_plan(plan_backwards_0);
  fftw_destroy_plan(plan_backwards_1);
}

void wavefunction::discrete_timestep_schroedinger(double dt_in, int steps, double direction, bool create_csv){
  double uu = -1.0;
  double dt = direction * dt_in;
  complex_d imag_i = complex_d(0.0, 1.0);

  complex_d* psi_0_temp = new complex_d[size_n *size_m];
  complex_d* psi_1_temp = new complex_d[size_n *size_m];

  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
        std::cout << "=";
      } else if(prog == pos){
        std::cout << ">";
      } else {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    for(int i = 0; i < size_n * size_m; i++){
      psi_0_temp[i] = psi_0[i];
      psi_1_temp[i] = psi_1[i];
    }

    double f = 2.0;
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        //(0,0) +cos(p1)+cos(p2)
        psi_0[index((i + 1) % size_n, j)]               += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index((i - 1 + size_n) % size_n, j)]      += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index(i, (j + 1) % size_m)]               += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index(i, (j - 1 + size_m) % size_m)]      += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];

        //(1,1) -cos(p1)-cos(p2)
        psi_1[index((i + 1) % size_n, j)]               += (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index((i - 1 + size_n) % size_n, j)]      += (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index(i, (j + 1) % size_m)]               += (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index(i, (j - 1 + size_m) % size_m)]      += (dt / (imag_i * f)) * psi_1_temp[index(i,j)];

        //subtract the rest from (0,0) and (1,1) [see discrete representation of the laplace operator]
        psi_0[index(i, j)]                              -= 1.0 * (dt / imag_i) * psi_0_temp[index(i,j)];
        psi_1[index(i, j)]                              -= 1.0 * (dt / imag_i) * psi_1_temp[index(i,j)];

        //apply the perturbation
        psi_1[index(i, j)]                              += uu * (dt / imag_i) * psi_0_temp[index(i,j)];
        psi_0[index(i, j)]                              += uu * (dt / imag_i) * psi_1_temp[index(i,j)];
      }
    }

    //apply the potential
    for(int i = 0; i < size_n * size_m; i++){
      if(potential[i] == 0){
        psi_0[i] = 0;
        psi_1[i] = 0;
      }
    }

    //renormalize the wavefunction
    this->renormalize();

    //eport the wavefunction to image
    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }
  }
  std::cout << std::endl;

}


void wavefunction::discrete_timestep_qwz(double dt_in, int steps, double direction, double uu, bool create_csv){
  double dt = direction * dt_in;
  complex_d imag_i = complex_d(0.0, 1.0);

  complex_d* psi_0_temp = new complex_d[size_n *size_m];
  complex_d* psi_1_temp = new complex_d[size_n *size_m];

  std::cout << "computing time development:" << std::endl;
  int barWidth = 70;
  double progress = 0.0;
  int pos = 0;
  for(int step_number = 0; step_number < steps; step_number++){
    //progress bar
    progress =  (1.0 * step_number) / (1.0 * steps);
    pos = barWidth * progress;
    std::cout << "[";
    for(int prog = 0; prog < barWidth; prog++){
      if(prog < pos){
        std::cout << "=";
      } else if(prog == pos){
        std::cout << ">";
      } else {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << "% \r";
    std::cout.flush();

    for(int i = 0; i < size_n * size_m; i++){
      psi_0_temp[i] = psi_0[i];
      psi_1_temp[i] = psi_1[i];
    }


    double f = 2.0;
    //timestep, apply D
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        psi_0[index((i + 1) % size_n, j)]               += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index((i - 1 + size_n) % size_n, j)]      += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index(i, (j + 1) % size_m)]               += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_0[index(i, (j - 1 + size_m) % size_m)]      += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];

        psi_1[index((i + 1) % size_n, j)]               -= (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index((i - 1 + size_n) % size_n, j)]      -= (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index(i, (j + 1) % size_m)]               -= (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_1[index(i, (j - 1 + size_m) % size_m)]      -= (dt / (imag_i * f)) * psi_1_temp[index(i,j)];

        psi_0[index((i + 1) % size_n, j)]               -= (dt / f) * psi_1_temp[index(i,j)];
        psi_0[index((i - 1 + size_n) % size_n, j)]      += (dt / f) * psi_1_temp[index(i,j)];
        psi_0[index(i, (j + 1) % size_m)]               -= (dt / (imag_i * f)) * psi_1_temp[index(i,j)];
        psi_0[index(i, (j - 1 + size_m) % size_m)]      += (dt / (imag_i * f)) * psi_1_temp[index(i,j)];

        psi_1[index((i + 1) % size_n, j)]               -= (dt / f) * psi_0_temp[index(i,j)];
        psi_1[index((i - 1 + size_n) % size_n, j)]      += (dt / f) * psi_0_temp[index(i,j)];
        psi_1[index(i, (j + 1) % size_m)]               += (dt / (imag_i * f)) * psi_0_temp[index(i,j)];
        psi_1[index(i, (j - 1 + size_m) % size_m)]      -= (dt / (imag_i * f)) * psi_0_temp[index(i,j)];

        psi_0[index(i, j)]                              -= 1.0 * (dt / imag_i) * psi_0_temp[index(i,j)];
        psi_1[index(i, j)]                              -= 1.0 * (dt / imag_i) * psi_1_temp[index(i,j)];

        psi_0[index(i, j)]                              += uu * (dt / imag_i) * psi_0_temp[index(i,j)];
        psi_1[index(i, j)]                              -= uu * (dt / imag_i) * psi_1_temp[index(i,j)];
      }
    }

    //apply the potential
    for(int i = 0; i < size_n * size_m; i++){
      if(potential[i] == 0){
        psi_0[i] = 0;
        psi_1[i] = 0;
      }
    }

    //renormalize the wavefunction
    this->renormalize();

    //eport the wavefunction to image
    if(create_csv){
      if(step_number % frame_rate == 0){
        this->to_csv(this->generate_path().c_str());
        frame_number++;
      }
    }
  }

  std::cout << std::endl;
}

void wavefunction::timestep(std::string potential, std::string method, double dt_in, int steps, double direction, double uu, bool create_csv, int fftwthreads){
  //potentials
  std::string qwz = "qwz";
  std::string schroedinger = "schroedinger";

  //methods
  std::string strang = "strang";
  std::string lie = "lie";
  std::string discrete = "discrete";

  //execute the corresponding method
  //to the sign of the direction: I messed up the sign of the matrix exponential and the easiest and safest way to fix this error is by using "-direction" instead of "direction"
  if(potential == qwz){
    if(method == strang){
      this->strang_timestep_qwz(dt_in, steps, -direction, uu, create_csv, fftwthreads);
    } else if(method == lie){
      this->lie_trotter_timestep_qwz(dt_in, steps, -direction, uu, create_csv, fftwthreads);
    } else if(method == discrete){
      this->discrete_timestep_qwz(dt_in, steps, direction, uu, create_csv);
    } else {
      std::cout << "Unknown method: " << method << std::endl;
      std::cout << "Available methods are: " << strang << ", " << lie << ", " << discrete << std::endl;
    }
  } else if(potential == schroedinger){
    if(method == strang){
      this->strang_timestep_schroedinger(dt_in, steps, direction, create_csv, fftwthreads);
    } else if(method == lie){
      this->lie_trotter_timestep_schroedinger(dt_in, steps, direction, create_csv, fftwthreads);
    } else if(method == discrete){
      this->discrete_timestep_schroedinger(dt_in, steps, direction, create_csv);
    } else {
      std::cout << "Unknown method: " << method << std::endl;
      std::cout << "Available methods are: " << strang << ", " << lie << ", " << discrete << std::endl;
    }
  } else {
      std::cout << "Unknown potential: " << potential << std::endl;
      std::cout << "Available potentials are: " << qwz << ", " << schroedinger << std::endl;
  }
}

void wavefunction::reduce_noise(double tol){
  for(int i = 0; i < size_n * size_m; i++){
    if(std::abs(psi_0[i]) + std::abs(psi_1[i]) < tol){
      psi_0[i] = complex_d(0.0, 0.0);
      psi_1[i] = complex_d(0.0, 0.0);
    }
  }
}

//normen:
double wavefunction::l2_norm() const {
  double norm = 0.0;
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      norm += std::pow(std::abs(psi_0[index(i,j)]), 2) + std::pow(std::abs(psi_1[index(i,j)]), 2);
    }
  }
  norm = std::sqrt(norm); 

  return norm;
}

double wavefunction::l2_norm(const wavefunction& wf) const {
  if(size_n != wf.return_n() || size_m != wf.return_m()){
    throw std::invalid_argument("Can not compute the norm diffrence, because the wavefunction have diffrent dimensions");
    return 0;
  }

  double norm = 0.0;
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      norm += std::pow(std::abs(psi_0[index(i,j)] - wf.return_psi_0(index(i,j))), 2) + std::pow(std::abs(psi_1[index(i,j)] - wf.return_psi_1(index(i,j))), 2);
    }
  }

  return std::sqrt(norm);
}

double wavefunction::l2_max() const {
  double max = 0.0;
  double temp;
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      temp = std::pow(std::abs(psi_0[index(i,j)]), 2) + std::pow(std::abs(psi_1[index(i,j)]), 2);
      if(max < temp){
        max = temp;
      }
    }
  }
  return sqrt(max);
}

void wavefunction::renormalize(){
  double norm = this->l2_norm();
  for(int i = 0; i < size_n * size_m; i++){
    psi_0[i] = psi_0[i] / norm;
    psi_1[i] = psi_1[i] / norm;
  }
}

//return functions:
int wavefunction::return_n() const {
  return size_n;
}

int wavefunction::return_m() const {
  return size_m;
}

complex_d wavefunction::return_psi_0(int i) const {
  return psi_0[i];
}

complex_d wavefunction::return_psi_1(int i) const {
  return psi_1[i];
}

void wavefunction::print() const {
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      std::cout << psi_0[index(i,j)] << " ";
    }
    std::cout << std::endl;
  }
}

void wavefunction::print_norm() const {
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      std::cout <<  std::sqrt(std::pow(std::abs(psi_0[index(i,j)]), 2.0) + std::pow(std::abs(psi_1[index(i,j)]),2)) << " ";
    }
    std::cout << std::endl;
  }
}

void wavefunction::to_csv(const char* name){
  std::ofstream read_data(name);
  std::stringstream strs;
  double temp;
  if(!read_data){
    std::cout << "ERROR: can't open output file: " << name << std::endl; 
  } else{
    std::string out = "x cord,y cord,z cord,scalar\n";
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        strs.str("");
        strs << i;
        out += strs.str();
        out +=  ",";

        strs.str("");
        strs << j;
        out += strs.str();
        out +=  ",";

        strs.str("");
        strs << 0;
        out += strs.str(); 
        out +=  ",";

        strs.str("");
        if(potential[index(i,j)] == 0){
          out += "Nan";
        } else {
          temp = std::sqrt(std::pow(std::abs(psi_0[index(i,j)]), 2.0) + std::pow(std::abs(psi_1[index(i,j)]),2));
          strs << temp;
          out += strs.str();
        }
        out += "\n";
      }
    }

    read_data << out;
    read_data.close();   
  }
}

//operators
wavefunction& wavefunction::operator=(const wavefunction& wf){
  if(this != &wf){

    size_n = wf.size_n;
    size_m = wf.size_m;

    delete[] potential;
    potential = new double[size_n * size_m];

    frame_rate = wf.frame_rate;
    frame_number = wf.frame_number;
    saving_name = wf.saving_name;
    saving_path = wf.saving_path;

    fftw_free(psi_0);
    fftw_free(psi_1);
    psi_0 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);
    psi_1 = (complex_d*) fftw_malloc(sizeof(complex_d) * size_n * size_m);

    for(int i = 0; i < size_n * size_m; i++){
      potential[i] = wf.potential[i];
      psi_0[i] = wf.psi_0[i];
      psi_1[i] = wf.psi_1[i];
    }
  }
  return *this;
}

//Zugriffsoperator fuer den Eintrag in Zeile i und Spalte j
complex_d wavefunction::operator()(int psi, int i, int j) const {
  if(psi == 0){
    return psi_0[index(i,j)]; 
  }
  else{
    return psi_1[index(i,j)]; 
  }
}

//Zugriffsoperator fuer den Eintrag in Zeile i und Spalte j
//Gibt Referenz auf diesen Eintrag zurueck
complex_d& wavefunction::operator()(int psi, int i, int j){
  if(psi == 0){
    return psi_0[index(i,j)]; 
  }
  else{
    return psi_1[index(i,j)]; 
  }
}
