#include <iostream>
#include <chrono>
#include <fstream>
#include <functional>

//Libary: matrix
#include "../matrix/matrix.h"

//Libary: settings manager
#include "../settings_manager/settings_manager.h"

//Libary: OMP
#include <omp.h>

//Libary: FFTW (Fast Fourie Transformation)
#include <fftw3.h>

//Libary: Wavefunction
#include "wavefunction/wavefunction.h"

//this project
#include "arithmetics.h"

matrix flip_matrix(matrix m){
  matrix out(m.m(), m.n());
  for(int i = 0; i < m.n(); i++){
    for(int j = 0; j < m.m(); j++){
      out(j,i) = m(i,j);
    }
  }
  
  return out;
}

matrix flip_horizontal(matrix m){
  matrix out(m.n(), m.m());
  for(int i = 0; i < m.n(); i++){
    for(int j = 0; j < m.m(); j++){
      out(m.n() - i - 1, j) = m(i,j);
    }
  }
  
  return out;
}

matrix flip_vertical(matrix m){
  matrix out(m.n(), m.m());
  for(int i = 0; i < m.n(); i++){
    for(int j = 0; j < m.m(); j++){
      out(i, m.m() - j - 1) = m(i,j);
    }
  }
  
  return out;
}
//help function needed to read the initial state from a file
int count_lines(const char* name){
  std::string line;
  std::ifstream file(name);
  int number_of_lines = 0;
  if(file.is_open()){
    //count the number of lines;
    while (getline(file, line)){
      number_of_lines++;
    }
    file.close();
  }
  return number_of_lines;
}

matrix create_potential(settings_manager settings, double phi){
  /* ---------------- Ussage ----------------
   * phi:       angle of the stripe
   * size_n:    vertical hight of the potential
   * a:         width of the stripe
   * d:         hight between the lower edge of the potential and the lower edge of stripe
   * width:     width of the horizontal part of the stripe (in percent of the total width)
   * per_node:  defines the numer of pixel per node (I don't know why I added the possibility to modify it, "1" is totaly fine and I'm way to lazy to remove it)
   */
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  double left_right = settings.get_double("left_right");

  //the hight of the field:
  //double h = size_n * 1.0;
  //the width of the field:
  double b = size_n * std::cos(phi) / std::sin(phi);
  double f = b * std::tan(phi);
  double slope = size_n / b;

  //The left_right region
  int size_m = b;
  left_right += (b - (size_m)) / 2;
  //add more nodes in the left right region
  //left_right += a;
  size_m += 2 * left_right;

  std::function<double(double)> line_1 = [slope, a, d, left_right](double x){
    if(x < left_right){
      return d - a;
    }
    double x_var = x - left_right;
    return (slope * x_var) - a + d;
  }; 
  std::function<double(double)> line_2 = [slope, d, left_right](double x){
    if(x < left_right){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d;
  }; 
  std::function<double(double)> line_3 = [slope, a, b, d, f, left_right](double x){
    if(x >= left_right + b){
      return d - a;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d - a - f;
  }; 
  std::function<double(double)> line_4 = [slope, b, d, f, left_right](double x){
    if(x >= left_right + b){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d - f;
  }; 

  matrix potential(size_n, size_m);
  for(int y_pos = 0; y_pos < size_n; y_pos++){
    for(int x_pos = 0; x_pos < size_m; x_pos++){
      if(y_pos < line_2(x_pos) && y_pos >= line_1(x_pos)){
        potential(y_pos, x_pos) = 0;
      } else if(y_pos < line_4(x_pos) && y_pos >= line_3(x_pos)){
        potential(y_pos, x_pos) = 0;
      } else {
        potential(y_pos, x_pos) = 255;
      }
    }
  }

  return potential;
}

matrix create_potential(settings_manager settings){
  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int d = settings.get_double("lower_domain");
  double a = settings.get_double("line_width");

  matrix potential(size_n, size_m);
  for(int x_pos = 0; x_pos < size_m; x_pos++){
    for(int y_pos = 0; y_pos < size_n; y_pos++){
      if(y_pos < d || y_pos >= d + a){
        potential(y_pos, x_pos) = 255;
      } else {
        potential(y_pos, x_pos) = 0;
      }
    }
  }

  return potential;
}

wavefunction place_psi(settings_manager settings, wavefunction wf){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  double a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  //place psi on the grid
  int start_x = d + a;
  int start_y = (wf.return_n() / 2) - (psi_size_m / 2);
  wf.reset();
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      wf(0, i + start_y, j + start_x) = complex_d(psi_start(2 * j, 2 * i), psi_start(2 * j, 2 * i + 1));  
      wf(1, i + start_y, j + start_x) = complex_d(psi_start(2 * j + 1, 2 * i), psi_start(2 * j + 1,  2 * i + 1));
    }
  }

  return wf;
}

wavefunction place_psi(settings_manager settings, wavefunction wf, double phi){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  int d = settings.get_double("lower_domain");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
 
  int size_m = size_n * std::cos(phi) / std::sin(phi);
  double left_right = settings.get_double("left_right");
  size_m += 2 * left_right;

  //place psi on the grid
  int start_x = size_m - psi_size_m;
  wf.reset();
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      wf(0, i + start_x, j + d) = complex_d(psi_start(2 * j, 2 * i), psi_start(2 * j, 2 * i + 1));  
      wf(1, i + start_x, j + d) = complex_d(psi_start(2 * j + 1, 2 * i), psi_start(2 * j + 1,  2 * i + 1));
    }
  }

  return wf;
}

wavefunction place_psi(settings_manager settings, wavefunction wf, matrix potential){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  //place psi on the grid
  int start_x = (wf.return_n() / 2);
  int start_y = 0;
  for(int i = wf.return_m() - 1; i >= 0; i--){
    if(potential(start_x, i) == 0){
      start_y = i + 1;
      break;
    }
  }

  wf.reset();
  start_x -=  (psi_size_m / 2);
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      wf(0, start_x + i, start_y + j) = complex_d(psi_start(2 * j, 2 * i), psi_start(2 * j, 2 * i + 1));  
      wf(1, start_x + i, start_y + j) = complex_d(psi_start(2 * j + 1, 2 * i), psi_start(2 * j + 1,  2 * i + 1));
    }
  }

  return wf;
}

wavefunction place_psi_rotated(settings_manager settings, wavefunction wf, double phi){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  double d = settings.get_double("lower_domain");
  double a = settings.get_double("line_width");
  double left_right = settings.get_double("left_right");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  //compute the line the which psi has to follow
  double h = size_n;
  double b = h * std::cos(phi) / std::sin(phi);
  double f = b * std::tan(phi);
  double slope = h / b;
  int size_m = b;
  left_right += (b - size_m) / 2;
  size_m += 2 * left_right;

  std::function<double(double)> line_2 = [slope, d, left_right, size_m](double x){
    if(x < left_right){
      return d;
    }
    double x_var = x - left_right;
    if( slope  * x_var + d > size_m ){ 
      return (slope * x_var) + d - size_m;
    }
      return (slope * x_var) + d;
  }; 

  //place psi on the grid
  wf.reset();
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      int x_pos = i * cos(phi) + j * sin(phi) + left_right;
      int y_pos = j * cos(phi) - i * sin(phi) + d;
      wf(0, x_pos, y_pos) = complex_d(psi_start(2 * j, 2 * i), psi_start(2 * j, 2 * i + 1));  
      wf(1, x_pos, y_pos) = complex_d(psi_start(2 * j + 1, 2 * i), psi_start(2 * j + 1,  2 * i + 1));
    }
  }

  return wf;
}

wavefunction clean_top(settings_manager settings, wavefunction wf, matrix potential){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  //place psi on the grid
  int start_x = (wf.return_n() / 2);
  int start_y = 0;
  for(int i = wf.return_m() - 1; i >= 0; i--){
    if(potential(start_x, i) == 0){
      start_y = i + 1;
      break;
    }
  }

  for(int i = 0; i < wf.return_n(); i++){
    for(int j = start_y; j < wf.return_m(); j++){
      wf(0, i, j) = complex_d(0.0, 0.0);
      wf(1, i, j) = complex_d(0.0, 0.0);
    }
  }

  return wf;
}

wavefunction clean_top(settings_manager settings, wavefunction wf){
  double a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  std::string saving_path_psi = settings.get_string("saving_path_psi");

  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  int start_y = (wf.return_n() / 2) -  (psi_size_m / 2);
  int start_x = d + a + 2;

  for(int i = start_y; i < wf.return_n(); i++){
    for(int j = 0; j < wf.return_m(); j++){
      wf(0, i, j) = complex_d(0,0);
      wf(1, i, j) = complex_d(0,0);
    }
  }
  for(int i = 0; i < start_y; i++){
    for(int j = start_x; j < wf.return_m(); j++){
      wf(0, i, j) = complex_d(0,0);
      wf(1, i, j) = complex_d(0,0);
    }
  }

  return wf;
}

wavefunction clean_bottom(settings_manager settings, wavefunction wf, double phi, int edge){
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  int d = settings.get_double("lower_domain");
  
  int size_m = size_n * std::cos(phi) / std::sin(phi);
  double left_right = settings.get_double("left_right");
  size_m += 2 * left_right;

  //place psi on the grid
  for(int i = 0; i < size_m; i++){
    for(int j = 0; j < edge; j++){
      wf(0, i, j + d) = complex_d(0.0, 0.0);
      wf(1, i, j + d) = complex_d(0.0, 0.0);
    }
  }

  return wf;
}

void save_wavefunction(settings_manager settings, wavefunction wf, const char* name, double phi){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());

  matrix save(psi_size_n * 2, psi_size_m * 2);

  int start_x = d;
  int start_y = 0;
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      save(2 * j, 2 * i) = wf(0, start_y + i, start_x + j).real();
      save(2 * j, 2 * i + 1) = wf(0, start_y + i, start_x + j).imag();
      save(2 * j + 1, 2 * i) = wf(1, start_y + i, start_x + j).real();
      save(2 * j + 1, 2 * i + 1) = wf(1, start_y + i, start_x + j).imag();
    }
  }

  save.to_file(name);
}

wavefunction reset_boundary(settings_manager settings, wavefunction wf, int edge){
  double a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  std::string saving_path_psi = settings.get_string("saving_path_psi");

  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());
  
  int start_y = (wf.return_n() / 2) -  (psi_size_m / 2);
  int start_x = d + a + edge;

  for(int i = 0; i < wf.return_n(); i++){
    for(int j = start_x; j < wf.return_m(); j++){
      wf(0, i, j) = complex_d(0,0);
      wf(1, i, j) = complex_d(0,0);
    }
  }

  return wf;
}

wavefunction reset_boundary(settings_manager settings, wavefunction wf, int edge, double phi){
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  double left_right = settings.get_double("left_right");
  double per_node = 1.0;

  //the hight of the field:
  double h = size_n; 
  //the width of the field:
  double b = h * std::cos(phi) / std::sin(phi);
  double f = b * std::tan(phi);
  double slope = h / b;

  //The left_right region
  int size_m = b / per_node;
  left_right += (b - (size_m)) / 2;
  //add more nodes in the left right region
  //left_right += a;
  size_m += 2 * left_right;

  std::function<double(double)> line_2 = [slope, d, left_right](double x){
    if(x < left_right){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d;
  }; 
  std::function<double(double)> line_4 = [slope, b, d, f, left_right](double x){
    if(x >= left_right + b){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d - f;
  }; 

  //
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      if(i < line_2(j) + edge && i >= line_2(j)){
        //wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        //wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
      } else if(i < line_4(j) + edge && i >= line_4(j)){
        //wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        //wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
      } else {
        wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
      }
    }
  }
  
  return wf;
}

wavefunction reset(settings_manager settings, wavefunction wf, double phi){
  //reading psi from file
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  
  int psi_size_n = count_lines(saving_path_psi.c_str()) / 2;
  int psi_size_m = settings.get_double("psi_size_m");
  matrix psi_start(2 * psi_size_n, 2 * psi_size_m, saving_path_psi.c_str());

  //save the parts of psi we want to keep
  complex_d*  save_0 = new complex_d[psi_size_n * psi_size_m];
  complex_d*  save_1 = new complex_d[psi_size_n * psi_size_m];
  int start_x = d;
  int start_y = 0;
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      save_0[(i * psi_size_n) + j] = wf(0, start_y + i, start_x + j);
      save_1[(i * psi_size_n) + j] = wf(1, start_y + i, start_x + j);
    }
  }

  //delete everything
  wf.reset();

  //restore what we want to keep
  for(int i = 0; i < psi_size_m; i++){
    for(int j = 0; j < psi_size_n; j++){
      wf(0, start_y + i, start_x + j) = save_0[(i * psi_size_n) + j]; 
      wf(1, start_y + i, start_x + j) = save_1[(i * psi_size_n) + j]; 
    }
  }

  //renormalize the output
  wf.renormalize();
  
  delete[] save_0;
  delete[] save_1;

  return wf;
}

double bulk_loss(settings_manager settings, wavefunction wf, int edge, double phi){
  int size_n = settings.get_double("size_n");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");
  double left_right = settings.get_double("left_right");
  double per_node = 1.0;

  //the hight of the field:
  double h = size_n; 
  //the width of the field:
  double b = h * std::cos(phi) / std::sin(phi);
  double f = b * std::tan(phi);
  double slope = h / b;

  //The left_right region
  int size_m = b / per_node;
  left_right += (b - (size_m)) / 2;
  size_m += 2 * left_right;

  std::function<double(double)> line_2 = [slope, d, left_right](double x){
    if(x < left_right){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d;
  }; 
  std::function<double(double)> line_4 = [slope, b, d, f, left_right](double x){
    if(x >= left_right + b){
      return d;
    }
    double x_var = x - left_right;
    return (slope * x_var) + d - f;
  }; 

  //wavefunction bulk_wf = wf;
  wf.renormalize();
  double norm = 0.0;
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      if(i < line_2(j) + edge && i >= line_2(j)){
        //bulk_wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        //bulk_wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
      } else if(i < line_4(j) + edge && i >= line_4(j)){
        //bulk_wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        //bulk_wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
      } else {
        //edge_wf(0, size_m - j - 1, i) = complex_d(0.0, 0.0);
        //edge_wf(1, size_m - j - 1, i) = complex_d(0.0, 0.0);
        norm += std::pow(std::abs(wf(0, size_m - j - 1, i)),2);
        norm += std::pow(std::abs(wf(1, size_m - j - 1, i)),2);
      }
    }
  }
  
  //return bulk_wf.l2_norm() / wf.l2_norm();
  return norm;
}

double bulk_loss(settings_manager settings, wavefunction wf, int edge){
  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");

  int start = a + d + edge;

  wf.renormalize();
  double norm = 0.0;
  for(int i = start; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      norm += std::pow(std::abs(wf(0,j,i)), 2);
      norm += std::pow(std::abs(wf(1,j,i)), 2);
    }
  }
  
  return norm;
}

double bulk_loss(settings_manager settings, wavefunction wf, int edge, matrix potential){
  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int a = settings.get_double("line_width");
  double d = settings.get_double("lower_domain");

  matrix saved = potential;
  for(int i = 0; i < saved.n(); i++){
    for(int j = 0; j < saved.m(); j++){
      if(potential(i,j) == 0){
        saved(i,j) = 0;
      } else {
        saved(i,j) = 1;
        for(int rel_i = -edge; rel_i <= edge; rel_i++){
          for(int rel_j = -edge; rel_j <= edge; rel_j++){
            if(potential((i + rel_i + saved.n()) % saved.n(), (j + rel_j + saved.m()) % saved.m()) == 0){
              saved(i,j) = 0;
              break;
            }
          }
        }
      }

    }
  }


  wf.renormalize();
  double norm = 0.0;
  for(int i = 0; i < saved.n(); i++){
    for(int j = 0; j < saved.m(); j++){
      if(saved(i,j) == 1){
        norm += std::pow(std::abs(wf(0,i,j)), 2);
        norm += std::pow(std::abs(wf(1,i,j)), 2);
      }
    }
  }
  
  return norm;
}


/*
 * starting by some eigenstate psi (given by psi_start) of the unperturbed potential, "optimize" tries to find
 * an eigenstate of the given potential with a perturbation
 */
void optimize(int steps, settings_manager settings, bool save_result){
  /* ---------------- Ussage ----------------
   * phi:         angle of the stripe
   * size_n:       vertical hight of the potential
   * width:       width of the stripe
   * iterations:  number of iterations done to find a stable 
   * steps:       number of time steps done in each iteration step
   * dt:          size of the single time time steps
   * frame_rate:  frame for the final result
   * saving_name: Name of for the PNG-image of the potential (if empty no image will be saved)
   */
  int size_n = settings.get_double("size_n");
  int iterations = settings.get_double("iterations");
  int frame_rate = settings.get_double("frame_rate");
  double dt = settings.get_double("time_step_size");
  double uu = settings.get_double("perturbation_v");

  std::string saving_path_psi = settings.get_string("saving_path_psi");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");
  std::string method = settings.get_string("method");

  //create the potential
  std::cout << "creating the test potential...";
  double phi = settings.get_double("angle_phi");
  matrix pot = create_potential(settings, phi);
  pot = flip_matrix(pot);
  pot = flip_horizontal(pot);
  std::cout << "DONE" << std::endl;

  //saving the potential as a PNG image
  if(saving_path_image!= ""){
    std::string temp = saving_path_image + ".png";
    pot.to_png(temp.c_str());
  }

  //Why ever the image is flipped in the first place... but okay, that is not realy complicated
  //matrix image = flip_matrix(potential);
  wavefunction wf(pot.n(), pot.m(), pot.return_data());
  wf = place_psi(settings, wf, phi);

  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());
  int fftwthreads = fftw_init_threads(); 
  int max_threads = 1;
  if(fftwthreads == 0){
    std::cout << "There is an error with fftw openMP" << std::endl;
  } else {
    max_threads = omp_get_max_threads();
    std::cout << "Maximal number of threads: " << max_threads << std::endl;
    std::cout << "Using " << max_threads << " threads for the computation" << std::endl;
  }

  std::string direct = "direct";
  std::string lie = "lie";
  std::string strang = "strang";
  for(int i = 0; i < steps; i++){
    std::cout << "Step " << i << "/" << steps << std::endl;
    wf.timestep(potential, method, dt, iterations, 1.0, uu, true, max_threads);
    wf = reset_boundary(settings, wf, 3, phi);
    wf = clean_bottom(settings, wf, phi, 3);
    wf.renormalize();
    wf.timestep(potential, method, dt, iterations, -1.0, uu, true, max_threads);
    wf = reset(settings, wf, phi);
    wf.renormalize();
  }
  std::cout << "Step " << steps << "/"  << steps << std::endl;
  wf.timestep(potential, method, dt, iterations, 1.0, uu, true, max_threads);

  if(save_result){
    std::string save = "optimize.dat";
    save_wavefunction(settings, wf, save.c_str(), phi);
  } else {
    std::string awnser;
    std::cout << "Do you want to save the result? (y/n) << ";
    std::cin >> awnser;
    std::string yes = "yes";
    std::string y = "y";
    if(awnser == yes || awnser == y){
      std::cout << "Please enter the saving path << ";
      std::string save;
      std::cin >> save;
      
      save_wavefunction(settings, wf, save.c_str(), phi);
    }
  }
}

void test_potential(settings_manager settings){
  const double PI = 3.141592653589793;

  int size_n = settings.get_double("size_n");
  int iterations = settings.get_double("iterations");
  int frame_rate = settings.get_double("frame_rate");
  double dt = settings.get_double("time_step_size");
  double phi = PI * settings.get_double("angle_phi");
  double d = settings.get_double("lower_domain");
  double uu = settings.get_double("perturbation_v");

  std::string saving_path_psi = settings.get_string("saving_path_psi");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");
  std::string method = settings.get_string("method");

  bool reduce = false;
  std::string temp = "true";
  if( settings.get_string("reduce_noise") == temp){
    reduce = true;
  }

  //create the potential
  std::cout << "creating the test potential...";
  matrix pot = create_potential(settings, phi);
  pot = flip_matrix(pot);
  pot = flip_horizontal(pot);
  std::cout << "DONE" << std::endl;

  //saving the potential as a PNG image
  if(saving_path_image!= ""){
    std::string temp = saving_path_image + ".png";
    pot.to_png(temp.c_str());
  }

  wavefunction wf(pot.n(), pot.m(), pot.return_data());
  wf = place_psi(settings, wf, phi);
  //wf = place_psi_rotated(settings, wf, phi);
  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());
  int fftwthreads = fftw_init_threads(); 
  int max_threads = 1;
  if(fftwthreads == 0){
    std::cout << "There is an error with fftw openMP" << std::endl;
  } else {
    max_threads = omp_get_max_threads();
    std::cout << "Maximal number of threads: " << max_threads << std::endl;
    std::cout << "Using " << max_threads << " threads for the computation" << std::endl;
  }

  if(reduce){
    int time_30 = 40.0 / dt;
    int iterations_part_2 = (iterations - time_30 < 0) ? 0 : (iterations - time_30);
    int iterations_part_1 = iterations - iterations_part_2;

    wf.timestep(potential, method, dt, iterations_part_1, 1.0, uu, true, max_threads);
    wf = clean_bottom(settings, wf, phi, 5);
    wf.reduce_noise(0.01);
    wf.renormalize();
    wf.timestep(potential, method, dt, iterations_part_2, 1.0, uu, true, max_threads);
  } else {
    wf.timestep(potential, method, dt, iterations, 1.0, uu, true, max_threads);
  }

  wf.renormalize();
  std::cout << "The loss into the bulk is: " << std::endl;
  std::cout << "Edge \t  Loss" << std::endl;  
  for(int i = 1; i < 10; i++){
    wavefunction temp1 = wf;
    temp1 = reset_boundary(settings, temp1, i, phi);
    std::cout << i << "\t | " << 1 - temp1.l2_norm() << std::endl;
  }

}

void test_eigen_vector(settings_manager settings){
  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int iterations = settings.get_double("iterations");
  int frame_rate = settings.get_double("frame_rate");
  int a = settings.get_double("line_width");
  double dt = settings.get_double("time_step_size");
  double d = settings.get_double("lower_domain");
  double uu = settings.get_double("perturbation_v");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");
  std::string method = settings.get_string("method");
 
  bool reduce = false;
  std::string temp = "true";
  if( settings.get_string("reduce_noise") == temp){
    reduce = true;
  }

  //create the potential
  std::cout << "Creating the test potential...";
  matrix pot = create_potential(settings);
  pot = flip_matrix(pot);
  std::cout << "DONE" << std::endl;

  //saving the potential as a PNG image
  if(saving_path_image!= ""){
    std::string temp = saving_path_image + ".png";
    pot.to_png(temp.c_str());
  }
  
  //create the wavefunction
  wavefunction wf(pot.n(), pot.m(), pot.return_data());
  wf = place_psi(settings, wf);
  
  /*
  wf.reset();
  wf(0, size_n / 2, d + a) = complex_d(1.0 / sqrt(2), 0.0);
  wf(0, size_n / 2, d + a) = complex_d(0.0, 1.0 / sqrt(2));
  */
 
  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());

  if(reduce){
    int time_30 = 30.0 / dt;
    int iterations_part_2 = (iterations - time_30 < 0) ? 0 : (iterations - time_30);
    int iterations_part_1 = iterations - iterations_part_2;

    wf.timestep(potential, method, dt, iterations_part_1, 1.0, uu, true);
    wf = clean_top(settings, wf);
    wf.reduce_noise(0.02);
    wf.renormalize();
    wf.timestep(potential, method, dt, iterations_part_2, 1.0, uu, true);
  } else {
    wf.timestep(potential, method, dt, iterations, 1.0, uu, true);
  }
  
  wf.renormalize();
  std::cout << "The loss into the bulk is: " << std::endl;
  std::cout << "Edge \t | Loss" << std::endl;  
  for(int i = 1; i < 10; i++){
    wavefunction temp1 = wf;
    temp1 = reset_boundary(settings, temp1, i);
    std::cout << i << "\t | "  <<  1 - temp1.l2_norm()  << std::endl;
  }

}

//creates an potential from an image and performs the time evolution
void from_image(settings_manager settings, const char* path){
  std::string saving_path_psi = settings.get_string("saving_path_psi");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");
  std::string method = settings.get_string("method");
  int iterations = settings.get_double("iterations");
  int frame_rate = settings.get_double("frame_rate");
  double dt = settings.get_double("time_step_size");
  double uu = settings.get_double("perturbation_v");

  bool reduce = false;
  std::string temp = "true";
  if( settings.get_string("reduce_noise") == temp){
    reduce = true;
  }

  matrix pot = image_to_matrix(path);
  pot = flip_matrix(pot);
  pot = flip_vertical(pot);

  wavefunction wf(pot.n(), pot.m(), pot.return_data());
  //pot = flip_vertical(pot);
  wf = place_psi(settings, wf, pot);

  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());
  int fftwthreads = fftw_init_threads(); 
  int max_threads = 1;
  if(fftwthreads == 0){
    std::cout << "There is an error with fftw openMP" << std::endl;
  } else {
    max_threads = omp_get_max_threads();
    std::cout << "Maximal number of threads: " << max_threads << std::endl;
    std::cout << "Using " << max_threads << " threads for the computation" << std::endl;
  }

  if(reduce){
    int time_30 = 30.0 / dt;
    int iterations_part_2 = (iterations - time_30 < 0) ? 0 : (iterations - time_30);
    int iterations_part_1 = iterations - iterations_part_2;

    wf.timestep(potential, method, dt, iterations_part_1, 1.0, uu, true, max_threads);
    wf = clean_top(settings, wf, pot);
    wf.reduce_noise(0.02);
    wf.renormalize();
    wf.timestep(potential, method, dt, iterations_part_2, 1.0, uu, true, max_threads);
  } else {
    wf.timestep(potential, method, dt, iterations, 1.0, uu, true, max_threads);
  }

  wf.renormalize();
  std::cout << "The loss into the bulk is: " << bulk_loss(settings, wf, 2, pot) << std::endl;


}

void stability_test(settings_manager settings){
  settings.set("size_n", 80);
  settings.set("size_m", 80);
  settings.set("psi_size_n", 40);
  settings.set("psi_size_m", 40);
  settings.set("saving_path_psi", "eigen_vector.dat");
  int iterations = settings.get_double("iterations");
  double dt = settings.get_double("time_step_size");
  double uu = settings.get_double("perturbation_v");
  int frame_rate = settings.get_double("frame_rate");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = "qwz";
  std::string method = settings.get_string("method");

  int size_n = 80;
  int size_m = 80;

  //create the potential
  std::cout << "Creating the test potential...";
  matrix image(size_n, size_m); //create_potential(settings);
  for(int  i = 0; i < size_n; i++){
    for(int j = 0; j < 20; j++){
      image( i, j ) = 0;
    }
    for(int j = 20; j < 60; j++){
      if(i < 40){
        image( i, j ) = 255;
      } else {
        image( i, j ) = 0;
      }
    }
    for(int j = 60; j < size_m; j++){
      image( i, j ) = 0;
    }
  }
  std::cout << "DONE" << std::endl;

  //saving the potential as a PNG image
  if(saving_path_image!= ""){
    std::string temp = saving_path_image + ".png";
    image.to_png(temp.c_str());
  }
  
  //create the wavefunction
  wavefunction wf(image.n(), image.m(), image.return_data());
  wf = place_psi(settings, wf);

  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());

  wf.timestep(potential, method, dt, iterations, 1.0, uu, true);
}

/*
 * computes the magnitude of the error for the Lie-Trotter splitting and
 * for the Strang splitting
 */

void error_test(settings_manager settings, int number_of_tests, const char* saving_name){
  const double PI = 3.141592653589793;

  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int iterations = settings.get_double("iterations");
  int frame_rate = settings.get_double("frame_rate");
  
  double phi = PI * settings.get_double("angle_phi");
  double d = settings.get_double("lower_domain");
  double uu = settings.get_double("perturbation_v");

  std::string saving_path_psi = settings.get_string("saving_path_psi");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");

  //create the potential
  std::string qwz = "qwz";
  std::cout << "creating the test potential...";
  matrix pot(size_n, size_m);
  if(potential == qwz){
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        pot(i, j) = 255;
      }
    }
  } else {
    pot = create_potential(settings);
    pot = flip_matrix(pot);
  }
  std::cout << "DONE" << std::endl;

  wavefunction wf(pot.n(), pot.m(), pot.return_data());

  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());
  int fftwthreads = fftw_init_threads(); 
  int max_threads = 1;
  if(fftwthreads == 0){
    std::cout << "There is an error with fftw openMP" << std::endl;
  } else {
    max_threads = omp_get_max_threads();
    std::cout << "Maximal number of threads: " << max_threads << std::endl;
    std::cout << "Using " << max_threads << " threads for the computation" << std::endl;
  }

  double dt = 0.001;
  iterations = 1.0 / dt;
  std::string li = "lie";
  std::string st = "strang";
  std::string di = "discrete";

  std::cout << "Compute the reference:" << std::endl;
  //wf.reset(size_n / 2, size_m / 2); 
  wf.reset();
  wf = place_psi(settings, wf);
  wf.timestep(potential, st, dt, iterations, 1.0, uu, false, max_threads);
  wavefunction reference_st = wf;
  reference_st.renormalize();

  //wf.reset(size_n / 2, size_m / 2); 
  wf.reset();
  wf = place_psi(settings, wf);
  wf.timestep(potential, li, dt, iterations, 1.0, uu, false, max_threads);
  wavefunction reference_li = wf;
  reference_li.renormalize();

  //wf.reset(size_n / 2, size_m / 2); 
  wf.reset();
  wf = place_psi(settings, wf);
  wf.timestep(potential, di, dt, iterations, 1.0, uu, false);
  wavefunction reference_di = wf;
  reference_di.renormalize();

  auto start = std::chrono::high_resolution_clock::now();  
  auto end = std::chrono::high_resolution_clock::now();  
  std::chrono::duration<double, std::milli> executen_time_ms = start - end; 

  dt = settings.get_double("time_step_size");
  iterations = 1.0 / dt;
  matrix norm_diff(4, number_of_tests);
  matrix time_diff(4, number_of_tests);
  for(int counter = 0; counter < number_of_tests; counter++){
    norm_diff(0, counter) = dt;
    time_diff(0, counter) = dt;

    std::cout << "Test the Lie-Trotter Splitting, with step with: " << dt << std::endl;
    //wf.reset(size_n / 2, size_m / 2); 
    wf.reset();
    wf = place_psi(settings, wf);
    start = std::chrono::high_resolution_clock::now(); 
    wf.timestep(potential, li, dt, iterations, 1.0, uu, false, max_threads);
    end = std::chrono::high_resolution_clock::now(); 
    executen_time_ms = end - start;
    wf.renormalize();
    norm_diff(1, counter) = wf.l2_norm(reference_li);
    time_diff(1, counter) = executen_time_ms.count() / 1000.0;
    
    std::cout << "Test the Strang-Trotter Splitting, with step with: " << dt << std::endl;
    //wf.reset(size_n / 2, size_m / 2); 
    wf.reset();
    wf = place_psi(settings, wf);
    start = std::chrono::high_resolution_clock::now(); 
    wf.timestep(potential, st, dt, iterations, 1.0, uu, false, max_threads);
    end = std::chrono::high_resolution_clock::now(); 
    executen_time_ms = end - start;
    wf.renormalize();
    norm_diff(2, counter) = wf.l2_norm(reference_st);
    time_diff(2, counter) = executen_time_ms.count() / 1000.0;

    std::cout << "Test the direct method, with step with: " << dt << std::endl;
    //wf.reset(size_n / 2, size_m / 2); 
    wf.reset();
    wf = place_psi(settings, wf);
    start = std::chrono::high_resolution_clock::now(); 
    wf.timestep(potential, di, dt, iterations, 1.0, uu, false, max_threads);
    end = std::chrono::high_resolution_clock::now(); 
    executen_time_ms = end - start;
    wf.renormalize();
    norm_diff(3, counter) = wf.l2_norm(reference_di);
    time_diff(3, counter) = executen_time_ms.count() / 1000.0;
    
    dt = dt / 2.0;
    iterations = iterations * 2.0;
  }

  std::cout << "Norm Difference:" << std::endl;
  std::cout << "dt\t| lie-trotter\t| strang\t| direct" << std::endl;
  std::cout << "--------+---------------+---------------+------------" << std::endl;
  for(int counter = 0; counter < number_of_tests; counter++){
    std::cout << norm_diff(0, counter) << "\t| " << norm_diff(1, counter) << "\t| " << norm_diff(2, counter) << "\t| " << norm_diff(3, counter) << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Run time:" << std::endl;
  std::cout << "dt\t| lie-trotter\t| strang\t| direct" << std::endl;
  std::cout << "--------+---------------+---------------+------------" << std::endl;
  for(int counter = 0; counter < number_of_tests; counter++){
    std::cout << time_diff(0, counter) << "\t| " << time_diff(1, counter) << "\t| " << time_diff(2, counter) << "\t| " << time_diff(3, counter) << std::endl;
  }

  std::string temp = "_error.dat";
  std::string saving_name_error = saving_name + temp;
  temp = "_time.dat";
  std::string saving_name_time  = saving_name + temp;
  norm_diff.to_file(saving_name_error.c_str());
  time_diff.to_file(saving_name_time.c_str());
}


void find_stable(settings_manager settings, int steps){
  int size_n = settings.get_double("size_n");
  int size_m = settings.get_double("size_m");
  int frame_rate = settings.get_double("frame_rate");
  int a = settings.get_double("line_width");
  double dt = settings.get_double("time_step_size");
  double d = settings.get_double("lower_domain");
  double uu = settings.get_double("perturbation_v");
  std::string saving_path_image = settings.get_string("saving_path_image");
  std::string saving_path_csv = settings.get_string("saving_path_csv");
  std::string saving_name_csv = settings.get_string("saving_name_csv");
  std::string potential = settings.get_string("potential");
  std::string method = settings.get_string("method");
 
  //create the potential
  std::cout << "Creating the test potential...";
  matrix pot = create_potential(settings);
  pot = flip_matrix(pot);
  std::cout << "DONE" << std::endl;

  //saving the potential as a PNG image
  if(saving_path_image!= ""){
    std::string temp = saving_path_image + ".png";
    pot.to_png(temp.c_str());
  }
  
  //create the wavefunction
  wavefunction wf(pot.n(), pot.m(), pot.return_data());
  wf = place_psi(settings, wf);
  wf.saving_options(frame_rate, saving_path_csv.c_str(), saving_name_csv.c_str());

  double left = 0;
  double right = 2 * dt;
  for(int i = 0; i < steps; i++){
    std::cout << "New dt = " << dt << ", left = " << left << ", right = " << right << std::endl;
    dt = ((right - left) / 2) + left;
    int iterations = 30.0 / dt;

    wf.timestep(potential, method, dt, iterations, 1.0, uu, false);
    wf = clean_top(settings, wf);
    wf.renormalize();
    wf.timestep(potential, method, dt, iterations, 1.0, uu, false);
 
    double loss =  bulk_loss(settings, wf, 4);
    std::cout << "The loss into the bulk is: " << loss << std::endl;

    if(loss < 0.2){
      left = (right - left) / 2.0 + left;
    } else {
      right = (right - left) / 2.0 + left;
    }
    if(left == right){
      break;
    }

  }
  std::cout << "Final dt = " << dt << std::endl;
}
