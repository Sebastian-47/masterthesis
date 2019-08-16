#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>


#include <Eigen/Dense>

#include "colour.h"
#include "matrix.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define index(i,j) ((i)*size_m+(j))

//Constructor
//defult constructor
matrix::matrix():			
  size_m(0),
  size_n(0),
  data(NULL)
{}

//Copy constructor: (Rule of 3: Part 1 of 3)
matrix::matrix(const matrix& m){		
  size_n = m.size_n;
  size_m = m.size_m;
  data = new double[m.size_n *  m.size_m];
  for( int i = 0 ; i < size_m * size_n; i++ ){
    data[i] = m.data[i]; 
  }
}

//Constructor für eine 'n' x 'm' matrix (initialisiert mit 0)
matrix::matrix(int n, int m){
  size_n = n;
  size_m = m;
  data = new double[n * m];
  for(int i = 0; i < n * m; i++){
    data[i] = 0;
  }
}

//Constructor um eine 'n' x 'm' Matrix zu erstellen un dmit den Datan aus dem double array M_data_ zu initialisieren
matrix::matrix(int n, int m, double* data_){	
  size_n = n;
  size_m = m;
  data = new double[n * m];
  for( int i = 0 ; i < n * m ; i++ ){
    data[i] = data_[i]; 
  }
}

//Constructer um die Matrix aus einer Datei heraus zu lesen
matrix::matrix(int size_n_, int size_m_, const char* name){
  size_n = size_n_;
  size_m = size_m_;
  data = new double[size_n * size_m];

  std::ifstream read_data;                          			
  read_data.open(name, std::ios_base::in); 
  if(!read_data){    /*Fange den Fall ab, dass die Datei nicht geöffnet werden kann*/
    std::cout << BOLDRED << "ERROR: " << RESET << "can't open input file: " << name << std::endl;
    for(int i = 0; i < size_n * size_m; i++){
      data[i] = 0;
    } 
  } else {
    std::string line;	
  
    int k = 0;    /*Zeile*/
    int j;        /*Spalte*/
    int length;   /*Anzahl der Zeichen in der aktuellen Zeile (der Input Datei)*/
    bool found_number = false;
    while(getline(read_data, line)){                             /*lese "data" Zeile für Zeile:*/
      if(k == size_n){                                           /*Fange den Fall ab, dass die Input Datein mehr zeilen hat als die zu erstellebde Matrix*/
        break;
      }
      j = 0; 
      length = line.length();
      std::string temp_str;
      for(int position = 0; position < length; position++){             /*Laufe die Zeile Zeichen für Zeichen durch*/
        if(line[position] == '.' && found_number){                      /*Ein dezimal Punkt wurde in der Zahl gefunden*/
          temp_str += '.';
          position++;
        }
        if(line[position] == ',' && 48 <= line[position] && line[position] <= 57  && found_number){   /*Ein dezimal Komma wurde innheralb einer Zahl gefunden, dies wird als dezimal Punkt interpretiert*/
          temp_str += '.';
          position++;
        }
        if(48 <= line[position] && line[position] <= 57){               /*Das aktuelle Zeichen ist eine Zahl*/
          if(j == size_m){                                              /*Fange den Fall ab, dass die Zeile mehr Zahlen enthät, als die Matrix Spalten hat*/
            std::cout << BOLDRED << "Warning:" << RESET << "while reading the file " << BLUE << name << RESET << " found in line " << k << " more than " << j-1 
                      << " numbers! Please check the input file" << std::endl;
            break;
          }
          found_number = true;
          temp_str = temp_str + line[position];
        }
        if(found_number){ /*Überprüfe ob eine Zahl endet, falls eine Zahl endet*/
          if(48 > line[position] || line[position] > 57){
            data[index(k,j)] = std::stod(temp_str); 		
            temp_str = "";
            found_number = false;
            j++;
          }
        }     
        if(line[position] == '-'){                                      /*Anfang einer negativen Fall wurde gefunden*/
          if(48 > line[position] && line[position] > 57){               /*Fange den Fall ab, dass ein "wildes Minus" in der Zeile steht*/
            std::cout << BOLDRED << "Warning:" << RESET << "while reading the file " << BLUE << name << RESET << " I found in line " << k 
                      << " a minus sign, which is not in front of a number, I'm just going to ignore this minus sign" << std::endl; 
          } else{
            temp_str = temp_str + line[position];
            found_number = true;
          }
        }
        if(position == length - 1 && found_number){
            data[index(k,j)] = std::stod(temp_str); 		
            temp_str = "";
            found_number = false;
            j++;
        }
      }
      if(j != size_m){                                              /*Fange den Fall ab, dass weniger als M_dim_m Zahlen in der Zeile stehen*/
        std::cout << BOLDRED << "Warning:" << RESET << "found only " << j << " numbers (instead of " << size_m << ")! Filling the rest with 0. " 
                  << "Please check the input file" << std::endl;
        for(int fill = j; fill < size_m; fill++){
          data[index(k,fill)] = 0;
        }
      }
      k++;
    }
    if(k != size_n){                                            /*Fange den Fall ab, dass die input Datei weniger als M_dim_n Zeilen hat*/
      std::cout << BOLDRED << "Warning:" << RESET << "found only " << k << " lines (instead of " << size_n << ")! Filling the rest with 0. " 
                << "Please check the input file" << std::endl;
      for(int fill_n = k; fill_n < size_n; fill_n++){
        for(int fill_m = 0; fill_m < size_m; fill_m++){
          data[index(fill_n, fill_m)] = 0;
        }
      }
    }
    read_data.close();
  }
}

//destructor: (Rule of 3: Part 2 of 3)
matrix::~matrix(){	
  delete[] data;
}

//retrun options
int matrix::n() const{		/*return the number of lines*/
  return size_n;
}

int matrix::m() const{		/*retruns the number of columns*/
  return size_m;
}

int matrix::length() const{	/*returns the length of the array data*/
  return size_m * size_n;
}

void matrix::print() const{				//Prints the matrix to the console
  for(int i = 0; i < size_n; i++){
    for(int j = 0; j < size_m; j++){
      std::cout << data[index(i, j)] << "	";
    }
    std::cout << std::endl;
  }	
}

double* matrix::return_data() const{
  return data;
}

void matrix::to_file(const char* name) const{ 	//Prints the matrix to the file 'name'
  std::ofstream read_data(name);
  std::stringstream strs;
  if(!read_data){
    std::cout << BOLDRED << "ERROR: " << RESET << "can't open output file: " << name << std::endl; 
  } else{
    std::string out = "";
    for(int i = 0; i < size_n; i++){
      for(int j = 0; j < size_m; j++){
        strs.str("");
	strs  <<  std::setprecision(9) << std::showpoint << std::fixed <<  data[index(i, j)];
	out += strs.str();
        out +=  "	";
      }
      out += "\n";
    }

    read_data << out;
    read_data.close();   
  }
}

void matrix::to_csv(const char* name) const { 	//Prints the matrix to the file 'name'
  std::ofstream read_data(name);
  std::stringstream strs;
  if(!read_data){
    std::cout << BOLDRED << "ERROR: " << RESET << "can't open output file: " << name << std::endl; 
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
	strs << data[index(j, i)];
	out += strs.str();
        out += "\n";
      }
    }

    read_data << out;
    read_data.close();   
  }
}

void matrix::to_png(const char* name) const{
  int comp = 1;
  int stride_in_bytes = size_m;
  unsigned char* image_data = new unsigned char[size_n * size_m];
  for(int i = 0; i < size_n * size_m; i++){
    image_data[i] = data[i];
  }     

  stbi_write_png(name, size_m, size_n, comp, image_data, stride_in_bytes);

  delete[] image_data;

}

matrix image_to_matrix(const char* name){
  int width;
  int hight;
  int channels;
  int number_of_bytes_per_pixel = 1;

  unsigned char* data = stbi_load(name, &width, &hight, &channels, number_of_bytes_per_pixel);
  double* matrix_data = new double[(width * hight * number_of_bytes_per_pixel)];

  for(int i = 0; i < width * hight * number_of_bytes_per_pixel; i++){
    matrix_data[i] = data[i];
  }
  
  matrix out(hight, width * number_of_bytes_per_pixel, matrix_data);

  delete[] matrix_data;
  delete[] data;
  
  return out;
}

//Operatoren
//Zugriffsoperator fuer den Eintrag in Zeile i und Spalte j
double matrix::operator()(int i, int j) const {
  return data[index(i,j)]; 
}

//Zugriffsoperator fuer den Eintrag in Zeile i und Spalte j
//Gibt Referenz auf diesen Eintrag zurueck
double& matrix::operator()(int i,int j){
  return data[index(i,j)]; 
}

//operators:
//Zuweise operator (Rule of 3: Part 3 of 3)
matrix& matrix::operator=(const matrix& m){	
  if(this != &m){
    int length = m.size_n * m.size_m;
    double* new_data = new double[length];

    for(int i = 0; i < length; i++){
      new_data[i] = m.data[i];        
    }
    delete[] data;
    data = new_data;
    size_n = m.size_n;
    size_m = m.size_m;
  }
  return *this;
}

bool matrix::operator==(const matrix& m){	/*compares two matrices, return true when all enterys are the same*/
	if(m.size_n != size_n){
		return false;
	}
	if(m.size_m != size_m){
		return false;
	}
	for(int i = 0; i < size_m; i++){
	for(int j = 0; j < size_m; j++){
		if(data[(i * size_m) + j] != m.data[(i * m.size_m) + j]){
			return false;
		}
	}}
	return true;
}


//arithmetics
matrix matrix::operator%(matrix m){	/*the matrix product*/
	if(size_m != m.n()){
		std::cout << BOLDRED << "ERROR: " << RESET << "the matrices have incompatible shapes" << std::endl;

		return *this;
	}	
	
	double* out_data = new double[size_n * m.m()];
	for(int i = 0; i < size_n * m.m(); i++){
		out_data[i] = 0;
	}
	
	for(int k = 0; k < m.m(); k++){
	for(int j = 0; j < size_n; j++){
	for(int i = 0; i < size_m; i++){
		out_data[(m.m() * j) + k] += data[(size_m * j) + i] * m(i, k);
	}}}
	matrix out(size_n, m.m(), out_data);
	delete[] out_data;
	return out;
}

matrix matrix::operator*(const matrix& m){	/*the tensor product*/
	/*yupp, i'm that lazy and just using that gcc has return value optimization https://en.wikipedia.org/wiki/Return_value_optimization#Compiler_support */
	double* temp_data = new double[size_n * m.size_n * size_m * m.size_m];

	for(int i1 = 0; i1 < size_n; i1++){
	for(int j1 = 0; j1 < size_m; j1++){
		for(int i2 = 0; i2 < m.size_n; i2++){
		for(int j2 = 0; j2 < m.size_m; j2++){
			temp_data[((i2 + (m.size_n*i1)) * (size_m * m.size_m)) + (j2 + (m.size_m*j1))] = data[(i1*size_m) + j1] * m.data[(i2*m.size_m) + j2];
		}}	
	}}

	matrix out(size_n * m.size_n, size_m * m.size_m, temp_data);
	delete[] temp_data;

	return out;
}

matrix matrix::operator*(double d){		/*multiply every entery with d*/
	int size = size_n * size_m;
	double* temp_data = new double[size];
	for(int i = 0; i < size; i++){
		temp_data[i] = data[i] * d;
	}
	matrix out(size_n, size_m, temp_data);
	delete[] temp_data;
	return out;

}

matrix matrix::operator/(double d){		/*divide every entery with d*/
	if(d == 0){
		std::cout << BOLDRED << "ERROR: " << RESET << "you fool just divided by zero " << std::endl;
	}
	return *this * (1 / d);
}

matrix matrix::operator+(const matrix& m){
	if(m.size_n != size_n || m.size_m != size_m){
		std::cout << BOLDRED << "ERROR: " << RESET << "can't add two matrizes with diffrent size" << std::endl;
		return *this;
	}
	int size = size_n * size_m;
	double* temp_data = new double[size];
	for(int i = 0; i < size; i++){
		temp_data[i] = data[i] + m.data[i];
	}
	matrix out(size_n, size_m, temp_data);
	delete[] temp_data;
	return out;
}

matrix matrix::operator-(const matrix& m){
	if(m.size_n != size_n || m.size_m != size_m){
		std::cout << BOLDRED << "ERROR: " << RESET << "can't subtract two matrizes with diffrent size" << std::endl;
		return *this;
	}
	int size = size_n * size_m;
	double* temp_data = new double[size];
	for(int i = 0; i < size; i++){
		temp_data[i] = data[i] - m.data[i];
	}
	matrix out(size_n, size_m, temp_data);
	delete[] temp_data;
	return out;
}
