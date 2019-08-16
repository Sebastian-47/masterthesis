#include <iostream>

#include <fstream>
#include <sstream>

//Libary: matrix
#include "settings_manager.h"

settings_manager::settings_manager(){
  number_of_settings = 0;
  saving_path = "";
  values_d = NULL;
  values_s = NULL;
  options = NULL;
}

settings_manager::settings_manager(const settings_manager& s){
  number_of_settings = s.number_of_settings;
  saving_path = s.saving_path;
  values_d = new double[number_of_settings];
  values_s = new std::string[number_of_settings];
  options = new std::string[number_of_settings];
  for(int i = 0; i < number_of_settings; i++){
    values_d[i] = s.values_d[i];
    values_s[i] = s.values_s[i];
    options[i] = s.options[i];
  }
}

settings_manager::settings_manager(const char* saving_path_){
  saving_path = std::string(saving_path_);

  std::string line;
  std::ifstream settings_file(saving_path.c_str());
  if(settings_file.is_open()){
    //count the number of lines;
    int number_of_lines = 0;
    while (getline(settings_file, line)){
      number_of_lines++;
    }
    settings_file.close();
    settings_file = std::ifstream(saving_path.c_str());
    number_of_settings = number_of_lines;
    
    values_d = new double[number_of_lines];
    values_s = new std::string[number_of_lines];
    options = new std::string[number_of_lines];
    int line_counter = 0;
    while (getline(settings_file, line)){
      std::string buffer;                 // Have a buffer string
      std::stringstream ss(line);         // Insert the string into a stream
      int string_pos = 0;
      
      char delim = ' ';
      while(getline(ss, buffer, delim)){
        if(string_pos == 0){
          options[line_counter] = buffer; 
        }
        if(string_pos == 1){
          try{
            values_d[line_counter] = std::stod(buffer);
            values_s[line_counter] = "double";
          } catch(const std::invalid_argument&){
            values_d[line_counter] = 0;
            values_s[line_counter] = buffer;
          }
        }
        if(string_pos == 2){
          std::cout << "Warning: There are to many arguments in Line " << line_counter + 1 << " from file " << saving_path << std::endl;
        }
        string_pos++;
      }
      line_counter++;
      }
    settings_file.close();
  } else {
    std::cout << "Error: Unable to open file" << std::endl; 
  }
}

//destructor:
settings_manager::~settings_manager(){
  delete[] values_d;
  delete[] values_s;
  delete[] options;
}

void settings_manager::load(){
  *this = settings_manager(saving_path.c_str());
}

void settings_manager::load(const char* saving_path_){
  *this = settings_manager(saving_path_);
}

void settings_manager::save() const{
  std::ofstream settings_file(saving_path.c_str());
  std::string temp = "double";
  std::stringstream ss;
  if(settings_file.is_open()){
    for(int i = 0; i < number_of_settings; i++){
      settings_file << options[i];
      settings_file << " ";
      if(values_s[i] == temp){
        ss.str("");
        ss << values_d[i]; 
        settings_file << ss.str();
      } else {
        settings_file << values_s[i];
      }
      settings_file << "\n";
    }
    settings_file.close();
  } else {
    std::cout << "Error: Unable to open file " << saving_path << std::endl;
  } 
}

void settings_manager::save(const char* saving_path_){
  std::string str = saving_path_;
  saving_path = str;
  this->save();
}

void settings_manager::print() const{
  std::cout << "+---------------------------------------------------------------+" << std::endl;
  std::cout << "|The current settings:                                          |" << std::endl;
  std::cout << "+-------+-------------------------------+-----------------------+" << std::endl;
  std::cout << "|Nr     | Option                        | Value                 |" << std::endl;
  std::cout << "+-------+-------------------------------+-----------------------+" << std::endl;
  for(int i = 0; i < number_of_settings; i++){
    std::cout << "|" << i + 1 << "\t| " << options[i];
    if(options[i].length() < 13){
      std::cout << "\t\t\t| ";
    } else if(options[i].length() < 18){
      std::cout << "\t\t| ";
    } else {
      std::cout << "\t| "; 
    }
    std::string temp = "double";
    if(values_s[i] == temp){
      std::cout << values_d[i] << "\t\t\t|";
    } else {
      std::cout << values_s[i];
      if(values_s[i].length() < 5){
        std::cout << "\t\t\t|"; 
      } else if(values_s[i].length() < 14){
        std::cout << "\t\t|"; 
      } else {
        std::cout << "\t|"; 
      }
    }
    std::cout << std::endl;
  }
  std::cout << "+-------+-------------------------------+-----------------------+" << std::endl;
}

void settings_manager::set(const char* option, double value){
  std::string str(option);
  bool found = false;
  for(int i = 0; i < number_of_settings; i++){
    if(options[i] == str){
      found = true;
      values_d[i] = value;
      break;
    }
  }
  if(!found){
    std::cout << "Error: There is no option " << str << std::endl;
  }
}

void settings_manager::set(const char* option, const char* value){
  std::string str(option);
  bool found = false;
  for(int i = 0; i < number_of_settings; i++){
    if(options[i] == str){
      found = true;
      values_s[i] = std::string(value);
      break;
    }
  }
  if(!found){
    std::cout << "Error: There is no option " << str << std::endl;
  }
}

double settings_manager::get_double(const char* option) const{
  std::string str(option);
  for(int i = 0; i < number_of_settings; i++){
    if(options[i] == str){
      return values_d[i];
    }
  }
  std::cout << "Error: There is no option " << str << std::endl;
}

std::string settings_manager::get_string(const char* option) const{
  std::string str(option);
  for(int i = 0; i < number_of_settings; i++){
    if(options[i] == str){
      return values_s[i];
    }
  }
  std::cout << "Error: There is no option " << str << std::endl;
}

void settings_manager::user_input(){
  this->print();
  std::cout << "Enter number of the option you want to change: (enter '0' to return)" << std::endl;
  int input;
  bool check = true;
  while(check){
    std::cout << "Choose the option you want to change: ";
    std::string user_input;
    std::cin >> user_input;

    check = false;
    try{
      input = stoi(user_input);
    } catch(const std::invalid_argument&){
      std::cout << "Please enter only a number" << std::endl;
      check = true;
      input = 0;
    }
    if(input > number_of_settings || input < 0){
      std::cout << "Please enter a number between 0 and" << number_of_settings << std::endl;
      check = true;
      input = 0;
    }
  }
  if(input != 0){
    std::string user_option;
    std::cout << "Enter the new value: ";
    std::cin >> user_option; 
    std::string str = "double";
    check = true;
    if(values_s[input - 1] == str){
      double user_double = values_d[input - 1];
      try{
        user_double = stod(user_option);
      } catch(const std::invalid_argument&){
        std::cout << "The option " << options[input - 1] << " only takes numbers as argument!" << std::endl;
        check = false;
      }
      values_d[input - 1] = user_double;
    } else {
      values_s[input - 1] = user_option;
    }
  }
  if(check){
    std::cout << "Changed succsessfully the option " << options[input - 1] << "!" << std::endl;
  } else {
    std::cout << "No option was changed!" << std::endl;
  }
}       

settings_manager& settings_manager::operator=(const settings_manager& s){
  if(this != &s){
    number_of_settings = s.number_of_settings;
    saving_path = s.saving_path;
    
    delete[] values_d;
    delete[] values_s;
    delete[] options;

    values_d = new double[number_of_settings];
    values_s = new std::string[number_of_settings];
    options = new std::string[number_of_settings];

    for(int i = 0; i < number_of_settings; i++){
      values_d[i] = s.values_d[i];
      values_s[i] = s.values_s[i];
      options[i] = s.options[i];
    }
  }
  return *this;
}
