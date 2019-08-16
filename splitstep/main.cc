#include <iostream>
#include <chrono>

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
#include "arithmetics.h"

void welcome(){
  std::cout << "   _____         __ _  __         _____  __                       " << std::endl; 
  std::cout << "  / ___/ ____   / /(_)/ /_       / ___/ / /_ ___   ____           " << std::endl;  
  std::cout << "  \\__ \\ / __ \\ / // // __/______ \\__ \\ / __// _ \\ / __ \\          " << std::endl; 
  std::cout << " ___/ // /_/ // // // /_ /_____/___/ // /_ /  __// /_/ /          " << std::endl; 
  std::cout << "/____// .___//_//_/ \\__/       /____/ \\__/ \\___// .___/           " << std::endl; 
  std::cout << "    _/_/        __                             /_/  _             " << std::endl; 
  std::cout << "   /  _/____   / /_ ___   ____ _ _____ ____ _ / /_ (_)____   ____ " << std::endl; 
  std::cout << "   / / / __ \\ / __// _ \\ / __ `// ___// __ `// __// // __ \\ / __ \\" << std::endl; 
  std::cout << " _/ / / / / // /_ /  __// /_/ // /   / /_/ // /_ / // /_/ // / / /" << std::endl; 
  std::cout << "/___//_/ /_/ \\__/ \\___/ \\__, //_/    \\__,_/ \\__//_/ \\____//_/ /_/ " << std::endl; 
  std::cout << "                       /____/                                     " << std::endl; 
}

void print_options(){
  std::cout << "Options:" << std::endl;
  std::cout << "0  >> exit" << std::endl;
  std::cout << "1  >> change options" << std::endl;
  std::cout << "2  >> save options" << std::endl;
  std::cout << "3  >> load options" << std::endl;
  std::cout << "4  >> show options" << std::endl;
  std::cout << "5  >> run optimization" << std::endl;
  std::cout << "6  >> run potential test (only one iteration of 5)" << std::endl;
  std::cout << "7  >> run eigen vector test" << std::endl;
  std::cout << "8  >> run potential test with an image as potential" << std::endl;
  std::cout << "9  >> run the stability test" << std::endl;
  std::cout << "10 >> run error test" << std::endl;
  std::cout << "11 >> find minimal stable step width" << std::endl;
}

void user_interface(){
  welcome();
  settings_manager settings("options.dat");
  settings.print();
  
  print_options();  

  int number_of_settings = 11;
  int input = 1;
  bool check = true;
  
  std::string exit = "exit";
  std::string help = "help";

  while(input != 0){
  while(check){
    std::cout << "Please enter the number of what you want to do: << ";
    std::string user_input;
    std::cin >> user_input;

    check = false;
    if(user_input == exit){
      input = 0;
    } else if(user_input == help){
      print_options();
      check = true;
    } else {
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
  }
  check = true;

  switch(input){
  case 1:{ //change some options
    settings.user_input();
    }
    break;
  case 2:{ //save the current options
    std::string saving_path;
    std::cin >> saving_path;
    settings.save(saving_path.c_str());
    }
    break;
  case 3:{ //load new options
    std::string loading_path;
    std::cin >> loading_path;
    settings.load(loading_path.c_str());
    }
    break;
  case 4:{ //print the current options
    settings.print();
    }
    break;
  case 5:{ //run optimization
    std::cout << "Run an optimization:" << std::endl;
    int steps;
    std::string steps_input;
    std::cout << "Enter the Number of optimation steps, that will be performed: << ";
    std::cin >> steps_input;
    try{
      steps = stoi(steps_input);
    } catch(const std::invalid_argument&){
      std::cout << "Please enter only a number" << std::endl;
      check = true;
      break;
    }
    std::cout << "This may take some time..." << std::endl << std::endl;
    optimize(steps, settings);
    }
    break;
  case 6:{ //run potential test
    std::cout << "Run an simple test, with an irrational cut potential:" << std::endl << std::endl;
    test_potential(settings);
  }
  break;
  case 7:{ //run eigenvector test
    std::cout << "Run an eigen vector test:" << std::endl << std::endl;
    test_eigen_vector(settings);
  }
  break;
  case 8:{ //run eigenvector test
    std::cout << "Run an potential test with an image as potential:" << std::endl << std::endl;
    std::string path;
    std::cout << "please enter the path of the imput image << ";
    std::cin >> path;
    from_image(settings, path.c_str());
  }
  break;
  case 9:{ //run stability test
    std::cout << "Run the stability test:" << std::endl;
    stability_test(settings);
  }
  break;
  case 10:{ //run error test
    int steps;
    std::string steps_input;
    std::cout << "Enter the Number tests, that will be performed: << ";
    std::cin >> steps_input;
    try{
      steps = stoi(steps_input);
    } catch(const std::invalid_argument&){
      std::cout << "Please enter only a number" << std::endl;
      check = true;
      break;
    }
    std::string saving_path;
    std::cout << "Enter the saving path for the result: << ";
    std::cin >> saving_path;
    error_test(settings, steps, saving_path.c_str());
  }
  break;
  case 11:{ //run stability test
    int steps;
    std::string steps_input;
    std::cout << "Enter the Number tests, that will be performed: << ";
    std::cin >> steps_input;
    try{
      steps = stoi(steps_input);
    } catch(const std::invalid_argument&){
      std::cout << "Please enter only a number" << std::endl;
      check = true;
      break;
    }

    find_stable(settings, steps);
  }
  break;
  default:
    std::cout << "Good bye" << std::endl;
  }
  std::cout << std::endl << "----------------------------------------------------" << std::endl;
  if(input != 0){
    print_options();
  }
  }
}

void with_arguments(int argc, const char *argv[]){
  welcome();
  std::cout << "Calles Splitstep with arguments!, further userinput is disabled" << std::endl;
  std::cout << "The settings used are:" << std::endl;
  
  settings_manager settings("options.dat");
  settings.print();
  
  int converted = 0;
  bool check = true;
  std::string input = argv[1];
  try{
    converted = stoi(input);
  } catch(const std::invalid_argument&){
    std::cout << "The first argument has to be a number, the available options are:" << std::endl;
  }

  if(converted < 5 || converted > 11){
    std::cout << "The inpute has to be a number between 5 and 10, the available options are:" << std::endl;
    std::cout << "5  >> run optimization" << std::endl;
    std::cout << "6  >> run potential test (only one iteration of 5)" << std::endl;
    std::cout << "7  >> run eigen vector test" << std::endl;
    std::cout << "8  >> run potential test with an image as potential" << std::endl;
    std::cout << "9  >> run the stability test" << std::endl;
    std::cout << "10 >> run error test" << std::endl;
    std::cout << "11 >> find minimal stable step width" << std::endl;
    check = false;
  }
  
  if(check){
    switch(converted){
      case 5:{ //run optimization
        std::cout << "Run an optimization:" << std::endl;
        int steps;
        if(argc != 3){
          std::cout << "The optimzation only takes 2 arguments!";
          break;
        }
        std::string steps_input = argv[2];
        try{
          steps = stoi(steps_input);
        } catch(const std::invalid_argument&){
          std::cout << "Please enter only a number" << std::endl;
          break;
        }
        optimize(steps, settings, true);
        }
        break;
      case 6:{ //run potential test
        if(argc != 2){
          std::cout << "The optimzation only takes 1 argument!";
          break;
        }
        std::cout << "Run an simple test, with an irrational cut potential:" << std::endl << std::endl;
        test_potential(settings);
      }
      break;
      case 7:{ //run eigenvector test
        if(argc != 2){
          std::cout << "The optimzation only takes 1 argument!";
          break;
        }
        std::cout << "Run an eigen vector test:" << std::endl << std::endl;
        test_eigen_vector(settings);
      }
      break;
      case 8:{ //Compute the time evolution from an image
        if(argc != 3){
          std::cout << "The optimzation only takes 2 argument!";
          break;
        }
        std::cout << "Run an potential test with an image as potential:" << std::endl << std::endl;
        from_image(settings, argv[2]);
      }
      break;
      case 9:{ //run stability test
        if(argc != 2){
          std::cout << "The optimzation only takes 1 argument!";
          break;
        }
        std::cout << "Run the stability test:" << std::endl;
        stability_test(settings);
      }
      break;
      case 10:{ //run error test
        int steps;
        if(argc != 3){
          std::cout << "The optimzation only takes 2 arguments!";
          break;
        }
        std::string steps_input = argv[2];
        try{
          steps = stoi(steps_input);
        } catch(const std::invalid_argument&){
          std::cout << "Please enter only a number" << std::endl;
          check = true;
          break;
        }
        error_test(settings, steps, steps_input.c_str());
      }
      case 11:{ //run stability test
        int steps;
        if(argc != 3){
          std::cout << "The optimzation only takes 2 arguments!";
          break;
        }
        std::string steps_input = argv[2];
        try{
          steps = stoi(steps_input);
        } catch(const std::invalid_argument&){
          std::cout << "Please enter only a number" << std::endl;
          check = true;
          break;
        }
        find_stable(settings, steps);
      }
    }
  }
}

int main(int argc, const char *argv[]){
  if(argc == 1){
    user_interface();
  } else {
    with_arguments(argc, argv);
  }
  return 0;
}
