#pragma once

//Libary: matrix
#include "../matrix/matrix.h"

//Libary: settings manager
#include "../settings_manager/settings_manager.h"

void optimize(int steps, settings_manager settings, bool save_result = false);

void test_potential(settings_manager settings);

void test_eigen_vector(settings_manager settings);

void from_image(settings_manager settings, const char* name);

void stability_test(settings_manager settings);

void error_test(settings_manager settings, int number_of_tests, const char* saving_name);

void find_stable(settings_manager settings, int steps);
