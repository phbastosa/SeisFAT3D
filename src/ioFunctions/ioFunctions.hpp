# ifndef IO_FUNCTIONS_CPP
# define IO_FUNCTIONS_CPP

# include <cmath>
# include <chrono>
# include <string>
# include <vector>
# include <sstream>
# include <fstream>
# include <iostream>
# include <algorithm>

# define PI 3.14159265359

bool str2bool(std::string s);

void import_binary_float(std::string path, float * array, int n);
void export_binary_float(std::string path, float * array, int n);

void import_text_file(std::string path, std::vector<std::string> &elements);

std::string catch_parameter(std::string target, std::string file);

std::vector<std::string> split(std::string s, char delimiter);

# endif