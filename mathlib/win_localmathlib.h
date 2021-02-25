#define NOMINMAX
#include <iostream>
#include <algorithm>
#include <Windows.h>
#include <shlwapi.h>
#pragma comment(lib,"shlwapi.lib")

std::string get_application_path(void);
std::string get_working_path(void);
