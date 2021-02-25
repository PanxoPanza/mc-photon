#include "win_localmathlib.h"

std::string get_application_path(void) {
	std::string s_path;

	char path[MAX_PATH];
	GetModuleFileNameA(NULL, path, MAX_PATH);				// Extract .exe path
	PathRemoveFileSpecA(path);											// remove .exe file
	s_path = path;																	// convert char to string
	std::replace(s_path.begin(), s_path.end(), '\\', '/');	// replace '\' by '/'

	return s_path;
}

std::string get_working_path(void){
	char cwd[1024];
	getcwd(cwd, sizeof(cwd));

	return std::string(cwd) + "/";
}
