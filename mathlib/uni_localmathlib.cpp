#include "uni_localmathlib.h"

std::string get_application_path(void) {
	std::string s_path;

	char path[PATH_MAX];
	int len = PATH_MAX;
	char szTmp[32];

	sprintf(szTmp, "/proc/%d/exe", getpid()); // "/proc/self/exe"
	int bytes = std::min((int)readlink(szTmp, path, len), len - 1);
	if (bytes >= 0) path[bytes] = '\0';

	s_path = dirname(path); // convert char to string (remove application name)

	return s_path;
}

std::string get_working_path(void){
	char cwd[1024];
	getcwd(cwd, sizeof(cwd));

	return std::string(cwd) + "/";
}
