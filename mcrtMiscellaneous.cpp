#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <string>
#include <math.h>
#include "MCRT_library.h"
//#include "localmathlib.h"
//#include "materials.h"

string CommentOut(string &str) {
	int i_com;
	string comment;

	if (str.empty()) {
		return "";
	}

	i_com = (int)str.find("//");
	if (i_com > -1) {
		comment = str.substr(i_com + 2, str.size() - (i_com + 2));
		str = str.substr(0, i_com);
	}
	return comment;
}


void GetArguments(string strArg, objArgsearch *strSearch, const int &nSearch,
	string strWhere) {

	vector<string> args;
	//bool *state = new bool[nSearch];
	objArgsearch *objlocal;

	int *iArray;
	double *dArray;
	bool *boolArray;
	string *strArray;
	cdouble *cdArray, cdval;
	Point3D *p3DArray, p3Dval;
	doublestr *dstrArray, dstrVal;

	bool bool_str;

	vstring token = Tokenize(strArg);
	if (token.size() > nSearch) {
		ErrorMsg(strWhere + ":  Number of inputs exceeds " + to_string(nSearch));
		return;
	}

	for (int i = 0; i < token.size(); i++) {
		args = EvalFunction(token.at(i));

		for (objlocal = strSearch; !objlocal->Name.empty(); objlocal++) {

			if (!token.at(i).compare(objlocal->Name)) {
				if (!objlocal->IsFound) { // Element has not been found yet

					objlocal->IsFound = true;
					switch (objlocal->Type)
					{
					case PA_DOUBLE:
						dArray = (double *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs &&
							args.size() >= objlocal->nArgs_min) {
							for (int i = 0; i < args.size(); i++) {
								*(dArray++) = stod(args.at(i));
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs));
						}
						break;

					case PA_INT:
						iArray = (int *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs &&
							args.size() >= objlocal->nArgs_min) {
							for (int i = 0; i < args.size(); i++) {
								*(iArray++) = (int)stod(args.at(i));
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs));
						}
						break;

					case PA_STRING:
						strArray = (string *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs &&
							args.size() >= objlocal->nArgs_min) {
							for (int i = 0; i < args.size(); i++) {
								*(strArray++) = args.at(i);
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs));
						}
						break;

					case PA_CDOUBLE:
						cdArray = (cdouble *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs*2 &&
							args.size() >= objlocal->nArgs_min*2) {
							for (int i = 0; i < args.size()/2; i++) {
								cdval = real(stod(args.at(i))) + II * real(stod(args.at(i + 1)));
								*(cdArray++) = cdval;
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs*2));
						}
						break;

					case PA_POINT3D:
						p3DArray = (Point3D *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs*3 &&
							args.size() >= objlocal->nArgs_min * 3) {
							for (int i = 0; i < args.size()/3; i++) {
								p3Dval.SetPoint3D(stod(args.at(i)),
									stod(args.at(i + 1)),
									stod(args.at(i + 2)));
								*(p3DArray++) = p3Dval;
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs*3));
						}
						break;

					case PA_DBLSTR:
						dstrArray = (doublestr *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs * 2 &&
							args.size() >= objlocal->nArgs_min * 2) {
							for (int i = 0; i < args.size() / 2; i++) {
								dstrVal.Val = stod(args.at(i));
								dstrVal.strfunc = args.at(i + 1);
								*(dstrArray++) = dstrVal;
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs*2));
						}
						break;

					case PA_BOOL:
						boolArray = (bool *)objlocal->Storage;
						if (args.size() <= objlocal->nArgs &&
							args.size() >= objlocal->nArgs_min) {
							for (int i = 0; i < args.size(); i++) {
								if (!args.at(i).compare("true"))  bool_str = true;
								if (!args.at(i).compare("false")) bool_str = false;
								*(boolArray++) = bool_str;
							}
						}
						else {
							ErrorMsg(strWhere + ": Number of arguments in " +
								objlocal->Name + " is: " + to_string(objlocal->nArgs));
						}
						break;

					default:
						break;
					};
				}
				else {
					ErrorMsg(strWhere + ": Only one instance of " + objlocal->Name + " is allowed");
				}
			}
		}
	}
}


string check_command_line(int argc, char* argv[]) {
	if (argc < 2) {
		return "";
	}
	else if (argc == 2) {
		return argv[1];
	}
	else {
		ErrorMsg("Only one setup file is allowed");
	}

	return 0;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status) {
	static thread_local long i1, i2, ma[56];   /* ma[0] is not used. */
	long        mj, mk;
	short       i, ii;

	if (Type == 0) {              /* set seed. */
		mj = MSEED - (Seed < 0 ? -Seed : Seed);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for (ii = 1; ii <= 4; ii++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		i1 = 0;
		i2 = 31;
	}
	else if (Type == 1) {       /* get a number. */
		if (++i1 == 56) i1 = 1;
		if (++i2 == 56) i2 = 1;
		mj = ma[i1] - ma[i2];
		if (mj < MZ) mj += MBIG;
		ma[i1] = mj;
		return (mj * FAC);
	}
	else if (Type == 2) {       /* get status. */
		for (i = 0; i < 55; i++) Status[i] = ma[i + 1];
		Status[55] = i1;
		Status[56] = i2;
	}
	else if (Type == 3) {       /* restore status. */
		for (i = 0; i < 55; i++)
			ma[i + 1] = Status[i];
		i1 = Status[55];
		i2 = Status[56];
	}
	else
		puts("Wrong parameter to RandomGen().");
	return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
