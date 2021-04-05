// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define _USE_MATH_DEFINES

#include "targetver.h"

#include <limits>

#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <time.h>  
#include <locale.h>
#include <Windows.h>
#include <direct.h>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

/* PARDISO prototype. */
extern "C" void pardisoinit(void   *, int    *, int *, int *, double *, int *);
extern "C" void pardiso(void   *, int    *, int *, int *, int *, int *,
	double *, int    *, int *, int *, int *, int *,
	int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec(int *, int *, double *, int *);
extern "C" void pardiso_printstats(int *, int *, double *, int *, int *, int *, double *, int *);


#include "Well.h"
#include "Props.h"
#include "Acoplamento_poco.h"
#include "props_res.h"
#include "reser_func_2.h"
#include "pressure_drop_wellbore.h"
#include "calc_wellbore.h"
#include "calc_wellbore_multi.h"
#include "unit_cv.h"
#include "coupling.h"


// TODO: reference additional headers your program requires here
