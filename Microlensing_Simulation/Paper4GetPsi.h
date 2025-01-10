#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <time.h>
#include <functional>
#include <algorithm>
#include <ctime>
#include <vector>
#include <thread>
#include <string>
#include <random>
#include "../spline.h"


using namespace std;


double* Preparation4CreatPhiKappaStar(string MacroType, double kappa, double gamma, double kappaStar, double MassPerSolarMass, double LensRedshift, int PrecisionFactor);
double* ReadMicroAndMinus(string MacroType, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Columns, long DisColumns, long DisRows, long EspBoost);
double* ReadCreatPsi(string MacroType, double kappa, double gamma, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Column, long DisColumns, long DisRows, double SourceY1, double SourceY2, long EspBoost);

double* SaddlePreparation4CreatPhiKappaStar(string MacroType, double kappa, double gamma, double kappaStar, double MassPerSolarMass, double LensRedshift, int PrecisionFactor);
double* SaddleReadMicroAndMinus(double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X2Set, long Columns, long DisColumns, long DisRows, long EspBoost);
double* SaddleReadCreatPsi(double kappa, double gamma, double kappa_star, long BoostLenRows, long MicroFileLenRows, long MicroFileLenColumns, long FileLenRows, double * FileX2Set, double * X1Set, double * X2Set, long Column, long DisColumns, long DisRows, double SourceY1, double SourceY2, long EspBoost);