//--------- DECLARATION OF ALL GLOBAL VARIABLES, FUNCTIONS AND CLASSES ---------
// from 09/01/2012 you have to add to your project constants.cpp
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
#ifndef _COULOMB_FUNCTIONS
#define _COULOMB_FUNCTIONS
extern double *DFC,*GC,*DGC,*FC,*SIGMA;
extern int *IEXP;
extern double FACL[51],DFACL[51];
void firstCoulombStart();
void Wittaker(double _eta,double mu,double X, double *W, double *WP, double EPS);
void GAMMA  (int Nn,double X,double Y,double *GR,double *GI);
void COULOM (int LL,double ETA1,double RO1,
   double Ff[],double FP[],double G[],double GP[],double SIGMA[],int IEXP[]);
#endif
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
#ifndef _ALL_V_P_H
#define _ALL_V_P_H
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
#define NDEFAULT 1001
#define SMD_MAX_LS_N 200
#define SMD_MAX_LS_P 200
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
/*#include "potential.h"
#include "bound.h"
#include "bounds.h"
#include "resonance.h"
#include "scatt.h"
#include "crsect.h"
#include "coul_semi.h" */
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
extern double
   hbar, c, mp, mn, mN, amu, h_mn, hc, large_val,
   h2_mn, e2, degTorad, radTodeg, small_add, small_e, zerro;
extern bool CondonShortley;
extern int NGauss;
extern double *xGauss,*wGauss;
extern char buf[];
extern char NucName[121][5];
extern bool doFirstCoulombJob;
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
bool printToFile(double *x,double *y,int n,char *fname);
double erfc_ap(double x);
double erf_ap(double x);
double erf(double x);
double getR(double r,int m1,int m2);
double getr(double R,int m1,int m2);
bool setGaussQuadrature(int N);
bool doGaussQuadrature(int N);
bool doLaguerreQuadrature(int N);
double Plm_cos_th(int L,int M,double x);
double Pl0_cos_th(int L,double x,double *p);
bool LegendrePL0ToArray(int L,double x,double *pl);
double d2LegandreP(int L,double x,double *pl,double *dpl);
complex Ylm(double cos_th, double fi, int l, int m, double *re, double *im);
double Ql0_cos_th(int L,double x);
bool LegendreQL0ToArray(int L,double x,double *ql);
void YlmToArray(double x, double phi, int l, complex *Y);
double Laguerre(double x,int n);
double D1Laguerre(double x,int n,double *L);
double erf(double x);
double fsign(double x,double y);
int isign(int x,int y);
bool chet(int n);
void imag_i_powL(int L,double *re,double *im);
void i_powL(int L,double *re,double *im);
complex I_powL(int L);
double Npow(double x,int N);
int _1pow(int n);
void calculateFACL();
void complexMultiplication(double re1,double im1,double re2,double im2,double *re,double *im);
double interpolate(double p,double f0,double f1);
double interpolate3Point(double p,double f_1,double f0,double f1);
double Interpolation(double x,double *X,double *F,int N);
double InterpolationCentrifugal(double x,double *X,double *F,int N,bool DoCheck,double Limit = 0);
double Interpolation2D(double x,double y,double *X,double *Y,double **F,int Nx,int Ny);
complex Interpolation2D_complex(double x,double y,double *X,double *Y,complex **F,int Nx,int Ny);
double Interpolation2D_Triangle(double x,double y,double *X,double *Y,double **F,int Nx,int Ny);
void SplineCoefficients(double *x, double *y, int n, double yp1, double ypn, double *y2);
double SplineInterpolation(double x, double *ya, double *xa, int n, double *y2a);
double NumericalDerivation
	(double x,double *X,double *F,int N,bool DoCheck = false,double InterpolLimit = 1.e150);
double _1stNumericalDerivation(int ix,int N,double *F,double h);
double _2ndNumericalDerivation(int ix,int N,double *F,double h);   
int getIndex(double x,double array[],int N);
double getPhiAngle(double x,double y);
char* getElementName(int z);
double getAngleBetweenTwoVectors(double r1,double x1,double y1,double z1,
	double r2,double x2,double y2,double z2,bool getCosine);
double getAngleBetweenTwoSphericalVectors(double r1,double th1,double phi1,
	double r2,double th2,double phi2,bool getCosine);

double j0(double x); 
double dj0(double x);
double j1(double x); 
double dj1(double x);
double j2(double x);
double n0(double x);
double dn0(double x);
double n1(double x); 
double dn1(double x);
double n2(double x);
double jl(int l,double x,double *jl_1);
double nl(int l,double x,double *nl_1);
double djl(int l,double x);
double dnl(int l,double x);
void jBessel(double *Re,double *Im,int l);
void nBessel(double *Re,double *Im,int l);
void c_jBessel(complex j[], int l);
void c_nBessel(complex n[], int l);
void c_hHankel(complex h[],int l,int plus_minus);
void GreenL(double pre, double pim, double r1, double r2, int l,double *Re, double *Im, bool DoRicatti);
void _GreenL(double pre,double pim,double r1, double r2,int l, double Re[],double Im[], bool DoRicatti);
double Clebsh_Gordon(int JX1,int JX2,int JX3,int MX1,int MX2);
double dWigner(int J,int M1,int M2,double angle);

char* getElementName(int z);
char* getNucleusName(int a,int z,bool formated);
char* getJp(int L,int s);
char* getJ(int _2J);
void AiBi(double z,double *ai,double *bi);
bool JacobyToJacoby(int mi,int mj,int mk,double Ecm,
                    double q0,double thq0,double fiq0,
                              double thp0,double fip0,
                    double *q1,double *thq1,double *fiq1,
                    double *p1,double *thp1,double *fip1);

void tridag(complex a[], complex b[], complex c[], complex r[], complex u[],long n);
void bandec(complex **a, long n, int m1, int m2, complex **al, long *indx, long *d);
void banbks(complex **a, long n, int m1, int m2, complex **al, long *indx, complex *b);
//------------------------------------------------------------------------------------------------
extern char *FILE_NAME_MASS;
bool readMasses();
double getMass(int Z,int A);
//------------------------------------------------------------------------------------------------
extern double FACTORIAL[];
double factorial(int n);
long double long_factorial(int n);
double dFact(int n);
double CG0(int II, int JJ, int KK);
double PHASEF(int N);
double YXFCT(int M,int N);
double DELR(int J1,int J2,int J3);
double RACAH(int J1,int J2,int J3,int J4,int J5,int J6);
double RACAH_2(double J1,double J2,double J3,double J4,double J5,double J6);
double S6J(double J1,double J2,double J3,double J4,double J5,double J6);
double WINEJ(int J1,int J2,int J3,int J4,int J5,int J6,int J7,int J8,int J9);
int e_ijk(int i,int j,int k);
void ReynalRevai(int n, int l1, int l2, int _n, int _l1, int _l2, int L,
                 double cos_phi, double sin_phi,
                 double *re, double *im);
double JacobiPolinomial(double z,int n,double l1,double l2);
double NK_lxly(int K,int n,int lx,int ly);
bool HyperHarmonics_LM(int K,int lx,int ly,int L,int M,
                       double cos_thx,double phix,
                       double cos_thy,double phiy,
                       double cos_alpha,double sin_alpha,
                       double *re,double *im);
bool HyperHarmonics_MxMy(int K,int lx,int ly,int mx,int my,
                         double cos_thx,double phix,
                         double cos_thy,double phiy,
                         double cos_alpha,double sin_alpha,
                         double *re,double *im);
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
double ClusterFormFactor(double q,int A,int Z,bool shell_model);
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
bool JacobyToJacoby_spherical
   (double sin_Phi, double cos_Phi,
    double  x0,double  thx0,double  fix0,double  y0,double  thy0,double  fiy0,
    double *x1,double *thx1,double *fix1,double *y1,double *thy1,double *fiy1);
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
extern __int64 IY64;
double Urand();
double Urand64();
double NormalDistrib();
void setRandomInitialValue();
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
void progonka_Collins(complex *A, complex *B,complex *X,int N);
void progonka(complex *C,complex *A,complex *D, complex *B,complex *X,int N);
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
#endif//  _ALL_V_P_H
