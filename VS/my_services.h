//----------------------------------------------------------------------
// file containce declaration of different service functions
//----------------------------------------------------------------------
// 1. global constants
extern const char NucName[][5];

// 2. methods and functions
void firstCoulombStart();

void doTrim(char *source, char * dest, char c);
const char* getElementName(int z);
double Npow(double x, int N);
double getPhiAngle(double x, double y);
double fsign(double x, double y);
int isign(int x, int y);
bool chet(int n);
double PHASEF(int N);
double factorial(int n);
long double long_factorial(int n);
int getIndex(double x, double array[], int N);
double Interpolation(double x, double *X, double *F, int N);
bool SLAU_GaussJordan(double **a, double *b, int N, bool DoInverse);