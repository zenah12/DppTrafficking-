#ifndef GLOBAL
#define GLOBAL //Global.hpp


//Define global fixed parameters, these I can use to define matrices etc but cannot be modified in main
#define numx (62)             // discretization in space, xdivs in Daniel's code
#define numsave (100)              // every how many steps to save profiles
#define L_P (62)                // define tissue length
#define numt (10000000)           // maximum number of time steps
#define store_n (5)         // prediod of storing points
#define Pi (3.141592653589793238462643383279502884197169399375)

#define T_th (350)

#define datN (14) //18-22 for WT, 16 for dally, 16 for LOPpent, 11 for small, 15 small pent

//define global parameters that can be modified in main and used throughput
extern double b;        //bleaching depth
extern double h;        //length of ROI
extern double d;        //initial position for bleaching; typically set to 0


//unknown parameters in equations
extern double kon;    // mu in analysis
extern double d0;
extern double koff;   //this is alpha_off
extern double ki;
extern double k2;
extern double kr;
extern double ko;
extern double k1;
extern double k;
extern double kN;
extern double pkr;
extern double lambda;
extern double koffD;
extern double konD;


// upper limit for paramater boundaries
extern double kon_max;    // mu in analysis
extern double koff_max;
extern double ki_max;
extern double k2_max;
extern double ko_max;
extern double kN_max;
extern double k_max;
extern double b_max;
extern double lambda_max;

// upper limit for paramater boundaries
extern double kon_min;    // mu in analysis
extern double koff_min;
extern double ki_min;
extern double k2_min;
extern double ko_min;
extern double kN_min;
extern double k_min;
extern double b_min;
extern double lambda_min;
//parameters in equations that are fixed
extern double r;
extern double a;
extern double lambdaexp;
extern double rho;
extern double phii;
extern int n_source;
extern double nu;

extern double phi;      //flux from source (how do we measure this?!)

extern double dt;      // time discretization
extern double dx;      // space discretization

extern int FRAP;            // run FRAP if = 1, dont't otherwise
extern int nano;            // run nanobody if = 1, don't otherwise

extern int i_max;           // maximum index for bleaching in FRAP
extern int i_min;
extern double T;
extern double R2;
extern double chi2;
extern double pacc;
extern double probacc;

extern double epsilon;

#endif
