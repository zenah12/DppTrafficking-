#ifndef GLOBAL
#define GLOBAL //Global.hpp

#define Pi (3.141592653589793238462643383279502884197169399375)

//Define global fixed parameters, these I can use to define matrices etc but cannot be modified in main
#define numx (62)                   // discretization in space, xdivs in Daniel's code
#define numsave (100)               // every how many steps to save profiles
#define L_P (62)                    // define tissue length
#define numt (10000000)             // maximum number of time steps
#define store_n (5)                 // prediod of storing points
#define T_th (2)                    // number of repeats in current job
#define datN (14)                   //length of experimental FRAP

//define global parameters that can be modified in main and used throughput
extern double b;                    //bleaching depth
extern double h;                    //length of bleaching ROI
extern double d;                    //initial position for bleaching; typically set to 0, i.e. mext to source

//unknown parameters in equations
extern double kon;      // dimensionless binding rate = kon/koff in model
extern double d0;       // dimensionless diffusion constant = D/(kr*a)
extern double koff;     // dimensionless unbinding rate = koff/kr in model
extern double ki;       // dimensionless immobilization rate = ki/kr in model
extern double k2;       // dimensionless degradation in immobile pool = k2/kr in model
extern double kr;       // recycling rate (kr in model)
extern double ko;       // dimensionless output rate (=ko/kr in model)
extern double k1;       // dimensionless degradation in mobile pool (=k1/kr in model)
extern double k;        // dimensionless rate of endocytosis (=k/kr in model)
extern double kN;       // dimensionless rate of effective endocytosis (=k/kr in model)
extern double pkr;      // recycling rate
extern double lambda;   // decay length
extern double koffD;    // dimenstional koff rate for sampling in dimensional space
extern double konD;     // dimenstional kon rate for sampling in dimensional space
extern double kl;       // dimensionless leakage rate  (=kl/kr in model)
extern double klD;      // dimensional leakage rate

// upper limit for paramater sampling intervals
extern double kon_max;
extern double koff_max;
extern double ki_max;
extern double k2_max;
extern double ko_max;
extern double kN_max;
extern double k_max;
extern double b_max;
extern double lambda_max;

// upper limit for paramater sampling intervals
extern double kon_min;
extern double koff_min;
extern double ki_min;
extern double k2_min;
extern double ko_min;
extern double kN_min;
extern double k_min;
extern double b_min;
extern double lambda_min;
//parameters in equations that are fixed
extern double a;                // cell diamter
extern double lambdaexp;        // decay length
extern double rho;              // extracellular fraction analytical
extern double rhoComp;          // extracellular fraction numerical
extern double phii;             // immobile fraction
extern int n_source;            // source size
extern double nu;               // production rate at source

extern double dt;           // time discretization
extern double dx;           // space discretization

extern int i_max;           // maximum index for bleaching in FRAP
extern int i_min;           // minimum index for bleaching in FRAP
extern double T;
extern double R2;
extern double chi2;

extern double epsilon;      // steady state check

#endif
