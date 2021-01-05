#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>

#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"
#include "Dynamics.hpp"
#include "Annealing_stats.hpp"


using namespace std;

//declare parameters
//FRAP parameters
double h=11.0/3.0;        //length of ROI
double d=0.0;        //initial position for bleaching; typically set to 0
double pkr=0.0022; //WT: 0.00137   Dally: 0.0011  LOP: 0.000754924    Dally small: 0.00112

////////////////////////////////////
// parameters in equations
// varried in annealing algorithm
double b=0.195;        //bleaching depth
double kon=0.93;    //this is mu in analysis
double koff=0.13;   //this is alpha_off
double konD = 0.01;
double koffD = 0.01;
double ki=0.0;
double k2=0.0;
double ko= 0.00022/pkr;//WT: 0.0002113   Dally: 0.00019 LOP: 0.0000869939   Dally small: 0.000155
double kN = 0.0;//0.15/pkr;  //0.184877/pkr;//63.2;
// errors from nanobody and FRAP experiments
double eko=0.000014/pkr;
double ekN= kN*0.2;//0.0067/pkr;
double eb=0.1;//0.01;
// max values
//  free parameters
double kon_max = 1.0;//1.5;
double koff_max = 1.0;
double ki_max = ko;
double k2_max = 0.0001/pkr;
double k_max = 10.0;
// parameteres contrained by nanoboy and FRAP
double ko_max = 0.00018/pkr;//+eko;
double kN_max = 0.0042/pkr;//kN;
double b_max = 0.3;//WT:0.27   Dally:0.22  LOP:0.12 Dally small: 0.24  LOPpent: 0.18
double lambda_max = 22.68;
// min values
//  free parameters
double kon_min = -1.0;//0.5;
double koff_min = -1.0;
double ki_min = 0.00001/pkr;//0.001;
double k2_min = 0.00001/pkr;//0.01;
double k_min = 10.0;
// parameteres contrained by nanoboy and FRAP
double ko_min = 0.0001/pkr;// - eko;
double kN_min = 0.0026/pkr;//kN - ekN;
double b_min = 0.2;//b - eb;
double lambda_min = 19.08;
// computed or fixed
double r=1.0;
double a=1.0;
double kr=1.0;
double lambdaexp = 29.0/2.6;//WT: 24.7   Dally: 14.7  LOP: 29.0
double d0=5.36;
double k1=ko-ki;
double k= 100.0;//kN*(1.0+1.0/kon);
////////////////////////////////////
double phi = 1.0;      //flux from source
int n_source = 10.0-2.0;
double nu = 0.1;
double dx = 0.0;
double dt = 0.0;
double epsilon = pow(10.0, -15.0);     // stability threshold

// other parameters
double phii = 0.0;
double rho = 0.0;

// program params
int i_min = numx/2 + n_source/2 + 1; //FRAP next to the source on the RHS
int i_max=i_min+4;
double meanDat = 0.0;
double SStot = 0.0;

int FRAP = 1;
double FRAP_t = pkr*4000.0;
// simulated annealing parameters
double T = 1.0;
double alpha = 0.9;
int niter = 10;//10;
double R2_acc = 0.0;//the R2 value of the presently accepted set
double R2 = 0.0;//the R2 value of the working set
double chi2_acc = 10.0;//the chi2 value of the presently accepted set
double chi2 = 0.0;//the chi2 value of the working set
double acc = 0.0; // auxiliary acceptance parameter
double pacc = 0.2;
double probacc = 0.1;//0.1;
int keep = 0;

double timeCheck = pkr*36.0;
double timePassed = 0.0;

int indexPrint = 0;
double par1 = 0.0;
double par2 = 0.0;

//Add objects
Dynamics myDynamics;
Annealing_Stats myAnnealing;
//////////////////////////////////////////////////////////////////////////////
////////////////////////// M A I N   P R O G R A M //////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {
    
	

    indexPrint = atoi(argv[1]);
    double randomize = 1.0*indexPrint;	
    srand((unsigned)time(NULL)+indexPrint);
    //int repeat = 10;
    par1 = atof(argv[2]);
    par2 = atof(argv[3]);
    //// koff FIXED ///
    koff_min = par1;
    koff_max = par1+1.0;
    //// kon FIXED ///
    kon_min = par2;
    kon_max = par2+1.0;
    
    // for repeat
    int count=0;
    int stabCheck = 0;
    int i=0;
    int saving_count=0;
    int ii = 0;
    int j=0;
    int all_store_i =0;
    int good_store_i=0;
    
    //for(repeat=1; repeat<10; repeat++){
    count=0;
    stabCheck = 0;
    i=0;
    saving_count=0;
    ii = 0;
    j=0;
    all_store_i =0;
    good_store_i=0;
    chi2_acc = 10000.0;
    R2_acc = 0.0;
    
    // initiate koff and kon values
    kon = kon_max + (kon_max-kon_min)*unifRand();
    koff = koff_max + (koff_max-koff_min)*unifRand();
    // initiate pool profiles
    Init_N(myDynamics.uL);
    Init_N(myDynamics.uSR);
    Init_N(myDynamics.uSL);
    Init_N(myDynamics.uSE);
    Init_N(myDynamics.uSi);
    // initiate pool recovery integral profile logs
    Init_FRAP_logs(myDynamics.iuL_F_log);
    Init_FRAP_logs(myDynamics.iuSR_F_log);
    Init_FRAP_logs(myDynamics.iuSL_F_log);
    Init_FRAP_logs(myDynamics.iuSE_F_log);
    Init_FRAP_logs(myDynamics.iuSi_F_log);
    Init_FRAP_logs(myDynamics.iuT_F_log);
    // Initiate par triplets and set initial values and intervals
    // randomly choose koff and kon within interval
    koffD=koff_min+(koff_max-koff_min)*unifRand(); // sample dimensional parameter
    koff = koffD / pkr; //nondimensionalise
    konD=kon_min+(kon_max-kon_min)*unifRand(); // sample dimensional parameter
    kon = konD / koffD;  // nondimensionalise
    // set up triplets for varying pars
//    Init_3(myAnnealing.kon_mat, kon, kon_min, kon_max);
//    Init_3(myAnnealing.koff_mat, koff, koff_min, koff_max);
    Init_3(myAnnealing.kon_mat, konD, kon_min, kon_max);
    Init_3(myAnnealing.koff_mat, koffD, koff_min, koff_max);
    Init_3(myAnnealing.ki_mat,ki, ki_min, ki_max);
    Init_3(myAnnealing.k2_mat,k2, k2_min, k2_max);
    Init_3(myAnnealing.b_mat,b, b_min, b_max);
    Init_3(myAnnealing.kN_mat,kN, kN_min, kN_max);
    Init_3(myAnnealing.k_mat,k, k_min, k_max);
    Init_3(myAnnealing.ko_mat,ko, ko_min, ko_max);
    Init_3(myAnnealing.lambda_mat, lambdaexp, lambda_min, lambda_max);

    // initiate paramater storing matrices
    Init_Params(myAnnealing.params);
    Init_Params(myAnnealing.params_accepted);
    // compute dx
    dx = L_P/(double numx);
    /////////////////////////////////////////////////////////////////////////
    ////////////////////     R E A D  D A T A        ////////////////////////
    /////////////////////////////////////////////////////////////////////////
    // WT data
//    double exp_t[] = {30.0*pkr, 210.0*pkr, 390.0*pkr, 570.0*pkr, 750.0*pkr, 930.0*pkr, 1170.0*pkr, 1410.0*pkr, 1650.0*pkr, 1890.0*pkr, 2130.0*pkr, 2370.0*pkr, 2610.0*pkr, 2850.0*pkr, 3090.0*pkr, 3330.0*pkr, 3570.0*pkr, 3810.0*pkr};// , 4050.0*pkr, 4290.0*pkr, 4530.0*pkr, 4770.0*pkr};
//    double exp_rec[] = {0.273142951, 0.382482307, 0.433335448, 0.460559701, 0.488227873, 0.499212588, 0.541135657, 0.549482215, 0.574959688, 0.585431538, 0.618031903, 0.62023424, 0.621802463, 0.639226199, 0.645531011, 0.650644305, 0.65666496, 0.648831298};//, 0.677532663, 0.697748071, 0.702280769, 0.712563424};
    // Dally data
//    double exp_rec[] = {0.229243397, 0.323138439, 0.382371345, 0.438307066, 0.466498852, 0.469757185, 0.504373765, 0.526524911, 0.532156465, 0.544214179, 0.555829365, 0.565758352, 0.571331295, 0.579842649, 0.574511395, 0.598389798};
  //  double exp_t[] = {30*pkr, 210*pkr, 390*pkr, 570*pkr, 750*pkr, 930*pkr, 1170*pkr, 1410*pkr, 1650*pkr, 1890*pkr, 2130*pkr, 2370*pkr, 2610*pkr, 2850*pkr, 3090*pkr, 3330*pkr};
    // choose dally
//    double exp_rec[] = {0.23618772, 0.333342435, 0.39347313, 0.436071374, 0.46215654, 0.468546102, 0.502312427, 0.516603037, 0.525190785, 0.538257273, 0.542418833, 0.556788516, 0.562778861, 0.561447447, 0.551781925, 0.573265464};
//    double exp_t[] = {30.0*pkr, 210.0*pkr, 390.0*pkr, 570.0*pkr, 750.0*pkr, 930.0*pkr, 1170.0*pkr, 1410.0*pkr, 1650.0*pkr, 1890.0*pkr, 2130.0*pkr, 2370.0*pkr, 2610.0*pkr, 2850.0*pkr, 3090.0*pkr, 3330.0*pkr};
    // LOP data
     //   double exp_rec[] = {0.124691347, 0.282336249, 0.323322006, 0.399188561, 0.445510692, 0.474124201, 0.49300545, 0.506695052, 0.51012764, 0.516407094, 0.540654177, 0.543153296, 0.562490053, 0.553854972, 0.59618856, 0.598377974};
     //   double exp_t[] = {30.0*pkr, 170.0*pkr, 310.0*pkr, 590.0*pkr, 870.0*pkr, 1150.0*pkr, 1430.0*pkr, 1710.0*pkr, 1990.0*pkr, 2270.0*pkr, 2550.0*pkr, 2830.0*pkr, 3110.0*pkr, 3390.0*pkr, 3670.0*pkr, 3950.0*pkr};
    // LOP data 2
   // double exp_rec[] = {0.16041108, 0.28196564, 0.34075043, 0.415454953, 0.452690447, 0.482196669, 0.501588328, 0.510187201, 0.517222695, 0.519805904, 0.541247679, 0.545851551, 0.556132941 ,0.55237859};
   // double exp_t[] = {30.0*pkr, 170.0*pkr, 310.0*pkr, 590.0*pkr, 870.0*pkr, 1150.0*pkr, 1430.0*pkr, 1710.0*pkr, 1990.0*pkr, 2270.0*pkr, 2550.0*pkr, 2830.0*pkr, 3110.0*pkr, 3390.0*pkr};
    // LOPpent data
   //double exp_rec[] = {0.176241239, 0.291968732, 0.329751083, 0.4006611, 0.44003465, 0.45946071, 0.452067893, 0.459441576, 0.473554682, 0.495809391, 0.505443284, 0.509782911, 0.505673187, 0.500783022, 0.509057943, 0.508724442};
    //double exp_t[] = {30*pkr, 170.0*pkr, 310.0*pkr, 590.0*pkr, 870.0*pkr, 1150.0*pkr, 1430.0*pkr, 1710.0*pkr, 1990.0*pkr, 2270.0*pkr, 2550.0*pkr, 2830.0*pkr, 3110.0*pkr, 3390.0*pkr, 3670.0*pkr, 3950.0*pkr};
    // LOP small
    //double exp_rec[] = {0.184990551, 0.321013138, 0.376725909, 0.429430162, 0.481954049, 0.475701949, 0.490627007, 0.529691258, 0.508415005, 0.511564161, 0.565276728};
    //double exp_t[] = {30.0*pkr, 170.0*pkr, 310.0*pkr, 610.0*pkr, 910.0*pkr, 1210.0*pkr, 1510.0*pkr, 2110.0*pkr, 2410.0*pkr, 2710.0*pkr, 3010.0*pkr};
    // LOP small Pent 85
//    double exp_t[] = {30*pkr, 170*pkr, 310*pkr, 590*pkr, 870*pkr, 1150*pkr, 1430*pkr, 1710*pkr, 1990*pkr, 2270*pkr, 2550*pkr, 2830*pkr, 3110*pkr, 3390*pkr, 3670*pkr};
//    double exp_rec[] = {0.08341685, 0.154930745, 0.176015133, 0.219514349, 0.247903028, 0.263381857, 0.276526047, 0.298455833, 0.312962969, 0.312853854, 0.326178169, 0.349010409, 0.354510522, 0.354130775, 0.379009065};
    // LOP small Pent 100
//    double exp_t[] = {30*pkr, 170*pkr, 310*pkr, 590*pkr, 870*pkr, 1150*pkr, 1430*pkr, 1710*pkr, 1990*pkr, 2270*pkr, 2550*pkr, 2830*pkr, 3110*pkr, 3390*pkr};
//    double exp_rec[] = {0.092084639, 0.168187861, 0.179855207, 0.241154013, 0.251013542, 0.25776181, 0.276770007, 0.298105232, 0.303927342, 0.306528559, 0.314648319, 0.335404981, 0.339572582, 0.338663899};
    // LOP small Pent OE
    double exp_t[] = {30*pkr, 170*pkr, 310*pkr, 590*pkr, 870*pkr, 1150*pkr, 1430*pkr, 1710*pkr, 1990*pkr, 2270*pkr, 2550*pkr, 2830*pkr, 3110*pkr, 3390*pkr};
    
    double exp_rec[] = {0.305001869, 0.370121651, 0.417454733, 0.445716055, 0.48489662, 0.497089568, 0.534394893, 0.541708851, 0.556015229, 0.555615714, 0.549831597, 0.543891193, 0.547563848, 0.564777858, 0.574995677, 0.576721155};
    // find data mean and SS_tot to be used later
    meanDat = myAnnealing.DataMean(exp_rec);
    SStot = myAnnealing.DataSS_Tot(exp_rec);
    /////////////////////////////////////////////////////////////////////////
    ////////////////////////     ABC ROUTINE        /////////////////////////
    /////////////////////////////////////////////////////////////////////////
    T = 1.0;
   // printf("numt/store_n+1: %.d\n",numt/store_n+1);
    while(T<T_th){
        // reset things
        R2 = 0.0;
        chi2 = 0.0;
        Init_FRAPrec(myAnnealing.FRAP_rec);
        // i_max = 0;
        h=11.0/3.0;
        // initiate pool profiles
        Init_N(myDynamics.uL);
        Init_N(myDynamics.uSR);
        Init_N(myDynamics.uSL);
        Init_N(myDynamics.uSE);
        Init_N(myDynamics.uSi);
        // initiate pool recovery integral profile logs
        Init_FRAP_logs(myDynamics.iuL_F_log);
        Init_FRAP_logs(myDynamics.iuSR_F_log);
        Init_FRAP_logs(myDynamics.iuSL_F_log);
        Init_FRAP_logs(myDynamics.iuSE_F_log);
        Init_FRAP_logs(myDynamics.iuSi_F_log);
        Init_FRAP_logs(myDynamics.iuT_F_log);
        // sample parameters (update parameter values to be used globaly)
//        kon = pow(10.0, myAnnealing.sample_par(myAnnealing.kon_mat));
//        koff = pow(10.0, myAnnealing.sample_par(myAnnealing.koff_mat));
        // comment below for dimensionaless sampling
        konD = pow(10.0, myAnnealing.sample_par(myAnnealing.kon_mat));
        koffD = pow(10.0, myAnnealing.sample_par(myAnnealing.koff_mat));
        koff = koffD/pkr;
        kon = konD/koffD;
        // end comment for dimensionless sampling
        k2 = myAnnealing.sample_par(myAnnealing.k2_mat);
//        ki_max = k2;
//        Init_3(myAnnealing.ki_mat,ki, ki_min, ki_max);
        ko = myAnnealing.sample_par(myAnnealing.ko_mat);
        ki_max = ko;
        ki = myAnnealing.sample_par(myAnnealing.ki_mat);
        kN = myAnnealing.sample_par(myAnnealing.kN_mat);
//        k = myAnnealing.sample_par(myAnnealing.k_mat);
        b = myAnnealing.sample_par(myAnnealing.b_mat);
        k= kN*(1.0+1.0/kon);
        // ki=(0.5 + 0.1*unifRand() - 0.1*unifRand())*k2;
        //ki = rand_Normal(0.5, 0.2)*k2;
        //while(ki>ko || ki<0) ki = rand_Normal(0.5, 0.2)*k2;
        k1 = ko - ki;
        //ki=2.3*k2+(2.5*k2-2.3*k2)*unifRand();
        lambdaexp = myAnnealing.sample_par(myAnnealing.lambda_mat)/2.6;
	// compute value for k,d0, rho and phi
        d0 = myAnnealing.f_d0();
        rho = myAnnealing.f_rho();
        phii = myAnnealing.f_phi();
        // get dt for current parameters
        dt = myDynamics.Get_dt();
        // Initiate pool profiles near SS with current paremeters
        Initiate_Profile_SS(myDynamics.uL, 0);
        Initiate_Profile_SS(myDynamics.uSR, 1);
        Initiate_Profile_SS(myDynamics.uSL, 2);
        Initiate_Profile_SS(myDynamics.uSE, 3);
        Initiate_Profile_SS(myDynamics.uSi, 4);
        // Initiate saved profiles
        Init_N(myDynamics.uL_save);
        Init_N(myDynamics.uSR_save);
        Init_N( myDynamics.uSL_save);
        Init_N(myDynamics.uSE_save);
        Init_N(myDynamics.uSi_save);
        /////////////////////////////////////////////////////////////////////////
        ////////////////////     SS  +   F R A P         ////////////////////////
        /////////////////////////////////////////////////////////////////////////
        // steady state simulation //
        count = 0;
        stabCheck = 0;
        timePassed =0.0;
        Initiate_source(myDynamics.source);
        while(stabCheck!=1){
            // update profiles
            myDynamics.UpdateProfiles(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.source);
            timePassed = timePassed + dt; //count*dt;
            //check stability if timePassed passed
            if(timeCheck<timePassed && timePassed>0.0){
                stabCheck = myDynamics.StabilityCheck(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.uL_save, myDynamics.uSR_save, myDynamics.uSL_save, myDynamics.uSE_save, myDynamics.uSi_save);
                // stability check and log update
                copyNtoN(myDynamics.uL, myDynamics.uL_save);
                copyNtoN(myDynamics.uSR, myDynamics.uSR_save);
                copyNtoN(myDynamics.uSL, myDynamics.uSL_save);
                copyNtoN(myDynamics.uSE, myDynamics.uSE_save);
                copyNtoN(myDynamics.uSi, myDynamics.uSi_save);
                // reset stability check counter to 0
                timePassed = 0.0;
            }
        }// while loop
        // FRAP simulation //
        // bleach and follow recovery
        for(i=0; i<5; i++) myDynamics.prebleach[i]=0.0;// initiate prebleach matrix
        // h = dx*(1.0*(i_max));   // update length of bleaching region
        myDynamics.PreBleach_Integrate(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.prebleach);//get prebleach value for each pool
        myDynamics.Bleach(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi);// bleach //
        // recovery loop //
        count=0;
        saving_count = 0;
        // save first integral before dynamics resume
        myDynamics.FRAP_Integrate(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.iuL_F_log, myDynamics.iuSR_F_log, myDynamics.iuSL_F_log, myDynamics.iuSE_F_log, myDynamics.iuSi_F_log,myDynamics.iuT_F_log, saving_count, count, myDynamics.prebleach);
        saving_count++;
        //printf("sim rec full first point:  %.5f\n", myDynamics.iuT_F_log[0][0]);
        while(count*dt<FRAP_t){
            myDynamics.UpdateProfiles(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.source);// update profiles
            count++;
          // store current profile
            if(count % store_n == 0){
                myDynamics.FRAP_Integrate(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.iuL_F_log, myDynamics.iuSR_F_log, myDynamics.iuSL_F_log, myDynamics.iuSE_F_log, myDynamics.iuSi_F_log, myDynamics.iuT_F_log, saving_count, count, myDynamics.prebleach);
                saving_count++;
            }
            //printf("count %.d\n",count);
            //printf("saving count %.d\n",saving_count);
        }// while
        //////////////////////////////////////////////////////////////////
        ///////////// STEADY STATE AND FRAP FINISHED ////////////////////
        /////////////////////////////////////////////////////////////////
        // print total recovery
       // for(j=0; j<saving_count; j++) printf("sim rec full :  %.5f\n", myDynamics.iuT_F_log[0][j]);
       // printf("\n");
        // use FRAP recovery data to construct a matrix comparable to experimental data
        myAnnealing.get_sims_rec(myDynamics.iuT_F_log, exp_t, exp_rec, myAnnealing.FRAP_rec);
        // for(j=0; j<datN; j++) printf("sim rec:  %.5f\n", myAnnealing.FRAP_rec[j]);
        // GET R^2 and chi^2 values between data and simulation
        R2 = myAnnealing.f_R2(exp_rec, myAnnealing.FRAP_rec, SStot);
        chi2 = myAnnealing.f_chi2(exp_rec, myAnnealing.FRAP_rec);
        // store present run
        myAnnealing.UpdateParams(myAnnealing.params, all_store_i);
        all_store_i++;
        T++;
        //printf("here: d0: %.10f kon: %.10f koff: %.10f k2%.10f  ki%.5f  b%.5f  k:%.5f  T: %.5f chi^2: %.5f    R^2: %.5f\n", d0, kon, koff, k2, ki, b, k, T, chi2, R2);
    }// while T
    Write_Params(myAnnealing.params, all_store_i, indexPrint, 0);
    //}//for repeat
	return 0;
}



