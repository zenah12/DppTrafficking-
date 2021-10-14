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
#include "Sampling_stats.hpp"


using namespace std;

//declare parameters (definitions in Global.hpp)
//FRAP parameters
double h=11.0/3.0;          //length of ROI normalized by cell diameter
double d=0.0;               //initial position for bleaching; typically set to 0
double pkr=0.00172;         //recycling rate measured by nanobody experiment
double kl = 0.0001/pkr;     //dimensioless leakage rate for simulations that incorporate loss of Dpp in extracellualr space

// parameters in equations
// all parameters are nondimensionalized by dividing by the recycling rate pkr
// varried in sampling algorithm
double b=0.0;                   //bleaching depth
double d0=0.0;                  // dimensionless diffusion coefficient
double kon=0.0;                 // fit parameter, sampled
double koff=0.0;                // fit parameter, sampled
double konD = 0.0;              // dimensional kon; to be nondimensionalized following sampling
double koffD = 0.0;             // dimensional koff; to be nondimensionalized following sampling
double ki=0.0;
double k2=0.0;
double k1=ko-ki;
double k= 0.0;
double ko= 0.00022/pkr;         // from nanobody experiment
double kN = 0.0;                // effective internalization; CI for sampling set by nanobody experiment
double lambdaexp = 0.0;          // will be adjusted following sampling

// ranges for sampling algorithm
// max values
double kon_max = 1.0;//1.5;
double koff_max = 1.0;
double ki_max = ko;
double k2_max = 0.0001/pkr;
double k_max = 10.0;
// parameteres constrained by nanobody, FRAP, steady state decay length
double ko_max = 0.00024/pkr;    //CI for sampling set by nanobody experiment
double kN_max = 0.014/pkr;      //CI for sampling set by nanobody experiment
double b_max = 0.16;            //CI for sampling set by FRAP
double lambda_max = 29.08;      //CI for sampling set by steady state decay length measurement
// min values
double kon_min = -1.0;          // lower limit for sampling interval set by input parameter; line 127
double koff_min = -1.0;         // lower limit for sampling interval set by input parameter; line 128
double ki_min = 0.00001/pkr;    // lower limit set by long FRAP experiment
double k2_min = 0.00001/pkr;    // lower limit set by long FRAP experiment
double k_min = 0.0;            // k is not sampled but computed from sampled kN, kon and koff (Section 2.2.1 in Supplementary)
// parameteres constrained by nanobody, FRAP, steady state decay length
double ko_min = 0.00013/pkr;
double kN_min = 0.012/pkr;
double b_min = 0.14;
double lambda_min = 27.16;

// computed or fixed
//double r=1.0;
double a=1.0;                           //cell dimater
int n_source = 10.0-2.0;
double nu = 0.1;                        // production rate at source cells, FRAP curves are normalized and so nu independent
double dx = 0.0;
double dt = 0.0;
double epsilon = pow(10.0, -15.0);     // stability threshold for steady state

// initiate other parameters
double phii = 0.0;                      // immobile fraction
double rho = 0.0;                       // extracellular fraction (analytical)
double rhoComp = 0.0;                   // extracellular fraction (numerical)

// program params
int i_min = numx/2 + n_source/2 + 1;    // start of FRAP ROI; next to the source on the RHS
int i_max=i_min+4;                      // end of FRAP ROI
double T = 1.0;

// statistics
double R2 = 0.0;                        //the R2 value of the working set
double chi2 = 0.0;                      //the chi2 value of the working set
double meanDat = 0.0;                   //mean of experimental data
double SStot = 0.0;                     //sum square of experimental data

double FRAP_t = pkr*4000.0;             //duration of FRAP experiment in seconds
double timeCheck = pkr*36.0;            //sampling interval for DFRAP and stability check
double timePassed = 0.0;

int indexPrint = 0;                     // inout paramater to randomize sampling seed
double par1 = 0.0;                      // used to speficy range of sampling for kon through HPC jobs
double par2 = 0.0;                      // used to speficy range of sampling for koff through HPC jobs

//Add objects
Dynamics myDynamics;
Sampling_stats mySampling;

// ******************************
// M A I N   P R O G R A M
// ******************************

int main(int argc, char ** argv) {
    
    // ******************************
    // R E A D  E X P E R I M E N T A L  D A T A  I N
    // ******************************
    
    // LOP data L = 144um
    double exp_rec[] = {0.16041108, 0.28196564, 0.34075043, 0.415454953, 0.452690447, 0.482196669, 0.501588328, 0.510187201, 0.517222695, 0.519805904, 0.541247679, 0.545851551, 0.556132941 ,0.55237859};
    double exp_t[] = {30.0*pkr, 170.0*pkr, 310.0*pkr, 590.0*pkr, 870.0*pkr, 1150.0*pkr, 1430.0*pkr, 1710.0*pkr, 1990.0*pkr, 2270.0*pkr, 2550.0*pkr, 2830.0*pkr, 3110.0*pkr, 3390.0*pkr};
    // data statistics to be used for R^2 computation
    meanDat = mySampling.DataMean(exp_rec);
    SStot = mySampling.DataSS_Tot(exp_rec);

    // ******************************
    // READ INPUT VARIABLES AND INITIATE ARRAYS
    // ******************************
    
    // randomize intial conditions
    indexPrint = atoi(argv[1]);
    double randomize = 1.0*indexPrint;	
    srand((unsigned)time(NULL)+indexPrint);

    par1 = atof(argv[2]);   // read in paramet that determines lower limi of koff sampling interval (in powers of 10)
    par2 = atof(argv[3]);   // read in paramet that determines lower limi of kon sampling interval (in powers of 10)
    
    // specify koff sampling interval for current run (in powers of 10)
    koff_min = par1;
    koff_max = par1+1.0;
    // specify kon sampling interval for current run (in powers of 10)
    kon_min = par2;
    kon_max = par2+1.0;
    
    // indices for run
    int count=0;
    int stabCheck = 0;
    int i=0;
    int saving_count=0;
    int ii = 0;
    int j=0;
    int all_store_i =0;
    int good_store_i=0;
    
    count=0;
    stabCheck = 0;
    i=0;
    saving_count=0;
    ii = 0;
    j=0;
    all_store_i =0;
    good_store_i=0;

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
    // Initiate par triplets (max, min, current) values; and set initial values and intervals
    // randomly choose koff and kon within interval
    koffD=koff_min+(koff_max-koff_min)*unifRand(); // sample dimensional parameter
    koff = koffD / pkr; //nondimensionalise
    konD=kon_min+(kon_max-kon_min)*unifRand(); // sample dimensional parameter
    kon = konD / koffD;  // nondimensionalise
    // set up triplets for varying pars
    Init_3(mySampling.kon_mat, konD, kon_min, kon_max);
    Init_3(mySampling.koff_mat, koffD, koff_min, koff_max);
    Init_3(mySampling.ki_mat,ki, ki_min, ki_max);
    Init_3(mySampling.k2_mat,k2, k2_min, k2_max);
    Init_3(mySampling.b_mat,b, b_min, b_max);
    Init_3(mySampling.kN_mat,kN, kN_min, kN_max);
    Init_3(mySampling.k_mat,k, k_min, k_max);
    Init_3(mySampling.ko_mat,ko, ko_min, ko_max);
    Init_3(mySampling.lambda_mat, lambdaexp, lambda_min, lambda_max);
    // initiate paramater storing matrices
    Init_Params(mySampling.params);
    Init_Params(mySampling.params_accepted);
    // compute dx
    dx = L_P/(double numx);
    
    // ******************************
    // ******************************
    // ******************************
    // ABC ROUTINE
    // ******************************
    // ******************************
    // ******************************
    
    T = 1.0;
    while(T<T_th){
        // reset things
        R2 = 0.0;
        chi2 = 0.0;
        Init_FRAPrec(mySampling.FRAP_rec);
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
        // comment below for dimensionaless sampling
        konD = pow(10.0, mySampling.sample_par(mySampling.kon_mat));
        koffD = pow(10.0, mySampling.sample_par(mySampling.koff_mat));
        koff = koffD/pkr;
        kon = konD/koffD;
        // end comment for dimensionless sampling
        k2 = mySampling.sample_par(mySampling.k2_mat);
        ko = mySampling.sample_par(mySampling.ko_mat);
        ki_max = ko;
        ki = mySampling.sample_par(mySampling.ki_mat);
        kN = mySampling.sample_par(mySampling.kN_mat);
        b = mySampling.sample_par(mySampling.b_mat);
        k= kN*(1.0+1.0/kon);
        k1 = ko - ki;
        lambdaexp = mySampling.sample_par(mySampling.lambda_mat)/2.6;
        // compute value for k,d0, rho and phi
        d0 = mySampling.f_d0();
        rho = mySampling.f_rho();
        phii = mySampling.f_phi();
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
        
        // ******************************
        // ******************************
        // STEADY STATE AND FRAP NUMERICS
        // ******************************
        // ******************************
        
        // ******************************
        // steady state simulation
        // ******************************
        
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
        // compute rho (extracellular fraction) for current parameter set
        rhoComp = myDynamics.getRhoNum(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi);
        
        // ******************************
        // FRAP simulation
        // ******************************
        
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
        while(count*dt<FRAP_t){
            myDynamics.UpdateProfiles(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.source);// update profiles
            count++;
          // store current profile
            if(count % store_n == 0){
                myDynamics.FRAP_Integrate(myDynamics.uL, myDynamics.uSR, myDynamics.uSL, myDynamics.uSE, myDynamics.uSi, myDynamics.iuL_F_log, myDynamics.iuSR_F_log, myDynamics.iuSL_F_log, myDynamics.iuSE_F_log, myDynamics.iuSi_F_log, myDynamics.iuT_F_log, saving_count, count, myDynamics.prebleach);
                saving_count++;
            }

        }// while
        // STEADY STATE AND FRAP FINISHED

        // use FRAP recovery data to construct a matrix comparable to experimental data
        mySampling.get_sims_rec(myDynamics.iuT_F_log, exp_t, exp_rec, mySampling.FRAP_rec);
        // Compute R^2 and chi^2 values comparing data and simulation
        R2 = mySampling.f_R2(exp_rec, mySampling.FRAP_rec, SStot);
        chi2 = mySampling.f_chi2(exp_rec, mySampling.FRAP_rec);
        // store all parameters and quality of fit for present run
        mySampling.UpdateParams(mySampling.params, all_store_i);
        all_store_i++;
        T++;
    }// while T
    
    // store all sampled parameters with quality of fit
    Write_Params(mySampling.params, all_store_i, indexPrint, 0);
	return 0;
}



