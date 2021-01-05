//
//  Annealing_stats.hpp
//  1D_simple
//
//  Created by Zena H on 23.05.17.
//  Copyright © 2017 Zena Hadjibvasiliou. All rights reserved.
//

//#ifndef Annealing_stats_hpp
//#define Annealing_stats_hpp

//#include <stdio.h>

//#endif /* Annealing_stats_hpp */

class Annealing_Stats
{
private:
    int x;
public:
    
    // matrix that holds all parameter values tested with R2, phi and rho value
    double params[14][T_th];
    // matrix that holds  parameter values accepted with R2, phi and rho value
    double params_accepted[14][T_th];
    // the 3 dimensional vectors below hold the value of the relevatn parameter and the upper and lower limits of its working interval.
    double kon_mat[3];
    double koff_mat[3];
    double ki_mat[3];
    double k2_mat[3];
    double b_mat[3];
    double kN_mat[3];
    double ko_mat[3];
    double k_mat[3];
    double lambda_mat[3];
    double FRAP_rec[datN];


    void UpdateParInterval(double parMat[3], double T_now, double par_min, double par_max);
    double Rsquared(double data[10][2], double sims[10][2]);
    double f_rho();
    double f_d0();
    double f_phi();
    double sample_par(double par[3]);
    void UpdateParams(double pars[14][T_th], int index);
    double DataMean(double mat[datN]);
    double DataSS_Tot(double mat[datN]);
    double f_R2(double data[datN], double sims[datN],  double ssTot_dat);
    double f_chi2(double data[datN], double sims[datN]);
    double Iterpolate(double FRAP_rec[2][numt/store_n+1], int i, double sims_rec[datN], double time);
    void get_sims_rec(double FRAP_rec[2][numt/store_n+1], double exp_tt[datN], double exp_r[datN], double sims_rec[datN]);
    int KeepSet(double S, double Si);
    
    void RevertParams(double pars[14][T_th], int index, double konmat[3], double koffmat[3], double kimat[3], double k2mat[3], double kN_mat[3], double komat[3], double bmat[3]);

};
