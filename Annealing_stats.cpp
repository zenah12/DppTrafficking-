//
//  Annealing_stats.cpp
//  1D_simple
//
//  Created by Zena H on 23.05.17.
//  Copyright © 2017 Zena Hadjibvasiliou. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "Dynamics.hpp"
#include "Annealing_stats.hpp"

// This function takes a vector with three entries holding the present value for a parameter, the max_i and min_i of its working intervals
// and update the working interval according to the present value of T
void Annealing_Stats::UpdateParInterval(double parMat[3], double T_now, double par_min, double par_max){
    double dist1=0.0;
    double dist2=0.0;
    
    dist1 = parMat[0] - par_min;
    dist2= par_max - parMat[0];
    
    // update interval
    parMat[1] = par_min + (1.0-T_now)*dist1;
    parMat[2] = par_max - (1.0-T_now)*dist2;
}

// give a triple pars holiding a parameter value and its lower and upper limit it samples a new values and updates the old value
double Annealing_Stats::sample_par(double pars[3]){
    double par = 0.0;


    par = pars[1] + (pars[2]-pars[1])*unifRand();
    pars[0] = par;
    //printf("T: %.4f par: %.5f min: %.25f   max:%.15f \n", T, pars[0],pars[1],pars[2]);


    return(par);
}

double Annealing_Stats::f_d0(){
    double value  = 0.0;
    
    value = (k*koff*kon*(4.0*k*ko*lambdaexp*lambdaexp + koff*(-1.0 + 4.0*ko*lambdaexp*lambdaexp)))/(4.0*(k + koff)*(k*ko + koff + ko*koff));
    return(value);
}

double Annealing_Stats::f_rho(){
    double value  = 0.0;
    
    value = (k2*(1.0 + ko)*koff*kon + exp(1.0/lambdaexp)*k2*(2.0*k*ko + (1.0 + ko)*koff*(2.0 + kon)))/((k*ki + k2*(1.0 + k + ko))*koff*kon + exp(1.0/lambdaexp)*(2.0*k*k2*ko + k*(k2 + ki)*koff*kon + k2*(1.0 + ko)*koff*(2.0 + kon)));
    return(value);
}

double Annealing_Stats::f_phi(){
    double value  = 0.0;
    
    value = ((1.0 + exp(1.0/lambdaexp))*k*ki*koff*kon)/((k*ki + k2*(1.0 + k + ko))*koff *kon + exp(1.0/lambdaexp)*(2.0*k*k2*ko + k*(k2 + ki)*koff*kon + k2*(1.0 + ko)*koff*(2.0 + kon)));
    return(value);
}

double Annealing_Stats::DataMean(double mat[datN]){
    double sum=0.0;
    int n=datN;
    int i=0;
    
    for(i=0; i<n; i++) sum = sum + mat[i];
    sum = sum/(n*1.0);
    return(sum);
}//

double Annealing_Stats::DataSS_Tot(double mat[datN]){
    double mean=0.0;
    double dif = 0.0;
    int n=datN;
    int i=0;
    
    mean = DataMean(mat);
    //printf("mean of data  = %.5f\n", mean);
    
    for(i=0; i<n; i++) dif = dif + (mat[i] - mean)*(mat[i]-mean);
    //printf("ssRTot  = %.5f\n", dif);

    return(dif);
}//


double Annealing_Stats::f_R2(double data[datN], double sims[datN],  double ssTot_dat){
    double ssRes = 0.0;
    double value=0.0;
    int n=datN;
    int i=0;
    
    for(i=0; i<n; i++) ssRes = ssRes + (data[i] - sims[i])*(data[i] - sims[i]);
    value = 1.0 - ssRes/ssTot_dat;
    //printf("R2: %.10f\n", value);
    return(value);
}//


double Annealing_Stats::f_chi2(double data[datN], double sims[datN]){
    double tot = 0.0;
    int n=datN;
    int i=0;
    
    for(i=0; i<n; i++) tot = tot + (data[i] - sims[i])*(data[i] - sims[i])/data[i];
    
    return(tot);
}//


void Annealing_Stats::UpdateParams(double pars[14][T_th], int index){
    pars[0][index] = d0;
    pars[1][index] = kon;
    pars[2][index] = koff;
    pars[3][index] = k;
    pars[4][index] = ki;
    pars[5][index] = k2;
    pars[6][index] = kN;
    pars[7][index] = b;
    pars[8][index] = phii;
    pars[9][index] = rho;
    pars[10][index] = R2;
    pars[11][index] = chi2;
    pars[12][index] = lambdaexp;
    pars[13][index] = ko;
}


double Annealing_Stats::Iterpolate(double FRAP_rec[2][numt/store_n+1], int i, double sims_rec[datN], double time){
    double value = 0.0;
    // iterpolate between values at i and i-1
    //printf("time: %.5f  i: %.d idt:%.5f y1: %.5f    y0: %.5f\n", time, i, i*1.0*dt, FRAP_rec[0][i], FRAP_rec[0][i-1]);
    
    value = FRAP_rec[0][i-1] + (time - FRAP_rec[1][i-1])*(FRAP_rec[0][i] - FRAP_rec[0][i-1])/(FRAP_rec[1][i] - FRAP_rec[1][i-1]);
    return(value);
}

void Annealing_Stats::get_sims_rec(double FRAP_rec[2][numt/store_n+1], double exp_tt[datN], double exp_r[datN], double sims_rec[datN]){
    int i = 0;
    int j = 0;
    
    for(i=0; i<datN; i++){
        // find index of simulated time that is >= exp_tt[i]
        while(exp_tt[i]>FRAP_rec[1][j]){
            j++;
        }// while
        if((FRAP_rec[1][j]-exp_tt[i])< 0.0000001){ // if same, then sims_rec[i] = FRAP_rec[1][index]
            sims_rec[i] = FRAP_rec[0][j];
        }// if
        else{ // if not the same, interpolate between index and index-1
            sims_rec[i] = Iterpolate(FRAP_rec, j, sims_rec, exp_tt[i]);
        }//else
     //   printf("time exp: %.4f  time sim: %.4f  rec sim:%.4f \n ", exp_tt[i], FRAP_rec[1][j],sims_rec[i]);
    }// for
}// get_sims_rec

// returns 0 if the working set is to be rejected and 1 if it is to be accepted
int Annealing_Stats::KeepSet(double Si, double S){
    double acc=0.0;
    
    acc = exp((Si-S)/T);

    // if acc is below 1 the new set can only be accepted with a probability
    if(acc<1.0 && (1.0-Si/S)<pacc*T && unifRand()<probacc){
        acc = 1.0;
    }
    //printf("T: %.5f acc: %.25f   S:%.15f  Si%.15f\n", T, acc, S, Si);
    if(acc>=1.0) return(1);
    else return(0);
}

void Annealing_Stats::RevertParams(double pars[14][T_th], int index, double konmat[3], double koffmat[3], double kimat[3], double k2mat[3], double kNmat[3], double komat[3], double bmat[3]){
    
    d0 = pars[0][index];
    kon = pars[1][index];
    konmat[0] = kon;
    koff =pars[2][index];
    koffmat[0] = koff;
    k = pars[3][index];
    ki = pars[4][index];
    kimat[0] = ki;
    k2 = pars[5][index];
    k2mat[0] = k2;
    kN = pars[6][index];
    kNmat[0]=kN;
    b = pars[7][index];
    bmat[0]=b;
    ko=pars[13][index];
}
























