#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "Dynamics.hpp"


///////////////////////////////////////////////////////////////////////
// Here all the functions used in the object Kinetics_SS are defined //
//////////////////////////////////////////////////////////////////////

// This object contains functions that solve for the steady state profile of the different pools.

//get input Delta for each cell from apical interactions
void Dynamics::UpdateProfiles(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double source[numx]){
    double uL_1[numx];
    double uSR_1[numx];
    double uSL_1[numx];
    double uSE_1[numx];
    double uSi_1[numx];
    
    int i=0;
    
    //initiate matrices
    Init_N(uL_1);
    Init_N(uSR_1);
    Init_N(uSL_1);
    Init_N(uSE_1);
    Init_N(uSi_1);
    
    //copy matrices
    copyNtoN(u_L, uL_1);
    copyNtoN(u_SR, uSR_1);
    copyNtoN(u_SL, uSL_1);
    copyNtoN(u_SE, uSE_1);
    copyNtoN(u_Si, uSi_1);
    
    // update profiles in middle
    for(i=1; i<(numx-1); i++){
        // upadate uL pool
        u_L[i] = uL_1[i] + dt*(d0*(uL_1[i-1]-2.0*uL_1[i]+uL_1[i+1])/(dx*dx) +koff*(uSR_1[i] + uSL_1[i+1]) - kon*koff*uL_1[i] + nu*source[i]);
        // update bound right and left pools
        u_SR[i] = uSR_1[i] +dt*(-(koff+k)*uSR_1[i]+0.5*kon*koff*uL_1[i]+0.5*uSE_1[i]);
        u_SL[i] = uSL_1[i] +dt*(-(koff+k)*uSL_1[i]+0.5*kon*koff*uL_1[i-1]+0.5*uSE_1[i]);
        // update endosomal pool
        u_SE[i] = uSE_1[i] + dt*(k*(uSL_1[i]+uSR_1[i]) - (1+ko)*uSE_1[i]);
        // update immobile pool
        u_Si[i] = uSi_1[i] + dt*(ki*uSE_1[i] - k2*uSi_1[i]);
    }
    
    // update edges; i = 0 and i = N-1 for u_L
    //reflective boundary conditions for u_L and u_SL
    u_L[0] = u_L[1];
    u_L[numx-1] = u_L[numx-2];
    u_SL[0] = u_SL[1];
    u_SL[numx-1] = u_SL[numx-2];
    //update u_SR, u_SE, s_Si
    u_SR[0] = uSR_1[0] +dt*(-(koff+k)*uSR_1[0]+0.5*kon*koff*uL_1[0]+0.5*uSE_1[0]);
    u_SR[numx-1] = uSR_1[numx-1] +dt*(-(koff+k)*uSR_1[numx-1]+0.5*kon*koff*uL_1[numx-1]+0.5*uSE_1[numx-1]);
    u_SE[0] = uSE_1[0] + dt*(k*(uSL_1[0]+uSR_1[0]) - (1+ko)*uSE_1[0]);
    u_SE[numx-1] = uSE_1[numx-1] + dt*(k*(uSL_1[numx-1]+uSR_1[numx-1]) - (1+ko)*uSE_1[numx-1]);
    u_Si[0] = uSi_1[0] + dt*(ki*uSE_1[0] - k2*uSi_1[0]);
    u_Si[numx-1] = uSi_1[numx-1] + dt*(ki*uSE_1[numx-1] - k2*uSi_1[numx-1]);
    
}//UpdateProfiles

double Dynamics::Get_dt(){
    double dtl=0.0;
    double dtsplus = 0.0;
    double dtse = 0.0;
    double min=10.0;
    
    // Courant-Levy-Friedrich criterion
    dtl = 2.0*pow(dx,2.0)/(4.0*d0+pow(dx,2.0)*(kon*koff+koff))*0.8;
    dtsplus = 2.0*dx*dx/(4.0/8.0+dx*dx*(kon*koff+1.0+koff+k))*0.8;
    dtse= 2.0*dx*dx/(4.0*k/8.0+dx*dx*(k+1.0+ko))*0.8;
    printf("dtl: %10f\n", dtl);
    
    // find min
    if(dtl<min) min = dtl;
    if(dtsplus<min) min = dtsplus;
    if(dtse<min) min = dtse;
    if(min>pkr*3600.0/100.0) min = pkr*3600.0/100.0; //makde sure we have enough resolution
  
    return(min);
}

// finds the first index where the position is > bleaching depth b
int Dynamics::get_i_max(){
    int value = 0;
    double b_now =0.0;
    
    while(b_now<h){
        value++;
        b_now = b_now+dx;
    }
    
    return(value+1);
}

void Dynamics::Bleach(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx]){
    int i=0;
    //printf("imax in bleach: %.d\n\n", i_max);
    // bleach all pools from i_min to i_max
    for(i=i_min; i<i_max; i++){
        u_L[i] = u_L[i]*b;
        u_SR[i] = u_SR[i]*b;
        u_SL[i] = u_SL[i]*b;
        u_SE[i] = u_SE[i]*b;
        u_Si[i] = u_Si[i]*b;
    }// for i
}

// this function integrates from position 0 to i_max in each pool
// and stores the integral and time for each pool.
void Dynamics::PreBleach_Integrate(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double prebl[5]){
    int i = 0;
    double prefactor=0.0;
    prefactor = 0.5*dx;
    
    for(i=i_min; i<i_max; i++){
        prebl[0] = prebl[0] + 2.0*u_L[i];
        prebl[1] = prebl[1] + 2.0*u_SR[i];
        prebl[2] = prebl[2] + 2.0*u_SL[i];
        prebl[3] = prebl[3] + 2.0*u_SE[i];
        prebl[4] = prebl[4] + 2.0*u_Si[i];
    }
    prebl[0] = prefactor*(prebl[0] - u_L[i_min] - u_L[i_max-1]);
    prebl[1] = prefactor*(prebl[1] - u_SR[i_min] - u_SR[i_max-1]);
    prebl[2] = prefactor*(prebl[2] - u_SL[i_min] - u_SL[i_max-1]);
    prebl[3] = prefactor*(prebl[3] - u_SE[i_min] - u_SE[i_max-1]);
    prebl[4] = prefactor*(prebl[4] - u_Si[i_min] - u_Si[i_max-1]);
    
}

// this function integrates from position 1 to i_max-1 in each pool and stores the integral and time for each pool and for the total.
void Dynamics::FRAP_Integrate(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx],double iu_L_log[2][numt/store_n+1], double iu_SR_log[2][numt/store_n+1], double iu_SL_log[2][numt/store_n+1], double iu_SE_log[2][numt/store_n+1], double iu_Si_log[2][numt/store_n+1],double iu_T_log[2][numt/store_n+1], int time_i, int time_TOTj, double prebl[5]){
    int i = 0;
    double prefactor=0.0;
    prefactor = 0.5*dx;

    iu_L_log[0][time_i] = 0.0;
    iu_SR_log[0][time_i] = 0.0;
    iu_SL_log[0][time_i] = 0.0;
    iu_SE_log[0][time_i] = 0.0;
    iu_Si_log[0][time_i] = 0.0;
    
    for(i=i_min; i<i_max; i++){
        iu_L_log[0][time_i] = iu_L_log[0][time_i] + 2.0*u_L[i];
        iu_SR_log[0][time_i] = iu_SR_log[0][time_i] + 2.0*u_SR[i];
        iu_SL_log[0][time_i] = iu_SL_log[0][time_i] + 2.0*u_SL[i];
        iu_SE_log[0][time_i] = iu_SE_log[0][time_i] + 2.0*u_SE[i];
        iu_Si_log[0][time_i] = iu_Si_log[0][time_i] + 2.0*u_Si[i];
    }
    // total pool
//    iu_T_log[0][time_i] = prefactor*(iu_L_log[0][time_i] + iu_SR_log[0][time_i] + iu_SL_log[0][time_i] + iu_SE_log[0][time_i] + iu_Si_log[0][time_i] - (u_L[1] + u_L[i_max-1] + u_SR[1] + u_SR[i_max-1] + u_SL[1] + u_SL[i_max-1] + u_SE[1] + u_SE[i_max-1] + u_Si[1] +u_Si[i_max-1]))/(prebl[0]+prebl[1]+prebl[2]+prebl[3]+prebl[4]);
    iu_T_log[0][time_i] = prefactor*(iu_L_log[0][time_i] + iu_SR_log[0][time_i] + iu_SL_log[0][time_i] + iu_SE_log[0][time_i] + iu_Si_log[0][time_i] - (u_L[i_min] + u_L[i_max-1] + u_SR[i_min] + u_SR[i_max-1] + u_SL[i_min] + u_SL[i_max-1] + u_SE[i_min] + u_SE[i_max-1] + u_Si[i_min] +u_Si[i_max-1]))/(prebl[0]+prebl[1]+prebl[2]+prebl[3]+prebl[4]);
    // other pools
    iu_L_log[0][time_i] =  prefactor*(iu_L_log[0][time_i] - u_L[i_min] - u_L[i_max-1])/prebl[0];
    iu_SR_log[0][time_i] = prefactor*(iu_SR_log[0][time_i] - u_SR[i_min] - u_SR[i_max-1])/prebl[1];
    iu_SL_log[0][time_i] = prefactor*(iu_SL_log[0][time_i] - u_SL[i_min] - u_SL[i_max-1])/prebl[2];
    iu_SE_log[0][time_i] = prefactor*(iu_SE_log[0][time_i] - u_SE[i_min] - u_SE[i_max-1])/prebl[3];
    iu_Si_log[0][time_i] = prefactor*(iu_Si_log[0][time_i] - u_Si[i_min] - u_Si[i_max-1])/prebl[4];
    // time entry
    iu_T_log[1][time_i] = time_TOTj*dt;
    iu_L_log[1][time_i] = time_TOTj*dt;
    iu_SR_log[1][time_i] = time_TOTj*dt;
    iu_SL_log[1][time_i] = time_TOTj*dt;
    iu_SE_log[1][time_i] = time_TOTj*dt;
    iu_Si_log[1][time_i] = time_TOTj*dt;
    
    // printf("time: %.9f  total pool rec:%.9f\n",iu_T_log[1][time_i], iu_T_log[0][time_i]);

}



//returns the mean value of the matrix
double meanMat(double mat[numx]){
    int i = 0;
    double sum = 0.0;
    
    for(i=0; i<numx; i++) sum = sum + mat[i];
    sum = sum/(double numx);
    return (sum);
}//meanMat

// computes the average difference between all positions of given matrix between the previous and current saving points
double mat_dif(double mat1[numx], double mat2[numx]){
    int i = 0;
    double sum = 0.0;
    
    for(i=0; i<numx; i++){
        sum = sum + sqrt(pow((mat1[i] - mat2[i]),2.0));
    }
    
    return(sum/(double numx));
}


int Dynamics::StabilityCheck(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double Su_L[numx], double Su_SR[numx], double Su_SL[numx], double Su_SE[numx], double Su_Si[numx]){
    
    double difL = 0.0;
    double difSR = 0.0;
    double difSL = 0.0;
    double difSE = 0.0;
    double difSi = 0.0;
    double maxDif = 100.0;
    double meanL = 0.0;
    double meanSR = 0.0;
    double meanSL = 0.0;
    double meanSE = 0.0;
    double meanSi = 0.0;
    
    difL = mat_dif(u_L, Su_L);
    difSR = mat_dif(u_SR, Su_SR);
    difSL = mat_dif(u_SL, Su_SL);
    difSE = mat_dif(u_SE, Su_SE);
    difSi = mat_dif(u_Si, Su_Si);

    
    meanL= meanMat(u_L);
    meanSR = meanMat(u_SR);
    meanSL = meanMat(u_SL);
    meanSE = meanMat(u_SE);
    meanSi = meanMat(u_Si);
    
    // find max difference ratio to mean
    maxDif = difL/meanL;
    if(difSR/meanSR>maxDif) maxDif = difSR/meanSR;
    if(difSL/meanSL>maxDif) maxDif = difSL/meanSL;
    if(difSE/meanSE>maxDif) maxDif = difSE/meanSE;
    if(difSi/meanSi>maxDif) maxDif = difSi/meanSi;
    
    //printf("maxDif: %.10f\n", maxDif);
    
    // check if max difference is above stability threshhold
    if(maxDif > epsilon) return(0);
    else return(1);
}//















