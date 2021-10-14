#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"


// ******************************
// Initiating matirces
// ******************************

//initiate profiles
void Init_N(double mat[numx]){
    int i=0;
    
    for (i=0; i<numx; i++) mat[i] = 0.0;
}//Init_N

//initiate profile logs
void Init_NT(double mat[numx][numt/store_n]){
    int i=0;
    int j=0;
    
    for (i=0; i<numx; i++)
        for(j=0; j<numt/store_n; j++) mat[i][j] = 0.0;
}//Init_N

void Init_FRAP_logs(double mat[2][numt/store_n+1]){
    int i=0;
    int j=0;
    
    for (i=0; i<2; i++)
        for(j=0; j<(numt/store_n+1); j++) mat[i][j] = 0.0;
}

//initiate 3vector
void Init_3(double mat[3], double par, double min, double max){
    int i=0;
    
    for (i=0; i<3; i++) mat[i] = 0.0;
    mat[0] = par;
    mat[1] = min;
    mat[2] = max;
}//Init_N

//initiate parans logs
void Init_Params(double mat[14][T_th]){
    int i=0;
    int j=0;
    
    for (i=0; i<14; i++)
        for(j=0; j<T_th; j++) mat[i][j] = 0.0;
}//Init_Params

//initiate
void Init_FRAPrec(double mat[datN]){
    int i=0;
    
    for (i=0; i<datN; i++) mat[i] = 0.0;
}//Init_FRAPrec

// initiate source
void Initiate_source(double source[numx]){
    int i = 0;
    
    for (i=0; i<numx; i++) {
        if(i<=(numx/2 + n_source/2) && i>=(numx/2 - n_source/2)) source[i] = 1.0;
        else source[i] = 0.0;
    }
}// Initiate source

void Initiate_Profile_SS(double mat[numx], int index){
    double value = 0.0;
    double con1=0.0;
    int i =0;
    
    double phi=1.0;
    
    con1=(k*ko + (1.0 + ko)*koff);
    value = phi*dx*exp((0.0*dx + 2.0*dx*(numx-1.0) + dx)/lambdaexp)/(exp(2.0*0.0/lambdaexp)-exp(2.0*dx*(numx-1.0)/lambdaexp))/(1.0-exp(dx/lambdaexp));
    switch (index) {
        case 0: //L
            for (i=0; i<numx; i++) mat[i] = value*exp(-i*dx/lambdaexp);
            break;
        case 1: //SR
            for (i=0; i<numx; i++) mat[i] = value*exp(-i*dx/lambdaexp)*koff*kon*(k*(1.0+exp(-1.0/lambdaexp))+2.0*con1)/4.0/(k+koff)/con1;

            break;
        case 2: //SL
            for (i=0; i<numx; i++) mat[i] = value*exp(-i*dx/lambdaexp)*koff*kon*(k*(1.0+exp(-1.0/lambdaexp))+2.0*con1*exp(-1.0/lambdaexp))/4.0/(k+koff)/con1;

            break;
        case 3: //SE
            for (i=0; i<numx; i++) mat[i] = value*exp(-i*dx/lambdaexp)*koff*kon*k*(1.0+exp(-1.0/lambdaexp))/2.0/con1;
            break;
        default: //SI
            for (i=0; i<numx; i++) mat[i] = (ki/k2)*value*exp(-i*dx/lambdaexp)*koff*kon*k*(1.0+exp(-1.0/lambdaexp))/2.0/con1;
            break;
    }
}
