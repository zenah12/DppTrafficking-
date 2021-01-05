#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"



void  WriteDifferentiationTimes(double mat[numx], int rep){
    char name[550];
    
    FILE *pfile=NULL;
    
    //sprintf(name, "Differentiation_N%.d_bN%.1f_bD%.1f_g%.3f_q%.2f_a%.3f_h%.1f_m1%.1f_b%.1f_Fm%.2f_Fse%.2f_Nth%.4f_mindx%.2f_Te%.2f_Tr%.2f_RS%.d_ain%.5f_bin%.5f_%.d.txt",N, betaN, betaD, Gamma, q, a, h, m1, b, F_mean, F_se, N_thresh, min_dx, t_ext, t_ret, armRetraction,alpha_b,  alpha_f,   rep);
    sprintf(name,"FileName");
    
    pfile = fopen(name, "w+");
    
    int i=0;
    for (i=0; i<numx; i++) fprintf(pfile, "%.5f \n", mat[i]);
    fprintf(pfile, "\n");
    
    fclose(pfile);
}//WritePosition


void  Write_Final_Profile(double mat[numx],  int time, int index){
    char name[250];
    int j=0;
    FILE *pfile=NULL;
    
    
    sprintf(name, "Profile_L%.4f_%.d_%.d.txt", 50.0, time, index);
    
    pfile = fopen(name, "w+");
    
    
    for(j=0; j<numx; j++){
        fprintf(pfile, "%.10f   %.10f\n", j*dx, mat[j]);
    }
    fclose(pfile);
}//WriteDelta

void  Write_Recovery(double mat[2][numt/store_n+1],  int index){
    char name[250];
    int j=0;
    FILE *pfile=NULL;
    
    
    sprintf(name, "Recovery_L%.4f_%.d.txt", 50.0, index);
    
    pfile = fopen(name, "w+");
    
    
    for(j=0; j<numt/store_n; j++){
        fprintf(pfile, "%.10f   %.10f\n", mat[1][j], mat[0][j]);
    }
    fclose(pfile);
}//WriteDelta

void  Write_Params(double mat[14][T_th],  int max, int repeat, int index){
    char name[350];
    int j=0;
    int i=0;
    FILE *pfile=NULL;
    
    
    sprintf(name, "ABC_Paramas_koffmax%.4f_koffmin%.4f_konmax%.4f_konmin%.4fk2min%.4f_k2max%.4f_kimin%.4f_kimax%.4f_kmin%.4f_kmax%.4flambda%.3f_RUN%.d_%.d.txt", koff_max, koff_min, kon_max, kon_min, k2_min, k2_max, ki_min, ki_max,k_min, k_max,  lambdaexp,repeat, 0);
    
    pfile = fopen(name, "w+");
    
    for(i=0; i<max; i++){
        for(j=0; j<14; j++){
            fprintf(pfile, "%.10f   ", mat[j][i]);
        }
        fprintf(pfile,"\n");
    }
    fclose(pfile);
}//WriteDelta
