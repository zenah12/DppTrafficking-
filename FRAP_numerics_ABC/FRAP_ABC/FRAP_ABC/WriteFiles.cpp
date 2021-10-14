#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"

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
}//Write_Params
