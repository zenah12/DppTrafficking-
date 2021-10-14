#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"

// ********************************************************************
// FUNTIONS
// This file holds generic functions that will be useful for our analysis
// ********************************************************************

//rand_Normal
double rand_Normal (double mean, double sigma){
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mean + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mean + sigma * (double) X1);
}//rand_Normal

//create random number from 0 to 1
double unifRand()
{
    return rand() / double(RAND_MAX);
}//unifRand

//random integer from rangeLow to rangeHight
int uniform_distribution(int rangeLow, int rangeHigh)
{
    int myRand = (int)rand();
    int range = rangeHigh - rangeLow + 1; //+1 makes it [rangeLow, rangeHigh], inclusive.
    int myRand_scaled = (myRand % range) + rangeLow;
    return myRand_scaled;
}//uniform_distribution

//absolute value of difference between two doubles
double Absolute(double a, double b){
    if(a>b) return(a-b);
    else return(b-a);
}//Absolute

double min2(double a, double b){
    if(a>b) return(b);
    else return(a);
}//min

double max2(double a, double b){
    if(a<b) return(b);
    else return(a);
}//max


void copyNtoN(double original[numx], double copy[numx]){
    int i=0;
    
    for (i=0; i<numx; i++)
        copy[i] = original[i];
    
}













