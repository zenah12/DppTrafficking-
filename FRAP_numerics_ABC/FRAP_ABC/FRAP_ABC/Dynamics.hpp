class Dynamics
{
private:
	int x;
public:

    // current profiles for pools
    double uL[numx];
    double uSR[numx];
    double uSL[numx];
    double uSE[numx];
    double uSi[numx];
    // source matrix
    double source[numx];
    // saved profiles for pools
    double uL_save[numx];
    double uSR_save[numx];
    double uSL_save[numx];
    double uSE_save[numx];
    double uSi_save[numx];
    // log for normalised pool integrals during FRAP
    double iuL_F_log[2][numt/store_n+1];
    double iuSR_F_log[2][numt/store_n+1];
    double iuSL_F_log[2][numt/store_n+1];
    double iuSE_F_log[2][numt/store_n+1];
    double iuSi_F_log[2][numt/store_n+1];
    double iuT_F_log[2][numt/store_n+1];
    
    // prebleach values
    double prebleach[5];

    void UpdateProfiles(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double source[numx]);
    double Get_dt();
    void Bleach(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx]);
    void PreBleach_Integrate(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double prebl[5]);
    void FRAP_Integrate(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx],
                        double iu_L_log[2][numt/store_n+1], double iu_SR_log[2][numt/store_n+1], double iu_SL_log[2][numt/store_n+1], double iu_SE_log[2][numt/store_n+1], double iu_Si_log[2][numt/store_n+1], double iu_T_log[2][numt/store_n+1], int time_i, int time_TOTj, double prebl[5]);
    int get_i_max();
    int StabilityCheck(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx], double Su_L[numx], double Su_SR[numx], double Su_SL[numx], double Su_SE[numx], double Su_Si[numx]);
    
    double getRhoNum(double u_L[numx], double u_SR[numx], double u_SL[numx], double u_SE[numx], double u_Si[numx]);
};


