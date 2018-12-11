#ifndef ISI_SAGE_H_
#define ISI_SAGE_H_
#include "sage.h"

using namespace Eigen;
using namespace std;

class IsiSage: protected Sage{
public:
    // construction for all params
    IsiSage(int spot, int cycle_idx, int max_iter_num, int max_path_num):
    Sage(spot, cycle_idx, max_iter_num, max_path_num){};
    // construction for required params  
    IsiSage(int sopt, int cycle_idx): Sage(sopt, cycle_idx){};
    void ParamsInit();
    void ParamsIterationUpdate();
    RowVectorXcd InitTau(const complex_mat data, const int iL);
    void InitDoppler(const RowVectorXcd tau_data, const int iL);
    void InitAoA(const RowVectorXcd tau_data, const int iL);
    void InitAoD(const RowVectorXcd tau_data, const int iL);
    void UpdateAlpha(const RowVectorXcd tau_data, const int iL);
    RowVectorXcd SignalConstruct(const int iL);
    RowVectorXcd UpdateTau(complex_mat data, const int iL);
    void UpdateAoA(RowVectorXcd tau_data, const int iL);
    void UpdateAoD(RowVectorXcd tau_data, const int iL);
    void UpdateDoppler(RowVectorXcd tau_data, const int iL);
    double ObjectiveFunction(RowVectorXcd tau_data, double doppler, MatrixX2cd C1,
    MatrixX2cd C2, const int iL);
    void run();
    void SaveResult();
    ~IsiSage(){};

};

#endif