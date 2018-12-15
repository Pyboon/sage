#ifndef SAGE_H_
#define SAGE_H_
#include "Eigen/Eigen"
#include "config.h"
#include <string> 
#include <vector> 

using namespace Eigen;
using namespace std;

class Sage{
protected:
    Config config;
    complex_mat CIR;
    int spot;
    int cycle_idx;
    int max_iter_num;
    int max_path_num;
    vector<complex_mat> TxArrayEphi;
    vector<complex_mat> TxArrayEtheta;
    vector<complex_mat> RxArrayEphi;
    vector<complex_mat> RxArrayEtheta;
    SageResult result;
    ResultTol tol;
public:
    // construction
    Sage(string path_cir, string path_antenna, string path_result, int max_iter_num = 10, int max_path_num = 10):
    max_iter_num(max_iter_num), max_path_num(max_path_num){
        config.path_cir = path_cir;
        config.path_antenna = path_antenna;
        config.path_result = path_result;
        tol.tau = 2;
        tol.doppler = 1;
        tol.phi_rx = 2;
        tol.theta_rx = 2;
        tol.phi_tx = 2;
        tol.theta_tx = 2;
    };     
    // init config
    virtual void ConfigInit(int spot, int cycle_idx, int Tx = 32, int Rx = 56, int pn_code = 127, int cyc_num = 4,
    double fc = 3.5e9, double bd = 100e6, double cycle_rate = 26.981) = 0;
    // init CIR
    virtual void InitCIR();
    // read CIR
    virtual void ReadCIR();
    void CalculatePDP();
    // read antenna response
    virtual void ReadAntennaResponse();
    // get response at one angle
    virtual MatrixX2cd getAntennaResponse(double theta, double phi, const int k = 1);
    // iteration convergence 
    virtual bool IterConvergence(SageResult last_result, SageResult cur_result);
    ~Sage(){};
};




#endif