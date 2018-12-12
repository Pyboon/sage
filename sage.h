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
public:
    // construction
    Sage(string path_cir, string path_antenna, string path_result, int max_iter_num = 10, int max_path_num = 10):
    max_iter_num(max_iter_num), max_path_num(max_path_num){
        config.path_cir = path_cir;
        config.path_antenna = path_antenna;
        config.path_result = path_result;
    };     
    // init config
    virtual void ConfigInit(int spot, int cycle_idx, int Tx = 32, int Rx = 56, int pn_code = 127, int cyc_num = 4,
    double fc = 3.5e9, double bd = 100e6, double cycle_rate = 26.981) = 0;
    // init CIR
    void InitCIR();
    // read CIR
    void ReadCIR();
    void CalculatePDP();
    // read antenna response
    void ReadAntennaResponse();
    // get response at one angle
    MatrixX2cd getAntennaResponse(double theta, double phi, const int k = 1);
    ~Sage(){};
};




#endif