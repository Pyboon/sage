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
    Sage(int spot, int cycle_idx, int max_iter_num = 10, int max_path_num = 10):
    spot(spot), cycle_idx(cycle_idx), max_iter_num(max_iter_num), 
    max_path_num(max_path_num){
    };     
    // init config
    void ConfigInit(string path_cir = "../../data/", string path_antenna = "../../antenna_bin", string path_result = "../result", 
    int Tx = 32, int Rx = 56, int pn_code = 127, int cyc_num = 4, double fc = 3.5e9, double bd = 100e6,
    double cycle_rate = 26.981);
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