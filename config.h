#ifndef _CONFIG_H_
#define _CONFIG_H_
#include"Eigen/Dense"

using namespace std;
using namespace Eigen;

typedef MatrixXcd complex_mat;
const double SPEED_OF_LIGHT = 3e8;
const double PI = 3.14159265;

class Config{
public:
    // antenna config
    int Rx;
    int Tx;
    int cycle_num;
    int channel_num;
    int PN_CODE;
    int sample_codes;
    string path_antenna;
    bool is_pg;
    string path_result;
    // channel config
    string path_cir;
    double fc;
    double lambda;
    double bandwidth;
    int subchannel_num;
    int subchannel_withPG_num; 
    double T;
    int Tsc;
    int Tt;
    int Tr;
    int Tcy;
    int switch_interval;
    double cycle_interval;
    double cycle_rate;
    double max_doppler;
    double min_theta_tx;
    double max_theta_tx;
    double min_phi_tx;
    double max_phi_tx;
    double min_theta_rx;
    double max_theta_rx;
    double min_phi_rx;
    double max_phi_rx;
    // quantization steps
    double doppler_step;
    double angle_step;
    int delay_step;
};

struct SageResult{
    VectorXi tau;
    VectorXd doppler;
    VectorXd theta_tx;
    VectorXd phi_tx;
    VectorXd theta_rx;
    VectorXd phi_rx;
    complex_mat alpha;
};
struct ResultTol{
    int tau;
    double doppler;
    double theta_rx;
    double phi_rx;
    double theta_tx;
    double phi_tx;
};

#endif