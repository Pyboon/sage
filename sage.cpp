#include<iostream>
#include<fstream>
#include<string>
#include<ctime>
#include "sage.h"
#include "config.h"
#include "util.h"

using namespace std;
using namespace Eigen;


void Sage::InitCIR(){
    int rows = config.sample_codes;
    int cols = config.subchannel_num * config.cycle_num;
    this->CIR = MatrixXcd::Zero(rows, cols);
    cout<<"CIR MAT size: "<<this->CIR.rows()<<" * "<<this->CIR.cols()<<endl;
}
void Sage::ReadCIR(){
    // generate ReCIR and ImCIR path
    string path = config.path_cir;
    char file_ReCIR[128] = "";
    char file_ImCIR[128] = "";
    sprintf(file_ReCIR, "%02d/fCIR_r_%02d_%04d", this->spot, this->spot, this->cycle_idx);
    sprintf(file_ImCIR, "%02d/fCIR_i_%02d_%04d", this->spot, this->spot, this->cycle_idx);
    string path_ReCIR = path + "/" + file_ReCIR;
    string path_ImCIR = path + "/" + file_ImCIR;
    cout<<path_ReCIR<<endl;
    cout<<path_ImCIR<<endl;
    ifstream fReCIR(path_ReCIR, ios::binary);
    ifstream fImCIR(path_ImCIR, ios::binary);
    if(!fReCIR || !fImCIR){
        cerr<<"ReCIR or ImCIR is not exist!"<<endl;
        abort();
    }
    if(!(CIR.rows())){
        InitCIR();
    }
    int rows = CIR.rows();
    int cols_one_cycle = config.subchannel_num;
    if(config.is_pg){
        cols_one_cycle = config.subchannel_withPG_num;
        Util::log("Protection-Guard in CIR!");
    }
    double real_point = 0.0;
    double imag_point = 0.0;
    // in matlab convert, fwrite is save in binary with column order
    // so here read from binary with columns order
    int real_col = 0;
    for(int k = 0; k<config.cycle_num; ++k){
        for(int j = 1; j<=cols_one_cycle; ++j){
                for(int i = 0; i<rows; ++i){
                    fReCIR.read((char*)&real_point, sizeof(double));
                    fImCIR.read((char*)&imag_point, sizeof(double));
                    if(config.is_pg && j%(config.Rx+1)==0){
                        continue; // only read but not write to CIR
                    }else{
                        CIR(i,real_col).real(real_point);
                        CIR(i,real_col).imag(imag_point);
                    }
                }
                if(!(config.is_pg && j%(config.Rx+1)==0)){
                    ++real_col;
                }
            }
    }
    
    fReCIR.close();
    fImCIR.close();
    // for(int i = 0; i<rows; ++i){
    //     for(int j = 0;j<1; ++j){
    //         cout<<i<<","<<j<<":"<<CIR(i,j)<<endl;
    //     }
    // }
}

void Sage::ReadAntennaResponse(){
    // init Antenna array 
    string path = config.path_antenna;    
    for(int i = 0; i < config.Tx; ++i){
        TxArrayEphi.emplace_back(MatrixXcd::Zero(181,91));
        TxArrayEtheta.emplace_back(MatrixXcd::Zero(181,91));
    }
    for(int i = 0; i < config.Rx; ++i){
        RxArrayEphi.emplace_back(MatrixXcd::Zero(181,91));
        RxArrayEtheta.emplace_back(MatrixXcd::Zero(181,91));
    }
    // read antenna response
    double temp_response = 0.0;
    ifstream fTxReEphi(path+"/TxReEphi", ios::binary);
    ifstream fTxImEphi(path+"/TxImEphi", ios::binary);
    ifstream fTxReEtheta(path+"/TxReEtheta", ios::binary);
    ifstream fTxImEtheta(path+"/TxImEtheta", ios::binary);
    if(!fTxReEphi|| !fTxImEphi || !fTxReEtheta || !fTxImEtheta){
        cerr<<"TxEphi or TxEtheta is not exist!"<<endl;
        abort();
    }
    for(int k = 0; k< config.Tx; ++k){
        for(int j = 0; j< 91; ++j){
            for(int i = 0; i< 181; ++i){
                fTxReEphi.read((char*)&temp_response, sizeof(double));
                TxArrayEphi[k](i,j).real(temp_response);
                fTxImEphi.read((char*)&temp_response, sizeof(double));
                TxArrayEphi[k](i,j).imag(temp_response);
                fTxReEtheta.read((char*)&temp_response, sizeof(double));
                TxArrayEtheta[k](i,j).real(temp_response);
                fTxImEtheta.read((char*)&temp_response, sizeof(double));
                TxArrayEtheta[k](i,j).imag(temp_response);
            }
        }
    }

    ifstream fRxReEphi(path+"/RxReEphi", ios::binary);
    ifstream fRxImEphi(path+"/RxImEphi", ios::binary);
    ifstream fRxReEtheta(path+"/RxReEtheta", ios::binary);
    ifstream fRxImEtheta(path+"/RxImEtheta", ios::binary);
    if(!fRxReEphi|| !fRxImEphi || !fRxReEtheta || !fRxImEtheta){
        cerr<<"RxEphi or RxEtheta is not exist!"<<endl;
        abort();
    }
    for(int k = 0; k< config.Rx; ++k){
        for(int j = 0; j< 91; ++j){
            for(int i = 0; i< 181; ++i){
                fRxReEphi.read((char*)&temp_response, sizeof(double));
                RxArrayEphi[k](i,j).real(temp_response);
                fRxImEphi.read((char*)&temp_response, sizeof(double));
                RxArrayEphi[k](i,j).imag(temp_response);
                fRxReEtheta.read((char*)&temp_response, sizeof(double));
                RxArrayEtheta[k](i,j).real(temp_response);
                fRxImEtheta.read((char*)&temp_response, sizeof(double));
                RxArrayEtheta[k](i,j).imag(temp_response);
            }
        }
    }
    // for(int i = 0;i<181; ++i){
    //     cout<<RxArrayEphi[1](i,0)<<endl;
    // }
}

MatrixX2cd Sage::getAntennaResponse(double theta, double phi, const int k){
    int theta_index = (int)((theta+90)/180 * 90);
    int phi_index = (int)((phi+180)/360 * 180);
    // k == 1 : tx antenna response 
    if(k==1){
        MatrixX2cd response(config.Tx, 2);
        for(int i = 0; i < config.Tx; ++i){
            response(i, 1) = TxArrayEphi[i](phi_index, theta_index);
            response(i, 0) = TxArrayEtheta[i](phi_index, theta_index);
        }
        return response;

    } else if(k ==2){ // k ==2 : rx antenna response
        MatrixX2cd response(config.Rx, 2);
        for(int i = 0; i < config.Rx; ++i){
            response(i, 1) = RxArrayEphi[i](phi_index, theta_index);
            response(i, 0) = RxArrayEtheta[i](phi_index, theta_index);
        }
        return response;
    } else {
        cerr<<"parameters error"<<endl;
    }
}

bool Sage::IterConvergence(SageResult last_result, SageResult cur_result){
    VectorXi tau_diff = last_result.tau - cur_result.tau;
    VectorXd doppler_diff = last_result.doppler - cur_result.doppler;
    VectorXd theta_rx_diff = last_result.theta_rx - cur_result.theta_rx;
    VectorXd phi_rx_diff = last_result.phi_rx - cur_result.phi_rx;
    VectorXd theta_tx_diff = last_result.theta_tx - cur_result.theta_tx;
    VectorXd phi_tx_diff = last_result.phi_tx - cur_result.phi_tx;
    bool is_convergence = true;
    for(int iL = 0; iL < max_path_num; ++iL){
        if(abs(tau_diff(iL)) > tol.tau || abs(doppler_diff(iL)) > tol.doppler || 
        abs(theta_rx_diff(iL)) > tol.theta_rx || abs(phi_rx_diff(iL)) > tol.phi_rx ||
        abs(theta_tx_diff(iL)) > tol.theta_tx || abs(phi_tx_diff(iL)) > tol.phi_tx){
            is_convergence = false;
            break;
        }
    }
    return is_convergence;
}