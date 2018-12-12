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
    int rows = this->CIR.rows();
    if(!rows){
        InitCIR();
    }
    int cols = this->CIR.cols();
    if(config.is_pg){
        cols = config.subchannel_withPG_num * config.cycle_num;
        Util::log("Protection-Guard in CIR!");
    }
    double real_point = 0.0;
    double imag_point = 0.0;
    // in matlab convert, fwrite is save in binary with column order
    // so here read from binary with columns order
    int real_col = 0;
    for(int j = 0; j<cols; ++j){
        for(int i = 0; i<rows; ++i){
            fReCIR.read((char*)&real_point, sizeof(double));
            fImCIR.read((char*)&imag_point, sizeof(double));
            if(config.is_pg && j>0 && j%config.Rx==0){
                continue; // only read but not write to CIR
            }else{
                CIR(i,real_col).real(real_point);
                CIR(i,real_col).imag(imag_point);
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
            response(i, 0) = TxArrayEphi[i](phi_index, theta_index);
            response(i, 1) = TxArrayEtheta[i](phi_index, theta_index);
        }
        return response;

    } else if(k ==2){ // k ==2 : rx antenna response
        MatrixX2cd response(config.Rx, 2);
        for(int i = 0; i < config.Rx; ++i){
            response(i, 0) = RxArrayEphi[i](phi_index, theta_index);
            response(i, 1) = RxArrayEtheta[i](phi_index, theta_index);
        }
        return response;
    } else {
        cerr<<"parameters error"<<endl;
    }
}