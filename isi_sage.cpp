#include "isi_sage.h"
#include "Eigen/KroneckerProduct"
#include "config.h"
#include "util.h"
#include <fstream>

using namespace std;
using namespace Eigen;

void IsiSage::ConfigInit(int spot, int cycle_idx, int Tx, int Rx, int pn_code, int cyc_num, 
double fc, double bd, double cycle_rate){
    Util::log("config init...");
    this->spot = spot;
    this->cycle_idx = cycle_idx;
    // to be modified config
    config.Tx = Tx;
    config.Rx = Rx;
    config.cycle_num = cyc_num;
    config.PN_CODE = pn_code;
    config.fc = fc;
    config.bandwidth = bd;
    config.cycle_rate = cycle_rate;
    config.is_pg = true;
    // fixed config
    config.sample_codes = config.PN_CODE * 2;
    config.lambda = SPEED_OF_LIGHT / config.fc;    
    config.subchannel_num = config.Tx * config.Rx;
    config.subchannel_withPG_num = config.Tx * (config.Rx + 1);
    // souder 的采样原理设置了一个保护间隔
    config.channel_num = config.Tx * (config.Rx + 1) * config.cycle_num;
    // 此处的T为chip rate 就是一个码片的时间，一个PN序列有PN_CODE长度的码片
    config.T = 1 / config.bandwidth;
    // 收端一个天线阵列的扫描时间，即接收一个PN序列的时间，但是由于需要满足采样定理，因此为2倍的PN序列长度
    // 以下有关时间的都以码片为单位
    config.Tsc = config.sample_codes;
    config.switch_interval = config.Tsc;
    config.Tr = config.Tsc + config.switch_interval;
    config.Tt = config.Tr * (config.Rx + 1);
    config.Tcy = (int)(1/config.cycle_rate/config.T);
    config.cycle_interval = 1/config.cycle_rate;
    config.max_doppler = ceil(config.cycle_rate / 2);
    config.max_phi_rx = 180;
    config.min_phi_rx = -180;
    config.max_theta_rx = 75;
    config.min_theta_rx = -75;
    config.max_phi_tx = 75;
    config.min_phi_tx = -75;
    config.max_theta_tx = 75;
    config.min_theta_tx = -75;
    config.doppler_step = 1; // hz
    config.angle_step = 2; // degree
    config.delay_step = 1; // one T period 
    Util::log("config init done");
}
void IsiSage::ParamsInit(){
    Util::log("sage params init");
    //init result array and CT and CR array
    result.tau = VectorXi::Zero(max_path_num);
    result.doppler = VectorXd::Zero(max_path_num);
    result.theta_rx = VectorXd::Zero(max_path_num);
    result.phi_rx = VectorXd::Zero(max_path_num);
    result.theta_tx = VectorXd::Zero(max_path_num);
    result.phi_tx = VectorXd::Zero(max_path_num);
    result.alpha = complex_mat::Zero(4, max_path_num);

    // init sage
    int iL = 0;
    complex_mat data = CIR;
    while(iL<max_path_num){
        Util::log("============ Initing the " + to_string(iL) + " path =============");
        //init tau
        RowVectorXcd tau_data = InitTau(data, iL);
        // init doppler frequency
        InitDoppler(tau_data, iL);
        // init angle of arrival 
        InitAoA(tau_data, iL);
        // init angle of departure
        InitAoD(tau_data, iL);
        //init alpha
        UpdateAlpha(tau_data, iL);
        //  construct signal 
        data.row(result.tau(iL)-1) -= SignalConstruct(iL);
        cout<<"the path "<<iL+1<<" init result:"<<" tau="<<result.tau(iL)<<" doppler="
        <<result.doppler(iL)<<" phi_rx="<<result.phi_rx(iL)<<" theta_rx="<<result.theta_rx(iL)
        <<" phi_tx="<<result.phi_tx(iL)<<" theta_tx="<<result.theta_tx(iL)<<" alpha="<<result.alpha(0, iL)<<endl;
        ++iL;
    }
}

RowVectorXcd IsiSage::InitTau(const complex_mat data, const int iL){
    Util::log("init tau...");
    // go through all CIR find tau index with max power 
    int min_tau = 1;
    int max_tau = config.sample_codes;
    int step = config.delay_step;
    int tau_size = (max_tau-min_tau)/step + 1;
    VectorXd tau_power_sum(tau_size);
    tau_power_sum.setZero(tau_size);
    for(int i_tau = 0; i_tau < tau_size; ++i_tau){
        tau_power_sum(i_tau) = data.row(i_tau*step).norm();
    }
    // get index of max tau_power_sum
    int max_index = 0;
    tau_power_sum.maxCoeff(&max_index);
    result.tau(iL) = min_tau + max_index * step;
    //return tau_power
    return data.row(max_index);
}

void IsiSage::InitDoppler(const RowVectorXcd tau_data, const int iL){
    Util::log("init doppler...");
    // go through doppler frequecy bound to find index with max value
    double min_doppler = -config.max_doppler;
    double max_doppler = config.max_doppler;
    double step = config.doppler_step; // step of doppler is 0.1 hz
    int doppler_size = (int)(max_doppler-min_doppler)/step + 1;
    VectorXd doppler_data(doppler_size);
    doppler_data.setZero(doppler_size);
    complex<double> cycle_data = 0;
    complex<double> temp = 0;
    for(double i_step = 0; i_step < doppler_size; ++i_step){
        for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
            for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
                cycle_data = 0;
                for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
                    double t = i_cyc * config.cycle_interval + (i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                    temp.real(cos(2*PI*(min_doppler+i_step*step)*t));
                    temp.imag(-sin(2*PI*(min_doppler+i_step*step)*t));
                    cycle_data += temp * tau_data(i_cyc * config.subchannel_num + i_tx * config.Rx + i_rx);
                }
                doppler_data(i_step) +=norm(cycle_data);
            }
        }
    }
    int max_index = 0;
    doppler_data.maxCoeff(&max_index);
    result.doppler(iL) = min_doppler + max_index * config.doppler_step;
}

void IsiSage::InitAoA(const RowVectorXcd tau_data, const int iL){
    Util::log("init AOA...");
    MatrixXcd Y(config.Rx, config.Tx);
    complex<double> cycle_data = 0;
    complex<double> temp = 0;
    for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
        for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
            cycle_data = 0;
            for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
                double t = (i_cyc * config.Tcy + i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                temp.real(cos(2*PI* result.doppler(iL) * t));
                temp.imag(-sin(2*PI* result.doppler(iL) * t));
                cycle_data += temp * tau_data(i_cyc * config.subchannel_num + i_tx * config.Rx + i_rx);
            }
            Y(i_rx, i_tx) = cycle_data;
        }
    }
    double min_theta_rx = config.min_theta_rx;
    double max_theta_rx = config.max_theta_rx;
    double min_phi_rx = config.min_phi_rx;
    double max_phi_rx = config.max_phi_rx;
    double step = config.angle_step;
    int theta_size = (max_theta_rx - min_theta_rx)/step + 1;
    int phi_size = (max_phi_rx - min_phi_rx)/step + 1;
    MatrixXd aoa_data(theta_size, phi_size);
    aoa_data.setZero(theta_size, phi_size);
    Util::log("generating aoa_data...");
    int bar_size = theta_size * phi_size;
    ProgressBar p_bar(bar_size, 50, '#', '-');
    for(int i_theta = 0; i_theta < theta_size; ++i_theta){
        for(int i_phi = 0; i_phi < phi_size; ++i_phi){
                MatrixX2cd C2 = getAntennaResponse(min_theta_rx+i_theta*step, min_phi_rx+i_phi*step, 2);
                VectorXcd c_2_1 = C2.col(0) / (C2.col(0).norm());
                VectorXcd c_2_2 = C2.col(1) / (C2.col(1).norm());
                RowVectorXcd c_2_1_H = c_2_1.adjoint();
                RowVectorXcd c_2_2_H = c_2_2.adjoint();
            for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
                VectorXcd y_m1 = Y.col(i_tx);
                complex<double> part1 = c_2_1_H * y_m1;
                complex<double> part2 = c_2_2_H * y_m1;
                double p1_norm = norm(part1);
                double p2_norm = norm(part2);
                RowVectorXcd y_m1_H = y_m1.adjoint();
                complex<double> part3 = y_m1_H * c_2_2 * c_2_1_H * y_m1 * c_2_2_H * c_2_1;
                aoa_data(i_theta, i_phi) += p1_norm + p2_norm - 2 * part3.real();
            }
            ++p_bar;
        }
        p_bar.display();
    }
    p_bar.done();
    int i,j;
    aoa_data.maxCoeff(&i, &j);
    result.theta_rx(iL) = i * config.angle_step + min_theta_rx;
    result.phi_rx(iL) = j * config.angle_step + min_phi_rx;
}

void IsiSage::InitAoD(const RowVectorXcd tau_data, const int iL){
    Util::log("init AOD");
    // X_l 
    MatrixXcd X(config.Rx, config.Tx);
    complex<double> cycle_data = 0;
    complex<double> temp = 0;
    for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
        for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
            cycle_data = 0;
            temp = 0;
            for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
                double t = (i_cyc * config.Tcy + i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                temp.real(cos(2*PI* result.doppler(iL) * t));
                temp.imag(sin(-2*PI* result.doppler(iL) * t));
                cycle_data += temp * tau_data(i_cyc * config.subchannel_num + i_tx * config.Rx + i_rx);
            }
            X(i_rx, i_tx) = cycle_data;
        }
    }
    // c_2_p
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    VectorXcd c_2_1 = C2.col(0);
    VectorXcd c_2_2 = C2.col(1);
    Vector4cd f;
    Matrix2Xcd f_temp(2, config.Tx);
    f_temp.row(0) = c_2_1.adjoint() * X;
    f_temp.row(1) = c_2_2.adjoint() * X;
    double min_theta_tx = config.min_theta_tx;
    double max_theta_tx = config.max_theta_tx;
    double min_phi_tx = config.min_phi_tx;
    double max_phi_tx = config.max_phi_tx;
    double step = config.angle_step;
    int theta_size = (max_theta_tx - min_theta_tx)/step;
    int phi_size = (max_phi_tx - min_phi_tx)/step;
    MatrixXd aod_data(theta_size, phi_size);
    for(int i_theta = 0; i_theta < theta_size; ++i_theta){
        for(int i_phi = 0; i_phi < phi_size; ++i_phi){
            MatrixX2cd C1 = getAntennaResponse(i_theta * step + min_theta_tx, i_phi * step + min_phi_tx, 1);
            VectorXcd c_1_1_conj = C1.col(0).conjugate();
            VectorXcd c_1_2_conj = C1.col(1).conjugate();
            f(0) = f_temp.row(0) * c_1_1_conj;
            f(1) = f_temp.row(0) * c_1_2_conj;
            f(2) = f_temp.row(1) * c_1_1_conj;
            f(3) = f_temp.row(1) * c_1_2_conj;
            Matrix2cd part1 = C2.adjoint() * C2;
            Matrix2cd part2 = C1.adjoint() * C1;
            MatrixXcd D = kroneckerProduct(part1, part2);
            complex<double> z = f.adjoint() * D.inverse() * f;
            aod_data(i_theta, i_phi) = norm(z);
        }
    }
    int i,j;
    aod_data.maxCoeff(&i, &j);
    result.theta_tx(iL) = i * config.angle_step + min_theta_tx;
    result.phi_tx(iL) = j * config.angle_step + min_phi_tx;
}

void IsiSage::UpdateAlpha(const RowVectorXcd tau_data, const int iL){
    Util::log("updating Alpha of the " + to_string(iL) + "th path");
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), result.phi_tx(iL), 1);
    Matrix2cd part1 = C2.adjoint() * C2;
    Matrix2cd part2 = C1.adjoint() * C1;
    MatrixXcd D = kroneckerProduct(part1, part2);
    // X_l 
    MatrixXcd X(config.Rx, config.Tx);
    complex<double> cycle_data = 0;
    complex<double> temp = 0;
    for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
        for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
            cycle_data = 0;
            temp = 0;
            for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
                double t = (i_cyc * config.Tcy + i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                temp.real(cos(2*PI* result.doppler(iL) *t));
                temp.imag(sin(2*PI* result.doppler(iL) *t));
                cycle_data += temp * tau_data(i_cyc * config.subchannel_num + i_tx * config.Rx + i_rx);
            }
            X(i_rx, i_tx) = cycle_data;
        }
    }
    VectorXcd c_2_1 = C2.col(0);
    VectorXcd c_2_2 = C2.col(1);
    Vector4cd f;
    VectorXcd c_1_1_conj = C1.col(0).conjugate();
    VectorXcd c_1_2_conj = C1.col(1).conjugate();
    f(0) = c_2_1.adjoint() * X * c_1_1_conj;
    f(1) = c_2_1.adjoint() * X * c_1_2_conj;
    f(2) = c_2_2.adjoint() * X * c_1_1_conj;
    f(3) = c_2_2.adjoint() * X * c_1_2_conj;
    result.alpha.col(iL) = D.inverse() * f / config.cycle_num;
}

RowVectorXcd IsiSage::SignalConstruct(const int iL){
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), result.phi_tx(iL), 1);
    RowVectorXcd signal(config.subchannel_num * config.cycle_num);
    complex<double> temp = 0;
    for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
        for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
            for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
                double t = (i_cyc * config.Tcy + i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                temp.real(cos(2*PI* result.doppler(iL) *t));
                temp.imag(sin(2*PI* result.doppler(iL) *t));
                complex<double> pattern_signal = 0;
                for(int p1 = 0; p1 < 2; ++p1){
                    for(int p2 = 0; p2 < 2; ++p2){
                        pattern_signal += result.alpha(2*p1+p2,iL) * C2(i_rx, p1) * C1(i_tx, p2);
                    }
                }
                signal(i_cyc*config.subchannel_num+i_tx*config.Rx+i_rx) = temp * pattern_signal;
            }
        }
    }
    return signal;
}

void IsiSage::ParamsIterationUpdate(){
    Util::log("SAGE iteration step for updating params");
    int n_step = 1;
    SageResult last_update_result; 
    while(n_step <= max_iter_num){
        Util::log("====== Starting the " + to_string(n_step) + "th iteration step ======");
        // save init result of all path
        last_update_result = result;
        complex_mat y;
        for(int i_path = 0; i_path < max_path_num; ++i_path){
            y = CIR;
            //========== Expectation Step ==================
            for(int j = 0; j<max_path_num; ++j){
                if(j!=i_path){
                    y.row(result.tau(j)-1) -= SignalConstruct(j);
                }
            }
            //========== Maximization Step =================
            // update tau 
            RowVectorXcd tau_data = UpdateTau(y, i_path);
            // update AoA
            UpdateAoA(tau_data, i_path);
            // update AoD
            UpdateAoD(tau_data, i_path);
            // update doppler
            UpdateDoppler(tau_data, i_path);
            // update alpha 
            UpdateAlpha(tau_data, i_path);

        }
        // judge coverage
        if(IterConvergence(last_update_result, result)){
            Util::log("Iteration has coveraged !");
            break;
        }

        ++n_step;
    }
}
double IsiSage::ObjectiveFunction(RowVectorXcd tau_data, double doppler, MatrixX2cd C1, MatrixX2cd C2, const int iL){
    // X 
    MatrixXcd X(config.Rx, config.Tx);
    complex<double> cycle_data = 0;
    complex<double> temp = 0;
    for(int i_tx = 0; i_tx < config.Tx; ++i_tx){
        for(int i_rx = 0; i_rx < config.Rx; ++i_rx){
            cycle_data = 0;
            temp = 0;
            for(int i_cyc = 0; i_cyc < config.cycle_num; ++i_cyc){
                double t = (i_cyc * config.Tcy + i_tx * config.Tt + i_rx * config.Tr) * config.T; 
                temp.real(cos(2*PI* doppler *t));
                temp.imag(-sin(2*PI* doppler *t));
                cycle_data += temp * tau_data(i_cyc * config.subchannel_num + i_tx * config.Rx + i_rx);
            }
            X(i_rx, i_tx) = cycle_data;
        }
    }
    // f 
    RowVectorXcd c_2_1_H = C2.col(0).adjoint();
    RowVectorXcd c_2_2_H = C2.col(1).adjoint();
    VectorXcd c_1_1_conj = C1.col(0).conjugate();
    VectorXcd c_1_2_conj = C1.col(1).conjugate();
    Vector4cd f;
    f(0) = c_2_1_H * X * c_1_1_conj;
    f(1) = c_2_1_H * X * c_1_2_conj;
    f(2) = c_2_2_H * X * c_1_1_conj;
    f(3) = c_2_2_H * X * c_1_2_conj;
    // D
    Matrix2cd part1 = C2.adjoint() * C2;
    Matrix2cd part2 = C1.adjoint() * C1;
    MatrixXcd D = kroneckerProduct(part1, part2);
    complex<double> z = f.adjoint() * D.inverse() * f;
    return sqrt(norm(z));

}
RowVectorXcd IsiSage::UpdateTau(complex_mat data, const int iL){
    Util::log("updating Tau of the " + to_string(iL) + "th path");
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), result.phi_tx(iL), 1);
    int min_tau = 1;
    int max_tau = config.sample_codes;
    int step = config.delay_step;
    int tau_size = (max_tau-min_tau)/step + 1;
    VectorXd tau_data(tau_size);
    for(int i_tau = 0; i_tau < tau_size; ++i_tau){
        tau_data(i_tau) = ObjectiveFunction(data.row(i_tau*step), result.doppler(iL), C1, C2, iL);
    }
    int max_index = 0;
    tau_data.maxCoeff(&max_index);
    result.tau(iL) = min_tau + max_index*step;
    return data.row(max_index);
}

void IsiSage::UpdateAoA(RowVectorXcd tau_data, const int iL){
    Util::log("updating AoA of the " + to_string(iL) + "th path");
    MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), result.phi_tx(iL), 1);
    double min_theta_rx = config.min_theta_rx;
    double max_theta_rx = config.max_theta_rx;
    double min_phi_rx = config.min_phi_rx;
    double max_phi_rx = config.max_phi_rx;
    double step = config.angle_step;
    int theta_size = (max_theta_rx - min_theta_rx)/step+1;
    int phi_size = (max_phi_rx - min_phi_rx)/step+1;
    VectorXd theta_data(theta_size);
    for(int i_theta = 0; i_theta < theta_size; ++i_theta){
        MatrixX2cd C2 = getAntennaResponse(i_theta * step + min_theta_rx, result.phi_rx(iL), 2);
        theta_data(i_theta) = ObjectiveFunction(tau_data, result.doppler(iL), C1, C2, iL);
    }
    int max_index = 0;
    theta_data.maxCoeff(&max_index);
    result.theta_rx(iL) = max_index * config.angle_step + min_theta_rx;

    VectorXd phi_data(phi_size);
    int i_phi = 0;
    for(i_phi = 0; i_phi < phi_size; ++i_phi){
        MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), i_phi * step + min_phi_rx, 2);
        phi_data(i_phi) = ObjectiveFunction(tau_data, result.doppler(iL), C1, C2, iL);
    }
    phi_data.maxCoeff(&max_index);
    result.phi_rx(iL) = max_index * config.angle_step + min_phi_rx;
}

void IsiSage::UpdateAoD(RowVectorXcd tau_data, const int iL){
    Util::log("updating AoD of the " + to_string(iL) + "th path");
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    double min_theta_tx = config.min_theta_tx;
    double max_theta_tx = config.max_theta_tx;
    double min_phi_tx = config.min_phi_tx;
    double max_phi_tx = config.max_phi_tx;
    double step = config.angle_step;
    int theta_size = (max_theta_tx - min_theta_tx)/step+1;
    int phi_size = (max_phi_tx - min_phi_tx)/step+1;
    VectorXd theta_data(theta_size);
    for(int i_theta = 0; i_theta < theta_size; ++i_theta){
        MatrixX2cd C1 = getAntennaResponse(i_theta * step + min_theta_tx, result.phi_tx(iL), 1);
        theta_data(i_theta) = ObjectiveFunction(tau_data, result.doppler(iL), C1, C2, iL);
    }
    int max_index = 0;
    theta_data.maxCoeff(&max_index);
    result.theta_tx(iL) = max_index * config.angle_step + min_theta_tx;

    VectorXd phi_data(phi_size);
    for(int i_phi = 0; i_phi < phi_size; ++i_phi){
        MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), i_phi* step + min_phi_tx, 1);
        phi_data(i_phi) = ObjectiveFunction(tau_data, result.doppler(iL), C1, C2, iL);
    }
    phi_data.maxCoeff(&max_index);
    result.phi_tx(iL) = max_index * config.angle_step + min_phi_tx;
}

void IsiSage::UpdateDoppler(RowVectorXcd tau_data, const int iL){
    Util::log("updating Doppler of the " + to_string(iL) + "th path");
    MatrixX2cd C1 = getAntennaResponse(result.theta_tx(iL), result.phi_tx(iL), 1);
    MatrixX2cd C2 = getAntennaResponse(result.theta_rx(iL), result.phi_rx(iL), 2);
    // go through doppler frequecy bound to find index with max value
    double min_doppler = -config.max_doppler;
    double max_doppler = config.max_doppler;
    double step = config.doppler_step; // step of doppler is 0.1 hz
    int doppler_size = (int)(max_doppler-min_doppler)/step+1;
    VectorXd doppler_data(doppler_size);
    for(double i_step = 0; i_step < doppler_size; ++i_step){
        doppler_data(i_step) = ObjectiveFunction(tau_data, min_doppler+i_step * step, C1, C2, iL);
    }
    int max_index = 0;
    doppler_data.maxCoeff(&max_index);
    result.doppler(iL) = min_doppler+max_index * config.doppler_step;
}

void IsiSage::SaveResult(const string save_path){
    Util::log("Saving Sage Result at " + save_path);
    ofstream f_result(save_path, ios::binary);
    if(!f_result){
        cerr<<"sage result save Error!!!";
        abort();
    }
    //save sage metainfo 
    f_result.write((char*)&spot, sizeof(int));
    f_result.write((char*)&cycle_idx, sizeof(int));
    f_result.write((char*)&max_path_num, sizeof(int));
    f_result.write((char*)&max_iter_num, sizeof(int));
    // save sage result data
    for(int iL = 0; iL < max_path_num; ++iL){
        //save tau of the iL path
        int tau = result.tau(iL);
        f_result.write((char*)&tau, sizeof(int));
        //save doppler of the iL path
        double doppler = result.doppler(iL);
        f_result.write((char*)&doppler, sizeof(double));
        //save AOA of the iL path
        double theta = result.theta_rx(iL);
        double phi = result.phi_rx(iL);
        f_result.write((char*)&theta, sizeof(double));
        f_result.write((char*)&phi, sizeof(double));
        //save AoD of the iL path
        theta = result.theta_tx(iL);
        phi = result.phi_tx(iL);
        f_result.write((char*)&theta, sizeof(double));
        f_result.write((char*)&phi, sizeof(double));
        //save alpha of the iL path
        double alpha_real = result.alpha(0, iL).real();
        double alpha_imag = result.alpha(0, iL).imag();
        for(int i = 0; i < 4; ++i){
            alpha_real = result.alpha(i, iL).real();
            f_result.write((char*)&alpha_real, sizeof(double));
            alpha_imag = result.alpha(i, iL).imag();
            f_result.write((char*)&alpha_imag, sizeof(double));
        }
    }
}

void IsiSage::run(){
    // antenna response and CIR init
    ReadCIR();
    ReadAntennaResponse();
    ParamsInit();
    ParamsIterationUpdate();
    string save_path = config.path_result + "/spot_" + to_string(spot) + "_" + to_string(cycle_idx);
    SaveResult(save_path);
}
