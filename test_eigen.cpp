#include<iostream>
#include"Eigen/Dense"
#include"sage.h"

using namespace std;
using namespace Eigen;
typedef MatrixXcd Mat;
int main(){
    complex<double> a;
    a.real(3);
    a.imag(4);
    cout<<norm(a)<<endl;
    Mat cir;
    cout<<cir.rows()<<cir.cols()<<endl;
    VectorXcd test = VectorXcd::Random(3,1);
    cout<<test<<endl;
    VectorXcd test_H = test.adjoint();
    system("pause");
    return 0;
}