#include<iostream>
#include<string>
#include"isi_sage.h"
#include"config.h"

using namespace std;
int main(){

    IsiSage test_Sage(9,1); 
    // config init
    string path_cir = "../data/";
    string path_antenna = "../antenna_bin/";
    string path_result  = "../result/";
    test_Sage.ConfigInit(path_cir, path_antenna, path_result);
    test_Sage.run();
    return 0;
}