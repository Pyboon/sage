#include<iostream>
#include<string>
#include"isi_sage.h"
#include"config.h"

using namespace std;
int main(){
    // config init
    string path_cir = "../../data/low_fre_MIMO/";
    string path_antenna = "../../antenna_bin/";
    string path_result  = "../../result/";
    IsiSage test_Sage(path_cir, path_antenna, path_result); 
    test_Sage.ConfigInit(11, 1);
    test_Sage.run();
    return 0;
}