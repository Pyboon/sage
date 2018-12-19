#include<iostream>
#include<string>
#include"isi_sage.h"
#include"config.h"

using namespace std;
int main(){
    // config init
    string path_cir = "~/code/data/low_fre_MIMO/cut_cir_bin/";
    string path_antenna = "~/code/data/antenna_bin/";
    string path_result  = "~/code/data/result/";
    IsiSage test_Sage(path_cir, path_antenna, path_result, 1, 1); 
    test_Sage.ConfigInit(11, 1);
    test_Sage.run();
    return 0;
}