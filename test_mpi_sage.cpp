// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// An intro MPI hello world program that uses MPI_Init, MPI_Comm_size,
// MPI_Comm_rank, MPI_Finalize, and MPI_Get_processor_name.
//
#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <utility>
#include <stdlib.h>
#include"isi_sage.h"
#include"config.h"
using namespace std;

vector<pair<int, int> > ParseFileNum(vector<int> spots){
    int cycle_start = 1;
    int cycle_end = 197;
    int interval = 4;
    int spots_size = spots.size();
    vector<pair<int, int> > result;
    for(int i_spot = 0; i_spot < spots_size; ++i_spot){
        for(int i_cycle = cycle_start; i_cycle<=cycle_end; i_cycle += interval){
            // 此处将所有spot和cycle拼接在一起，然后均分给多个MPI计算
            result.push_back(make_pair(spots[i_spot], i_cycle));
        }
    }
    return result;
}

int main(int argc, char** argv) {
    vector<int> spots = {11, 14, 17};
    vector<pair<int, int> > file_nums = ParseFileNum(spots);
    int file_size = file_nums.size();
    // config init
    string home_path = getenv("HOME");
    string path_cir = home_path + "/code/sage/data/low_fre_MIMO/cut_cir_bin/";
    string path_antenna = home_path + "/code/sage/data/antenna_bin/";
    string path_result  = home_path + "/code/sage/data/result/";
    IsiSage test_Sage(path_cir, path_antenna, path_result, 20, 50); 

    double start_time, end_time;
    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();

    // Get the number of processes
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int local_size = floor(file_size / (double)comm_size);
    int left = file_size - local_size*comm_size;
    int start = 0;
    int end = 0;
    if(rank<left){
        start = (local_size+1) * rank;
        end = (local_size+1) * (rank+1);
    } else {
        start = (local_size+1)*left + (rank-left)*local_size;
        end = (local_size+1)*left + (rank-left+1)*local_size;
    }
    int spot = 0;
    int cycle_idx = 0;
    int count = 0;
    for(int i = start; i < end; ++i){
        spot = file_nums[i].first;
        cycle_idx = file_nums[i].second;
        printf("running spot_%d_cyc_%d from rank %d \n", spot, cycle_idx, rank);
        test_Sage.ConfigInit(spot, cycle_idx);
        test_Sage.run();
        ++count;
    }
    end_time = MPI_Wtime();
    printf("[%s:rank-%d]Calculation Finished, total cycle:%d, elapsed time:%lfs \n", 
    processor_name, rank, count, end_time - start_time);
    
    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
}
