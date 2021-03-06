#ifndef _UTIL_H_
#define _UTIL_H_

#include<iostream>
#include<iomanip> 
#include<ctime>
#include<cstdio>
#include<string>
#include<chrono>
#include<unistd.h>
#define MAX_LENGTH 1024 
using namespace std;

class Util{
public:
    static void log(const string message){
        char hostname[MAX_LENGTH];
        gethostname(hostname, MAX_LENGTH);
        string cur_time = GetCurTime();
        printf("[%s][%s] INFO:%s\n", hostname, cur_time.c_str(), message.c_str());
    }
    static string  GetCurTime(){
        time_t nowtime;  
        nowtime = time(NULL); 
        char tmp[64];   
        strftime(tmp,sizeof(tmp),"%Y-%m-%d %H:%M:%S",localtime(&nowtime));   
        return tmp;
    }
};

class ProgressBar {
private:
    unsigned int ticks = 0;

    const unsigned int total_ticks;
    const unsigned int bar_width;
    const char complete_char = '=';
    const char incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

public:
    ProgressBar(unsigned int total, unsigned int width, char complete, char incomplete) :
            total_ticks {total}, bar_width {width}, complete_char {complete}, incomplete_char {incomplete} {}

    ProgressBar(unsigned int total, unsigned int width) : total_ticks {total}, bar_width {width} {}

    unsigned int operator++() { return ++ticks; }

    void display() const
    {
        float progress = (float) ticks / total_ticks;
        unsigned int pos = (int) (bar_width * progress);

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now-start_time).count();

        std::cout << "[";

        for (unsigned int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << complete_char;
            else if (i == pos) std::cout << ">";
            else std::cout << incomplete_char;
        }
        std::cout << "] " << int(progress * 100.0) << "% "<<setiosflags(ios::fixed)<<setprecision(3)
        <<float(time_elapsed) / 1000.0 << "s\r";
        std::cout.flush();
    }

    void done() const
    {
        display();
        std::cout << std::endl;
    }
};
#endif