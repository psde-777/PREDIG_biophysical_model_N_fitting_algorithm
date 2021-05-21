#ifndef DPLIST_HXX_
#define DPLIST_HXX_
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <list>
#include <algorithm>
#include <functional>
#include <math.h>
#include <sys/time.h>
#include <tuple>
#include <unordered_map>



class DPList {//Structure of the DP distribution at a time t

public:
    double real_time;//real_time at which this is taken
    int t0;//step number
    std::vector<double> DP;//this was an array

    DPList(){
        this->real_time = 0;
        this->t0 = 0;
        //DP = new int[1];
        DP.push_back(0);
    }

    DPList(double real_time, int t0, int length_fibril){
        this->real_time = real_time;
        this->t0 = t0;
//        DP = new int[length_fibril];
        for(int i=0; i<length_fibril;i++)
            DP.push_back(0);
        //    DP[i] = 0;

    }

    DPList(int length_fibril){
        this->real_time = 0;
        this->t0 = 0;
//        DP = new int[length_fibril];
        for(int i = 0;i<length_fibril;i++)
            DP.push_back(0);
  //          DP[i] = 0;
    }


    ~DPList(){
     //   delete[] DP;
    }

};


#endif