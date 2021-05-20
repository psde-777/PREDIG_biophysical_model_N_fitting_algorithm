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
#include <iomanip>
#include <string>
#include <string.h>
#include <list>
#include <algorithm>
#include <functional>
#include <math.h>
#include <sys/time.h>

using namespace std;

double average_time(int N_simus, string input_file){

    double averaged_time = 0.;


    for(int j=1;j<=N_simus;j++){
        cout << j <<  endl;
        ifstream file1(input_file + to_string(j) + ".txt");
        string line;
        line.clear();
        double time = 0.;
        int DP = 0;
        double abundance = 0.;
        if(file1)
        {

            int i = 0;
            while(getline(file1,line))
            {
                istringstream in(line);
                in >> time >> DP >> abundance;
                i++;
    //            cout << "time = " << time << endl;
            }
            averaged_time += time;
        }
        else{
            cout << "File " << input_file << " does not exist. Exiting." << endl;
            exit(1);
        }

    }
    averaged_time/=N_simus;
    return averaged_time;
}


int main(int argc, char *argv[]){

    int N_simus = stoi(argv[1]);
    string input_file = argv[2];//"DP_distrib/DP_distrib_1.txt";
    string output_file = argv[3];
    int N_bins = stoi(argv[4]);
    int length_fibril = stoi(argv[5]);
    double time_fraction_first  = stod(argv[6]);
    double time_fraction_second = 1.-time_fraction_first;

    double averaged_time = average_time(N_simus,input_file);

    cout << "averaged_time = " << averaged_time << endl;

    int len_array = N_bins*(length_fibril+1);

    cout << "len_array = " << len_array << endl;


    vector<double> time_bins(N_bins*(length_fibril+1));
    vector<int> DP_bins(N_bins*(length_fibril+1));
    vector<double> fraction_bins(N_bins*(length_fibril+1));
    vector<int> vals_in_bin(N_bins*(length_fibril+1));

/*    double *time_bins = new double[N_bins*(length_fibril+1)];
    int *DP_bins = new int[N_bins*(length_fibril+1)];
    double *fraction_bins = new double[N_bins*(length_fibril+1)];
    int *vals_in_bin = new int[N_bins*(length_fibril+1)];*/




    double dt1 = time_fraction_first * averaged_time / double(time_fraction_second * N_bins -1);
    double dt2 = time_fraction_second * averaged_time / double(time_fraction_first * N_bins - 1);

    for(int i=0;i<time_fraction_second*len_array;i++){
        time_bins[i] = dt1*(1 + int(i/(length_fibril+1)));
        DP_bins[i] = i%(length_fibril+1);
        fraction_bins[i] = 0.;
        vals_in_bin[i] = 0;
    }
    for(int i=time_fraction_second*len_array;i<len_array;i++){
        time_bins[i] = time_fraction_first*averaged_time + dt2*(1 + int((i-time_fraction_second*len_array)/(length_fibril+1)));
        DP_bins[i] = i%(length_fibril+1);
        fraction_bins[i] = 0.;
        vals_in_bin[i] = 0;
    }
        double time = 0.;
        int DP = 0;
        double abundance = 0.;
        int k = 0;
        int index = 0;


//====== Run through data again ======

    
    for(int j=1;j<=N_simus;j++){
        cout << j <<  endl;
        ifstream file1(input_file + to_string(j) + ".txt");
        string line;
        line.clear();
        time = 0.;
        DP = 0;
        abundance = 0.;
        k = 0;
        index = 0;

        if(file1)
        {

            int i = 0;
            while(getline(file1,line))
            {
                istringstream in(line);
                in >> time >> DP >> abundance;
                if(time <= averaged_time){
                    if(time < time_fraction_first*averaged_time){
                        k = int(time/dt1);
                        index = k*(length_fibril+1)+DP;
                        if(DP_bins[index] == DP and index < len_array){
                            fraction_bins[index] += abundance;
                            vals_in_bin[index] += 1;
                        }
                        else{
                            cout << "whoops" << endl;
                            exit(1);
                        }
                    }
                    else{
                        k = time_fraction_second*N_bins + int((time - time_fraction_first*averaged_time)/dt2);
                        index = k*(length_fibril+1)+DP;
                        if(DP_bins[index] == DP and index < len_array){
                            fraction_bins[index] += abundance;
                            vals_in_bin[index] += 1;
                        }
                        else{
                            cout << "index = " << k*(length_fibril+1)+DP << "; k = " << k << "; N_bins = " << N_bins << endl;
                            cout << "time = " << time << "; DP in array = " << time_bins[k*(length_fibril+1)+DP] << endl;
                            cout << "DP = " << DP << "; DP in array = " << DP_bins[k*(length_fibril+1)+DP] << endl;
                            exit(1);
                        }
                    }                    
                }

                i++;
    //            cout << "time = " << time << endl;
            }
        }
        else{
            cout << "File " << input_file << " does not exist. Exiting." << endl;
            exit(1);
        }

    }






    ofstream output(output_file);
    for(int i=0;i<length_fibril;i++){
        output << 0 << "\t" << i << "\t" << 0 << endl;
    }
    output << 0 << "\t" << length_fibril << "\t" << 1 << endl;
    for(int i=0;i<len_array;i++){
        if(i!=0){
            if(time_bins[i] != time_bins[i-1]){
                output << endl;
            }
        }
        if(vals_in_bin[i] != 0){
            output << time_bins[i] << "\t" << DP_bins[i] << "\t" << fraction_bins[i]/double(vals_in_bin[i]) << endl;
        }
        else{
            output << time_bins[i] << "\t" << DP_bins[i] << "\t" << 0 << endl;
        }

    }    

    output.close();




    return 0;
}
