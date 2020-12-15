import numpy as np
import sys
import os.path


if len(sys.argv) != 4:
    print(len(sys.argv));
    print("Please provide the number of files, as well as the name of the input and output file");
    exit();

i = 1;
end = False;



print(sys.argv[2]+str(i)+".txt");

count = 0
countFound = 0
lastCount = 0#number of lines in previous file
max_size = 0
N_files = int(sys.argv[1]);#Number of files
#N_files = 20



my_file = sys.argv[2] + str(1) + ".txt";
#if my_file.is_file(my_file):
if os.path.isfile(my_file):
    countFound += 1
    data = np.loadtxt(sys.argv[2]+str(1)+".txt");
    Ncolums = data.shape[1]
    #data = np.loadtxt(sys.argv[2]+str(1)+".txt");
    max_size = data[:,0].size;
    min_size = max_size;
    final_time = data[-1,0]
    for i in range(2,N_files+1):
        my_file = sys.argv[2] + str(i) + ".txt";
        if os.path.isfile(my_file):
            countFound += 1
            count = 0;
            data = np.loadtxt(sys.argv[2]+str(i)+".txt");
            count = data[:,0].size;
            max_size = max(count,max_size);
            min_size = min(count,min_size);
            final_time += data[-1,0]
            print(count)
            count = lastCount

    print("Max_size: ", max_size)
    final_time /= countFound
    print("final time: " + str(final_time))
    if final_time > 100:
        final_time = 100
    print("final time after if: " + str(final_time))
    outFile = np.full((max_size,Ncolums),0.)


    x = np.linspace(0,final_time,max_size)
    print(x[-1])

#    outFile[:,0] = x[:];
    for i in range(N_files):
        filename = sys.argv[2]+str(i+1)+".txt";    
        my_file = sys.argv[2] + str(i+1) + ".txt";
        if os.path.isfile(my_file):
            data = np.loadtxt(filename);

            for j in range(data[1,:].size):
                outFile[:,j] += np.interp(x,data[:,0],data[:,j])


    #        outFile[1:,:]/=N_files


    outFile[:,:]/=countFound

    np.savetxt(sys.argv[3], outFile[:,:], delimiter="    ")
else:
    outFile = np.full((2,2),0.)
    outFile[0,1] = -9999
    outFile[1,0] = 70
    outFile[1,1] = -9999
    np.savetxt(sys.argv[3], outFile[:,:], delimiter="    ") 
           
