import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) != 5:
    print(len(sys.argv));
    print("Please provide the number of files, as well as the name of the input and output file, and the fibril length");
    exit();

i = 1;
end = False;



print(sys.argv[2]+str(i)+".txt");

count = 0
lastCount = 0#number of lines in previous file
max_size = 0
N_files = int(sys.argv[1]);#Number of files

fibril_length = int(sys.argv[4])
#N_files = 20
data = np.loadtxt(sys.argv[2]+str(1)+".txt");
Ncolums = data.shape[1]
data = np.loadtxt(sys.argv[2]+str(1)+".txt");
max_size = data[:,0].size;
min_size = max_size;



final_time = 0
final_time += data[-1,0]
for i in range(2,N_files+1):
    try:
        count = 0;
        data = np.loadtxt(sys.argv[2]+str(i)+".txt");
        count = data[:,0].size;
        max_size = max(count,max_size);
        min_size = min(count,min_size);
        final_time += data[-1,0]
        print(count)
        count = lastCount
    except FileNotFoundError:
        print("Last file");
        end = True;

print("Max_size: ", max_size)
final_time /= N_files
print("Final time:", final_time)
if (max_size/(fibril_length+1))/int(max_size/(fibril_length+1)) != 1:
	print("Something funny is going on with the DP_files... stopping")
	print((max_size/(fibril_length+1)))
	exit()

outFile = np.full((max_size,Ncolums),0.)

x = np.linspace(0,final_time,int(max_size/(fibril_length+1)))

for i in range(fibril_length+1):
    outFile[i::fibril_length+1,0] = x[:]


for i in range(N_files):
    filename = sys.argv[2]+str(i+1)+".txt";    
    data = np.loadtxt(filename);
    
    for j in range(data[1,:].size):
        for k in range(fibril_length+1):
            outFile[k::fibril_length+1,j] += np.interp(x,data[k::fibril_length+1,0],data[k::fibril_length+1,j])

#np.savetxt("raw_" + sys.argv[3], outFile[:,:], delimiter="    ")
outFile[:,:]/=N_files




#for i in range(N_files):
#    filename = sys.argv[2]+str(i+1)+".txt";    
#    data = np.loadtxt(filename)
#    plt.plot(data[:,0],data[:,1])
#plt.plot(outFile[:,0], outFile[:,1])
#plt.legend()
#plt.show()

#outFile[:,:]/=N_files
np.savetxt(sys.argv[3], outFile, delimiter=" ")
           
           
           
