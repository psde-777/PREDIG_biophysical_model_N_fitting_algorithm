
import sys
import numpy as np

def find_min(x,y):
    if x < y:
        return x
    elif x == y:
        return x
    else:
        return y


filepath = sys.argv[1]
filename = sys.argv[2]
N_generations = sys.argv[3]



N_families = 10
N_generations = 20
min_vars = np.full(N_families,0.)
family = 0
generation = 0
globalmin = 0
globalfam = 0
globalgen = 0
for i in range(1,N_families+1):
    for j in range(1,N_generations+1):
        try:
            data = np.loadtxt(filepath + str(i) + "/" + filename + str(j) + ".txt")[:,1]
            minimum = data[0]
            if i == 1 and j == 1:
                globalmin = minimum
                globalfam = 1
                globalgen = 1
            family = i
            generation = j
            for k in range(1, data.size):
                oldmin = minimum
                minimum = find_min(minimum,data[k])
                oldglobal = globalmin
                globalmin = find_min(globalmin,data[k])
                if globalmin < oldglobal:
                    globalfam = i
                    globalgen = j
                #if minimum < oldmin:
            #print("Family " + str(family) + "; generation " + str(generation) + ": lowest var = " + str(minimum))
            #min_vars[i] = minimum
        except IOError:
            pass
print("Smallest value in this run: Family " + str(globalfam) + ", generation " + str(globalgen) + "; var = " + str(globalmin))
if globalfam > 0:
    data = np.loadtxt(filepath + str(globalfam) + "/" + filename + str(N_generations) + ".txt")
    mymin = data[0,1]
    for i in range(1,N_generations):
         mymin = find_min(mymin,data[i,1])
    print("Variance of final generation: " + str(mymin))

