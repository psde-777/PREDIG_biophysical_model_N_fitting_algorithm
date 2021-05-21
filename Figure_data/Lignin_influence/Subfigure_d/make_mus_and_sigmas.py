import numpy as np
with open("mus.txt", "w") as f:
    with open("sigmas.txt", "w") as g:
        for i in np.linspace(0.02,1,50):
            f.write("%1.8f" % i + "\n")
            g.write("%1.8f" % i + "\n")
                
