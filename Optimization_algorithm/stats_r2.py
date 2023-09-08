import numpy as np
import os
from scipy.stats import linregress
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error


# Function to convert scientific notation to float
def sci_to_float(s):
    try:
        return float(s)
    except ValueError:
        # If the conversion fails, try to handle scientific notation
        parts = s.split("e")
        base = float(parts[0])
        exponent = float(parts[1].replace("+", ""))
        return base * 10**exponent


# Read keywords from keywords.txt file
with open('keywords.txt', 'r') as f:
    keywords = [line.strip() for line in f]


# Check for all keywords
for keyword in keywords:
    #both glucose and xylose data
    if os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_glc.txt") and os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_xyl.txt"): #both glucose and xylose data
        #read glucose data
        exp_data_glc = np.loadtxt("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_glc.txt", delimiter="\t", usecols=(0, 1))
        exp_time_glc = exp_data_glc[:, 0]
        exp_release_glc = exp_data_glc[:, 1]
        #read xylose data
        exp_data_xyl = np.loadtxt("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_xyl.txt", delimiter="\t", usecols=(0, 1))
        exp_time_xyl = exp_data_xyl[:, 0]
        exp_release_xyl = exp_data_xyl[:, 1]

        #read model output
        # Load data from the second file, specifying only the first two columns (Time and Glucose)
        model_data = np.loadtxt("BEST_FIT/best_Run/Output/saccharification/average_saccharification_" + keyword + ".txt", delimiter="\t", usecols=(0,1,2,3), converters={1: sci_to_float})
        model_time = model_data[:, 0]
        model_glucose = model_data[:, 1]
        model_xylose = model_data[:, 3]

        #Interpolate model data to match the time points of the experimental data
        interpolated_model_glucose = np.interp(exp_time_glc, model_time, model_glucose)
        interpolated_model_xylose = np.interp(exp_time_xyl, model_time, model_xylose)

        #Calculate stats
        r_squared_glc = r2_score(exp_release_glc, interpolated_model_glucose)
        r_squared_xyl = r2_score(exp_release_xyl, interpolated_model_xylose)
        rmse_glc = np.sqrt(mean_squared_error(exp_release_glc, interpolated_model_glucose))
        rmse_xyl = np.sqrt(mean_squared_error(exp_release_xyl, interpolated_model_xylose))

        # Print R-squared value
        print(f"sample name: {keyword}")
        print("-------------------------------")
        print(f"R-squared Glucose: {r_squared_glc:.4f}")
        print(f"R-squared Xylose: {r_squared_xyl:.4f}")
        print("-------------------------------")
        print(f"R-squared Glucose: {r_squared_glc:.4f}", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))
        print(f"R-squared Xylose: {r_squared_xyl:.4f}", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))


    #Only glucose data
    if os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_glc.txt") and not os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_xyl.txt"): #Only glucose data
        #read glucose data
        exp_data_glc = np.loadtxt("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_glc.txt", delimiter="\t", usecols=(0, 1))
        exp_time_glc = exp_data_glc[:, 0]
        exp_release_glc = exp_data_glc[:, 1]

        #read model output
        # Load data from the second file, specifying only the first two columns (Time and Glucose)
        model_data = np.loadtxt("BEST_FIT/best_Run/Output/saccharification/average_saccharification_" + keyword + ".txt", delimiter="\t", usecols=(0,1,2,3), converters={1: sci_to_float})
        model_time = model_data[:, 0]
        model_glucose = model_data[:, 1]

        #Interpolate model data to match the time points of the experimental data
        interpolated_model_glucose = np.interp(exp_time_glc, model_time, model_glucose)


        #Calculate stats
        r_squared_glc = r2_score(exp_release_glc, interpolated_model_glucose)
        rmse_glc = np.sqrt(mean_squared_error(exp_release_glc, interpolated_model_glucose))

        # Print R-squared value
        print(f"sample name: {keyword}")
        print("-------------------------------")
        print(f"R-squared Glucose: {r_squared_glc:.4f}")
        print("R-squared Xylose: N.A.")
        print("-------------------------------")
        print(f"R-squared Glucose: {r_squared_glc:.4f}", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))
        print("R-squared Xylose: N.A.", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))



    #only xylose data
    if os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_xyl.txt") and not os.path.exists("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_glc.txt"): #only xylose data
        #read xylose data
        exp_data_xyl = np.loadtxt("BEST_FIT/best_Run/Output/expe_data/expe_saccharification_" + keyword + "_xyl.txt", delimiter="\t", usecols=(0, 1))
        exp_time_xyl = exp_data_xyl[:, 0]
        exp_release_xyl = exp_data_xyl[:, 1]

        #read model output
        # Load data from the second file, specifying only the first two columns (Time and Glucose)
        model_data = np.loadtxt("BEST_FIT/best_Run/Output/saccharification/average_saccharification_" + keyword + ".txt", delimiter="\t", usecols=(0,1,2,3), converters={1: sci_to_float})
        model_time = model_data[:, 0]
        model_xylose = model_data[:, 3]

        #Interpolate model data to match the time points of the experimental data
        interpolated_model_xylose = np.interp(exp_time_xyl, model_time, model_xylose)

        #Calculate stats
        r_squared_xyl = r2_score(exp_release_xyl, interpolated_model_xylose)
        rmse_xyl = np.sqrt(mean_squared_error(exp_release_xyl, interpolated_model_xylose))

        # Print R-squared value
        print(f"sample name: {keyword}")
        print("-------------------------------")
        print("R-squared Glucose: N.A.")
        print(f"R-squared Xylose: {r_squared_xyl:.4f}")
        print("-------------------------------")
        print("R-squared Glucose: N.A.", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))
        print(f"R-squared Xylose: {r_squared_xyl:.4f}", file=open('BEST_FIT/fit-stats_'+ keyword +'.txt', 'a'))
