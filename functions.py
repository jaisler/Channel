import math
import numpy as np
from scipy import signal
import pandas as pd
import csv

def CalculatePSD(u, dt, params):
    """ Calculate PSD """
    # f contains the frequency components
    # S is the PSD
    pf = []
    pS = []
    fs = []
    #print(math.ceil(params['nfiles']/2))
    # Obtain the time step for each file and adjust parameters
    for i in range(params['nfilepsd']):
        dt[i] *= params['L']/params['Ub']
        fs.append(int(1/dt[i]))
        u[i] *= params['Ub']

    # Calculare PSD
    nameWin = 'hanning'
    for i in range(params['nfilepsd']):
        nx = int(len(u[i])/4)  # /6 number of windows - 3 or 6
        (f, S) = signal.welch(u[i], fs[i], scaling='density', 
                                    window=signal.get_window(nameWin, nx),
                                    nperseg=nx, noverlap=round(nx/2))
        pf.append(f)
        pS.append(S) # 10.0*np.log10(S) 
        del f, S
    return (pf, pS)

def CalculateTwoPointsCorrelation(u, params):
    """ Calculate two-points correlation function """

    # denominator - RMS
    k = 0
    tn = 0
    den = 0
    while (k < len(u)):
        meanz = 0
        for i in range(params['npcorr']): 
            meanz += u[i + k] * u[i + k]
        meanz /= (params['npcorr'])
        den += meanz
        k += params['npcorr']
        tn += 1
    den /= tn

    # numerator
    Rii = np.zeros(params['npcorr']) # correlation
    for l in range(params['npcorr']):
        k = 0
        tn = 0
        while (k < len(u)): # time
            for i in range(params['npcorr']): # space
                if (i + l < params['npcorr']):
                    Rii[l] += u[i + k] * u[i + l + k] 
                else:
                    Rii[l] += u[i + k] * u[i + l + 2 * (params['npcorr'] - 1 - i - l) + k]
            k += params['npcorr']
            tn += 1
        Rii[l] /= tn
        Rii[l] /= (params['npcorr'] * den) 

    return Rii


def CalculateOnePointCorrelation(u, params):
    """ Calculate one-point correlation function """

    # denominator - RMS
    k = 0
    tn = 0
    den = 0
    while (k < len(u)):
        meanz = 0
        for i in range(params['npcorr']): 
            meanz += u[i] * u[i]
        meanz /= (params['npcorr'])
        den += meanz
        k += params['npcorr']
        tn += 1
    den /= tn

    # This is the average of the square of the deviations 
    # of the variable from its mean value at that specific time: 
    tn = int(len(u) / params['npcorr'])
    Rii = np.zeros(tn) # correlation
    meanz = np.zeros(tn) 
    for l in range(tn):
        mean = 0  # Obtain the average in z 
        for i in range(params['npcorr']):
            mean += u[i + l * params['npcorr']]
        meanz[l] /= params['npcorr'] 

    i = 0 
    for l in range(tn):
        Rii[l] = meanz[0] * meanz[l] / den


    return Rii

def WriteVectorsToCSV(filename, headers, vector, matrix):
    """
    Writes two vectors to a CSV file with headers.

    Args:
    filename (str): The name of the file to write to.
    headers (list of str): The headers for the CSV file.
    vector (list): The vector (column).
    matrix (list): The matrix (columns).
    """

    # Open a file in write mode
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the header
        writer.writerow(headers)

        # Write the data
        for v1, v2, v3, v4 in zip(vector, matrix[0], matrix[1], matrix[2]):
            writer.writerow([v1,v2,v3,v4])

def LoadCSV(filename, delimiter=',', header='infer', skiprows=None, usecols=None):
    """
    Load data from a CSV file into a DataFrame.

    Args:
    filename (str): The path to the CSV file to load.
    delimiter (str, optional): The character used to separate fields. Defaults to ','.
    header (int, list of int, or str, optional): Row(s) to use as column names, and the start of data.
        Defaults to 'infer' which means pandas tries to guess the header location.
    skiprows (list-like, int or callable, optional): Line numbers to skip (0-indexed) or number of lines to skip (int) at the start of the file.
    usecols (list-like or callable, optional): Return a subset of the columns by specifying column names or indices.

    Returns:
    pandas.DataFrame: A DataFrame containing the loaded data.
    """
    return pd.read_csv(
        filename,
        delimiter=delimiter,
        header=header,
        skiprows=skiprows,
        usecols=usecols
    )