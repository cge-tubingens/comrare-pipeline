import subprocess
import sys
import os
import pandas as pd
import numpy as np

def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):

    """
    From GenoTools
    """
    
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res = subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')

    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')
    
def load_imiss_file(imiss_path:str)->pd.DataFrame:

    """
    Function to parse .imiss file to a pandas.DataFrame and log transform F_MISS column
    """

    # check path
    check_imiss = os.path.exists(imiss_path)

    if not check_imiss:
        raise FileNotFoundError("The imiss file cannot be found.")
    
    with open(imiss_path, 'r') as file:
        line_num=0
        for line in file:
            if line_num==0:
                frame = pd.DataFrame(columns=line.rstrip().split())
            else:
                frame.loc[line_num-1] = line.rstrip().split()
            line_num+=1

    frame['N_MISS'] = frame['N_MISS'].astype(float)
    frame['N_GENO'] = frame['N_GENO'].astype(float)
    frame['F_MISS'] = frame['F_MISS'].astype(float)

    frame['logF_MISS'] = np.log10(frame['F_MISS'])

    return frame

def load_het_file(het_path:str)->pd.DataFrame:

    """
    Function to parse .het file to a pandas.DataFrame and compute meanHet
    """

    # check path
    check_het = os.path.exists(het_path)

    if not check_het:
        raise FileNotFoundError("The het file cannot be found.")
    
    with open(het_path, 'r') as file:
        line_num=0
        for line in file:
            if line_num==0:
                frame = pd.DataFrame(columns=line.rstrip().split())
            else:
                frame.loc[line_num-1] = line.rstrip().split()
            line_num+=1

    frame['O(HOM)']= frame['O(HOM)'].astype(float)
    frame['E(HOM)']= frame['E(HOM)'].astype(float)
    frame['N(NM)'] = frame['N(NM)'].astype(float)
    frame['F']     = frame['F'].astype(float)

    frame['meanHet'] = (frame['N(NM)']- frame['O(HOM)'])/frame['N(NM)']

    return frame
