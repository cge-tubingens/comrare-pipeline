"""
Python module to perform sample quality control
"""

import os
import json

import pandas as pd
import numpy as np

from classes.Helpers import shell_do, load_imiss_file, load_het_file

class SampleQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_path:str) -> None:

        # check if paths are set
        if input_path is None or output_path is None:
            raise ValueError("Both input_path and output_path must be set upon initialization.")

        # Check path validity
        bed_path = os.path.join(input_path, input_name + '.bed')
        fam_path = os.path.join(input_path, input_name + '.fam')
        bim_path = os.path.join(input_path, input_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not os.path.exists(input_path) or not os.path.exists(output_path):
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not bed_check:
            raise FileNotFoundError(".bed file not found")
        if not fam_check:
            raise FileNotFoundError(".fam file not found")
        if not bim_check:
            raise FileNotFoundError(".bim file not found")

        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name

        with open(config_path, 'r') as file:
            self.config_dict = json.load(file)

        self.results_dir = os.path.join(output_path, 'sample_qc_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def run_ld_prune(self, ld_region_path:str)->dict:

        """
        Prunes samples based on Linkage Disequilibrium

        Parameters:

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'metrics': Metrics associated with the pruning, such as 'ld_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_path= self.output_path
        output_name= self.output_name

        maf = self.config_dict['maf']
        geno= self.config_dict['geno']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']


        # Check type of maf
        if not isinstance(maf, float):
             raise TypeError("maf should be of type float.")

        # Check type of geno
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")

        # Check type of mind
        if not isinstance(mind, float):
            raise TypeError("mind should be of type float.")
        
        # Check type of hwe
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.")
        
        # Check if maf is in range
        if maf < 0.05 or maf > 0.1:
            raise ValueError("maf should be between 0.05 and 0.1")
        
        # Check if geno is in range
        if geno < 0.05 or geno > 0.1:
            raise ValueError("geno should be between 0.05 and 0.1")
        
        # Check if mind is in range
        if mind < 0.1 or mind > 0.15:
            raise ValueError("mind should be between 0.1 and 0.15")
        
        # Check if hwe is in range
        if hwe < 0.00000001 or hwe > 0.001:
            raise ValueError("hwe should be between 0.00000001 and 0.001")

        step = "ld_prune"

        # generates prune.in and prune.out
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {ld_region_path} --range --indep-pairwise 50 5 0.2 --out {os.path.join(result_path, output_name+'_1')}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(result_path, output_name+'_1.prune.in')} --make-bed --out {os.path.join(output_path, output_name+'_1')}"

        print('command_one', plink_cmd1)
        print('command_two', plink_cmd2)

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        #listOfFiles = [f'{ld_temp}.log', f'{output_path}.log']
        #concat_logs(step, out_path, listOfFiles)

        process_complete = True

        outfiles_dict = {
            'plink_out': output_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict
    
    def run_heterozygosity_rate(self)->dict:

        """
        Function to identify individuals with elevated missing data rates or outlying heterozygosity rate
        """

        output_path= self.output_path
        output_name= self.output_name
        result_path= self.results_dir

        step = "heterozygosity_rate"

        # 
        plink_cmd1 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --missing --out {os.path.join(result_path, output_name)}"

        # 
        plink_cmd2 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --het --autosome --extract {os.path.join(result_path, output_name+'_1.prune.in')} --out {os.path.join(result_path, output_name)}"

        print('command_one', plink_cmd1)
        print('command_two', plink_cmd2)

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        df_het = load_het_file(
            het_path=os.path.join(result_path, output_name+'.het')
        )
        df_imiss = load_imiss_file(
            imiss_path=os.path.join(result_path, output_name+'.imiss')
        )

        # Compute the lower 2 standard deviation bound
        meanHet_lower = df_het['meanHet'].mean() - 2*df_het['meanHet'].std()

        # Compute the upper 2 standard deviation bound
        meanHet_upper = df_het['meanHet'].mean() + 2*df_het['meanHet'].std()

        mask = ((df_imiss['F_MISS']>=0.04) | (df_het['meanHet'] < meanHet_lower) | (df_het['meanHet'] > meanHet_upper))

        df = df_imiss[mask].reset_index(drop=True)
        df = df.iloc[:,0:2].copy()

        del df_het
        del df_imiss

        path_df = os.path.join(result_path, output_name+'.fail-imisshet-qc.txt')

        df.to_csv(path_or_buf=path_df, sep='\t', index=False, header=False)

        del df

        # Creation of cleaned binary file
        plink_cmd3 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --remove {path_df} --make-bed --out {os.path.join(result_path, output_name+'_2')}"

        shell_do(plink_cmd3, log=True)

        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def run_sex_check(self)->dict:

        """
        Function to identify individuals with discordant sex information
        """

        output_name= self.output_name
        result_path= self.results_dir

        sex_check = self.config_dict['sex_check']

        # check type sex_check
        if not isinstance(sex_check, list):
            raise TypeError("sex_check should be a list")
        if len(sex_check)!=2:
            raise ValueError("sex_check must be a list of length 2")
        
        for value in sex_check:
            if not isinstance(value, float):
                raise TypeError("sex_check values should be float")
            if 0 > value or value > 1:
                raise ValueError("sex_check values must be between 0 and 1")
        
        if sum(sex_check) != 1:
            raise ValueError("sex_check values should add to 1")
        
        step = "sex_check"

        # 
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --check-sex {sex_check[0]} {sex_check[1]} --keep-allele-order --extract {os.path.join(result_path, output_name+'_1.prune.in')} --out {os.path.join(result_path, output_name)}"

        shell_do(plink_cmd1, log=True)

        df = pd.read_csv(
            os.path.join(result_path, output_name+'.sexcheck'),
            sep='\s+'
        )

        df_probs = df[df['STATUS']=='PROBLEM'].reset_index(drop=True)
        df_probs.to_csv(
            os.path.join(result_path, output_name+'.sexprobs'),
            index=False,
            sep='\t'
        )

        df_probs = df_probs.iloc[:,0:2].copy()
        df_probs.to_csv(
            os.path.join(result_path, output_name+'.fail-sexcheck-qc.txt'),
            index=False,
            header=False,
            sep=' '
        )

        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --keep-allele-order --remove {os.path.join(result_path, output_name+'.fail-sexcheck-qc.txt')} --make-bed --out {os.path.join(result_path, output_name+'_3')}"

        shell_do(plink_cmd2, log=True)

        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict

    def run_relatedness_prune(self)->dict:

        """
        Function to identify duplicated or related individuals
        """

        result_path= self.results_dir
        output_path= self.output_path
        output_name= self.output_name

        step = "duplicate_relative_prune"

        to_remove = pd.DataFrame(columns=['FID', 'IID'])

        # Run genome
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --extract {os.path.join(result_path, output_name+'_1.prune.in')} --genome --out {os.path.join(output_path, output_name)}"

        # Generate new .imiss file
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --missing --out {os.path.join(output_path, output_name)}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # Load .imiss file
        df_imiss = pd.read_csv(
            os.path.join(output_path, output_name+'.imiss'),
            sep='\s+'
        )
        # Load .genome file
        df_genome = pd.read_csv(
            os.path.join(output_path, output_name+'.genome'),
            sep='\s+'
        )

        # Isolate duplicates or related samples
        df_dup = df_genome[df_genome['PI_HAT']>0.185].reset_index(drop=True)

        df_1 = pd.merge(
            df_dup[['FID1', 'IID1']], 
            df_imiss[['FID', 'IID', 'F_MISS']], 
            left_on=['FID1', 'IID1'],
            right_on=['FID', 'IID']
        ).drop(columns=['FID', 'IID'], inplace=False)

        df_2 = pd.merge(
            df_dup[['FID2', 'IID2']], 
            df_imiss[['FID', 'IID', 'F_MISS']], 
            left_on=['FID2', 'IID2'],
            right_on=['FID', 'IID']
        ).drop(columns=['FID', 'IID'], inplace=False)

        for k in range(len(df_dup)):

            if df_1.iloc[k,2]>df_2.iloc[k,2]:
                print('sample', df_1.iloc[k,0:2].to_list())
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()
            elif df_1.iloc[k,2]<df_2.iloc[k,2]:
                to_remove.loc[k] = df_2.iloc[k,0:2].to_list()
            else:
                to_remove.loc[k] = df_1.iloc[k,0:2].to_list()

        to_remove = to_remove.drop_duplicates(keep='last')
        to_remove.to_csv(
            os.path.join(result_path, output_name+'.fail-IBD1-qc.txt'),
            index=False,
            header=False
        )

        # Create cleaned binary files
        plink_cmd3 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --remove {os.path.join(result_path, output_name+'.fail-IBD1-qc.txt')} --make-bed --out ${os.path.join(result_path, output_name+'_4')}"

        shell_do(plink_cmd3, log=True)

        process_complete = True

        outfiles_dict = {
            'plink_out': result_path
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'output': outfiles_dict
        }

        return out_dict





