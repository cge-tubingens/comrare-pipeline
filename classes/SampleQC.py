"""
Python module to perform sample quality control
"""

import os
import json

from classes.Helpers import shell_do

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
    
    def run_heterozygosity_rate(self):

        """
        Function to identify individuals with elevated missing data rates or outlying heterozygosity rate
        """

        input_path = self.input_path
        input_name = self.input_name
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
