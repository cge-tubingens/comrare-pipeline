"""
Python module to perform variant quality control
"""

import os
import json
import shutil

import pandas as pd

from classes.Helpers import shell_do


class VariantQC:

    def __init__(self, input_path:str, input_name:str, output_path:str, output_name:str, config_path:str, dependables_path:str) -> None:

        # check if paths are set
        if input_path is None or output_path is None or dependables_path is None:
            raise ValueError("values for input_path, output_path and dependables_path must be set upon initialization.")

        # Check path validity
        bed_path = os.path.join(input_path, input_name + '.bed')
        fam_path = os.path.join(input_path, input_name + '.fam')
        bim_path = os.path.join(input_path, input_name + '.bim')

        bed_check = os.path.exists(bed_path)
        fam_check = os.path.exists(fam_path)
        bim_check = os.path.exists(bim_path)

        if not os.path.exists(input_path) or not os.path.exists(output_path):
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not os.path.exists(dependables_path):
            raise FileNotFoundError("dependables_oath is not a valid path")
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
        self.dependables = dependables_path

        with open(config_path, 'r') as file:
            self.config_dict = json.load(file)

        self.results_dir = os.path.join(output_path, 'sample_qc_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def missing_data_rate(self)->dict:

        """
        Function to identify all markers with an excessive missing rate
        """

        result_path      = self.results_dir
        output_name      = self.output_name

        chr = self.config_dict['chr']

        # check type for chr

        step = 'high_rate_missing_data' 

        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'.clean')} --keep-allele-order --missing --filter-males --chr {chr} --out {os.path.join(result_path, output_name+'.clean_m_only')}"

        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'.clean')} --akeep-allele-order --missing --not-chr {chr} --out {os.path.join(result_path, output_name+'.clean_not_y')}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        df_males = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_m_only.lmiss'),
            sep="\s+"
        )
        df_males = df_males[df_males['F_MISS']>0.2].reset_index(drop=True)
        males_snp = df_males['SNP'].to_list()

        with open(os.path.join(result_path, output_name+'.clean_m_only-fail-lmiss-qc.txt'), 'rw') as file:
            for snp in males_snp:
                file.write('%s\n' % snp)

        df_females = pd.read_csv(
            os.path.join(result_path, output_name+'.clean_not_y.lmiss'),
            sep="\s+"
        )
        df_females = df_females[df_females['F_MISS']>0.2].reset_index(drop=True)
        females_snp = df_females['SNP'].to_list()

        with open(os.path.join(result_path, output_name+'.clean_not_y-fail-lmiss-qc.txt'), 'rw') as file:
            for snp in females_snp:
                file.write('%s\n' % snp)

        fail_snp = males_snp + females_snp
        with open(os.path.join(result_path, output_name+'.clean-fail-lmiss-qc.txt'), 'rw'):
            for snp in fail_snp:
                file.write('%s\n' % snp)

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

    def different_genotype_call_rate(self)->dict:

        """
        Funtion to identify test markers for different genotype call rates between cases and controls
        """

        result_path      = self.results_dir
        output_name      = self.output_name

        step = 'different_genotype_case_control'

        plink_cmd = f"plink --bfile {os.path.join(result_path, output_name+'.clean')} --keep-allele-order --test-missing --out {os.path.join(result_path, output_name+'.clean')}"

        shell_do(plink_cmd, log=True)

        missing_file = os.path.join(result_path, output_name+'.clean.missing')

        # Open the input file for reading
        with open(missing_file, 'r') as infile:
            # Open the output file for writing
            with open(missing_file + ".clean-fail-diffmiss-qc.txt", 'w') as outfile:
                # Iterate over each line in the input file
                for line in infile:
                    # Remove leading whitespace from the line
                    line = line.strip()
                    # Split the line into fields based on whitespace
                    fields = line.split()
                    # Check if the first field is not 'CHR'
                    if fields[0] != 'CHR':
                        # Check if the fifth field is less than 0.00001
                        if float(fields[4]) < 0.00001:
                            # Write the second field to the output file
                            outfile.write(fields[1] + '\n')

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

    def remove_markers(self)->dict:

        """
        Function to remove markers failing quality control
        """

        result_path      = self.results_dir
        output_name      = self.output_name

        maf = self.config_dict['maf']
        geno= self.config_dict['geno']
        mind= self.config_dict['mind']
        hwe = self.config_dict['hwe']

        step = "remove_markers"

        fail_lmiss = os.path.join(result_path, output_name+'.clean-fail-lmiss-qc.txt')
        fail_diffmiss = os.path.join(result_path, output_name+'.clean-fail-diffmiss-qc.txt')

        output_file = os.path.join(result_path, output_name+'.clean-fail-markers-qc.txt')

        # Concatenate files
        with open(output_file, 'wb') as outfile:
            for file_path in [fail_lmiss, fail_diffmiss]:
                with open(file_path, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

        # Sort and remove duplicates
        with open(output_file, 'r') as file:
            lines = sorted(set(file.readlines()))

        # Write sorted and unique content back to the output file
        with open(output_file, 'w') as file:
            file.writelines(lines)

        plink_cmd = f"--bfile {os.path.join(result_path, output_name+'.clean')} --keep-allele-order --exclude {os.path.join(result_path, output_name+'.clean-fail-markers-qc.txt')} --maf {maf} --mind {mind} --hwe {hwe} --geno {geno} --make-bed --out {os.path.join(result_path, output_name+'.clean.final')}"

        shell_do(plink_cmd, log=True)

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
