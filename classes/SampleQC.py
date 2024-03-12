"""
Python module to perform sample quality control
"""

import os
import json

import pandas as pd

from classes.Helpers import shell_do

class SampleQC:

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

        # open config file
        with open(config_path, 'r') as file:
            self.config_dict = json.load(file)

        # create results folder
        self.results_dir = os.path.join(output_path, 'sample_qc_results')
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

    def run_ld_prune(self, ld_region_file:str)->dict:

        """
        Funtion to prunes samples based on Linkage Disequilibrium

        Parameters:
        - ld_region_file: string
            file name with regions with high Linkage Distribution

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path       = self.input_path
        input_name       = self.input_name
        result_path      = self.results_dir
        output_path      = self.output_path
        output_name      = self.output_name
        dependables_path = self.dependables

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
        
        # check existence of high LD regions file
        high_ld_regions_file = os.path.join(dependables_path, ld_region_file)
        if not os.path.exists(high_ld_regions_file):
            raise FileNotFoundError("File with high LD region was not found")

        step = "ld_prune"

        # generates prune.in and prune.out
        plink_cmd1 = f"plink --bfile {os.path.join(input_path, input_name)} --maf {maf} --geno {geno} --mind {mind} --hwe {hwe} --exclude {high_ld_regions_file} --range --indep-pairwise 50 5 0.2 --out {os.path.join(result_path, output_name+'_1')}"

        # prune and creates a filtered binary file
        plink_cmd2 = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --extract {os.path.join(result_path, output_name+'_1.prune.in')} --make-bed --out {os.path.join(output_path, output_name+'_1')}"

        # execute Plink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # report
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
        Function to identify individuals with elevated missing data rates or outlying heterozygosity rate.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        output_path= self.output_path
        output_name= self.output_name
        result_path= self.results_dir

        step = "heterozygosity_rate"

        # 
        plink_cmd1 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --missing --out {os.path.join(result_path, output_name)}"

        # 
        plink_cmd2 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --het --autosome --extract {os.path.join(result_path, output_name+'_1.prune.in')} --out {os.path.join(result_path, output_name)}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # load .het and .imiss files
        df_het = pd.read_csv(
            os.path.join(result_path, output_name+'.het'),
            sep="\s+"
        )
        df_imiss = pd.read_csv(
            os.path.join(result_path, output_name+'.imiss'),
            sep="\s+"
        )

        # compute Het mean
        df_het['meanHet'] = (df_het['N(NM)']-df_het['O(HOM)'])/df_het['N(NM)']
    
        # compute the lower 2 standard deviation bound
        meanHet_lower = df_het['meanHet'].mean() - 2*df_het['meanHet'].std()

        # compute the upper 2 standard deviation bound
        meanHet_upper = df_het['meanHet'].mean() + 2*df_het['meanHet'].std()

        # filter samples
        mask = ((df_imiss['F_MISS']>=0.04) | (df_het['meanHet'] < meanHet_lower) | (df_het['meanHet'] > meanHet_upper))

        df = df_imiss[mask].reset_index(drop=True)
        df = df.iloc[:,0:2].copy()

        del df_het
        del df_imiss

        # save samples that failed imiss-het QC
        path_df = os.path.join(result_path, output_name+'.fail-imisshet-qc.txt')
        df.to_csv(
            path_or_buf =path_df, 
            sep         ='\t', 
            index       =False, 
            header      =False
        )

        del df

        # create cleaned binary file
        plink_cmd3 = f"plink --bfile {os.path.join(output_path, output_name+'_1')} --keep-allele-order --remove {path_df} --make-bed --out {os.path.join(result_path, output_name+'_2')}"

        # execute PLink command
        shell_do(plink_cmd3, log=True)

        # report
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
        Function to identify individuals with discordant sex information.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
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

        # execute PLink command
        shell_do(plink_cmd1, log=True)

        # load file with sex analysis
        df = pd.read_csv(
            os.path.join(result_path, output_name+'.sexcheck'),
            sep='\s+'
        )

        # filter problematic samples and save file
        df_probs = df[df['STATUS']=='PROBLEM'].reset_index(drop=True)
        df_probs.to_csv(
            os.path.join(result_path, output_name+'.sexprobs'),
            index =False,
            sep   ='\t'
        )

        # save IDs of samples who failed sex check QC
        df_probs = df_probs.iloc[:,0:2].copy()
        df_probs.to_csv(
            os.path.join(result_path, output_name+'.fail-sexcheck-qc.txt'),
            index  =False,
            header =False,
            sep    =' '
        )

        # 
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_2')} --keep-allele-order --remove {os.path.join(result_path, output_name+'.fail-sexcheck-qc.txt')} --make-bed --out {os.path.join(result_path, output_name+'_3')}"

        # execute PLink command
        shell_do(plink_cmd2, log=True)

        # report
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
        Function to identify duplicated or related individuals.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        result_path= self.results_dir
        output_path= self.output_path
        output_name= self.output_name

        step = "duplicate_relative_prune"

        to_remove = pd.DataFrame(columns=['FID', 'IID'])

        # run genome
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --extract {os.path.join(result_path, output_name+'_1.prune.in')} --genome --out {os.path.join(output_path, output_name)}"

        # Generate new .imiss file
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --missing --out {os.path.join(output_path, output_name)}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # load .imiss file
        df_imiss = pd.read_csv(
            os.path.join(output_path, output_name+'.imiss'),
            sep='\s+'
        )
        # load .genome file
        df_genome = pd.read_csv(
            os.path.join(output_path, output_name+'.genome'),
            sep='\s+'
        )

        # isolate duplicates or related samples
        df_dup = df_genome[df_genome['PI_HAT']>0.185].reset_index(drop=True)

        df_1 = pd.merge(
            df_dup[['FID1', 'IID1']], 
            df_imiss[['FID', 'IID', 'F_MISS']], 
            left_on =['FID1', 'IID1'],
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
            header=False,
            sep=" "
        )

        # create cleaned binary files
        plink_cmd3 = f"plink --bfile {os.path.join(result_path, output_name+'_3')} --keep-allele-order --remove {os.path.join(result_path, output_name+'.fail-IBD1-qc.txt')} --make-bed --out {os.path.join(result_path, output_name+'_4')}"

        # execute PLink command
        shell_do(plink_cmd3, log=True)

        # report
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

    def delete_failing_QC(self)->None:

        """
        Function to remove samples that failed quality control.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        input_path = self.input_path
        input_name = self.input_name
        result_path= self.results_dir
        output_name= self.output_name

        step = "delete_sample_failed_QC"


        # load files with samples who failed one or several QC steps
        df_sex = pd.read_csv(
            os.path.join(result_path, output_name+'.fail-sexcheck-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        df_imiss = pd.read_csv(
            os.path.join(result_path, output_name+'.fail-imisshet-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        df_ibd1 = pd.read_csv(
            os.path.join(result_path, output_name+'.fail-IBD1-qc.txt'),
            sep      =' ',
            index_col=False,
            header   =None
        )

        # concatenate all failings samples
        df = pd.concat([df_sex, df_imiss, df_ibd1], axis=0)

        # sort values by the first column
        df.sort_values(by=df.columns[0], inplace=True)

        # drop duplicates
        df = df.drop_duplicates(keep='first')

        # save file with samples who failed QC
        df.to_csv(
            os.path.join(result_path, output_name+'.fail-qc_1-inds.txt'),
            sep   =' ',
            header=False,
            index =False
        )

        # 
        plink_cmd = f"plink --bfile {os.path.join(input_path, input_name)} --keep-allele-order --remove {os.path.join(result_path, output_name+'.fail-qc_1-inds.txt')} --make-bed --out {os.path.join(result_path, output_name+'.pre_ind_clean')}"

        # execute PLink command
        shell_do(plink_cmd, log=True)

        # report
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

    def divergent_ancestry_step_one(self, hapmap3r2_no_at_cg_snps_file:str, hapmap3r2_founders_name:str)->dict:

        """
        Function to identify subject with divergent ancestry.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('ld_prune').
            * 'output': Dictionary containing paths to the generated output files.
        """

        result_path= self.results_dir
        output_name= self.output_name
        dependables_path = self.dependables

        maf = self.config_dict['maf']

        # path to auxiliary files
        hapmap3r2_txt = os.path.join(dependables_path, hapmap3r2_no_at_cg_snps_file)
        hapmap_bed = os.path.join(dependables_path, hapmap3r2_founders_name+'.bed')
        hapmap_bim = os.path.join(dependables_path, hapmap3r2_founders_name+'.bim')
        hapmap_fam = os.path.join(dependables_path, hapmap3r2_founders_name+'.fam')

        # check if hapmap3r2 exists
        if not os.path.exists(hapmap3r2_txt):
            raise FileNotFoundError("hapmap3r2_no_at_cg_snps_file not found")
        
        # check if hapmap3r2 founders bed, bim, fam files exist
        if not os.path.exists(hapmap_bed):
            raise FileNotFoundError("hapmap_founders.bed file not found")
        if not os.path.exists(hapmap_bim):
            raise FileNotFoundError("hapmap_founders.bim file not found")
        if not os.path.exists(hapmap_fam):
            raise FileNotFoundError("hapmap_founders.fam file not found")

        step = "identify_samples_divergent_ancestry"

        #
        plink_cmd1 = f"plink --bfile {os.path.join(result_path, output_name+'.pre_ind_clean')} --autosome --keep-allele-order --maf {maf} --extract {hapmap3r2_txt} --make-bed --out {os.path.join(result_path, output_name+'.hapmap-snps')}"

        #
        plink_cmd2 = f"plink --bfile {os.path.join(result_path, output_name+'.hapmap-snps')} --autosome --keep-allele-order --maf {maf} --bmerge {hapmap_bed} {hapmap_bim} {hapmap_fam} --extract {os.path.join(result_path, output_name+'_1.prune.in')} --make-bed --out {os.path.join(result_path, output_name+'.hapmap3r2.pruned')}"

        # execute PLink commands
        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd, log=True)

        # path to missnp file
        missnp_file = os.path.join(result_path, output_name+'.hapmap3r2.pruned-merge.missnp')

        if os.path.exists(missnp_file):

            #
            plink_cmd3 = f"plink --bfile {os.path.join(result_path, output_name+'.hapmap-snps')} --autosome --keep-allele-order --maf {maf} --flip {missnp_file} --make-bed --out {os.path.join(result_path, output_name+'.hapmap-snps1')}"

            # 
            plink_cmd4 = f"plink --bfile {os.path.join(result_path, output_name+'.hapmap-snps1')} --autosome --keep-allele-order --maf {maf} --bmerge {hapmap_bed} {hapmap_bim} {hapmap_fam} --extract {os.path.join(result_path, output_name+'_1.prune.in')} --make-bed --out {os.path.join(result_path, output_name+'.hapmap3r2.pruned')}"

            # execute PLink commands
            cmds = [plink_cmd3, plink_cmd4]
            for cmd in cmds:
                shell_do(cmd, log=True)

            # rename .bim file
            bash_cmd1 = f"cp {os.path.join(result_path, output_name+'.hapmap3r2.pruned.bim')} {os.path.join(result_path, output_name+'.hapmap3r2.pruned.pedsnp')}"

            # rename .fam file
            bash_cmd2 = f"cp {os.path.join(result_path, output_name+'.hapmap3r2.pruned.fam')} {os.path.join(result_path, output_name+'.hapmap3r2.pruned.pedind')}"

            # execute bash commands
            bashes = [bash_cmd1, bash_cmd2]
            for cmd in bashes:
                os.system(cmd)

        else:

            # rename .bim file
            bash_cmd1 = f"cp {os.path.join(result_path, output_name+'.hapmap3r2.pruned.bim')} {os.path.join(result_path, output_name+'.hapmap3r2.pruned.pedsnp')}"

            # rename .fam file
            bash_cmd2 = f"cp {os.path.join(result_path, output_name+'.hapmap3r2.pruned.fam')} {os.path.join(result_path, output_name+'.hapmap3r2.pruned.pedind')}"

            # execute bash commands
            bashes = [bash_cmd1, bash_cmd2]
            for cmd in bashes:
                os.system(cmd)

        # report
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
    

    def run_pca_analysis(self)->dict:

        return None
