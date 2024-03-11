from classes.SampleQC import SampleQC

INPUT_PATH = '/mnt/0A2AAC152AABFBB7/PipeLine/data/inputData'
INPUT_NAME = 'test_1'
OUTPUT_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/outputData'
OUTPUT_NAME= 'results'
CONFIG_PATH= '/mnt/0A2AAC152AABFBB7/comrare-pipeline/config.JSON'
DEPEND_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/auxiliarData'

sample_QC = SampleQC(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_path=CONFIG_PATH,
    dependables_path=DEPEND_PATH
)

sample_QC.run_ld_prune(ld_region_file='high-LD-regions.txt')

sample_QC.run_heterozygosity_rate()

sample_QC.run_sex_check()

sample_QC.run_relatedness_prune()

sample_QC.delete_failing_QC()

sample_QC.divergent_ancestry_step_one(
    hapmap3r2_no_at_cg_snps_file='hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt', 
    hapmap3r2_founders_name='hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps'
)
