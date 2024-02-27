from classes.SampleQC import SampleQC

INPUT_PATH = '/mnt/0A2AAC152AABFBB7/PipeLine/data/inputData'
INPUT_NAME = 'test_1'
OUTPUT_PATH= '/mnt/0A2AAC152AABFBB7/PipeLine/data/outputData'
OUTPUT_NAME= 'results'
CONFIG_PATH= '/mnt/0A2AAC152AABFBB7/comrare-pipeline/config.JSON'

REGION_PATH = '/mnt/0A2AAC152AABFBB7/PipeLine/data/high-LD-regions_GRCh37.txt'

sample_QC = SampleQC(
    input_path=INPUT_PATH,
    input_name=INPUT_NAME,
    output_path=OUTPUT_PATH,
    output_name=OUTPUT_NAME,
    config_path=CONFIG_PATH
)

sample_QC.run_ld_prune(ld_region_path=REGION_PATH)

sample_QC.run_heterozygosity_rate()
