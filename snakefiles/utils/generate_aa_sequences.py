import sys
sys.path.append('/home/oystein/HATK')
from IMGT2Seq.__main__ import HATK_IMGT2Seq
from NomenCleaner.__main__ import HATK_NomenCleaner
from bMarkerGenerator.__main__ import HATK_bMarkerGenertor

hg = 19
imgt = 3631
imgt_dir = "/home/oystein/IMGTHLA/"
output_dir = "/home/oystein/test/IMGT2Seq_test"
hped = "/home/oystein/hla_imputation_pipeout/2025.09.16/CookHLA/cookhla_output.MHC.HLA_IMPUTATION_OUT.hped"
hat = f"{output_dir}/HLA_ALLELE_TABLE.imgt3631.hat"

HATK_IMGT2Seq(imgt, hg, output_dir, imgt_dir, False, False, False, False, False, False)
HATK_NomenCleaner(hped, hat, output_dir, False, False, False, False, False, False)

hibag_chped = "/home/oystein/test/IMGT2Seq_test.imgt3631.2field.chped"

HATK_bMarkerGenertor(hibag_chped, f"{output_dir}/HIBAG", 19, "/home/oystein/test/IMGT2Seq_test/HLA_DICTIONARY_AA.hg19.imgt3631", "/home/oystein/test/IMGT2Seq_test/HLA_DICTIONARY_SNPS.hg19.imgt3631", "/home/oystein/hla_imputation_pipeout/2025.09.16/extract_hla")