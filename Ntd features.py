from Bio import SeqIO  
import numpy as np
from Bio.Seq import Seq
from Bio.Data import IUPACData
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC_skew
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import GC123
from scipy.stats import entropy
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
MyIndex = {
"GCA": 0.586, "GCC": 0.122, "GCG": 0.424, "GCT": 1.000,
"AGA": 0.004, "AGG": 0.002, "CGA": 0.004, "CGC": 0.356,
"CGG": 0.004, "CGT": 1.000, "AAC": 1.000, "AAT": 0.051,
"GAC": 1.000, "GAT": 0.434, "TGC": 1.000, "TGT": 0.500,
"CAA": 0.124, "CAG": 1.000, "GAA": 1.000, "GAG": 0.259,
"GGA": 0.010, "GGC": 0.724, "GGG": 0.019, "GGT": 1.000,
"CAC": 1.000, "CAT": 0.291, "ATA": 0.003, "ATC": 1.000,
"ATT": 0.185, "CTA": 0.007, "CTC": 0.037, "CTG": 1.000,
"CTT": 0.042, "TTA": 0.020, "TTG": 0.020, "AAA": 1.000,
"AAG": 0.253, "ATG": 1.000, "TTC": 1.000, "TTT": 0.296,
"CCA": 0.135, "CCC": 0.012, "CCG": 1.000, "CCT": 0.070,
"AGC": 0.410, "AGT": 0.085, "TCA": 0.077, "TCC": 0.744,
"TCG": 0.017, "TCT": 1.000, "ACA": 0.076, "ACC": 1.000,
"ACG": 0.099, "ACT": 0.965, "TGG": 1.000, "TAC": 1.000,
"TAT": 0.239, "GTA": 0.495, "GTC": 0.066, "GTG": 0.221,
"GTT": 1.000
}
import zlib
import sys

with open('ntd_table.txt', 'a') as f:
  header= 'Name\tLength\tGC Content\tGC Skew_Stdev\tGC Skew_median\tGC Skew_mean\tGC Position 1\tGC Position 2\tGC Postion 3\tMelting temperature\tMolecular weight\tEntropy\tZlib compression ratio\tCAI\tAssigned'.split('\t')
  f.write('\t'.join(header)+'\n')
  for seq in SeqIO.parse("pea.assigned.cds.fa","fasta"):
    geneID=str(seq.id)
    seq=str(seq.seq)
    gc=GC(seq)
    gcs = GC_skew(seq)# a list, a gc skew value for each window
    gcs_std = np.std(gcs)
    gcs_median = np.median(gcs)
    gcs_mean = np.mean(gcs)
    gcp = GC123(seq)
    mt = MeltingTemp.Tm_NN(seq)
    mw = molecular_weight(seq)
    labels = seq.replace('A','1').replace('T','2').replace('G','3').replace('C', '4')
    labels = list(map(int, labels))
    en = entropy(labels)
    s = str.encode(seq)
    after = sys.getsizeof(zlib.compress(s))
    before = sys.getsizeof(s)
    CAI = CodonAdaptationIndex()
    CAI.set_cai_index(MyIndex)
    x = len(seq)/3
    if (x - int(x)) == 0:
      CAI_value=CAI.cai_for_gene(seq)
    else:
      CAI_value = "0"

    print('\t'.join(map(str, [geneID, len(seq), gc] + [gcs_std, gcs_median, gcs_mean, gcp[1], gcp[2],gcp[3], mt, mw, en, after/before, CAI_value, 1])), file=f)
 
 

