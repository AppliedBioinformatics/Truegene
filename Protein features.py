from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio.SeqUtils import molecular_weight
from Bio.Data import IUPACData
letters = IUPACData.protein_letters
a = 2.9
b = 3.9
with open('protein_table.txt', 'a') as f:
  header = 'Name\tLength (aa)\tHas start codon\tAmino acid percentage\tIsoelectric point\tGravy value\tMolecular weight\tInstability index\tSecondary structure fraction helix\tSecondary structure fraction turn\tSecondary structure fraction sheet\tflexibility\tAromaticity\tMolar extinction co-efficient reduced\tMolar extinction co-efficient oxidized\tAliphatic_index\tcharge\tpolar\tacidic\tbasic\ttiny\tnpolar\taliphaticity\tAssigned'.split('\t')
  before = header[:3]
  after = header[4:]
  new = []
  for l in letters:
    new.append('Amino acid percentage %s'%l)
  
  f.write('\t'.join(before + new + after) + '\n')

  for seq in SeqIO.parse('pea.assigned.pep.fa', 'fasta'):
    prot = PA(str(seq.seq))
    seq = seq.upper()
    gen=str(seq.seq)
    arg = gen.count('R')
    his = gen.count('H')
    lys = gen.count('K')
    asp = gen.count('D')
    glu = gen.count('E')
    ser = gen.count('S')
    thr = gen.count('T')
    asn = gen.count('N')
    glm = gen.count('Q')
    gly = gen.count('G')
    ala = gen.count('A')
    cys = gen.count('C')
    val = gen.count('V')
    leu = gen.count('L')
    ile = gen.count('I')
    met = gen.count('M')
    trp = gen.count('W')
    phe = gen.count('F')
    pro = gen.count('P')
    charge = ((arg+lys+asp+glu)/len(seq))
    polar = ((ser+thr+asn+glm)/len(seq))
    acidic = ((asp+glu)/len(seq))
    basic = ((arg+his+lys)/len(seq))
    tiny = ((gly+ala+cys+ser)/len(seq))
    npolar = ((ala+val+leu+gly+ile+met+trp+phe+pro)/len(seq))
    aliphaticity = ((ala+val+leu+ile)/len(seq))
    length = float(len(seq))
    alanine_per = (gen.count('A') / length )
    valine_per = (gen.count('V') / length )
    isoleucine_per = (gen.count('I') / length )
    leucine_per = (gen.count('L') / length )
    # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    aliphatic_index = (100 * (alanine_per + a * valine_per + b * (isoleucine_per + leucine_per )))
    
    percentages = prot.get_amino_acids_percent()
    aminos = []

    for l in letters:
        aminos.append(percentages[l])
    has_start = 0
    if seq.seq.startswith('M'):
        has_start = 1
    print('\t'.join(map(str, [seq.id, len(seq), has_start] + aminos +  [prot.isoelectric_point(), prot.gravy(), prot.molecular_weight(), prot.instability_index()] + list(prot.secondary_structure_fraction()) + [sum(prot.flexibility()), prot.aromaticity()] + list(prot.molar_extinction_coefficient()) + [aliphatic_index, charge, polar, acidic, basic,tiny,npolar,aliphaticity, 1])), file=f)




