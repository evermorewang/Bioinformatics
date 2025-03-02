#this is to calculate a theoretical minimum Us in an mRNA encoding the designed pepitide sequence
#using human A-syn as an example: uniprot P37840 https://www.uniprot.org/uniprotkb/P37840/feature-viewer
#Epitope 114-125, “EDMPVDPDNEAY” https://www.iedb.org/epitope/1094247
#a-syn seq:MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA
#epitope seq:EDMPVDPDNEAY
# Codon table with min U codons for each amino acid
codon_table = {
    'A': 'GCC',  # Ala: 0 U (GCU, GCC, GCA, GCG)
    'R': 'CGC',  # Arg: 0 U (CGU, CGC, CGA, CGG, AGA, AGG)
    'N': 'AAC',  # Asn: 0 U (AAU, AAC)
    'D': 'GAC',  # Asp: 0 U (GAU, GAC)
    'C': 'UGC',  # Cys: 1 U (UGU, UGC)
    'Q': 'CAG',  # Gln: 0 U (CAA, CAG)
    'E': 'GAG',  # Glu: 0 U (GAA, GAG)
    'G': 'GGC',  # Gly: 0 U (GGU, GGC, GGA, GGG)
    'H': 'CAC',  # His: 0 U (CAU, CAC)
    'I': 'AUC',  # Ile: 1 U (AUU, AUC, AUA)
    'L': 'CUC',  # Leu: 0 U (UUA, UUG, CUU, CUC, CUA, CUG)
    'K': 'AAG',  # Lys: 0 U (AAA, AAG)
    'M': 'AUG',  # Met: 1 U (AUG, start codon only)
    'F': 'UUC',  # Phe: 2 U (UUU, UUC)
    'P': 'CCC',  # Pro: 0 U (CCU, CCC, CCA, CCG)
    'S': 'AGC',  # Ser: 0 U (UCU, UCC, UCA, UCG, AGU, AGC)
    'T': 'ACC',  # Thr: 0 U (ACU, ACC, ACA, ACG)
    'W': 'UGG',  # Trp: 1 U (UGG only)
    'Y': 'UAC',  # Tyr: 1 U (UAU, UAC)
    'V': 'GUC',  # Val: 0 U (GUU, GUC, GUA, GUG)
    '*': 'UAG'   # Stop: 1 U (UAA, UAG, UGA)
}

# Function to count U's in a codon
def count_u(codon):
    return codon.count('U')

# Function to translate ORF to amino acids
def translate_orf(orf):
    aa_seq = []
    # for i in range(0, len(orf), 3):
    #     codon = orf[i:i+3]
    #     for aa, min_codon in codon_table.items():
    #         if codon == 'AUG' and i == 0:  # Start codon
    #             aa_seq.append('M')
    #             break
    #         elif codon in ['UAA', 'UAG', 'UGA']:  # Stop codons
    #             aa_seq.append('*')
    #             break
    #         elif codon == min_codon:
    #             aa_seq.append(aa)
    #             break

    for l_i in orf:
        for k,v in codon_table.items():
            if l_i==k:
                print(v)
                aa_seq.append(v)
                break
    print(aa_seq)
    return ''.join(aa_seq)

# Function to calculate min U
def calc_min_u(orf):
    aa_seq = translate_orf(orf)
    min_u_seq = ''
    for aa in aa_seq:
        min_u_seq += codon_table[aa]
    min_u_count = count_u(min_u_seq)
    return min_u_seq, min_u_count

# Your ORF
orf = "MEDMPVDPDNEAYF*"  # Spaces for readability, remove in practice
orf = orf.replace(" ", "")  # Clean input

aa_seq=translate_orf(orf)
print("Original U content %d minimum."% count_u(aa_seq))