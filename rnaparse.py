"""
Simple Script to parse 723 RNA sequences from one fasta file
File provided by RNA Strand Database, BETA lab, University of British Columbia

Dataset found in: "rna.fasta"
Data selected by selecting 16s Riobsomal RNA, and returning 723 sequences

For analysis of secondary structures, analysis can be performed on this dataset
here:
http://www.rnasoft.ca/strand/search_results.php?select%5B%5D=16S+Ribosomal+RNA&org=&source%5B%5D=Any+source&source_id=&select2=Any+length&first1=&last1=&exp_proven=Any&select4=Any+number&first2=&last2=&select5=Any&select_duplicate=All+molecules&seq=&abstractshapetype=ABSTRACT_SHAPE_5&abstractshape=&Submit=Perform+search&select6=Any+number&first3=&last3=&select7=Any+number&first4=&last4=&select8=Any+number&first5=&last5=&select12=Any+number&first7=&last7=&select13=Any+number&first8=&last8=&motif=&select19=Any+number&first10=&last10=&select20=Any+number&first11=&last11=&select25=Any+number&first13=&last13=&select26=Any+number&first14=&last14=&select27=absolute&select28=Any+number&first15=&last15=&select33=Any+number&first17=&last17=&select34=Any+number&first18=&last18=&select35=average&select36=absolute&select37=Any+number&first19=&last19=&select38=Any+number&first20=&last20=&select69=Any+number&first31=&last31=&select71=bands&select72=Any+number&first32=&last32=&select73=base+pairs&select74=Any+number&first33=&last33=&select75=Any&select76=Any+number&first34=&last34=&select55=any&select_ncbp_context_2=any&select_ncbp_context_4=any&select_ncbp_context_1=any&select_ncbp_context_5=any&select56=Any+number&first26=&last26=&start=0&limit=10
"""

from Bio import SeqIO
from r_fold import RNA_Fold
import local_align
fasta_seq = SeqIO.parse(open("rna.fasta"), "fasta")
fold_scores = []

def main_script():
  for fasta in fasta_seq:
    name, sequence = fasta.id, str(fasta.seq)
    r = RNA_Fold(sequence)
    s = r.fold(sequence)
    fold_scores.append((name, s))

main_script()
