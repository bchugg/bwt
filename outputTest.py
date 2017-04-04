from local_align import localAlignRNA as localAlignRNA
from kLocalRNA import kLocalFold
from Bio import SeqIO

def main():
    seq = open("testseq.txt", 'r')
    seq = seq.read().strip()
    print(kLocalFold(seq, 1))
    test = localAlignRNA(seq)
    test.computeAlignment()
    printAlignment(test.align_info)


def printAlignment(alignment):
		"""
		Print single alignment
		"""
		print("SCORE: ", alignment[2])
		print('Match', alignment[0][0], '-', alignment[0][1])
		print('With', alignment[1][0], '-', alignment[1][1])
		print("Alignment: ")
		print(alignment[3][0][0:60])
		print(alignment[3][1][0:60])
		print('\n')


main()
