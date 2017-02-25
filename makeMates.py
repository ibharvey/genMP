from Bio import SeqIO
from numpy import random

seqFileName = 'aDir/andSubdir/toPath.fa'

def makeMates(seqFileName,readLength=101,matePairLength=2000,numberOfPairs=200000,pairLengthSTD=100):
    refSeq = SeqIO.read(seqFileName,'fasta')
    seq1s = []
    seq2s = []
    for _ in range(numberOfPairs):
        tempStart = len(refSeq)
        tempStd = 1
        # Generate a random DNA molecule from the template
        while tempStart + matePairLength + tempStd > len(refSeq):
            tempStart = random.randint(0,len(refSeq))
            tempStd = random.normal(0,pairLengthSTD)
        seq1s.append(refSeq[tempStart:101])
        temp2End = tempStart + matePairLength + tempStd
        seq2s.append(refSeq[tempEnd-101:tempEnd].reverse_complement())
