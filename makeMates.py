from Bio import SeqIO
from numpy import random

seqFileName = 'aDir/andSubdir/toPath.fa'

def makeMates(seqFileName,readLength=101,matePairLength=2000,numberOfPairs=200000,pairLengthSTD=100,output1='output1.fastq',output2='output2.fastq'):
    refSeq = SeqIO.read(seqFileName,'fasta')
    seq1s = []
    seq2s = []
    for i in range(numberOfPairs):
        if i%1000==0:
            print i
        tempStd = int(random.normal(0,pairLengthSTD))
        tempStart = random.randint(0,len(refSeq)-matePairLength-tempStd)
        tempEnd = tempStart + matePairLength + tempStd
        seq1s.append(refSeq[tempStart:tempStart+readLength])
        # If we have no phred quality scores for these sequences, let's say it's 40
        if not seq1s[-1].letter_annotations.has_key('phred_quality'):
            seq1s[-1].letter_annotations['phred_quality'] = [40 for _ in range(len(seq1s[-1]))]
        # Biopython uses description to assign sequence names?
        seq1s[-1].description = '{0} {1}-{2} ({3}) length={4}'.format(refSeq.description, tempStart, tempEnd, 1, readLength)
        seq2s.append(refSeq[tempEnd-readLength:tempEnd].reverse_complement())
        if not seq2s[-1].letter_annotations.has_key('phred_quality'):
            seq2s[-1].letter_annotations['phred_quality'] = [40 for _ in range(len(seq2s[-1]))]
        seq2s[-1].description = '{0} {1}-{2} ({3}) length={4}'.format(refSeq.description, tempStart, tempEnd, 2, readLength)
        seq2s[-1].id = seq1s[-1].id
    SeqIO.write(seq1s, output1, 'fastq')
    SeqIO.write(seq2s, output2, 'fastq')


