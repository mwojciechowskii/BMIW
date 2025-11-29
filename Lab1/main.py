from typing import List
import gzip
from sequence import Sequence
import argparse
import re

def readFile(file: str) -> List["Sequence"]:
    
    seq = []
    if file.endswith(".gz"):
        openF = gzip.open
    else:
        openF = open
    with openF(file, 'rt') as opFile:
        lines = []
        for line in opFile:
            lines.append(line.rstrip())
            if len(lines) == 4:
                sequence = Sequence(lines[0], lines[1], lines[2], lines[3])
                seq.append(sequence)
                lines = []
    return seq

def removeExt(filename: str) -> str:
    return re.sub(r'\.(fastq|fq)\.gz$|\.fastq$|\.fq$', '', filename)

def removePairSuffix(filename: str) -> str:
    return re.sub(r'[_.]?(R?[12])$', '', filename)

def argParser():

    parser = argparse.ArgumentParser(description="Process multiple files")
    parser.add_argument("files", nargs='*', help="List of files to process")
    parser.add_argument("--paired", action="store_true", help="Treat files as paired-end and aggregate results")
    args = parser.parse_args()
    return args

def Analyze(fileSequences: dict[str, List]):

    for filename, seq in fileSequences.items(): 
        print(f"FILE: {filename}")
        print(f"Total num of read: {Sequence.CntReads(seq)}")
        print(f"Sum of nucleotides: {Sequence.SumFileNT(seq)}")
        readLen = Sequence.ReadLen(seq)
        print(f"Mean read length: {Sequence.MeanReadLen(seq):.2f}")
        print(f"Mean GC count: {Sequence.calcMeanGC(seq):.2f}%")
        Sequence.WriteMeanReadHist(readLen, filename)

def main():

    args = argParser()
    if not args.files:
        args.files = ["IL_1.fastq.gz"]
    FilesSequences = {}
    for file in args.files:
        seq = readFile(file)
        filename = removeExt(file)
        if args.paired:
            filename = removePairSuffix(filename)
        if filename in FilesSequences:
            FilesSequences[filename].extend(seq)
        else:
            FilesSequences[filename] = seq
    Analyze(FilesSequences)
        
if __name__ == "__main__":
    main()
