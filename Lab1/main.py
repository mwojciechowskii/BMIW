from typing import List
from sequence import Sequence
import argparse

def readFile(file: str) -> List["Sequence"]:
    
    sequences = []
    with open(file, 'r') as opFile:
        lines = []
        for line in opFile:
            lines.append(line.rstrip())
            if len(lines) == 4:
                seq = Sequence(lines[0], lines[1], lines[2], lines[3])
                sequences.append(seq)
                lines = []
    return sequences

def argParser():

    parser = argparse.ArgumentParser(description="Process multiple files")
    parser.add_argument("files", nargs='*', help="List of files to process")
    args = parser.parse_args()
    return args

def Analyze(fileSequences: Dict[str, List]):

    


    
    pass

def main():

    args = argParser()
    if not args.files:
        args.files = ["file1.fq", "file2.fq"]
    FilesSequences = {}
    for file in args.files:
        sequences = readFile(file)
        FilesSequences[file] = sequences

        
        
   



    sequences = readFile("test.fq")
    print(Sequence.calcMeanGC(sequences))
    print(Sequence.SumFileNT(sequences))
    print(Sequence.CntReads(sequences))
    print(Sequence.MeanReadLen(sequences))

if __name__ == "__main__":
    main()
