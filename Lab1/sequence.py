from dataclasses import dataclass
from typing import Dict, Iterable
from statistics import mean
import matplotlib.pyplot as plt

@dataclass
class Sequence:
    seq_ID: str
    seq: str
    opt: str
    quality: str
    GC_content: float = 0.0

    def __post_init__(self):
        self.NtCounts: Dict[str, int] = {'A': 0, 'T': 0, 'U': 0, 'C': 0, 'G': 0}
        self.CntNT()
        self.calcGC()

    def __repr__(self):
        return f"Sequence(seq_ID={self.seq_ID!r}, seq={self.seq!r}, opt={self.opt!r}, quality={self.quality!r})"

    def CntNT(self): 

        for i in self.seq:
            if i in self.NtCounts:
                self.NtCounts[i] += 1
        return self

    def calcGC(self) -> "Sequence":
        
        total: int = sum(self.NtCounts.values())
        if total == 0:
            self.GC_content = 0.0
            return self
        GC: int = self.NtCounts['C'] + self.NtCounts['G']
        self.GC_content = (GC/total) * 100 
        return self

    @staticmethod
    def calcMeanGC(sequences: Iterable["Sequence"]) -> float:

        gcList = []
        for i in sequences:
            if i.GC_content is None:
                i.calcGC()
            gcList.append(i.GC_content)

        return mean(gcList) if gcList else 0.0

    @staticmethod
    def SumFileNT(sequences: Iterable["Sequence"]):
        
        mySum = 0
        for i in sequences:
            mySum += sum(i.NtCounts.values())

        return mySum

    @staticmethod
    def CntReads(sequences: Iterable["Sequence"]):
        return len(list(sequences))

    @staticmethod
    def ReadLen(sequences: Iterable["Sequence"]):
        return [len(i.seq) for i in sequences]

    @staticmethod
    def MeanReadLen(sequences: Iterable["Sequence"]):
        readLen = Sequence.ReadLen(sequences)
        return mean(readLen)
    
    @staticmethod
    def WriteMeanReadHist(readLen: list[int], outfile: str):

        plt.figure(figsize=(10, 6))
        plt.hist(readLen, bins=20, color='blue', edgecolor='black', alpha=0.7)

        plt.title('Read Length Distribution')
        plt.xlabel('Read Length')
        plt.ylabel('Frequency')
        plt.savefig(outfile+'.png', format='png')

        plt.close()
        print(f"Plot saved for: {outfile}")
