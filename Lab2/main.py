from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import argparse
from typing import Dict, Iterator, Tuple, List
import statistics
import itertools

""" Postawilem sobie wyzwanie, zeby uzywac tu jak najwiecej iteratorow"""

unknownFunc = ('hypothetical protein', 'putative protein', 'uncharacterized protein', 'predicted protein')

def argParser():
    parser = argparse.ArgumentParser(description="Process multiple files")
    parser.add_argument("files", nargs='*', help="List of files to process")
    args = parser.parse_args()
    return args

def readFile(file: str) -> Iterator[SeqRecord]:
    return SeqIO.parse(file, "genbank")

def extractCDS(fileRecrods: Tuple[str, Iterator[SeqRecord]]) -> Tuple[str, List[SeqFeature]]:

    file, records = fileRecrods
    cds: List[SeqFeature] = []
    for rec in records:
        cds.extend([f for f in rec.features if f.type == "CDS"])
    return file, cds

def extractProteins(cdsRecords: Iterator[Tuple[str, List[SeqFeature]]]) -> Iterator[Tuple[str, List[str]]]:

    for filename, cdsList in cdsRecords:
        proteins: List[str] = []
        for cds in cdsList:
            products = cds.qualifiers.get('product', [])
            proteins.extend(products)
        yield filename, proteins

def isUnknown(protein: str) -> bool:
    return protein.lower() in unknownFunc

def cntProteinTypes(proteins: Iterator[Tuple[str, List[str]]]) -> Iterator[Tuple[str, int, int]]:

    for filename, proteinList in proteins:
        unknownCnt = sum(1 for p in proteinList if isUnknown(p))
        knownCnt = len(proteinList) - unknownCnt
        yield filename, knownCnt, unknownCnt

def geneMedian(cdsRecords: Iterator[Tuple[str, List[SeqFeature]]]) -> Dict[str, int]:

    medians = {}
    for filename, cdsList in cdsRecords:
        medians[filename]= statistics.median([len(gene) for gene in cdsList])
    return medians

def main():
    args = argParser()
    if not args.files:
        args.files = ["assemblyLight.gbff", "assemblyFull.gbff"]

    genFiles = {file: readFile(file) for file in args.files}
    cds = (extractCDS(fr) for fr in genFiles.items())
    #copying iterator 
    cds1, cds2= itertools.tee(cds, 2)
    medians = geneMedian(cds1)
    for key, value in medians.items():
        print(f"Gene length in median in file {key} is: {value}")

    proteins = extractProteins(cds2)
    results = list(cntProteinTypes(proteins))
    for file, knownNo, unknownNo in results:
        print(f"Amount of genes with known functionality for {file} is: {knownNo}, unknown functionality: {unknownNo}")

if __name__ == "__main__":
    main()

