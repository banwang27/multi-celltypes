#!/usr/bin/env python

# Prints all non-indel positions in MAF file with coordinates and variants from each species / individual 

# Usage python callSNPsFromMAF.py mafFile outFile

import sys

mafFile, outFile = sys.argv[1:]

outFile = open(outFile, 'w')

def parseChrom(chrom):
    fields = chrom.split(".")
    newChrom = ".".join(fields[1:])
    return(newChrom)


with open(mafFile, 'r') as maf:
    getPair = False

    for line in maf:

        if getPair == False:

            if line.startswith("#") or line.startswith(" ") or line.startswith("a"):
                continue

            elif line.startswith("s"):
                lineID, refChrom, refStart, alnLen, refStrand, refLen, refSeq = line.strip().split()
                getPair = True

        else:
            lineID, altChrom, altStart, alnLen, altStrand, altLen, altSeq = line.strip().split()

            refChrom = parseChrom(refChrom)
            altChrom = parseChrom(altChrom) 

            refGaps = 0
            altGaps = 0

            for i in range(0, len(altSeq)):

                if refSeq[i] == "-":
                    refGaps += 1
                elif altSeq[i] == "-":
                    altGaps += 1

                if refSeq[i] == "-" or altSeq[i] == "-": 
                    ref = refSeq[i].upper()
                    alt = altSeq[i].upper()
                    refPos = int(refStart) - refGaps + i + 1

                    if altStrand == "+":
                        altPos = int(altStart) - altGaps + i + 1
                    else:
                        altPos = int(altLen) - (int(altStart) - altGaps + i)

                    outList = [refChrom, str(refPos), refStrand, ref, alt, altStrand, str(altPos), altChrom]
                    outPut = "\t".join(outList) + '\n'
                    outFile.write(outPut)

            getPair = False

outFile.close()
