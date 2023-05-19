"""
splitSpeciesReads.py

Designed to separate reads from the two species in hybrid data.
Adapted from GetGeneASEbyReads.py by Peter Combs.

NOTE: Input bam files MUST be sorted by mate-pair.

NOTE: Human is always the reference file regardless of which SNP file is used.

Last edit: 07.27.2016

"""

from __future__ import print_function
from pysam import Samfile
from argparse import ArgumentParser, FileType
from collections import defaultdict, Counter
from os import path
import time

def get_phase(read, snps):
    phase = None
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos + 1 in snps and read.query_qualities[read_pos] >= 30:
            if phase == None:
                try:
                    # 1 if alternate, -1 if reference
                    phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0 # This SNP isn't in the dataset
            else:
                try:
                    new_phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0
                if new_phase != phase:
                    return 0 # read seems misphased
    return phase

def get_snps(snpfile):
    snps = defaultdict(dict)
    if path.exists(path.join(path.dirname(snpfile), 'true_hets.tsv')):
        print("using true hets")
        true_hets = {tuple(line.strip().split()):True
                     for line in open(path.join(path.dirname(snpfile), 'true_hets.tsv'))
                    }
    else:
        true_hets = defaultdict(lambda x: True)
    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')

    return snps

def split_reads(reads, ref_reads, alt_reads):
    read_results = Counter()

    ref_bam = Samfile(ref_reads, 'wb', template=reads)
    alt_bam = Samfile(alt_reads, 'wb', template=reads)

    prev_phase = None
    prev_read = None
    prev_qname = None
    test = 0

    for read in reads.fetch(until_eof=True):
        test += 1
        chrom = read.reference_name
        snps_on_chrom = snp_dict[chrom]
        phase = get_phase(read, snps_on_chrom)
        read_qname = read.qname
        if read_qname == prev_qname:
            read_results["tot_read"]+=1
            phase_set = set([phase, prev_phase])
            phase_set.discard(None)
            if len(phase_set) == 1:
                read_phase = phase_set.pop()
                if read_phase == -1:
                    ref_bam.write(read)
                    ref_bam.write(prev_read)
                    read_results["ref_read"]+=1
                elif read_phase == 1:
                    alt_bam.write(read)
                    alt_bam.write(prev_read)
                    read_results["alt_read"]+=1
                elif read_phase == 0:
                    read_results['misphased_read']+=1
            elif len(phase_set)==0:
                read_results["no_snps_read"]+=1
            else:
                read_results['misphased_read']+=1

        prev_read = read
        prev_phase = phase
        prev_qname = read_qname

    return(read_results)
    


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('snp_file', help="Tab delimited text file of known SNPs with 4 columns: chrom, 0-position, 1-position, ref_allele|alt_allele")
    parser.add_argument('reads', type=Samfile, help="BAM file that needs splitting")
    parser.add_argument('ref_reads', help="Output BAM file of reference reads.")
    parser.add_argument('alt_reads', help="Output BAM file of alternate reads.")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    print(time.strftime(("%b %d ") + time.strftime("%I:%M:%S")),"Loading SNPs...")
    args = parse_args()
    snp_dict = get_snps(args.snp_file)

    print(time.strftime(("%b %d ") + time.strftime("%I:%M:%S")),"Done with SNPs. Splitting reads...")
    read_results = split_reads(args.reads, args.ref_reads, args.alt_reads)

    print(time.strftime(("%b %d ") + time.strftime("%I:%M:%S")),"Finished!")

    print()
    print("RUN STATISTICS:")
    print("  Total input reads:", str(read_results['tot_read']))
    print("  Reads with no SNPs:", str(read_results['no_snps_read']), 
        "(" + "%.2f" % ((read_results['no_snps_read']/float(read_results['tot_read']))*100) + "%)")
    print("  Reads with misphased SNPs:", str(read_results['misphased_read']), 
        "(" + "%.2f" % ((read_results['misphased_read']/float(read_results['tot_read']))*100) + "%)")
    print("  Reference reads:", str(read_results['ref_read']), 
        "(" + "%.2f" % ((read_results['ref_read']/float(read_results['tot_read']))*100) + "%)")
    print("  Alternate reads:", str(read_results['alt_read']), 
        "(" + "%.2f" % ((read_results['alt_read']/float(read_results['tot_read']))*100) + "%)")

