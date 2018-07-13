
import glob
import os
import pdb
import pysam

class Coverage:
    depths = [1, 5, 10, 15, 20, 25, 30, 50, 100]
    def __init__(self, line):
        line_to_unpack = line.strip('\n').split('\t')
        chrom, pos, mean, median = line_to_unpack[0:4]
        cov = line_to_unpack[4:]
        self.chrom = chrom
        self.pos = int(pos)
        self.mean = float(mean)
        self.median = float(median)

        self.cov = dict(zip(self.depths, map(float, cov)))

    def __getitem__(self, key):
        return self.cov[key]

def open_coverage(folder):
    """ open folder of coverage files, one per chromosome
    """
    path = os.path.join(folder, '*.chr*.coverage.txt.gz')
    coverage = {}
    for path in glob.glob(path):
        basename = os.path.basename(path)
        chrom = [ x for x in basename.split('.') if x.startswith('chr') ][0]
        chrom = chrom.strip('chr')
        coverage[chrom] = pysam.TabixFile(path)

    return coverage

def good_coverage(vcf, chrom, start, end, threshold=0.75):
    """ checks if a site is well covered in the reference cohort
    """
    successes = 0
    failures = 0
    for line in vcf.fetch(chrom, start-1, end):
        cov = Coverage(line)
        if(cov[30] > threshold):
            successes += 1
        else:
            failures += 1
    return successes, failures
