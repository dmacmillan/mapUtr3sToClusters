import argparse
import os
import sys
import pysam
import gzip
import copy
import math
from Kleat import *
from scipy import stats

def parseT2G(fpath):
    results = {}
    with gzip.open(fpath, 'rb') as f:
        for line in f:
            line = line.strip().split('\t')
            results[line[0]] = line[1]
    return results

def getUtr3Length(utr3, site):
    if utr3.strand == '+':
        return site - utr3.start
    else:
        return utr3.end - site

def filterUtr3s(utr3s, t2g):
    genes = {}
    for utr3 in utr3s:
        gene = t2g[utr3.name]
        if gene not in genes:
            genes[gene] = [utr3]
        else:
            genes[gene].append(utr3)
    keep = {}
    for gene in genes:
        genes[gene] = sorted(genes[gene], key = lambda x: (x.start, x.end))
        if genes[gene][0].strand == '+':
            if len(set([x.start for x in genes[gene]])) == 1:
                keep[gene] = genes[gene]
        else:
            if len(set([x.end for x in genes[gene]])) == 1:
                keep[gene] = genes[gene]
        store = True
        for i in xrange(1, len(genes[gene])):
            current = genes[gene][i]
            previous = genes[gene][i-1]
            if (current.start == previous.start):
                continue
            if (current.start < previous.end):
                store = False
                break
        if store:
            keep[gene] = genes[gene]
    return keep

def groupKleatCleavage(kleat, result={}):
    for entry in kleat:
        if entry.chromosome not in result:
            result[entry.chromosome] = {entry.gene: [entry.cleavage_site]}
        elif entry.gene not in result[entry.chromosome]:
            result[entry.chromosome][entry.gene] = [entry.cleavage_site]
        else:
            result[entry.chromosome][entry.gene].append(entry.cleavage_site)
    return result
            
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Given a KLEAT file return 3\'UTR coordinates for each cleavage site, otherwise return None')
    
    parser.add_argument('-k1', '--kleats_1', nargs='+', help='Set of KLEAT files for the first dataset')
    parser.add_argument('-k2', '--kleats_2', nargs='+', help='Set of KLEAT files for the second dataset')
    parser.add_argument('-t2g', '--transcripts_to_genes', default='/home/dmacmillan/annotations/ensembl/originals/ucsc.ensemblToGeneName.original.gz', help='A mapping of transcripts to genes in TSV format')
    parser.add_argument('-t', '--tissue', help='Name of the tissue or group being analyzed')
    parser.add_argument('-n', '--names', nargs=2, help='Assign names to each dataset k1, k2')
    parser.add_argument('-mul', '--max_utr3_length', type=int, default=8000, help='Maximum 3\'UTR length to consider, default: 8000')
    parser.add_argument('-ml', '--min_utr3_length', type=int, default=20, help='Minimum 3\'UTR length to consider, default: 20')
    parser.add_argument('-u', '--utr3s', default='/home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.sorted.gz', help='A bed file containing 3\'UTR coordinates, the name column should be transcript ids. Default is: /home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.sorted.gz')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output files to, default is current directory')
    
    args = parser.parse_args()

    if not args.kleats_1 or not args.kleats_2:
        sys.exit('Must specify -k1 and -k2')
    
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    
    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())
    if not args.names:
        names = ('A', 'B')
    else:
        names = args.names
    # Store the results of running the script
    results = {}
    
    if args.transcripts_to_genes:
        t2g = args.transcripts_to_genes
        t2g = parseT2G(t2g)
    else:
        sys.exit('Must specify -t2g')

    utr3s = filterUtr3s(utr3s, t2g)
    
    statistics_file_path = os.path.join(args.outdir, '{}_statistics'.format(args.tissue))
    lengths_file_path = os.path.join(args.outdir, '{}_lengths'.format(args.tissue))
    lengths_file = open(lengths_file_path, 'w')
    lengths_file_header = 'TISSUE\tGENE'
    for name in names:
        lengths_file_header += '\t{}_UTR3_LENGTHS'.format(name.upper())
    lengths_file_header += '\n'
    lengths_file.write(lengths_file_header)

    k1 = {}
    k2 = {}
    for kleat in args.kleats_1:
        groupKleatCleavage(Kleat.KleatFile(kleat), k1)
    for kleat in args.kleats_2:
        groupKleatCleavage(Kleat.KleatFile(kleat), k2)

    avgs = {'k1': [], 'k2': []}
    for chrom in k1:
        if chrom not in k2:
            continue
        for gene in k1[chrom]:
            if gene not in k2[chrom]:
                continue
            try:
                utr3s[gene]
            except KeyError:
                continue
            lens_k1 = []
            lens_k2 = []
            for utr3 in utr3s[gene]:
                if (utr3.end - utr3.start) > args.max_utr3_length:
                    continue
                elif (utr3.end - utr3.start) < args.min_utr3_length:
                    continue
                for site in k1[chrom][gene]:
                    if (site > utr3.end) or (site < utr3.start):
                        continue
                    length = getUtr3Length(utr3, site)
                    #if length >= args.min_utr3_length:
                    lens_k1.append(length)
                for site in k2[chrom][gene]:
                    if (site > utr3.end) or (site < utr3.start):
                        continue
                    length = getUtr3Length(utr3, site)
                    #if length >= args.min_utr3_length:
                    lens_k2.append(length)
            if not lens_k1 or not lens_k2:
                continue
            avg_k1 = sum(lens_k1)/len(lens_k1)
            avg_k2 = sum(lens_k2)/len(lens_k2)
            avgs['k1'].append(avg_k1)
            avgs['k2'].append(avg_k2)
            lengths_file.write('{}\t{}\t'.format(args.tissue, gene))
            lengths_file.write((',').join([str(x) for x in sorted(lens_k1)]))
            lengths_file.write('\t')
            lengths_file.write((',').join([str(x) for x in sorted(lens_k2)]))
            lengths_file.write('\n')
                
    lengths_file.close()
    with open(statistics_file_path, 'w') as o:
        ttest = stats.ttest_ind(avgs['k1'], avgs['k2'])
        o.write('t-statistic\t{}\n'.format(ttest[0]))
        o.write('p-value\t{}'.format(ttest[1]))
