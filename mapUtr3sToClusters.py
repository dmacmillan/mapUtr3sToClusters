import argparse
import os
import sys
import pysam
import gzip
import copy
import math
import subprocess as sub
import scipy.stats as stat
import numpy as np

def parseT2G(fpath):
    results = {}
    with gzip.open(fpath, 'rb') as f:
        for line in f:
            line = line.strip().split('\t')
            results[line[0]] = line[1]
    return results

def genScatterPdf(coords_file, tissue, names, window, extend, mul, ss, rs, skew, kurt):
    if ss:
        ss = 'skipSame'
    if rs:
        rs = 'removeSymmetric'
    directory = os.path.dirname(os.path.abspath(coords_file))
    outfile = os.path.join(directory, 'plotScatter.r')
    title = 'utr3_lengths_{}_{}_{}'.format(tissue, names[0], names[1])
    params = ['window={}'.format(window)]
    params.append('maxUtr3Length={}'.format(mul))
    params.append('extend={}'.format(extend))
    if ss:
        params.append(ss)
    if rs:
        params.append(rs)
    params = ('_').join(params)
    pdf = title + '_' + params
    main = title
    rscript = 'library(LSD)\n' \
        'library(e1071)\n' \
        'setwd("{}")\n' \
        'coords = read.table("coordinates")\n' \
        'pdf("{}.pdf")\n' \
        'png("{}.png", res=150, width=1000, height=1000)\n' \
        'par(oma=c(0,0,2,0))\n' \
        'heatscatter(coords[,1], coords[,2], log="xy", main="{}", ' \
        'xlab="{}", ylab="{}")\n' \
        'abline(0,1,col=rgb(1,0,0,0.5))\n' \
        'dev.off()\n' \
        'dev.off()\n'.format(directory, pdf, pdf, main, names[0], names[1])
    dens_title = 'density_{}_{}_{}'.format(tissue, names[0], names[1])
    main = dens_title
    pdf = dens_title + '_' + params 
    density = 'density = unlist(read.table("distances"))\n' \
        'skew = round(skewness(density),3)\n' \
        'kurt = round(kurtosis(density),3)\n' \
        'pdf("{}.pdf")\n' \
        'png("{}.png", res=150, width=1000, height=1000)\n' \
        'plot(density(density), main="{}")\n' \
        'mtext(paste0("skewness: ", skew), side=1, line=3, adj=0, col="blue")\n' \
        'mtext(paste0("kurtosis: ", kurt), side=1, line=4, adj=0, col="blue")\n' \
        'abline(v=0, col=rgb(1,0,0,0.5))\n' \
        'dev.off()\n' \
        'dev.off()\n'.format(pdf, pdf, main)
    with open(outfile, 'w') as f:
        f.write(rscript + density)
    return rscript, outfile
    
# Where line equation is: (ax + by + c = 0)
# coord is a tuple (x0, y0)
def distanceToLine(coord, a=-1, b=1, c=0):
    factor = 1
    if coord[0] > coord[1]:
        factor = -1
    num = abs( (a*coord[0]) + (b*coord[1]) + c )
    denom = math.sqrt( (a**2) + (b**2) )
    return (factor*num) / denom

def getUtr3Length(utr3, site):
    if utr3.strand == '+':
        return site - utr3.start
    else:
        return utr3.end - site

def getCoords(list_of_lengths, skip_same=True, window=15):
    coords = []
    skipped = 0
    lol = list_of_lengths
    for i in xrange(len(lol)-1):
        for j in xrange(i+1, len(lol)):
            for k in xrange(len(lol[i])):
                for l in xrange(len(lol[j])):
                    coord = (lol[i][k], lol[j][l])
                    if abs(coord[0] - coord[1]) <= window:
                        if skip_same:
                            skipped += 1
                            continue
                    coords.append(coord)
    return (coords, skipped)

def removeSymmetric(coords, window=15):
    indices = set()
    length = len(coords)
    for i in xrange(length-1):
        for j in xrange(i+1, length):
            x0, x1, y0, y1 = coords[i][0], coords[i][1], coords[j][0], coords[j][1]
            if (abs(x0 - y1) <= window) and (abs(x1 - y0) <= window):
                indices.add(i)
                indices.add(j)
    result = [x for i,x in enumerate(coords) if i not in indices]
    return (result, len(indices))

def getMultiStopGenes(utr3s, t2g):
    genes = {}
    for utr3 in utr3s:
        gene = t2g[utr3.name]
        if gene not in genes:
            if utr3.strand == '+':
                genes[gene] = set((utr3.start,))
            else:
                genes[gene] = set((utr3.end,))
        else:
            if utr3.strand == '+':
                genes[gene].add(utr3.start)
            else:
                genes[gene].add(utr3.end)
    return set([x for x in genes if (len(genes[x]) > 1)])

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
            #elif current.start == previous.end:
            #    print gene
        if store:
            keep[gene] = genes[gene]
    return keep

#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#import seaborn as sns

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Given a KLEAT file return 3\'UTR coordinates for each cleavage site, otherwise return None')
    
    parser.add_argument('clusters', nargs='+', help='Clustered KLEAT bed file for specific dataset')
    parser.add_argument('-t2g', '--transcripts_to_genes', help='A mapping of transcripts to genes in TSV format')
    parser.add_argument('-t', '--tissue', help='Name of the tissue or group being analyzed')
    parser.add_argument('-n', '--names', nargs='+', help='Assign names to each dataset, with respect to clusters order')
    #parser.add_argument('--xlim', type=float, default=30000, help='The limit of the x-axis to display. Default is 30000')
    #parser.add_argument('--ylim', type=float, default=30000, help='The limit of the y-axis to display. Default is 30000')
    parser.add_argument('-w', '--window', type=int, default=15, help='The minimum number of bases apart for the cleavage sites to be considered different. Default: 15')
    parser.add_argument('-e', '--extend', type=int, default=100, help='Extend the search for cleavage calls downstream of the 3\'UTR end by this much. Default is 100')
    parser.add_argument('-mul', '--max_utr3_length', type=int, default=8000, help='Maximum 3\'UTR length to consider, default: 8000')
    parser.add_argument('-ss', '--skip_same', action='store_true', help='Do not yield coordinates for cleavage sites that are the same (within the window parameter of distance between them)')
    parser.add_argument('-rs', '--remove_symmetric', action='store_true', help='Remove sets of coordinates that are mirror images of each other across the 45 degree line within +/- the window parameter\'s distance of each other. I.E. if the window is 15 then (10,12) and (19,8) would be removed because 10 is within +/- 15 of 8 and 12 is within +/- 15 of 19')
    parser.add_argument('-u', '--utr3s', default='/home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.sorted.gz', help='A bed file containing 3\'UTR coordinates, the name column should be transcript ids. Default is: /home/dmacmillan/annotations/ensembl/ensembl.fixed.sorted.utr3s_only.sorted.gz')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output files to, default is current directory')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    
    utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())
    datasets = [pysam.TabixFile(x, parser=pysam.asBed()) for x in args.clusters]
    if not args.names:
        names = [chr(97+i) for i in range(len(datasets))]
    else:
        names = args.names
    # Store the results of running the script
    results = {}
    
    if args.transcripts_to_genes:
        t2g = args.transcripts_to_genes
        t2g = parseT2G(t2g)
    else:
        sys.exit('Must specify -t2g')

    #print 'Num 3\'UTRs: {}'.format(len([x for x in utr3s]))
    #utr3s = pysam.tabix_iterator(open(args.utr3s), parser=pysam.asBed())
    utr3s = filterUtr3s(utr3s, t2g)
    #print 'Num 3\'UTRs: {}'.format(sum([len(utr3s[gene]) for gene in utr3s]))
    #sys.exit()
    
    distances = []
    scatter_file_path = os.path.join(args.outdir, 'coordinates')
    scatter_file = open(scatter_file_path, 'w')
    
    total_symmetric = 0
    total_zero = 0
    for gene in utr3s:
        utr3 = max(utr3s[gene], key = lambda x: x.end)
        if utr3.end - utr3.start > args.max_utr3_length:
            continue
        dataset_lengths = []
        for idx, dataset in enumerate(datasets):
            if utr3.strand == '+':
                try:
                    sites = dataset.fetch(utr3.contig, utr3.start, utr3.end + args.extend)
                except ValueError:
                    continue
            else:
                try:
                    sites = dataset.fetch(utr3.contig, utr3.start - args.extend, utr3.end)
                except ValueError:
                    continue
            lengths = []
            for site in sites:
                length = getUtr3Length(utr3, site.end)
                lengths.append(length)
            dataset_lengths.append(lengths)
        try:
            if len(filter(None, dataset_lengths)) < 2:
                continue
        except TypeError:
            continue
        coords, num_zero = getCoords(dataset_lengths, skip_same=args.skip_same, window=args.window)
        total_zero += num_zero
        if args.remove_symmetric:
            coords, num_symmetric = removeSymmetric(coords, window=args.window)
            total_symmetric += num_symmetric
        for coord in coords:
            scatter_file.write('{}\t{}\n'.format(coord[0], coord[1]))
        distances += [distanceToLine(x) for x in coords]
    
    scatter_file.close()
    
    _max = round(max(distances), 2)
    _min = round(min(distances), 2)
    median = round(np.percentile(distances, 50), 2)
    mean = round(sum(distances)/len(distances))
    mode = round(
        max(set(distances), key = distances.count), 2
    )
    stdev = round(np.std(distances), 2)
    skew = round(stat.skew(distances), 2)
    kurt = round(stat.kurtosis(distances), 2)
    total_pos = sum(i > 0 for i in distances)
    total_neg = len(distances) - total_pos

    genScatterPdf(scatter_file_path, 
                  args.tissue,
                  names,
                  args.window,
                  args.extend,
                  args.max_utr3_length, 
                  args.skip_same, 
                  args.remove_symmetric,
                  skew,
                  kurt
    )

    with open(os.path.join(args.outdir, 'statistics'), 'w') as o:
        o.write('min:\t{}\n'.format(_min))
        o.write('max:\t{}\n'.format(_max))
        o.write('mean:\t{}\n'.format(mean))
        o.write('median:\t{}\n'.format(median))
        o.write('mode:\t{}\n'.format(mode))
        o.write('stdev:\t{}\n'.format(stdev))
        o.write('skew:\t{}\n'.format(skew))
        o.write('total_pos:\t{}\n'.format(total_pos))
        o.write('total_neg:\t{}\n'.format(total_neg))
        o.write('total_zero:\t{}\n'.format(total_zero))
        o.write('total_symmetric:\t{}\n'.format(total_symmetric))

    with open(os.path.join(args.outdir, 'distances'), 'w') as d:
        d.write(('\n').join([str(x) for x in distances]))
