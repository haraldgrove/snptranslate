#! /usr/bin/env python

import sys
import gzip
import time

t = time.time()

if len(sys.argv) < 2:
    print('Usage: vcf2geno.py <infile> <outfile>')
    sys.exit(0)

trans = {'A':'1','C':'2','G':'3','T':'4','1':'1','2':'2','3':'3','4':'4'}

pedigree = {}
try:
    with open(sys.argv[3],'r') as fin:
        for line in fin:
            l = line.strip().split()
            if len(l) < 1: continue
            pedigree[l[0]] = {'father':l[1],'mother':l[2]}
except:
    pass

markers = {} # Markerlist specified by user.
header = []
try:
    with open(sys.argv[3],'r') as fin:
        for line in fin:
            l = line.strip().split()
            if len(l) < 1: continue
            markers[l[1]] = {'chrom':l[0],'pos':l[3]}
            header.append(l[1])
except:
    pass

genos = []

if sys.argv[1][-3:] == '.gz':
    op = gzip.open
else:
    op = open

reading = False
count = 0
genos = {}
snps = {}
with op(sys.argv[1],'r') as fin:
    for line in fin:
        l = line.strip().split()
        if l[0] == '[Data]':
            fin.next()
            reading = True
            continue
        if not reading: continue
        try:
            marker,name,a1,a2,gc = l[:5]
        except:
            count += 1
            continue
        if len(markers) > 0 and marker not in markers: continue # Skip markers not in output
        if len(markers) == 0 and marker not in snps:
            snps[marker] = 1
            header.append(marker)
        if name not in genos:
            genos[name] = {}
        genos[name][marker] = a1,a2
sys.stdout.write('Skipped %d lines.\n' % count)

with open(sys.argv[2],'w') as fout:
    fout.write('#\t%s\n' % '\t'.join(header))
    for sample in genos:
        fout.write('%s\t%s\t%s' % (sample,'0','0'))
        for mark in header:
            try:
                a1,a2 = genos[sample][mark]
            except KeyError:
                a1,a2 = '0','0'
            fout.write('\t%s\t%s' % (a1,a2))
        fout.write('\n')

sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
