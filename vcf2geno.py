import sys
import gzip

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

header = '#'
genos = []

if sys.argv[1][-3:] == '.gz':
    op = gzip.open
else:
    op = open

with op(sys.argv[1],'r') as fin:
    for line in fin:
        if line.startswith('##'): continue
        if line.startswith('#'):
            l = line.strip('#').strip().split()
            if l[8] == 'FORMAT':
                samples = l[9:]
            else:
                raise Exception('Missing field "FORMAT", aborting\n')
            for i in xrange(len(samples)):
                if len(pedigree) > 0 and samples[i] not in pedigree: continue
                genos.append('%s\t0\t0' % samples[i])
            continue
        l = line.strip().split()
        if len(l) < 9: continue
        header += '\t%s' % l[2]
        m = trans.get(l[3],'0')+trans.get(l[4],'0')
        co = 0
        for i,gt in enumerate(l[9:]):
            if len(pedigree) > 0 and samples[i] not in pedigree: continue
            a1,a2 = m[int(gt[0])],m[int(gt[2])]
            genos[co] += '\t%s\t%s' % (a1,a2)
            co += 1

with open(sys.argv[2],'w') as fout:
    fout.write('%s\n' % header)
    for sample in genos:
        fout.write('%s\n' % sample)
