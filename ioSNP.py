#!/usr/bin/env python

# Version 1.0

# from __future__ import division, print_function
import sys
import argparse
import numpy as np
import time
import os

def readPedigree(pedfile):
    """ Columns: name,father,mother,family,sex,phenotype """
    with open(pedfile,'r') as fin:
        ped = {'pedlist':[]}
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            name,father,mother,family,sex,phenotype = '0','0','0','0','3','-9'
            if name == 'pedlist':
                sys.stderr.write('"pedlist" is not a legal samplename\n')
                sys.exit(1)
            if len(l) > 0: name = l[0]
            if len(l) > 1: father = l[1]
            if len(l) > 2: mother = l[2]
            if len(l) > 3: family = l[3]
            if len(l) > 4: sex = l[4]
            if len(l) > 5: phenotype = l[5]
            if name == '0': continue
            if name not in ped:
                ped[name] = {'father':father,
                             'mother':mother,
                             'family':family,
                             'sex':sex,
                             'phe':phenotype,
                             'pos':count,
                             'children':[]}
                count += 1
                ped['pedlist'].append(name)
            else:
                pass
                #if family not in ped[name]['family']:
                #    ped[name]['family'].append(family)
                #    ped['pedlist'].append(name)
    updatePed(ped)
    return ped

def updatePed(ped):
    # Assign children to parents and set sex of parents
    for n in ped['pedlist']:
        father,mother = ped[n]['father'],ped[n]['mother']
        if father in ped:
            ped[father]['sex'] = '1'
            ped[father]['children'].append(n)
        if mother in ped:
            ped[mother]['sex'] = '0'
            ped[mother]['children'].append(n)

def readMarkers(markerfile,ch):
    """
    Columns options:
      name,position,allele1,allele2 [,chromosome] (BEAGLE)
      chromosome,name,gendist,position (PLINK)
      name
    """
    def trans(s):
        if s in ['A','1']: return '1'
        if s in ['C','2']: return '2'
        if s in ['G','3']: return '3'
        if s in ['T','4']: return '4'
        return '0'

    with open(markerfile,'r') as fin:
        mark = {'marklist':[]}
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if len(l) == 0: continue
            name,position,a1,a2,chrom,rank,alias = '0',0,None,None,'0',0,None
            if name == 'marklist':
                sys.stderr.write('"marklist" is not a legal markername\n')
                sys.exit(1)
            if len(l) == 7: # Plink MAP, with three more columns showing reference and alternative alleles and an alias
                chrom,name,gendist,position,a1,a2,alias = l[0],l[1],l[2],l[3],l[4],l[5],l[6]
            elif len(l) == 6: # Plink MAP, with two more columns showing major and minor alleles
                chrom,name,gendist,position,a1,a2 = l[0],l[1],l[2],l[3],l[4],l[5]
            elif len(l) == 5:
                name,position,a1,a2,chrom = l[0],int(l[1]),l[2],l[3],l[4]
                rank = count
            elif len(l) == 4:
                chrom,name,gendist,position = l[0],l[1],l[2],l[3]
            elif len(l) == 1:
                name = l[0]
            if name not in mark and (chrom == ch or not ch):
                mark[name] = {'chrom':chrom,
                              'pos':position,
                              'a1':trans(a1),
                              'a1x':0,
                              'a2':trans(a2),
                              'a2x':0,
                              'rank':count,
                              'alias': alias}
                count += 1
                mark['marklist'].append(name)
    return mark

def main():
    # Gather a list of available conversion formats
    pathname = os.path.dirname(sys.argv[0])
    filelist = os.listdir(pathname+'/iomods')
    fl = []
    for fi in filelist:
        if fi.startswith('io') and fi.endswith('.py'):
            fl.append(fi[2:-3])
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Converts between the following SNP-formats:\n %s' % ' '.join(fl))
    parser.add_argument('-i','--input', help='Input file')
    parser.add_argument('-n','--informat',help='Input file format',default='Geno')
    parser.add_argument('-p','--pedigree',help='Pedigree file')
    parser.add_argument('-m','--markers',help='Marker file')
    parser.add_argument('-c','--chrom',help='Chromosome')
    parser.add_argument('-o','--output',help='Output file')
    parser.add_argument('--output2',help='Extra output file')
    parser.add_argument('-u','--outformat',help='Output file format',default='Geno')
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    args = parser.parse_args()
    if args.input == 'None': args.input = None
    if args.pedigree == 'None': args.pedigree = None
    if args.markers == 'None': args.markers = None
    # Reads pedigree and marker information
    if args.pedigree:
        rped = readPedigree(args.pedigree)
    else:
        rped = None
    if args.markers:
        rmark = readMarkers(args.markers,args.chrom)
    else:
        rmark = None
    # Connecting to the input/output classes
    try:
        __import__('iomods.io'+args.informat) # No idea why this is necessary
    except ImportError:
        parser.error('Unknown module: "%s"\n' % (str(args.informat)) )
    try:
        __import__('iomods.io'+args.outformat) # No idea why this is necessary
    except KeyError: #except ImportError:
        parser.error('Unknown module: "%s"\n' % (str(args.outformat)) )
    try:
        inform = sys.modules['iomods.io'+args.informat]
        gen1 = inform.Geno(rped,rmark)
    except ImportError:
        parser.error('Unsupported inputformat: "%s"\n' % str(args.informat) )
    try:
        outform = sys.modules['iomods.io'+args.outformat]
        gen2 = outform.Geno(rped,rmark)
    except ImportError:
        parser.error('Unsupported outputformat: "%s"\n' % str(args.outformat) )
    # Reads genotypes and writes them to file
    genos = gen1.readFile(args.input,rmark)
    if args.output2: gen2.writeFile(args.output,genos,args.output2)
    else: gen2.writeFile(args.output,genos)


if __name__ == '__main__':
    t = time.time()
    main()
    sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
