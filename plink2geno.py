#!/usr/bin/env python

# Version 1.0

# from __future__ import division, print_function
import sys
import argparse

def readMarkers(markerfile):
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
            elif len(l) == 4:
                chrom,name,gendist,position = l[0],l[1],l[2],l[3]
            if name not in mark:
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

def convertFile(args):
    trans = {'A':'1','C':'2','G':'3','T':'4','0':'0','1':'1','2':'2','3':'3','4':'4'}
    mark = readMarkers(args.input+'.map')
    if args.output: fout = open(args.output,'w')
    else: fout = sys.stdout
    fout.write('#\t%s\n' % ('\t'.join(mark['marklist'])))
    with open(args.input+'.ped','r') as fin:
        for line in fin:
            l = line.strip().split()
            if len(l) < 1: continue
            animal,father,mother,genos = l[1],l[2],l[3],l[6:]
            fout.write('%s\t%s\t%s\t%s\n' % (animal,father,mother,'\t'.join([trans[g] for g in genos])))
    if args.output: fout.close()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Converts from Plink to Genos')
    parser.add_argument('-i','--input', help='Input')
    parser.add_argument('-o','--output',help='Output file')
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    args = parser.parse_args()
    convertFile(args)

if __name__ == '__main__':
    main()
