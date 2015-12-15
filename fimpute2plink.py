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
        if s in ['A','1']: return 'A'
        if s in ['C','2']: return 'C'
        if s in ['G','3']: return 'G'
        if s in ['T','4']: return 'T'
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
            if len(l) >= 7: # Plink MAP, with three more columns showing reference and alternative alleles and an alias
                chrom,name,gendist,position,a1,a2,alias = l[0],l[1],l[2],l[3],l[4],l[5],l[6]
            elif len(l) == 6: # Plink MAP, with two more columns showing major and minor alleles
                chrom,name,gendist,position,a1,a2 = l[0],l[1],l[2],l[3],l[4],l[5]
            else:
                raise Exception('Map file requires columns 5 and 6 to be marker alleles\n')
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

def readPedigree(pedigreefile):
    with open(pedigreefile,'r') as fin:
        ped = {'pedlist':[]}
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if len(l) == 0: continue
            if len(l) >= 3: 
                name,father,mother = l[0],l[1],l[2]
            if name not in ped:
                ped[name] = {'father':father,
                              'mother': mother}
                count += 1
                ped['pedlist'].append(name)
    return ped

def convertFile(args):
    def trans(a,m1,m2):
        if a == '0': return m1+m1
        if a == '1': return m1+m2
        if a == '2': return m2+m2
        if a == '3': return m1+m2
        if a == '4': return m2+m1
        return '00'
    mark = readMarkers(args.mapfile)
    ped = readPedigree(args.pedfile)
    with open(args.infile,'r') as fin, open(args.output,'w') as fout:
        #fout.write('#\t%s\n' % ('\t'.join([m for m in mark['marklist'] if mark[m]['chrom']==args.chrom])))
        fin.next()
        for line in fin:
            l = line.strip().split()
            if len(l) < 1: continue
            animal,chip,geno = l
            father,mother = '0','0'
            family = father
            sex,pheno = '0','-9'
            fout.write('%s\t%s\t%s\t%s\t%s\t%s' % (family,animal,father,mother,sex,pheno))
            for i,m in enumerate(mark['marklist']):
                if args.chrom != '0' and mark[m]['chrom'] != args.chrom: continue
                g = trans(geno[i],mark[m]['a1'],mark[m]['a2'])
                fout.write('\t%s\t%s' % (g[0],g[1]))
            fout.write('\n')

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Converts from FImpute to Genos')
    parser.add_argument('-i','--infile', help='FImpute file')
    parser.add_argument('-o','--output',help='Output file')
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    parser.add_argument('-m','--mapfile',help='Map file')
    parser.add_argument('-p','--pedfile',help='Pedigree')
    parser.add_argument('-c','--chrom',help='Chromosome',default='0')
    args = parser.parse_args()
    convertFile(args)

if __name__ == '__main__':
    main()
