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
        ped = {'pedlist':[],'familylist':[]}
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
                if sex == '0': sex = '2'
                ped[name] = {'father':father,
                             'mother':mother,
                             'family':[family],
                             'sex':sex,
                             'phe':phenotype,
                             'pos':count,
                             'children':[]}
                count += 1
                ped['pedlist'].append(name)
                if family not in ped['familylist']: ped['familylist'].append(family)
            else:
                if family not in ped[name]['family']:
                    ped[name]['family'].append(family)
                    #ped['pedlist'].append(name)
                    if family not in ped['familylist']: ped['familylist'].append(family)
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
            ped[mother]['sex'] = '2'
            ped[mother]['children'].append(n)

def convertFile(infile,outfile,rped):
    with open(infile,'r') as fin, open(outfile,'w') as fout:
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split(None,3)
            if len(l) < 1: continue
            name,father,mother,geno = l[0],l[1],l[2],l[3]
            if rped is None:
                fid,sex,pheno = '0','3','-9'
            else:
                fid = rped[name]['family']
                sex = rped[name]['sex']
                pheno = rped[name]['phe']
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (fid,name,father,mother,sex,pheno,geno))

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Converts from CIGENE to Plink format')
    parser.add_argument('-i','--inputfile', help='Input file')
    parser.add_argument('-p','--pedigree',help='Pedigree file')
    parser.add_argument('-o','--output',help='Output file')
    args = parser.parse_args()
    # Reads pedigree and marker information
    if args.pedigree is not None and args.pedigree != 'None':
        rped = readPedigree(args.pedigree)
    else:
        rped = None
    convertFile(args.inputfile,args.output,rped)


if __name__ == '__main__':
    t = time.time()
    main()
    sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
