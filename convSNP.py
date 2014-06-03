#!/usr/bin/env python

# Version 1.0

# from __future__ import division, print_function
import sys
import argparse
import numpy as np
import time

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
                                  'family':[family],
                                  'sex':sex,
                                  'phe':phenotype,
                                  'pos':count,
                                  'children':[]}
                count += 1
                ped['pedlist'].append(name)
            else:
                if family not in ped[name]['family']:
                    ped[name]['family'].append(family)
                    ped['pedlist'].append(name)
    return ped

def updatePed():
    # Assign children to parents and set sex of parents
    for n in ped:
        father,mother = ped[n]['father'],ped[n]['mother']
        if father in ped:
            ped[father]['sex'] = '1'
            ped[father]['children'].append(n)
        if mother in ped:
            ped[mother]['sex'] = '0'
            ped[mother]['children'].append(n)

def readMarkers(markerfile):
    """
        Columns options:
          name,position,allele1,allele2,chromosome
          chromosome,rank,name,position
          name
    """
    with open(markerfile,'r') as fin:
        mark = {'marklist':[]}
        count = 0
        known = True
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if len(l) == 0: continue
            name,position,a1,a2,chrom,rank = '0',0,'0','0','0',0
            if name == 'marklist':
                sys.stderr.write('"marklist" is not a legal markername\n')
                sys.exit(1)
            if len(l) == 5 or len(l) == 4 and len(l[2]) == len(l[3]) == 1:
                name,position,a1,a2 = l[0],int(l[1]),l[2],l[3]
                if len(l) == 5:
                    chrom = l[4]
                rank = count
            elif len(l) == 4:
                chrom,rank,name,position = l[0],count,l[2],l[3]
            elif len(l) == 1:
                name = l[0]
            if name not in mark:
                mark[name] = {'chrom':chrom,
                                   'pos':position,
                                   'a1':a1,
                                   'a1x':0,
                                   'a2':a2,
                                   'a2x':0,
                                   'rank':rank}
                if a1 == '0' or a2 == '0': known = False
                count += 1
                mark['marklist'].append(name)
    return mark

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('-i','--genotypes', dest='genosfile',help='Input file')
    parser.add_argument('-n','--informat',dest='inform',help='Input file format',default='Geno')
    parser.add_argument('-p','--pedigree',dest='pedigreefile',help='Pedigree file')
    parser.add_argument('-m','--markers',dest='markerfile',help='Marker file')
    parser.add_argument('-o','--output',dest='outputfile',help='Output file')
    parser.add_argument('-u','--outformat',dest='outform',help='Output file format',default='Geno')
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    args = parser.parse_args()
    if args.genosfile == 'None': args.genosfile = None
    if args.pedigreefile == 'None': args.pedigreefile = None
    if args.markerefile == 'None': args.markerfile = None
    # Reads pedigree and marker information
    if args.pedigreefile:
        rped = readPedigree(args.pedigreefile)
    else:
        rped = {'pedlist':[]}
    if args.markerfile:
        rmark = readMarkers(args.markerfile)
    else:
        rmark = {'marklist':[]}
    # Connecting to the input/output classes
    try:
        __import__('convmods.io'+args.inform) # No idea why this is necessary
    except ImportError:
        parser.error('Unknown module: "%s"\n' % (str(args.inform)) )
    try:
        __import__('convmods.io'+args.outform) # No idea why this is necessary
    except ImportError:
        parser.error('Unknown module: "%s"\n' % (str(args.outform)) )
    try:
        inform = sys.modules['convmods.io'+args.inform]
        gen1 = inform.Geno(rped,rmark)
    except ImportError:
        parser.error('Unsupported inputformat: "%s"\n' % str(args.inform) )
    try:
        outform = sys.modules['convmods.io'+args.outform]
        gen2 = outform.Geno(rped,rmark)
    except ImportError:
        parser.error('Unsupported outputformat: "%s"\n' % str(args.outform) )
    # Changes default parameters, such as male/female coding
    gen2.updateped()
    gen2.updatemark()
    # Reads genotypes and writes them to file
    with open(args.genosfile,'r') as fin, open(args.outputfile,'w') as fout:
        for line in fin:
            gen2.write(fout,gen1.translate(line))


if __name__ == '__main__':
    t = time.time()
    main()
    sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
