#!/usr/bin/env python

# plotphase, version 1.0, 2014-07-24

from __future__ import division, print_function
import sys
import argparse
import numpy as np

class SNP(object):

    def __init__(self):
        self.ped = {}
        self.pedlist = []
        self.mark = {}
        self.marklist = []

    def readGenos(self,genofile):
        self.gen = np.zeros((len(self.ped)*2,len(self.mark)))
        self.gen[:] = np.nan
        with open(genofile,'r') as fin:
            icol = 0
            for line in fin:
                if line.startswith('#'): continue 
                l = line.strip().split()
                if len(l) < 1: continue
                irow = 0
                for irow,a in enumerate(l[5:]):
                    self.gen[irow,icol] = int(a)
                icol += 1

    def readPedigree(self,pedfile):
        with open(pedfile,'r') as fin:
            count = 0
            fin.next() # Assuming there is always a header line in the file
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) < 1: continue
                family,name,miss,father,mother,sex,phenotype = l
                if name == '0': continue
                if name not in self.ped:
                    self.ped[name] = {'father':father,
                                      'mother':mother,
                                      'family':[family],
                                      'sex':sex,
                                      'phe':phenotype,
                                      'rank':count,
                                      'children':[]}
                    count += 2
                    self.pedlist.append(name)
                else:
                    sys.stderr.write('Duplicated sample %s\n' % name)
                    sys.exit(1)

    def updatePed(self):
        # Assign children to parents and set sex (Plink style) of parents
        self.fatherlist = []
        self.motherlist = []
        for n in self.ped:
            father,mother = self.ped[n]['father'],self.ped[n]['mother']
            if father in self.ped:
                self.ped[father]['sex'] = '1'
                self.ped[father]['children'].append(n)
                if father not in self.fatherlist: self.fatherlist.append(father)
            if mother in self.ped:
                self.ped[mother]['sex'] = '2'
                self.ped[mother]['children'].append(n)
                if mother not in self.motherlist: self.motherlist.append(mother)

    def readMarkers(self,markerfile):
        with open(markerfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                chrom,name,position,a1,a2 = l[0],l[1],l[2],l[3],l[4]
                if name not in self.mark:
                    self.mark[name] = {'chrom':chrom,
                                       'pos':position,
                                       'a1':a1,
                                       'a2':a2,
                                       'rank':count}
                    count += 1
                    self.marklist.append(name)

#*****************************************************************************************************

    def writeGenos(self,pedfile,mapfile,outformat):
        def trans1(a,m):
            if a == 0: return m[0]
            if a == 1: return m[1]
            return '0'
        def trans2(a,m):
            if a == 0: return '0'
            if a == 1: return '1'
            return '?'

        mlist = [self.mark[m]['a1']+self.mark[m]['a2'] for m in self.marklist]
        fout = open(pedfile,'w')
        fmap = open(mapfile,'w')
        if outformat.lower() == 'plink':
            for sample in self.pedlist:
                irow = self.ped[sample]['rank']
                p = self.ped[sample]
                father,mother,family,sex,phe = p['father'],p['mother'],p['family'][0],p['sex'],p['phe']
                fout.write('%s\t%s\t%s\t%s\t%s\t%s' % (family,sample,father,mother,sex,phe) )
                hap1 = self.gen[irow,:]
                hap2 = self.gen[irow+1,:]
                for h1,h2,m in zip(hap1,hap2,mlist):
                    fout.write('\t%s\t%s' % (trans1(h1,m),trans1(h2,m)))
                fout.write('\n')
            for mark in self.marklist:
                m = self.mark[mark]
                fmap.write('%s\t%s\t%s\t%s\n' % (m['chrom'],mark,'0',m['pos']))
        elif outformat.lower() == 'ldhat':
            fout.write('%d\t%d\t%d\n' % (len(self.pedlist)*2,len(self.marklist),1) )
            for sample in self.pedlist:
                irow = self.ped[sample]['rank']
                hap1 = self.gen[irow,:]
                hap2 = self.gen[irow+1,:]
                fout.write('>%s_1\n%s\n' % (sample,''.join([trans2(h) for h in hap1])) )
                fout.write('>%s_2\n%s\n' % (sample,''.join([trans2(h) for h in hap2])) )
            fmap.write('%d\t%s\t%s\n' % (len(self.marklist),self.mark[self.marklist[-1]]['pos'],'C') )
            for mark in self.marklist:
                fmap.write('%s\n' % self.mark[self.marklist[-1]]['pos'])
        fmap.close()
        fout.close()

def convert2Haploview(args):
    samplelist = []
    with open(args.sample,'r') as fin:
        fin.next()
        fin.next()
        for line in fin:
            l = line.strip().split()
            samplelist.append(['0',l[1]])
            samplelist.append(['0',l[1]])
    with open(args.haps,'r') as fin:
        for line in fin:
            l = line.strip().split()
            chr,marker,pos,a1,a2,alleles = l[0],l[1],l[2],l[3],l[4],l[5:] # 1 m1 198 3 1 0 0 0 1 0 0
            a = [a1,a2]
            for i,al in enumerate(alleles):
                samplelist[i].append(a[int(al)])
    with open(args.outfile,'w') as fout:
        for el in samplelist:
            fout.write('%s\n' % '\t'.join(el))

def main():
    parser = argparse.ArgumentParser(description='Converts Shapeit haplotype files to Plink or LDhat.')
    parser.add_argument('sample',help='Sample-file')
    parser.add_argument('haps',help='Haps-file')
    parser.add_argument('-o','--outfile',help='Output Haploview file')
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    args = parser.parse_args()
    convert2Haploview(args)

if __name__ == '__main__':
    main()
