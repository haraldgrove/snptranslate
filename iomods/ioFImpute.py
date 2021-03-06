#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys
import numpy as np

class Geno(object):
    def __init__(self, ped=None, mark=None):
        if len(mark['marklist']) == 0:
            raise Exception('DMU format requires a marker file')
        self.ped = ped
        self.mark = mark
        self.header = False # If a header line has been written

    def updatePed(self,input=None):
        """
        Collects necessary (additional) pedigree information from the input file
        """
        count = 0
        ped = {'pedlist':[]}
        with open(input,'r') as fin:
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split(None,1)
                if len(l) < 1 or l[0] == 'ID': continue
                animal = l[0]
                ped['pedlist'].append(animal)
                if animal not in ped:
                    ped[animal] = {'father':'0',
                                    'mother':'0',
                                    'family':'F0',
                                    'sex':'3',
                                    'phe':'-9',
                                    'rank':count,
                                    'children':[]}
                count += 1
        return ped

    def updateMark(self,input=None):
        """
        Gathers information about markers, if missing
        """
        mark = {'marklist':[]}
        return mark

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        def trans(a,m1,m2):
            if a == '0': return m1+m1
            if a == '1': return m1+m2
            if a == '2': return m2+m2
            if a == '3': return m1+m2
            if a == '4': return m2+m1
            return '00'
        mark = self.updateMark(input)
        marklist = self.mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            for line in fin:
                l = line.strip().split()
                if len(l) < 1 or l[0] == 'ID': continue
                animal,chip,geno = l
                ra = ped[animal]['rank']
                for i,m in enumerate(marklist):
                    rm = self.mark[m]['rank']
                    try: a = geno[i]
                    except:
                        print('Tried accessing element %d of genotype string length %d' % (i,len(geno)))
                        sys.exit(1)
                    gen[ra,rm] = int(trans(a,self.mark[m]['a1'],self.mark[m]['a2']))
        results = {'ped':ped,'mark':self.mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes out FImpute format, does not keep phasing information
        Results contain the genotypes, markerlist and pedigree as were found in the input file
        Any external markerlist/pedigree are stored in self.mark/self.ped
        """
        def trans(a,m1,m2):
            if a[0] != a[1]: return '1'
            if a[0] == m1: return '0'
            if a[0] == m2: return '2'
            if a =='00': return '5'
            raise Exception ('Non-matching genotypes %s for marker alleles %s,%s\n' % (a,m1,m2))
            #sys.exit(1)
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('ID\tChip\tCall...\n')
            for animal in pedlist:
                ra = results['ped'][animal]['rank']
                fout.write('%s\t1\t' % (animal))
                for i,m in enumerate(marklist):
                    try:
                        rm = results['mark'][m]['rank']
                        a = gen[ra,rm]
                        if a == 0:
                            a = '00'
                        else:
                            a = '%.0f' % a
                    except KeyError:
                        if output2: continue
                        a = '00'
                        # This could be handled with more feedback to user
                    if len(a) < 2:
                        raise Exception('Unknwon allele "%s"\n' % a)
                    try: fout.write('%s' % trans(a,mark[m]['a1'],mark[m]['a2']))
                    except:
                        sys.stderr.write('Non-matching genotypes marker: %s, %s for marker alleles %s,%s\n' % (m,a,mark[m]['a1'],mark[m]['a2']))
                fout.write('\n')
        if not output2: return
        with open(output2,'w') as fout:
            for m in marklist:
                if m not in results['mark']: continue
                ch = mark[m]['chrom']
                pos = mark[m]['pos']
                fout.write('%s\t%s\t%s\t%s\n' % (ch,m,0,pos))

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
