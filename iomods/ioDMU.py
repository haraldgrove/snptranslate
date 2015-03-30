#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys

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
                if len(l) < 1: continue
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
        count = 0
        mark = {'marklist':[]}
        with open(input,'r') as fin:
            line = fin.next()
            if not line.startswith('#'): raise Exception('Missing marker information')
            mark['marklist'] = line.strip('#').strip().split()
        for name in mark['marklist']:
            if name not in mark:
                mark[name] = {'chrom':'0',
                              'pos':count,
                              'a1':None,
                              'a1x':0,
                              'a2':None,
                              'a2x':0,
                              'rank':count}
                count += 1
        return mark

    def readFile(self,input):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        def trans(a,m1,m2):
            if a == '0': return m1+m1
            if a == '1': return m1+m2
            if a == '2': return m2+m2
            return '00'
        mark = self.updateMark(input)
        marklist = mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                l = line.strip().split()
                if len(l) == len(marklist)+1: animal,geno = l[0],l[1:]
                else:
                    raise Exception('Found %d genotypes for %s, expected %d markers\n' % (len(l)-1,l[0],len(marklist)) )
                ra = ped[animal]['rank']
                for i,m in enumerate(marklist):
                    rm = mark[m]['rank']
                    a = geno[i]
                    gen[ra,rm] = int(trans(a,self.mark[m]['a1'],self.mark[m]['a2']))
        results = {'ped':ped,'mark':mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes out DMU format
        Results contain the genotypes, markerlist and pedigree as were found in the input file
        Any external markerlist/pedigree are stored in self.mark/self.ped
        """
        def trans(a,m1,m2):
            if a[0] != a[1]: return '1'
            if a[0] == m1: return '0'
            if a[0] == m2: return '2'
            return '-1'
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('#\t%s\n' % '\t'.join(marklist))
            for animal in pedlist:
                ra = results['ped'][animal]['rank']
                fout.write('%s' % (animal))
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
                    fout.write('\t%s' % trans(a,mark[m]['a1'],mark[m]['a2']))
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
