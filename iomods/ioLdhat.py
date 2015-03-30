#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2014
# This is GNU GPL Software: http://www.gnu.org/

# Description:

import sys
import numpy as np
import gzip

class Geno(object):
    def __init__(self, ped=None, mark=None):
        self.ped = ped
        self.mark = mark
        self.trans = {'1':'1','2':'2','3':'3','4':'4',
                      'A':'1','C':'2','G':'3','T':'4',
                      'DEL':'5','D':'5','5':'5',
                      'INS':'6','I':'6','6':'6'}

    def updatePed(self,input=None):
        """
        Collects necessary (additional) pedigree information from the input file
        """
        count = 0
        ped = {'pedlist':[]}
        with open(input,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    marklist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) == len(marklist)+1 or len(l) == len(marklist)*2+1: animal,father,mother = l[0],'0','0'
                elif len(l) == len(marklist)+3 or len(l) == len(marklist)*2+3: animal,father,mother = l[0],l[1],l[2]
                elif len(l) == 0: pass
                else:
                    raise Exception('Incorrect number of columns')
                ped['pedlist'].append(animal)
                if animal not in ped:
                    ped[animal] = {'father':father,
                                    'mother':mother,
                                    'family':'F0',
                                    'sex':'3',
                                    'phe':'-9',
                                    'rank':count,
                                    'children':[]}
                count += 1
        for animal in ped['pedlist']:
            father,mother = ped[animal]['father'],ped[animal]['mother']
            if father in ped:
                ped[father]['sex'] = 1
            if mother in ped:
                ped[mother]['sex'] = 0
        return ped

    def updateMark(self,input=None):
        """
        Gathers information about markers, if missing
        """
        if not input: raise Exception('Missing marker information')
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
        mark = self.updateMark(input)
        marklist = mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                l = line.strip().split()
                ra = ped[animal]['rank']
                for i,m in enumerate(marklist):
                    rm = mark[m]['rank']
                    a = geno[i]
                    gen[ra,rm] = int(self.trans.get(a[0],'0')+self.trans.get(a[1],'0'))
        results = {'ped':ped,'mark':mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes out new official CIGENE genotype format
        Results contain the genotypes, markerlist and pedigree as were found in the input file
        Any external markerlist/pedigree are stored in self.mark/self.ped
        """
        def trans(a):
            if a in 'A1': return '0'
            if a in 'C2': return '1'
            if a in 'G3': return '2'
            if a in 'T4': return '3'
            return '?'
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('%d\t%d\t%d\n' % (len(pedlist)*2,len(marklist),1))
            for animal in pedlist:
                if len(animal)>30: raise Exception('Sample names have to be less than 30 characters')
                ra = results['ped'][animal]['rank']
                hap1 = ''
                hap2 = ''
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
                    hap1 += trans(a[0])
                    hap2 += trans(a[1])
                fout.write('>%s_1\n' % animal)
                fout.write('\n'.join([hap1[h:min(h+80,len(marklist))] for h in xrange(0,len(marklist),80)]))
                fout.write('\n')
                fout.write('>%s_2\n' % animal)
                fout.write('\n'.join([hap2[h:min(h+80,len(marklist))] for h in xrange(0,len(marklist),80)]))
                fout.write('\n')
        if not output2: return
        with open(output2,'w') as fout:
            fout.write('%d\t%.3f\t%s\n' % (len(marklist),int(mark[marklist[-1]]['pos'])/1000.0,'C'))
            for m in marklist:
                if m not in results['mark']: continue
                pos = int(mark[m]['pos'])/1000.0
                fout.write('%.3f\n' % (pos))
