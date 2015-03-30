#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2012
# This is GNU GPL Software: http://www.gnu.org/


import sys
import numpy as np
import gzip

class Geno(object):

    def __init__(self, ped=None, mark=None):
        self.ped = ped
        self.mark = mark
        self.trans = {'1':'A','2':'C','3':'G','4':'T',
                      'A':'A','C':'C','G':'G','T':'T',
                      'DEL':'D','D':'D','5':'D',
                      'INS':'I','I':'I','6':'I'}

    def updatePed(self,input=None):
        """
        Collects necessary (additional) pedigree information from the input file
        """
        count = 0
        ped = {'pedlist':[]}
        with open(input,'r') as fin:
            for line in fin:
                if line.startswith('I'):
                    l = line.strip().split()
                    for animal in l[2::2]:
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
        if not input: raise Exception('Missing marker information')
        count = 0
        mark = {'marklist':[]}
        with open(input,'r') as fin:
            for line in fin:
                l = line.strip().split(None,2)
                if l[0] != 'M': continue
                name = l[1]
                if name not in mark:
                    mark[name] = {'chrom':'0',
                                  'pos':count,
                                  'a1':None,
                                  'a1x':0,
                                  'a2':None,
                                  'a2x':0,
                                  'rank':count}
                    count += 1
                    mark['marklist'].append(name)
        return mark

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        trans = {'1':1,'2':2,'3':3,'4':4,
                 'A':1,'C':2,'G':3,'T':4,
                 'DEL':5,'D':5,'5':5,'INS':6,'I':6,'6':6}
        mark = self.updateMark(input)
        marklist = mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            for line in fin:
                l = line.strip().split()
                if l[0] != 'M': continue
                rm = mark[l[1]]['rank']
                for i,p in enumerate(pedlist):
                    ra = ped[p]['rank']
                    a1,a2 = l[2+i*2],l[2+i*2+1]
                    gen[ra,rm] = trans.get(a1,0)*10+trans.get(a2,0)
        results = {'ped':ped,'mark':mark,'gen':gen}
        return results

    def writeFile(self,output,results):
        """
        Writes out BEAGEL v3
        """
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('I\tid\t%s\n' % '\t'.join([p+"\t"+p for p in pedlist]))
            for i,m in enumerate(marklist):
                fout.write('M\t%s' % m)
                try:
                    rm = results['mark'][m]['rank']
                except KeyError:
                    sys.stderr.write('ERROR: Missing marker "%s" in genotype file\n' % m)
                    rm = None
                for animal in pedlist:
                    try:
                        ra = results['ped'][animal]['rank']
                    except KeyError:
                        sys.stderr.write('ERROR: Missing sample "%s" in genotype file\n' % animal)
                        ra = None
                    if ra and rm:
                        a = gen[ra,rm]
                    else:
                        a = 0
                    if a == 0:
                        a = '00'
                    else:
                        a = '%.0f' % a
                    fout.write('\t%s\t%s' % (self.trans.get(a[0],'0'),self.trans.get(a[1],'0')))
                fout.write('\n')
        fout.close()
