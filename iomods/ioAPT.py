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
                if line.startswith('#'): continue
                l = line.strip().split()
                if l[0] == 'probeset_id':
                    ped['pedlist'] = l[1:]
                for animal in l[1:]:
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
                if line.startswith('#'): continue
                name = line.strip().split()[0]
                if name == 'probeset_id': continue
                mark['marklist'].append(name)
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

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        mark = self.updateMark(input)
        marklist = mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if l[0] == 'probeset_id':
                    names = l[1:]
                    print('Found %d names\n' % len(names))
                    continue
                marker,genos = l[0],l[1:]
                #marker = mark[prb]
                try:
                    m1,m2 = emark[marker]['a1'],emark[marker]['a2']
                except KeyError:
                    #print('Marker %s not found' % marker)
                    continue
                icol = mark[marker]['rank']
                for ind,name in enumerate(names):
                    irow = ped[name]['rank']
                    el = genos[ind]
                    if el == '0':
                        gen[irow,icol] = int(m1+m1)
                    elif el == '1':
                        gen[irow,icol] = int(m1+m2)
                    elif el == '2':
                        gen[irow,icol] = int(m2+m2)
        results = {'ped':ped,'mark':mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes out new official CIGENE genotype format
        Results contain the genotypes, markerlist and pedigree as were found in the input file
        Any external markerlist/pedigree are stored in self.mark/self.ped
        """
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
