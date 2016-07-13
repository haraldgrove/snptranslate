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

    def isgzip(inputfile):
        """
        Determines if input file ends in .gz
        """
        if inputfile.lower().endswith(('.gz')):
            return True
    
    def updatePedMark(self,isgzip,input):
        """
        Collects necessary (additional) pedigree information from the input file
        """
        pcount = 0
        mcount = 0
        ped = {'pedlist':[]}
        mark = {'marklist':[]}
        csample = 1
        cmark = 0
        
        if isgzip(input):
            op = gzip.open
        else:
            op = open
        with op(input,'r') as fin:
            for line in fin:
                if line.startswith('[Header]'):
                    for i in xrange(0,8): fin.next()
                    header = fin.next().strip().split('\t')
                    csample = header.index('Sample ID')
                    cmark = header.index('SNP Name')
                    if csample == -1:
                        raise Exception('Header line found, but column with "Sample ID" is missing.')
                    if cmark == -1:
                        raise Exception('Header line found, but columne with "SNP Name" is missing.')
                    continue
                l = line.strip().split('\t')
                try:
                    animal = l[csample]
                    marker = l[cmark]
                except IndexError:
                    continue
                if animal not in ped:
                    ped['pedlist'].append(animal)
                    ped[animal] = {'father':'0',
                                    'mother':'0',
                                    'family':'F0',
                                    'sex':'3',
                                    'phe':'-9',
                                    'rank':pcount,
                                    'children':[]}
                    pcount += 1
                if marker not in mark:
                    mark['marklist'].append(marker)
                    mark[marker] = {'chrom':'0',
                                    'pos':mcount,
                                    'a1':None,
                                    'a1x':0,
                                    'a2':None,
                                    'a2x':0,
                                    'rank':mcount}
                    mcount += 1
        return ped,mark

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        ped,mark = self.updatePedMark(input)
        marklist = mark['marklist']
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        gcscore = np.zeros((len(pedlist),len(marklist)))
        csample = 1
        cmark = 0
        if isgzip(input):
            op = gzip.open
        else:
            op = open
        with op(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                if line.startswith('[Header]'):
                    for i in xrange(0,8): fin.next()
                    header = fin.next().strip().split('\t')
                    csample = header.index('Sample ID')
                    cmark = header.index('SNP Name')
                    if csample == -1:
                        raise Exception('Header line found, but column with "Sample ID" is missing.')
                    if cmark == -1:
                        raise Exception('Header line found, but columne with "SNP Name" is missing.')
                    continue
                l = line.strip().split('\t')
                try:
                    animal = l[csample]
                    marker = l[cmark]
                except IndexError:
                    continue
                ra = ped[animal]['rank']
                a1,a2 = l[2],l[3]
                try: gc = l[4]
                except IndexError: gc = -1
                #try:
                #    if float(gc) > 0 and float(gc) < 0.2:
                #        continue
                #except ValueError:
                #    pass
                rm = mark[marker]['rank']
                gen[ra,rm] = int(self.trans.get(a1,'0')+self.trans.get(a2,'0'))
                try: gcscore[ra,rm] = gc
                except ValueError: gcscore[ra,rm] = -1
        results = {'ped':ped,'mark':mark,'gen':gen,'gc':gcscore}
        return results

    def writeFile(self,output,results):
        """
        Writes out new official CIGENE genotype format
        """
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            for animal in pedlist:
                ra = results['ped'][animal]['rank']
                for i,m in enumerate(marklist):
                    rm = results['mark'][m]['rank']
                    a = gen[ra,rm]
                    if a == 0:
                        a = '00'
                    else:
                        a = '%.0f' % a
                    try: gc = results['gc'][ra,rm]
                    except KeyError: gc = -1
                    fout.write('%s\t%s\t%s\t%s\t%d\n' % (m,animal,a[0],a[1],gc))
    
