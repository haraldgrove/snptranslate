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

    def updatePedMark(self,input=None):
        """
        Collects necessary (additional) pedigree information from the input file
        """
        pcount = 0
        mcount = 0
        ped = {'pedlist':[]}
        mark = {'marklist':[]}
        firstline = True
        with open(input,'r') as fin:
            for line in fin:
                if line.startswith('[Header]'):
                    for i in xrange(0,8):
                        fin.next()
                    continue
                l = line.strip().split('\t')
                if firstline:
                    for animal in l:
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
                        else:
                            raise Exception('Duplicate samples')
                    firstline = False
                    continue
                marker = l[0]
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
        firstline = True
        with open(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                if line.startswith('[Header]'):
                    for i in xrange(0,9): fin.next()
                    continue
                elif firstline:
                    firstline = False
                    continue
                l = line.strip().split('\t')
                marker = l[0]
                rm = mark[marker]['rank']
                for i,g in enumerate(l[1:]):
                    ra = ped[pedlist[i]]['rank']
                    try:
                        a,gc = g.split('|')
                    except:
                        a = g
                        gc = -1
                    try:
                        if float(gc) > 0 and float(gc) < 0.0:
                            continue
                    except ValueError:
                        pass
                    gen[ra,rm] = int(self.trans.get(a[0],'0')+self.trans.get(a[1],'0'))
                    gcscore[ra,rm] = gc
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
            fout.write('%s\n' % '\t'.join(pedlist))
            for i,m in enumerate(marklist):
                rm = results['mark'][m]['rank']
                fout.write(m)
                for animal in pedlist:
                    ra = results['ped'][animal]['rank']
                    a = gen[ra,rm]
                    if a == 0:
                        a = '00'
                    else:
                        a = '%.0f' % a
                    try: gc = results['gc'][ra,rm]
                    except KeyError: gc = -1
                    fout.write('\t%s\t%d\n' % (a,gc))
    
