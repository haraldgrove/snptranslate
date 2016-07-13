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
                animal,father,mother = l[0],l[1],l[2]
                ped['pedlist'].append(animal)
                if animal not in ped:
                    ped[animal] = {'father':father,
                                    'mother':mother,
                                    'family':'F0',
                                    'sex':'3',
                                    'phe':'-9',
                                    'rank':count,
                                    'children':[]}
                else:
                    pass
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
        count = 0
        mark = {'marklist':[]}
        with open(input,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mnames = line.strip('#').strip().split()
                    mark['marklist'] = mnames
                    for name in mnames:
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
        raise Exception('Missing marker information in file')

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
            for line in fin:
                pass
                #Do something to the input line
        results = {'ped':ped,'mark':mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes Geno
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
                out1,out2 = '',''
                ra = results['ped'][animal]['rank']
                family = ped[animal]['family']
                out1 = '%s\t%s' % (family,animal)
                out2 = '%s\t%s' % (family,animal)
                for i,m in enumerate(marklist):
                    try:
                        rm = results['mark'][m]['rank']
                        a = gen[ra,rm]
                    except KeyError:
                        a = 0
                    if a == 0:
                        a = '00'
                    else:
                        a = '%.0f' % a
                    try:
                        out1 += '\t%s' % a[0]
                        out2 += '\t%s' % a[1]
                    except:
                        sys.stderr.write('ERROR: Unexpected allele: [%s,%s,%s]\n' % (animal,m,a))
                        sys.exit(1)
                fout.write('%s\n%s\n' % (out1,out2))
        if not output2: return
        with open(output2,'w') as fout:
            for m in marklist:
                if m not in results['mark']: continue
                pos = mark[m]['pos']
                fout.write('%s\t%s\n' % (m,pos)) 
    
