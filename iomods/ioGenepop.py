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
        family = -1
        pops = False
        ped = {'pedlist':[]}
        with open(input,'r') as fin:
            for line in fin:
                l = line.strip().split()
                if l[0].lower() == "pop":
                    pops = True
                    family += 1
                if pops and not l[0].lower() == "pop":
                    animal = l[0]
                    ped['pedlist'].append(animal)
                    if animal not in ped:
                        ped[animal] = {'father':'0',
                                        'mother':'0',
                                        'family':'F' + str(family),
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
        count = 0
        mark = {'marklist':[]}
        print input
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
        #raise Exception('Missing marker information in file')

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        if not self.mark: raise Exception('Missing marker file (.MAP)')
        mark = self.mark
        marklist = mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        pops = False
        count = -1
        with open(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                l = line.strip().split()
                if l[0].lower() == "pop":
                    pops = True
                    count = count + 1
                if pops and not l[0].lower() == "pop":
                    if len(l) == len(marklist)+2: animal,geno,family = l[0],l[2:],'F'+str(count)
                    elif len(l) == 0: pass
                    else:
                        raise Exception('Found %d genotypes for %s, expected %d markers\n' % (len(l)-3,l[0],len(marklist)) )
                    ra = ped[animal]['rank']
                    if len(geno) == len(marklist):
                        for i,m in enumerate(marklist):
                            rm = mark[m]['rank']
                            a = geno[i]
                            if len(a) == 4:
                                gen[ra,rm] = int(self.trans.get(a[1],'0')+self.trans.get(a[3],'0'))
                            else:
                                raise Exception('Unknown combination of alleles: %s' % a )
        results = {'ped':ped,'mark':self.mark,'gen':gen}
        return results

    def writeFile(self,output,results,output2=None):
        """
        Writes Plink
        """
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        try: familylist = ped['familylist']
        except KeyError: familylist = ['0']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('Comment line\n')
            for m in marklist:
                fout.write('%s\n' % m)
            fout.write('pop\n')
            for fam in familylist:
                for animal in pedlist:
                    try: ra = results['ped'][animal]['rank']
                    except KeyError: continue
                    if fam not in ped[animal]['family']: continue
                    fout.write('%s' % (animal))
                    for i,m in enumerate(marklist):
                        try: rm = results['mark'][m]['rank']
                        except KeyError:
                            #sys.stderr.write('Missing: %s\n' % m)
                            continue
                        a = gen[ra,rm]
                        if a == 0:
                            a = '00'
                        else:
                            a = '%.0f' % a
                        fout.write('\t0%s0%s' % (a[0],a[1]))
                    fout.write('\n')
        if not output2: return
        with open(output2,'w') as fout:
            for m in marklist:
                if m not in results['mark']: continue
                ch = mark[m]['chrom']
                pos = mark[m]['pos']
                fout.write('%s\t%s\t%s\t%s\n' % (ch,m,0,pos)) 
    
