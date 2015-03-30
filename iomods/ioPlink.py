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
                l = line.strip().split()
                family,animal,father,mother,sex,ph = l[0],l[1],l[2],l[3],l[4],l[5]
                ped['pedlist'].append(animal)
                if animal not in ped:
                    ped[animal] = {'father':father,
                                    'mother':mother,
                                    'family':family,
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
        raise Exception('Missing marker information')

    def readFile(self,input,emark=None):
        """
        Reads whole file into memory and returns a dictionary with a numpy array
        """
        if not self.mark: raise Exception('Missing marker file (.MAP)')
        marklist = self.mark['marklist']
        ped = self.updatePed(input)
        pedlist = ped['pedlist']
        gen = np.zeros((len(pedlist),len(marklist)))
        with open(input,'r') as fin:
            for line in fin:
                if line.strip().startswith('#'): continue
                l = line.strip().split()
                if len(l) == len(marklist)+1 or len(l) == len(marklist)*2+1: animal,geno = l[0],l[1:]
                elif len(l) == len(marklist)+6 or len(l) == len(marklist)*2+6: family,animal,father,mother,sex,ph,geno = l[0],l[1],l[2],l[3],l[4],l[5],l[6:]
                elif len(l) == 0: pass
                else:
                    raise Exception('Found %d genotypes for %s, expected %d markers\n' % (len(l)-3,l[0],len(marklist)) )
                ra = ped[animal]['rank']
                if len(geno) == len(marklist):
                    for i,m in enumerate(marklist):
                        rm = mark[m]['rank']
                        a = geno[i]
                        if len(a) == 2:
                            gen[ra,rm] = int(self.trans.get(a[0],'0')+self.trans.get(a[1],'0'))
                        elif len(a) == 1:
                            gen[ra,rm] = int(self.trans.get(a,'0')+self.trans.get(a,'0'))
                        elif len(a) == 4: # Includes either 'DEL' or 'INS'
                            if a[:3] == 'DEL': gen[ra,rm] =  int('5'+self.trans.get(a[3],'0'))
                            elif a[1:] == 'DEL': gen[ra,rm] = int(self.trans.get(a[0],'0')+'5')
                            elif a[:3] == 'INS': gen[ra,rm] = int('6'+self.trans.get(a[3],'0'))
                            elif a[1:] == 'INS': gen[ra,rm] = int(self.trans.get(a[0],'0')+'6')
                            else:
                                raise Exception('Error in importing, unknown allele: %s' % a )
                        elif len(a) == 6:
                            if a == 'DELDEL': gen[ra,rm] =  int('5'+'5')
                            elif a == 'DELINS': gen[ra,rm] =  int('5'+'6')
                            elif a == 'INSDEL': gen[ra,rm] =  int('6'+'5')
                            elif a == 'INSINS': gen[ra,rm] =  int('6'+'6')
                            else:
                                raise Exception('Error in importing, unknown allele: %s' % a )
                        else:
                            raise Exception('Unknown combination of alleles: %s' % a )
                elif len(geno) == 2*len(marklist):
                    for i,m in enumerate(marklist):
                        rm = self.mark[m]['rank']
                        a1 = geno[i*2]
                        a2 = geno[i*2+1]
                        gen[ra,rm] = int(self.trans.get(a1,'0')+self.trans.get(a2,'0'))
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
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            for animal in pedlist:
                try: ra = results['ped'][animal]['rank']
                except KeyError: continue
                father = ped[animal]['father']
                mother = ped[animal]['mother']
                family = ped[animal]['family']
                sex = ped[animal]['sex']
                ph = ped[animal]['phe']
                if sex == '0': sex = '2'
                fout.write('%s\t%s\t%s\t%s\t%s\t%s' % (family,animal,father,mother,sex,ph))
                for i,m in enumerate(marklist):
                    try: rm = results['mark'][m]['rank']
                    except KeyError: continue
                    a = gen[ra,rm]
                    if a == 0:
                        a = '00'
                    else:
                        a = '%.0f' % a
                    fout.write('\t%s\t%s' % (a[0],a[1]))
                fout.write('\n')
        if not output2: return
        with open(output2,'w') as fout:
            for m in marklist:
                if m not in results['mark']: continue
                ch = mark[m]['chrom']
                pos = mark[m]['pos']
                fout.write('%s\t%s\t%s\t%s\n' % (ch,m,0,pos)) 
    
