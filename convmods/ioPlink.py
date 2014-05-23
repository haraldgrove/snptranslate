#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys

class Geno(object):
    def __init__(self, ped=None, mark=None):
        self.ped = ped
        self.mark = mark
        self.header = False # If a header line has been written

    def updateped(self):
        trans = {'1':'1','0':'2'}
        lped = self.ped['pedlist']
        for name in lped:
            self.ped[name]['sex'] = trans.get(self.ped[name]['sex'],'3')

    def updatemark(self):
        """ Changes markeralleles from 'ACGT' to '1234' """
        trans = {'A':'1','1':'1','C':'2','2':'2','G':'3','3':'3','T':'4','4':'4'}
        lmark = self.mark['marklist']
        for mark in lmark:
            self.mark[mark]['a1'] = trans.get(self.mark[mark]['a1'],'0')
            self.mark[mark]['a2'] = trans.get(self.mark[mark]['a2'],'0')

    def translate(self,line):
        """ convert a string to a Geno object """
        
        def trans(g,m):
            if g[0] != g[1]: return '1'
            if g[0] == m[0]: return '0'
            if g[0] == m[1]: return '2'
            return 'nan'

        l = line.strip().split()
        animal,genos = l[1],l[6:]
        if len(self.mark['marklist']) > 0:
            lmark = self.mark['marklist']
            rgen = [animal]+[trans(genos[i*2:i*2+2],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
        else:
            rgen = [animal]+genos
        return rgen # [sample,a1,a2,a3,....,an]

    def write(self,fout,line):
        """ Converts a line to a genotype-line and writes it """
        
        def trans(a,m):
            if a == '0': return m[0]+'\t'+m[0]
            if a == '1': return m[0]+'\t'+m[1]
            if a == '2': return m[1]+'\t'+m[1]
            return '0\t0'

        if not line: return
        animal = line[0]
        if animal in self.ped:
            father,mother = self.ped[animal]['father'],self.ped[animal]['mother']
            sex = self.ped[animal]['sex']
            phe = self.ped[animal]['phe']
            family = self.ped[animal]['family'][0]
        else:
            father,mother,sex,phe,family = '0','0','3','-9','0'
        if len(self.mark['marklist']) > 0:
            lmark = self.mark['marklist']
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (family,animal,father,mother,sex,phe,
                    '\t'.join([trans(line[i+1],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)])))
        else:
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (family,animal,father,mother,sex,phe,'\t'.join(line[1:])))
                

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
