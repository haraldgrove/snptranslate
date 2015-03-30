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
        self.header = False

    def updateped(self):
        # Changes code for female from 2 to 0
        for sample in self.ped:
            if self.ped[sample]['sex'] == 2:
                self.ped[sample]['sex'] = 0

    def updatemark(self):
        trans = {'A':'1','1':'1','C':'2','2':'2','G':'3','3':'3','T':'4','4':'4'}
        lmark = self.mark['marklist']
        for mark in lmark:
            self.mark[mark]['a1'] = trans.get(self.mark[mark]['a1'],'0')
            self.mark[mark]['a2'] = trans.get(self.mark[mark]['a2'],'0')

    def translate(self,line):
        """ convert a string to a Geno object """
        if line.strip().startswith('#'):
            return None
        elif len(self.mark['marklist']) > 0:
            lmark = self.mark['marklist']
            l = line.strip().split()
            if len(l) == len(lmark)+1:
                return l,{}
            elif len(l) == len(lmark)*2+1:
                animal,geno = l[0],l[1:]
                return [animal]+[geno[i*2]+geno[i*2+1] for i,name in enumerate(lmark)],{}
            elif len(l) == len(lmark)*2+3:
                animal,geno = l[0],l[3:]
                return [animal]+[geno[i*2]+geno[i*2+1] for i,name in enumerate(lmark)],{'father':l[1],'mother':l[2]}
            elif len(l) == len(lmark)+3:
                return l[0]+l[3:],{'father':l[1],'mother':l[2]}
            elif len(l) == 0:
                return None
            else:
                sys.stdout.write('Found %d genotypes for %s, expected %d markers\n' % (len(l)-3,l[0],len(self.mark)) )
                return None
        else:
            animal,geno = l[0],l[3:]
            return [animal]+[geno[i*2]+geno[i*2+1] for i,name in enumerate(lmark)],{'father':l[1],'mother':l[2]}

    def write(self,fout,line,extra):
        """ Converts a line to a genotype-line and writes it """
        if not line: return
        animal = line[0]
        try:
            father = self.ped[animal]['father']
        except KeyError:
            if 'father' in extra: father = extra['father']
            else: father = '0'
        try:
            mother = self.ped[animal]['mother']
        except KeyError:
            if 'mother' in extra: mother = extra['mother']
            else: mother = '0'
        if not self.header:
            if len(self.mark['marklist']) > 0:
                fout.write('#\t%s\n' % ('\t'.join([m for m in self.mark['marklist']])))
            else:
                fout.write('#\t%s\n' % ('\t'.join(['m'+str(i) for i in enumerate(line[1:])])))
            self.header = True
        fout.write('%s\t%s\t%s\t%s\n' % (animal,father,mother,'\t'.join(line[1:])))
                

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
