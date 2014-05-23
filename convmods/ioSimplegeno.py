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

    def translate(self,line):
        """ convert a string to a Geno object """
        
        def trans(g,m):
            if g[0] != g[1]: return '1'
            if g[0] == m[0]: return '0'
            if g[0] == m[1]: return '2'
            return 'nan'

        lmark = self.mark['marklist']
        l = line.strip().split()
        animal,genos = l[0],l[1:]
        rgen = [trans(genos[i*2:i*2+2],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
        return rgen # [sample,a1,a2,a3,....,an]

    def write(self,fout,line):
        """ Converts a line to a genotype-line and writes it """
        
        def trans(a,m):
            if a == '0': return m[0]+'\t'+m[0]
            if a == '1': return m[0]+'\t'+m[1]
            if a == '2': return m[1]+'\t'+m[1]
            return '0\t0'

        lmark = self.mark['marklist']
        animal = line[0]
        father,mother = self.ped[animal]['father'],self.ped[animal]['mother']
        fout.write('%s\t%s\t%s\t%s\n' % (animal,father,mother,
                                         '\t'.join([trans(line[i+1],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)])))
                

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
