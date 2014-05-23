#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys

class Geno(object):
    def __init__(self, ped=None, mark=None):
        if len(mark['marklist']) == 0:
            sys.stderr.write('DMU format requires a marker file\n')
            sys.exit(1)
        self.ped = ped
        self.mark = mark
        self.header = False # If a header line has been written

    def updateped(self):
        pass

    def updatemark(self):
        """ Changes markeralleles from 'ACGT' to '1234' """
        trans = {'A':'1','1':'1','C':'2','2':'2','G':'3','3':'3','T':'4','4':'4'}
        lmark = self.mark['marklist']
        for mark in lmark:
            self.mark[mark]['a1'] = trans.get(self.mark[mark]['a1'],'0')
            self.mark[mark]['a2'] = trans.get(self.mark[mark]['a2'],'0')

    def translate(self,line):
        """ convert a string to a Geno object """
        if line.startswith('#'): return None
        return line.strip().split()

    def write(self,fout,line):
        """ Converts a line to a genotype-line and writes it """
        if not self.header:
            fout.write('#\t%s\n' % '\t'.join([m for m in self.mark['marklist']]))
            self.header = True
        try:
            fout.write('%s\n' % ('\t'.join([l for l in line])))
        except TypeError:
            return

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
