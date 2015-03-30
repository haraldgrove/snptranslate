#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# IO-functions for CRIMAP file format

import sys
import ioGeno
import libPed
import libMark

class Geno(ioGeno.Geno):

    def importFile(self):

        def convert(s):
            if len(s) == 1: return s,s
            elif len(s) == 2: return s[0],s[1]

        """ Merges duplicate samples """
        self.ped = libPed.Ped()
        head,firstline = False,True
        for line in self.fgeno:
            l = line.strip()
            if firstline:
                l = l.strip(',').split(',')
                self.mark = libMark.Mark(l)
                firstline = False
            else:
                l = l.split(',')
                animal,genolist,genos = l[0],l[1:],[]
                for el in genolist:
                    if len(el) == 1: el = el+el
                    if len(el) == 0:
                        sys.stderr.write('Error in importFile: %s\n' % (animal) )
                        sys.exit(1)
                    genos.append(el[0])
                    genos.append(el[1])
                self.addGenotype(animal,genos)
                self.ped.addAnimal(animal,'0','0','F0','3')
        self.ped.updateSex()
        

    def __str__(self):
        pass
        
    def fastPrint(self,fout):
        pass 
        
    def printData(self, fout, ped = [],markers = []):
        """ Use when modifications have been done to the markerlist
        """
        pass

def writeData(outfile, cm, ped, markobj, chrom, warn = False, outfile2 = None, opt = ''):
    pass

def trans(c):
    if c in 'A1': return '1'
    elif c in 'C2': return '2'
    elif c in 'G3': return '3'
    elif c in 'T4': return '4'
    elif c in '5': return '5'
    else: return '0'
 
def main():
    print "Not a standalone program, exiting."
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
