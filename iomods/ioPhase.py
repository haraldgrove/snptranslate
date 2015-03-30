#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2012
# This is GNU GPL Software: http://www.gnu.org/


import sys
import ioGeno
import libPed
import libMark

class Geno(ioGeno.Geno):

    def importFile(self):
        """ Merges duplicate samples """
        print "This option is not available at this point in time"
        sys.exit(1)
        self.ped = libPed.Ped()
        self.mark = libMark.Mark()
        self.ped.addAnimal(animal,dam,sire,family,sex)
        self.addGenotype(animal,geno.split())
        self.ped.updateSex()
        

    def __str__(self):
        pass
        
    def fastPrint(self,fout):
        pass 
        
    def printData(self, fout, ped = [],markers = []):
        """ Use when modifications have been done to the markerlist
        """
        pass

def writeData(outfile, cm, ped, markobj, chrom, warn = False, outfile2=None, opt=''):
    sep = '\t'
    fout = open(outfile,'w')
    markers = markobj.getMarkers(chrom)
    fout.write(str(len(ped))+'\n')
    fout.write(str(len(markers))+'\n')
    fout.write('P')
    for marker in markers:
        minfo = markobj[marker]
        fout.write(sep+str(minfo[3]))
    fout.write('\n'+'S'*len(markers)+'\n')
    for family in ped.getFamilies():
        for animal in ped.getFamilyMembers(family):
            fout.write(animal+'\n')
            hap1,hap2 = '',''
            for marker in markers:
                try:
                    try: a1,a2 = cm.getHapAlleles(animal,marker)
                    except AttributeError: a1,a2 = cm[animal,marker]
                except KeyError: a1,a2 = '0','0'
                if a1 in ['0','9']: a1 = '?'
                if a2 in ['0','9']: a2 = '?'
                hap1 += trans(a1)
                hap2 += trans(a2)
            fout.write(hap1+'\n')
            fout.write(hap2+'\n')
    fout.close()

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
