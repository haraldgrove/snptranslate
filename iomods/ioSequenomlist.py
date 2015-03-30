#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2013
# This is GNU GPL Software: http://www.gnu.org/


import sys
import ioGeno
import libPed
import libMark

class Geno(ioGeno.Geno):

    def importFile(self):
        """ Merges duplicate samples """
        self.ped = libPed.Ped()
        self.mark = libMark.Mark()
        marks = {}
        genos = {}
        reading = False
        firstline = True
        self.gcdict = {}
        for line in self.fgeno:
            l = line.strip().split(',')
            if line.startswith('#'): continue
            if len(l) == 6: (animal,marker,a,t1,gc,t3) = l[0],l[1],l[2],l[3],l[4],l[5]
            elif len(l) < 4: continue
            else:
                sys.stderr.write('Unexpected number of elements: %s \n' % line.strip())
                sys.exit(1)
            if animal == 'Sample': continue
            if a == 'NA': a1,a2 = '0','0'
            elif len(a) == 1: a1,a2 = a,a
            elif a == 'DEL': a1,a2 = 'D','D'
            elif a == 'INS': a1,a2 = 'I','I'
            elif '.' in a:
                try: a1,a2 = a.split('.')
                except ValueError:
                    sys.stderr.write('Unknown allele %s\n' % a)
                    sys.exit(1)
                if a2 == 'DEL': a2 = 'D'
                if a1 == 'DEL': a1 = 'D'
                if a2 == 'INS': a2 = 'I'
                if a1 == 'INS': a1 = 'I'
            elif len(a) == 2: a1,a2 = a[0],a[1]
            else:
                sys.stderr.write('Unknown allele %s\n' % a)
                sys.exit(1)
            if animal not in genos: genos[animal],self.gcdict[animal] = {},{}
            if marker not in marks: marks[marker] = 1
            genos[animal][marker] = [a1,a2]
            self.gcdict[animal][marker] = gc
        for mark in marks:
            self.mark.addMarker(mark,'99')
        for animal in genos:
            gen = []
            for mark in marks:
                gen += genos[animal][mark]
            self.ped.addAnimal(animal,'0','0','F0')
            self.addGenotype(animal,gen)
        

    def __str__(self):
        pass
        
    def fastPrint(self,fout):
        pass 
        
    def printData(self, fout, ped = [],markers = []):
        """ Use when modifications have been done to the markerlist
        """
        pass

def writeData(outfile, cm, ped, markobj, chrom, warn = False, outfile2=None, opt=''):
    print "This option is not available at this point in time"
    sys.exit(1)
    fout = open(outfile,'w')
    for animal in ped.getAnimals():
        fout.write(animal)
        for marker in markobj.getMarkers():
            try: a1,a2 = cm[animal,marker]
            except: a1,a2 = '0','0'
            fout.write('\t%s\t%s' % (trans(a1),trans(a2)) )
        fout.write('\n')
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
