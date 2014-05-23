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
        # Could check that male==1 and female==0
        pass

    def updatemark(self):
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

        if line.strip().startswith('#'):
            return None
        elif len(self.mark['marklist']) > 0:
            lmark = self.mark['marklist']
            l = line.strip().split()
            if len(l) == len(lmark)+1:
                animal,geno = l[0],l[1:]
                rgen = [animal]+[trans(geno[i],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
            elif len(l) == len(lmark)*2+1:
                animal,geno = l[0],l[1:]
                rgen = [animal]+[trans(geno[i*2:i*2+2],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
            elif len(l) == len(lmark)*2+3:
                animal,geno = l[0],l[3:]
                rgen = [animal]+[trans(geno[i*2:i*2+2],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
            elif len(l) == len(lmark)+3:
                animal,geno = l[0],l[3:]
                rgen = [animal]+[trans(geno[i],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)]
            elif len(l) == 0:
                return None
            else:
                sys.stdout.write('Found %d genotypes for %s, expected %d markers\n' % (len(l)-3,l[0],len(self.mark)) )
                return None
                #sys.exit(1)
            return rgen
        else:
            l = line.strip().split()
            if len(l) < 4:
                return None
            elif len(l[3]) == 2:
                animal,geno = l[0],l[3:]
                rgen = [animal]+[el[0]+el[1] for el in geno]
            elif len(l[3]) == 1:
                animal,geno = l[0],l[3:]
                rgen = [animal]+geno
            else:
                return None
            return rgen

    def write(self,fout,line):
        """ Converts a line to a genotype-line and writes it """
        
        def trans(a,m):
            if a == '0': return m[0]+m[0]
            if a == '1': return m[0]+m[1]
            if a == '2': return m[1]+m[1]
            return '00'

        if not line: return
        animal = line[0]
        try: father = self.ped[animal]['father']
        except KeyError: father = '0'
        try: mother = self.ped[animal]['mother']
        except KeyError: mother = '0'
        if len(self.mark['marklist']) > 0:
            if not self.header:
                fout.write('#\t%s\n' % ('\t'.join([m for m in self.mark['marklist']])))
                self.header = True
            lmark = self.mark['marklist']
            fout.write('%s\t%s\t%s\t%s\n' % (animal,father,mother,
                                         '\t'.join([trans(line[i+1],self.mark[name]['a1']+self.mark[name]['a2']) for i,name in enumerate(lmark)])))
        else:
            fout.write('%s\t%s\t%s\t%s\n' % (animal,father,mother,'\t'.join(line[1:])))
                

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
