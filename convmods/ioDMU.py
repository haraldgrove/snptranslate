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
        pass

    def translate(self,line):
        """ convert a string to a Geno object """
        def trans(self,a,mark):
            m1,m2 = self.mark[mark]['a1'],self.mark[mark]['a2']
            if a == '0': return m1+m1
            if a == '1': return m1+m2
            if a == '2': return m2+m2
            return '00'
        if line.startswith('#'): return None
        l = line.strip().split()
        return l[0]+[trans(e) for e in l[1:]]

    def writeFile(self,output,results,output2=None):
        """
        Writes out DMU format
        Results contain the genotypes, markerlist and pedigree as were found in the input file
        Any external markerlist/pedigree are stored in self.mark/self.ped
        """
        def trans(a,m1,m2):
            if a[0] != a[1]: return '1'
            if a[0] == m1: return '0'
            if a[0] == m2: return '2'
            return '-1'
        gen = results['gen']
        if self.ped: ped = self.ped
        else: ped = results['ped']
        pedlist = ped['pedlist']
        if self.mark: mark = self.mark
        else: mark = results['mark']
        marklist = mark['marklist']
        with open(output,'w') as fout:
            fout.write('#\t%s\n' % '\t'.join(marklist))
            for animal in pedlist:
                ra = results['ped'][animal]['rank']
                fout.write('%s' % (animal))
                for i,m in enumerate(marklist):
                    try:
                        rm = results['mark'][m]['rank']
                        a = gen[ra,rm]
                        if a == 0:
                            a = '00'
                        else:
                            a = '%.0f' % a
                    except KeyError:
                        if output2: continue
                        a = '00'
                        # This could be handled with more feedback to user
                    if len(a) < 2:
                        raise Exception('Unknwon allele "%s"\n' % a)
                    fout.write('\t%s' % trans(a,mark[m]['a1'],mark[m]['a2']))
                fout.write('\n')
        if not output2: return
        with open(output2,'w') as fout:
            for m in marklist:
                if m not in results['mark']: continue
                ch = mark[m]['chrom']
                pos = mark[m]['pos']
                fout.write('%s\t%s\t%s\t%s\n' % (ch,m,0,pos))

def main():
    print('Not a standalone program, exiting.')
    
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
