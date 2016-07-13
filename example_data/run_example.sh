#!/bin/bash
echo "Converting FinalReport_toy.txt"
../ioSNP.py -i FinalReport_toy.txt -n Genomelist -o FinalReport_toy.ped -u Plink
echo "Converting FinalReport_toy.txt.gz"
../ioSNP.py -i FinalReport_toy.txt.gz -n Genomelist -o FinalReport_toy.gz.ped -u Plink
