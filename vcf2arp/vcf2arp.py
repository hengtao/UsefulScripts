#!/usr/bin/env python3

import re,time,os,sys
import traceback
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def vcf2arp(vcfile, groupfile, arpfile, header):
    groupdict = dict()
    haplodict = dict()
    i = 0
    genotype = ""
    samplepos = dict()
    with open(groupfile, 'r') as g, open(vcfile, 'r') as v, open(arpfile, 'w+') as o:
        for line in g:
            if line.strip():
                groupcontent = re.split("\s+", line.strip())
                samplename = groupcontent[0]
                groupname = groupcontent[1]
                if groupname not in groupdict:
                    groupdict[groupname] = list()
                    groupdict[groupname].append(samplename)
                else:
                    groupdict[groupname].append(samplename)
        for line in v:
            if line.startswith("##"):
                pass
            elif line.startswith("#CHROM"):
                linecontent = line.strip().split("\t")
                for samplename in linecontent[9:]:
                    haplodict[samplename] = ''
                    samplepos[i] = samplename
                    i += 1
            else:
                linecontent = line.strip().split("\t")
                ref_allele = linecontent[3]
                alt_allele = linecontent[4]
                genotypelist =  linecontent[9:]
                for i in range(len(genotypelist)):
                    if genotypelist[i].startswith("0/0"):
                        genotype = ref_allele*2
                    elif genotypelist[i].startswith("1/1"):
                        genotype = alt_allele * 2
                    elif genotypelist[i].startswith("0/1"):
                        genotype = ref_allele + alt_allele
                    else:
                        genotype = "??"
                    haplodict[samplepos[i]] += genotype
        groupnum = len(groupdict)
        o.write('''
[Profile]
    Title={header}
    NbSamples={groupnum}
    GenotypicData=0
    DataType=DNA
    LocusSeparator=NONE
[Data]
    [[Samples]]\n'''.format(header = header, groupnum = groupnum))
        for groupname in groupdict:
            o.write("\t\tSampleName=\"{groupname}\"\n".format(groupname = groupname))
            samplesize = len(groupdict[groupname])
            o.write("\t\tSampleSize=\"{samplesize}\"\n".format(samplesize = samplesize))
            o.write("\t\tSampleData= {\n")
            for samplename in groupdict[groupname]:
                if samplename in haplodict:
                    o.write("\t\t\t{samplename}\t1\t{haplotype}\n".format(samplename = samplename, haplotype = haplodict[samplename]))
            o.write("\t\t}\n")
        o.write("""
[[Structure]]
    StructureName="A group of {groupnum} populations analyzed for DNA"
    NbGroups=1\n""".format(groupnum = groupnum))
        o.write('\tGroup= {\n')
        for groupname in groupdict:
            o.write("\t\t\"{groupname}\"\n".format(groupname = groupname))
        o.write("\t}\n")

def main(args):
    vcfile = args.vcfile
    groupfile = args.groupfile
    arpfile = args.arpfile
    header = args.header
    vcf2arp(vcfile, groupfile, arpfile, header)

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class=HelpFormatter, description='''
Usage:

python3 {scriptpath}/vcf2arp.py <args> <args>....


'''.format(scriptpath=scriptpath))

    parse.add_argument('-vcfile', '--vcfile', required=True, dest="vcfile", help="file in vcf format",
                       type=str, nargs='?')
    parse.add_argument('-groupfile', '--groupfile', required=True, dest="groupfile", help="group file define indivs affilitation, ",
                       type=str, nargs='?')
    parse.add_argument('-arpfile', '--arpfile', required=True, dest="arpfile",
                       help="Output file in arp format", type=str, nargs='?')
    parse.add_argument('-header', '--header', required=True, dest="header",
                       help="header of profile to describe the arp", type=str, nargs='?')


    args = parse.parse_args()

    main(args)


