#!/usr/bin/env python3

import re,time,os,sys
import argparse
import subprocess
import configparser
import glob
import gzip
from concurrent.futures import ThreadPoolExecutor as ThreadPool
import logging
from Bio import SeqIO

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def snpseq(seq, value, start, end,  vcfdict):
    seq = [i for i in seq]
    if value[0] in vcfdict:
        for posit in range(start, end):
            if str(posit) in vcfdict[value[0]]:
                seq[posit - start] = vcfdict[value[0]][str(posit)]
    seq = "".join(seq)
    return seq

def RetrieveSeq(args,logger):
    gff3 = args.gff3
    vcf = args.vcf
    genome = args.genome
    genes = args.genenames
    genelist = genes.split(",")
    outdir = args.outdir
    gene = dict()
    cds = dict()
    cdslist = list()
    intron = dict()
    utr3 = dict()
    utr3list = list()
    utr5 = dict()
    utr5list = list()
    vcfdict = dict()

    snpdict = {
               "AT":"[A/T]", "AC":"[A/C]", "AG":"[A/G]",
               "TA":"[T/A]", "TC":"[T/C]", "TG":"[T/G]",
               "CA":"[C/A]", "CT":"[C/T]", "CG":"[C/G]",
               "GA":"[G/A]", "GT":"[G/T]", "GC":"[G/C]"
               }

    logger.info("Processing genome fa file...")
    record_dict = SeqIO.index(genome, "fasta")
    logger.info("Processing gff3 file...")
    flag = 0
    dictname = ''
    with open(vcf, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                linecontent = line.strip().split("\t")
                chrom = linecontent[0]
                posit = str(linecontent[1])
                refer = linecontent[3]
                alter = linecontent[4]
                varia = refer + alter
                if len(alter) > 1:
                    continue
                else:
                    if chrom not in vcfdict:
                        vcfdict[chrom] = dict()
                        vcfdict[chrom][posit] = snpdict[varia]
                    else:
                        vcfdict[chrom][posit] = snpdict[varia]

    with open(gff3, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                linecontent = line.strip().split()
                chr = linecontent[0]
                feature = linecontent[2]
                start = linecontent[3]
                end = linecontent[4]
                direction = linecontent[6]
                name = linecontent[8].split(";")[0].replace("ID=","")
                if feature == "gene":
                    for genename in genelist:
                        if genename == name:
                            if dictname != "":
                                if utr5list:
                                    utr5[dictname] = []
                                    utr5[dictname].extend(utr5list)
                                if utr3list:
                                    utr3[dictname] = []
                                    utr3[dictname].extend(utr3list)
                                if cdslist:
                                    cds[dictname] = []
                                    cds[dictname].extend(cdslist)
                            dictname = genename
                            gene[dictname] = (chr, start, end, direction)
                            flag = 1
                            cdslist = list()
                            utr3list = list()
                            utr5list = list()
                            break
                        else:
                            flag = 0
                elif feature == "three_prime_UTR" and flag == 1:
                    if dictname not in utr3:
                        #utr5[dictname] = (chr, start, end, direction)
                        utr3list.append([chr, start, end])
                    else:
                        for anylist in utr3list:
                            if int(start) < int(anylist[1]) and int(end) > int(anylist[2]):
                                utr3list.append([chr, start, end])
                                utr3list.remove(anylist)
                            elif int(start) > int(anylist[1]) and int(start) < int(anylist[2]) and int(end) > int(anylist[2]):
                                utr3list.append([chr, anylist[1], end])
                                utr3list.remove(anylist)
                            elif int(start) < int(anylist[1]) and int(anylist[2]) < int(end) and int(end) < int(anylist[2]):
                                utr3list.append([chr, start, anylist[1]])
                                utr3list.remove(anylist)
                            elif int(start) > int(anylist[2]) or int(end) < int(anylist[1]):
                                utr3list.append([chr, start, end])
                            else:
                                pass
                    #utr3[dictname] = (chr, start, end, direction)
                elif feature == "five_prime_UTR" and flag == 1:
                    if dictname not in utr5:
                        #utr5[dictname] = (chr, start, end, direction)
                        utr5list.append([chr, start, end])
                    else:
                        for anylist in utr5list:
                            if int(start) < int(anylist[1]) and int(end) > int(anylist[2]):
                                utr5list.append([chr, start, end])
                                utr5list.remove(anylist)
                            elif int(start) > int(anylist[1]) and int(start) < int(anylist[2]) and int(end) > int(anylist[2]):
                                utr5list.append([chr, anylist[1], end])
                                utr5list.remove(anylist)
                            elif int(start) < int(anylist[1]) and int(anylist[2]) < int(end) and int(end) < int(anylist[2]):
                                utr5list.append([chr, start, anylist[1]])
                                utr5list.remove(anylist)
                            elif int(start) > int(anylist[2]) or int(end) < int(anylist[1]):
                                utr5list.append([chr, start, end])
                            else:
                                pass
                        #utr5[dictname].extend([(chr, start, end, direction)])
                elif feature == "CDS" and flag == 1:
                    if dictname not in cds:
                        #cds[dictname] = (chr, start, end, direction)
                        cdslist.append([chr, start, end])
                    else:
                        for anylist in cdslist:
                            if int(start) < int(anylist[1]) and int(end) > int(anylist[2]):
                                cdslist.append([chr, start, end])
                                cdslist.remove(anylist)
                            elif int(start) > int(anylist[1]) and int(start) < int(anylist[2]) and int(end) > int(anylist[2]):
                                cdslist.append([chr, anylist[1], end])
                                cdslist.remove(anylist)
                            elif int(start) < int(anylist[1]) and int(anylist[2]) < int(end) and int(end) < int(anylist[2]):
                                cdslist.append([chr, start, anylist[1]])
                                cdslist.remove(anylist)
                            elif int(start) > int(anylist[2]) or int(end) < int(anylist[1]):
                                cdslist.append([chr, start, end])
                            else:
                                pass
    utr5[dictname] = []
    utr3[dictname] = []
    cds[dictname] = []
    utr5[dictname].extend(utr5list)
    utr3[dictname].extend(utr3list)
    print(dictname)
    print(cdslist)
    cds[dictname].extend(cdslist)
    print(cds[dictname])
    tmplist = list()
    for key, value in gene.items():
        outfile = outdir + "/" + key + ".fa"
        with open(outfile, 'w') as o:
            o.write(">{key}\n".format(key = key))
            if value[3] == "+":
                seq = record_dict[value[0]][int(value[1]) - 5001:int(value[1]) - 1].seq
                start = int(value[1]) - 5001
                end = int(value[1]) - 1
                sequence = snpseq(seq, value, start, end,  vcfdict)
                o.write("##UPSTREAM##{up5k}".format(up5k = sequence))
                #o.write("##UPSTREAM##{up5k}".format(up5k = record_dict[value[0]][int(value[1]) - 5001:int(value[1]) - 1].seq))
                if key in utr5:
                    if len(utr5[key]) == 1:
                        print(utr5[key])
                        seq = record_dict[value[0]][int(utr5[key][0][1]):int(utr5[key][0][2])].seq
                        start = int(utr5[key][0][1])
                        end = int(utr5[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##UTR5##{utr5}".format(utr5 = sequence))
                        #o.write("##UTR5##{utr5}".format(utr5 =  record_dict[value[0]][int(utr5[key][1]):int(utr5[key[2]])].seq))
                    else:
                        for anyutr5 in sorted(utr5[key]):
                            tmplist.append(anyutr5[1])
                            tmplist.append(anyutr5[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##UTR5##{utr5}".format(utr5 = sequence))
                            #o.write("##UTR5##{utr5}".format(utr5=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq
                                start = int(tmplist[suffix + 1]) + 1
                                end = int(tmplist[suffix + 2]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron = sequence))
                                #o.write("##INTRON##{intron}".format(intron = record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq))
                        tmplist = list()
                if key in cds:
                    if len(cds[key]) == 1:
                        seq = record_dict[value[0]][int(cds[key][0][1]):int(cds[key][0][2])].seq
                        start = int(cds[key][0][1])
                        end = int(cds[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##CDS##{cds}".format(cds = sequence))
                        #o.write("##CDS##{cds}".format(cds = record_dict[value[0]][int(cds[key][1]):int(cds[key[2]])].seq))
                    else:
                        for anycds in sorted(cds[key]):
                            tmplist.append(anycds[1])
                            tmplist.append(anycds[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##CDS##{cds}".format(cds = sequence))
                            #o.write("##CDS##{cds}".format(cds=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq
                                start = int(tmplist[suffix + 1]) + 1
                                end = int(tmplist[suffix + 2]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron = sequence))
                                #o.write("##INTRON##{intron}".format(intron=record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq))
                        tmplist = list()
                if key in utr3:
                    if len(utr3[key]) == 1:
                        seq = record_dict[value[0]][int(utr3[key][0][1]):int(utr3[key][0][2])].seq
                        start = int(utr3[key][0][1])
                        end = int(utr3[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##UTR3##{utr3}".format(utr3=sequence))
                        #o.write("##UTR3##{utr3}".format(utr3 =  record_dict[value[0]][int(utr3[key][1]):int(utr3[key[2]])].seq))
                    else:
                        for anyutr3 in sorted(utr3[key]):
                            tmplist.append(anyutr3[1])
                            tmplist.append(anyutr3[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##UTR3##{utr3}".format(utr3=sequence))
                            #o.write("##UTR3##{utr3}".format(utr3=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq
                                start = int(tmplist[suffix + 1]) + 1
                                end = int(tmplist[suffix + 2]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron=sequence))
                                #o.write("##INTRON##{intron}".format(intron = record_dict[value[0]][int(tmplist[suffix + 1]) + 1:int(tmplist[suffix + 2]) - 1].seq))
                        tmplist = list()
                seq = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5000].seq
                start = int(value[2]) + 1
                end = int(value[2]) + 5000
                sequence = snpseq(seq, value, start, end, vcfdict)
                o.write("##DOWNSTREAM##{down5k}".format(down5k=sequence))
                #o.write("##DOWNSTREAM##{down5k}".format(down5k=record_dict[value[0]][int(value[2]) + 1:int(value[2] + 5000)].seq))

            else:
                seq = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5001].seq.reverse_complement()
                start = int(value[2]) + 1
                end = int(value[2]) + 5000
                sequence = snpseq(seq, value, start, end, vcfdict)
                o.write("##UPSTREAM##{up5k}".format(up5k=sequence))
                #o.write("##UPSTREAM{up5k}##".format(up5k = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5000].seq.reverse_complement()))
                if key in utr5:
                    if len(utr5[key]) == 1:
                        seq = record_dict[value[0]][int(utr5[key][0][1]):int(utr5[key][0][2])].seq.reverse_complement()
                        start = int(utr5[key][0][1])
                        end = int(utr5[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##UTR5##{utr5}".format(utr5=sequence))
                        #o.write("##UTR5##{utr5}".format(utr5 =  record_dict[value[0]][int(utr5[key][1]):int(utr5[key[2]])].seq.reverse_complement()))
                    else:
                        for anyutr5 in sorted(utr5[key], reverse = True):
                            tmplist.append(anyutr5[1])
                            tmplist.append(anyutr5[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##UTR5##{utr5}".format(utr5=sequence))
                            #o.write("##UTR5##{utr5}".format(utr5=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()
                                start = int(tmplist[suffix + 3]) + 1
                                end = int(tmplist[suffix]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron=sequence))
                                #o.write("##INTRON##{intron}".format(intron = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()))
                        tmplist = list()
                if key in cds:
                    if len(cds[key]) == 1:
                        seq =  record_dict[value[0]][int(cds[key][0][1]):int(cds[key][0][2])].seq.reverse_complement()
                        start = int(cds[key][0][1])
                        end = int(cds[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##CDS##{cds}".format(cds=sequence))
                        #o.write("##CDS##{cds}".format(cds =  record_dict[value[0]][int(cds[key][1]):int(cds[key[2]])].seq.reverse_complement()))
                    else:
                        for anycds in sorted(cds[key], reverse = True):
                            tmplist.append(anycds[1])
                            tmplist.append(anycds[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##CDS##{cds}".format(cds=sequence))
                            #o.write("##CDS##{cds}".format(cds=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()
                                start = int(tmplist[suffix + 3]) + 1
                                end = int(tmplist[suffix]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron=sequence))
                                #o.write("##INTRON##{intron}".format(intron = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()))
                        tmplist = list()
                if key in utr3:
                    if len(utr3[key]) == 1:
                        seq = record_dict[value[0]][int(utr3[key][0][1]):int(utr3[key][0][2])].seq.reverse_complement()
                        start = int(utr3[key][0][1])
                        end = int(utr3[key][0][2])
                        sequence = snpseq(seq, value, start, end, vcfdict)
                        o.write("##UTR3##{utr3}".format(utr3=sequence))
                        #o.write("##UTR3##{utr3}".format(utr3 =  record_dict[value[0]][int(utr3[key][1]):int(utr3[key[2]])].seq.reverse_complement()))
                    else:
                        for anyutr3 in sorted(utr3[key], reverse = True):
                            tmplist.append(anyutr3[1])
                            tmplist.append(anyutr3[2])
                        for suffix in range(0,len(tmplist),2):
                            seq = record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()
                            start = int(tmplist[suffix])
                            end = int(tmplist[suffix + 1])
                            sequence = snpseq(seq, value, start, end, vcfdict)
                            o.write("##UTR3##{utr3}".format(utr3=sequence))
                            #o.write("##UTR3##{utr3}".format(utr3=record_dict[value[0]][int(tmplist[suffix]):int(tmplist[suffix + 1])].seq.reverse_complement()))
                            if suffix + 3 < len(tmplist):
                                seq = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()
                                start = int(tmplist[suffix + 3]) + 1
                                end = int(tmplist[suffix]) - 1
                                sequence = snpseq(seq, value, start, end, vcfdict)
                                o.write("##INTRON##{intron}".format(intron=sequence))
                                #o.write("##INTRON##{intron}".format(intron = record_dict[value[0]][int(tmplist[suffix + 3]) + 1:int(tmplist[suffix]) - 1].seq.reverse_complement()))
                        tmplist = list()
                #print(int(value[1]))
                seq = record_dict[value[0]][int(value[1]) - 5002:int(value[1]) - 1].seq.reverse_complement()
                #seq = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5000].seq.reverse_complement()
                #print(seq)
                start = int(value[1]) - 5002
                end = int(value[1]) - 1
                sequence = snpseq(seq, value, start, end, vcfdict)
                o.write("##DOWNSTREAM##{down5k}".format(down5k=sequence))
                #o.write("##DOWNSTREAM##{down5k}".format(down5k=record_dict[value[0]][int(value[1]) - 1:int(value[1]) - 5001].seq.reverse_complement()))

def main(args):
    logfile = args.logfile

    logger = logging.getLogger(__name__)
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    handler = logging.FileHandler(logfile)
    handler.formatter = formatter
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.formatter = formatter
    logger.addHandler(handler)
    logger.addHandler(console_handler)
    logger.setLevel(level=logging.INFO)

    RetrieveSeq(args, logger)

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
This script is used for calc tag counts of all samples.
Usage:

python3 {scriptpath}/retrieve_geneseq_from_gtf_vcf_genome.py -gff3 <gff3> -vcf <vcf> -genome <genome.fa> -outdir <outdir>


'''.format(scriptpath = scriptpath))

    parse.add_argument('-gff3', '--gff3', required = True, dest = "gff3", help = "genome annotation gff3 file", type = str, nargs = '?')
    parse.add_argument('-vcf', '--vcf', required = True, dest = "vcf", help = "snp vcf file", type = str, nargs = '?')
    parse.add_argument('-genome', '--genome', required=True, dest="genome", help="genome fa file", type=str, nargs='?')
    parse.add_argument('-genenames', '--genenames', required=True, dest="genenames", help="comma delimited gene names", type=str, nargs='?')
    parse.add_argument('-outdir', '--outdir', required=True, dest="outdir", help="output dir", type=str, nargs='?')
    parse.add_argument('-logfile', '--logfile', required = True, dest = "logfile", help = "Log file to record procedures of processing of this script", metavar = "Log file to record procedures of processing of this script", type = str, nargs = '?')


    args = parse.parse_args()

    main(args)
