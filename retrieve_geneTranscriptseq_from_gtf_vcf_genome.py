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
import string

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def replaceseq(sequence, replacedict):
    tmplist = []
    seq = [i for i in sequence]
    for i in seq:
        if i in replacedict:
            tmplist.append(replacedict[i])
        else:
            tmplist.append(i)
    tmp = "".join(tmplist)
    return tmp

def snpseq(seq, value, start, end,  vcfdict, strand, snpdict):
    #maptrans = str.maketrans('ATCG', 'TAGC')
    maptrans = string.maketrans('ATCG', 'TAGC')
    seq = [i for i in seq]
    if value[0] in vcfdict:
        for posit in range(start, end):
            if str(posit) in vcfdict[value[0]]:
                varia = vcfdict[value[0]][str(posit)]
                if strand == "+":
                    if len(str(posit))%3 == 0:
                        seq[posit - start] = "[" + str(posit) + "::" + snpdict[varia]
                    elif len(str(posit))%3 == 1:
                        seq[posit - start] = "[" + str(posit) + ":" + snpdict[varia]
                    else:
                        seq[posit - start] = "[" + str(posit) + ":::" + snpdict[varia]
                else:
                    varia = varia.translate(maptrans)
                    if len(str(posit))%3 == 0:
                        seq[end - posit - 1] = "[" + str(posit) + "::" + snpdict[varia]
                    elif len(str(posit))%3 == 1:
                        seq[end - posit - 1] = "[" + str(posit) + ":" + snpdict[varia]
                    else:
                        seq[end - posit - 1] = "[" + str(posit) + ":::" + snpdict[varia]
                #seq[posit - start] = vcfdict[value[0]][str(posit)]
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
    mRNA = dict()
    mRNAdict = dict()
    cds = dict()
    cdslist = list()
    intron = dict()
    utr3 = dict()
    utr3list = list()
    utr5 = dict()
    utr5list = list()
    vcfdict = dict()
    replace = {
                "B":"A/T]", "D":"A/C]", "E":"A/G]",
                "F":"T/A]", "H":"T/C]", "I":"T/G]",
                "J":"C/A]", "K":"C/T]", "L":"C/G]",
                "M":"G/A]", "N":"G/T]", "O":"G/C]"
    }
    snpdict = {
               "AT":"B", "AC":"D", "AG":"E",
               "TA":"F", "TC":"H", "TG":"I",
               "CA":"J", "CT":"K", "CG":"L",
               "GA":"M", "GT":"N", "GC":"O"
               }

    logger.info("Processing genome fa file...")
    record_dict = SeqIO.index(genome, "fasta")
    logger.info("Processing gff3 file...")
    flag = 0
    dictname = ''
    mRNAname = ''
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
                        vcfdict[chrom][posit] = varia
                        #vcfdict[chrom][posit] = snpdict[varia]
                    else:
                        vcfdict[chrom][posit] = varia
                        #vcfdict[chrom][posit] = snpdict[varia]

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
                            dictname = genename
                            gene[dictname] = (chr, start, end, direction)
                            mRNA[dictname] = []
                            utr5[dictname] = dict()
                            utr3[dictname] = dict()
                            cds[dictname] = dict()
                            flag = 1
                            break
                        else:
                            flag = 0
                elif feature == "mRNA" and flag ==1:
                    mRNAname = name
                    mRNA[dictname].append(mRNAname)
                    mRNAdict[mRNAname] = [chr, start, end, direction]
                    #mRNAdict[mRNAname]   
                    utr5[dictname][mRNAname] = []
                    utr3[dictname][mRNAname] = []
                    cds[dictname][mRNAname] = []
                elif feature == "three_prime_UTR" and flag == 1:
                    utr3[dictname][mRNAname].append([chr, start, end])
                elif feature == "five_prime_UTR" and flag == 1:
                    utr5[dictname][mRNAname].append([chr, start, end])
                elif feature == "CDS" and flag == 1:
                    cds[dictname][mRNAname].append([chr, start, end])
    tmplist = list()
    for key, value in gene.items():
        outfile = outdir + "/" + key + ".fa"
        with open(outfile, 'w') as o:
            #o.write(">{key}\n".format(key = key))
            if value[3] == "+":
                strand = "+"
                #seq = record_dict[value[0]][int(value[1]) - 5001:int(value[1]) - 1].seq
                #start = int(value[1]) - 5001
                #end = int(value[1]) - 1
                #sequence = snpseq(seq, value, start, end,  vcfdict)
                #o.write("##UPSTREAM##{up5k}".format(up5k = sequence))
                for mRNAname in mRNA[key]:
                    o.write(">{mRNAname}\n".format(mRNAname = mRNAname))
                    seq = record_dict[value[0]][int(mRNAdict[mRNAname][1]) -5001:int(mRNAdict[mRNAname][1])-1].seq
                    start = int(mRNAdict[mRNAname][1]) -5001
                    end = int(mRNAdict[mRNAname][1])-1 
                    sequence = snpseq(seq, value, start, end,  vcfdict, strand, snpdict)
                    sequence = replaceseq(sequence, replace)
                    o.write("##UPSTREAM##{up5k}\n".format(up5k = sequence))
                    if utr5[key][mRNAname]:
                        print(utr5[key][mRNAname])
                        if len(utr5[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(utr5[key][mRNAname][0][1]) -1:int(utr5[key][mRNAname][0][2])].seq
                            start = int(utr5[key][mRNAname][0][1])-1
                            end = int(utr5[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            sequence = replaceseq(sequence, replace)
                            o.write("##UTR5##{utr5}\n".format(utr5 = sequence))
                        #o.write("##UTR5##{utr5}".format(utr5 =  record_dict[value[0]][int(utr5[key][1]):int(utr5[key[2]])].seq))
                        else:
                            for anyutr5 in sorted(utr5[key][mRNAname]):
                                tmplist.append(anyutr5[1])
                                tmplist.append(anyutr5[2])
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                o.write("##UTR5##{utr5}\n".format(utr5 = sequence))
                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 1]):int(tmplist[suffix + 2])].seq
                                    start = int(tmplist[suffix + 1])
                                    end = int(tmplist[suffix + 2])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron = sequence))
                            tmplist = list()

                    if cds[key][mRNAname]:
                        if len(cds[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(cds[key][mRNAname][0][1])-1:int(cds[key][mRNAname][0][2])].seq
                            start = int(cds[key][mRNAname][0][1])-1
                            end = int(cds[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            seqtmp = []
                            if len(sequence) % 3 == 0:
                                for i in range(0, len(sequence), 3):
                                    seqtmp.append(sequence[i:i+3])
                            elif len(sequence) % 3 == 1:
                                for i in range(0, len(sequence) - 1, 3):
                                    seqtmp.append(sequence[i:i+3])
                                seqtmp.append(sequence[-1])
                            elif len(sequence) % 3 == 2:
                                for i in range(0, len(sequence) - 2, 3):
                                    seqtmp.append(sequence[i:i+3])
                                seqtmp.append(sequence[-2:])
                            sequence = ' '.join(seqtmp)
                            sequence = replaceseq(sequence, replace)
                            o.write("##CDS##{cds}\n".format(cds = sequence))
                        else:
                            for anycds in sorted(cds[key][mRNAname]):
                                tmplist.append(anycds[1])
                                tmplist.append(anycds[2])
                            remainder = 3
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                seqtmp = []
                                if remainder == 3:
                                    if len(sequence) % 3 == 0:
                                        for i in range(0, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif len(sequence) % 3 == 1:
                                        for i in range(0, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-1])
                                    elif len(sequence) % 3 == 2:
                                        for i in range(0, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                elif remainder == 2:
                                    seqtmp.append(sequence[0:2])
                                    if (len(sequence)-2) % 3 == 0:
                                        for i in range(2, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-2) % 3 == 1:
                                        for i in range(2, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-1])
                                    elif (len(sequence)-2) % 3 == 2:
                                        for i in range(2, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                elif remainder == 1:
                                    seqtmp.append(sequence[0:1])
                                    if (len(sequence)-1) % 3 == 0:
                                        for i in range(1, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-1) % 3 == 1:
                                        for i in range(1, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-1) % 3 == 2:
                                        for i in range(1, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                remainder = 3-(len(sequence)-remainder)%3
                                sequence = ' '.join(seqtmp)
                                sequence = replaceseq(sequence, replace)
                                o.write("##CDS##{cds}\n".format(cds = sequence))
                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 1]):int(tmplist[suffix + 2])].seq
                                    start = int(tmplist[suffix + 1])
                                    end = int(tmplist[suffix + 2])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron = sequence))
                            tmplist = list()
                    if utr3[key][mRNAname]:
                        if len(utr3[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(utr3[key][mRNAname][0][1])-1:int(utr3[key][mRNAname][0][2])].seq
                            start = int(utr3[key][mRNAname][0][1])-1
                            end = int(utr3[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            sequence = replaceseq(sequence, replace)
                            o.write("##UTR3##{utr3}\n".format(utr3=sequence))
                        else:
                            for anyutr3 in sorted(utr3[key][mRNAname]):
                                tmplist.append(anyutr3[1])
                                tmplist.append(anyutr3[2])
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                sequence = replaceseq(sequence, replace)
                                o.write("##UTR3##{utr3}\n".format(utr3=sequence))
                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 1]):int(tmplist[suffix + 2])].seq
                                    #start = int(tmplist[suffix + 1]) + 1
                                    #end = int(tmplist[suffix + 2]) - 1
                                    start = int(tmplist[suffix + 1])
                                    end = int(tmplist[suffix + 2])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron=sequence))
                            tmplist = list()
                    seq = record_dict[value[0]][int(mRNAdict[mRNAname][2]) + 1:int(mRNAdict[mRNAname][2])+5001].seq
                    start = int(mRNAdict[mRNAname][2]) + 1
                    end = int(mRNAdict[mRNAname][2])+5001
                    sequence = snpseq(seq, value, start, end,  vcfdict, strand, snpdict)
                    sequence = replaceseq(sequence, replace)
                    o.write("##DOWNSTREAM##{down5k}\n--------------------------\n".format(down5k=sequence))
                    #o.write("##UPSTREAM##{up5k}\n--------------------------\n".format(up5k = sequence))
                #seq = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5001].seq
                #start = int(value[2]) + 1
                #end = int(value[2]) + 5001
                #sequence = snpseq(seq, value, start, end, vcfdict)
                #sequence = replaceseq(sequence, replace)
                #o.write("##DOWNSTREAM##{down5k}".format(down5k=sequence))

            else:
                #seq = record_dict[value[0]][int(value[2]) + 1:int(value[2]) + 5001].seq.reverse_complement()
                #start = int(value[2]) + 1
                #end = int(value[2]) + 5001
                #sequence = snpseq(seq, value, start, end, vcfdict)
                #sequence = replaceseq(sequence, replace)
                #o.write("##UPSTREAM##{up5k}".format(up5k=sequence))
                strand = "-"
                for mRNAname in mRNA[key]:
                    o.write(">{mRNAname}\n".format(mRNAname = mRNAname))
                    seq = record_dict[value[0]][int(mRNAdict[mRNAname][2]) + 1:int(mRNAdict[mRNAname][2]) + 5001].seq.reverse_complement()
                    start = int(mRNAdict[mRNAname][2]) + 1
                    end = int(mRNAdict[mRNAname][2]) + 5001
                    sequence = snpseq(seq, value, start, end,  vcfdict, strand, snpdict)
                    sequence = replaceseq(sequence, replace)
                    o.write("##UPSTREAM##{up5k}\n".format(up5k = sequence))
                    if utr5[key][mRNAname]:
                        if len(utr5[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(utr5[key][mRNAname][0][1])-1:int(utr5[key][mRNAname][0][2])].seq.reverse_complement()
                            start = int(utr5[key][mRNAname][0][1])-1
                            end = int(utr5[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            sequence = replaceseq(sequence, replace)
                            o.write("##UTR5##{utr5}\n".format(utr5=sequence))
                        else:
                            for anyutr5 in sorted(utr5[key][mRNAname], reverse = True):
                                tmplist.append(anyutr5[1])
                                tmplist.append(anyutr5[2])
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq.reverse_complement()
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                sequence = replaceseq(sequence, replace)
                                o.write("##UTR5##{utr5}\n".format(utr5=sequence))
                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 3]):int(tmplist[suffix])].seq.reverse_complement()
                                    start = int(tmplist[suffix + 3])
                                    end = int(tmplist[suffix])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron=sequence))
                            tmplist = list()
                    if cds[key][mRNAname]:
                        if len(cds[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(cds[key][mRNAname][0][1])-1:int(cds[key][mRNAname][0][2])].seq.reverse_complement()
                            start = int(cds[key][mRNAname][0][1])-1
                            end = int(cds[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            seqtmp = []
                            if len(sequence) % 3 == 0:
                                for i in range(0, len(sequence), 3):
                                    seqtmp.append(sequence[i:i+3])
                            elif len(sequence) % 3 == 1:
                                for i in range(0, len(sequence) - 1, 3):
                                    seqtmp.append(sequence[i:i+3])
                                seqtmp.append(sequence[-1])
                            elif len(sequence) % 3 == 2:
                                for i in range(0, len(sequence) - 2, 3):
                                    seqtmp.append(sequence[i:i+3])
                                seqtmp.append(sequence[-2:])
                            sequence = ' '.join(seqtmp)
                            sequence = replaceseq(sequence, replace)
                            o.write("##CDS##{cds}\n".format(cds=sequence))
                        else:
                            for anycds in sorted(cds[key][mRNAname], reverse = True):
                                tmplist.append(anycds[1])
                                tmplist.append(anycds[2])
                            remainder = 3
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq.reverse_complement()
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                seqtmp = []
                                if remainder == 3:
                                    if len(sequence) % 3 == 0:
                                        for i in range(0, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif len(sequence) % 3 == 1:
                                        for i in range(0, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-1])
                                    elif len(sequence) % 3 == 2:
                                        for i in range(0, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                elif remainder == 2:
                                    seqtmp.append(sequence[0:2])
                                    if (len(sequence)-2) % 3 == 0:
                                        for i in range(2, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-2) % 3 == 1:
                                        for i in range(2, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-1])
                                    elif (len(sequence)-2) % 3 == 2:
                                        for i in range(2, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                elif remainder == 1:
                                    seqtmp.append(sequence[0:1])
                                    if (len(sequence)-1) % 3 == 0:
                                        for i in range(1, len(sequence), 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-1) % 3 == 1:
                                        for i in range(1, len(sequence) - 1, 3):
                                            seqtmp.append(sequence[i:i+3])
                                    elif (len(sequence)-1) % 3 == 2:
                                        for i in range(1, len(sequence) - 2, 3):
                                            seqtmp.append(sequence[i:i+3])
                                        seqtmp.append(sequence[-2:])
                                remainder = 3-(len(sequence)-remainder)%3
                                sequence = ' '.join(seqtmp)
                                sequence = replaceseq(sequence, replace)
                                o.write("##CDS##{cds}\n".format(cds=sequence))

                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 3]):int(tmplist[suffix])].seq.reverse_complement()
                                    start = int(tmplist[suffix + 3])
                                    end = int(tmplist[suffix])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron=sequence))
                            tmplist = list()
                    if utr3[key][mRNAname]:
                        if len(utr3[key][mRNAname]) == 1:
                            seq = record_dict[value[0]][int(utr3[key][mRNAname][0][1])-1:int(utr3[key][mRNAname][0][2])].seq.reverse_complement()
                            start = int(utr3[key][mRNAname][0][1])-1
                            end = int(utr3[key][mRNAname][0][2])
                            sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                            sequence = replaceseq(sequence, replace)
                            o.write("##UTR3##{utr3}\n".format(utr3=sequence))
                        else:
                            for anyutr3 in sorted(utr3[key][mRNAname], reverse = True):
                                tmplist.append(anyutr3[1])
                                tmplist.append(anyutr3[2])
                            for suffix in range(0,len(tmplist),2):
                                seq = record_dict[value[0]][int(tmplist[suffix])-1:int(tmplist[suffix + 1])].seq.reverse_complement()
                                start = int(tmplist[suffix])-1
                                end = int(tmplist[suffix + 1])
                                sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                sequence = replaceseq(sequence, replace)
                                o.write("##UTR3##{utr3}\n".format(utr3=sequence))
                                if suffix + 3 < len(tmplist):
                                    seq = record_dict[value[0]][int(tmplist[suffix + 3]):int(tmplist[suffix])].seq.reverse_complement()
                                    start = int(tmplist[suffix + 3])
                                    end = int(tmplist[suffix])
                                    sequence = snpseq(seq, value, start, end, vcfdict, strand, snpdict)
                                    sequence = replaceseq(sequence, replace)
                                    o.write("##INTRON##{intron}\n".format(intron=sequence))
                            tmplist = list()
                    seq = record_dict[value[0]][int(mRNAdict[mRNAname][1]) - 5002:int(mRNAdict[mRNAname][1])-1].seq.reverse_complement()
                    start = int(mRNAdict[mRNAname][1]) - 5002
                    end = int(mRNAdict[mRNAname][1])-1
                    sequence = snpseq(seq, value, start, end,  vcfdict, strand, snpdict)
                    sequence = replaceseq(sequence, replace)
                    o.write("##DOWNSTREAM##{down5k}\n--------------------------\n".format(down5k=sequence))
                    #o.write("##UPSTREAM##{up5k}\n--------------------------\n".format(up5k = sequence))
                #seq = record_dict[value[0]][int(value[1]) - 5002:int(value[1]) - 1].seq.reverse_complement()
                #start = int(value[1]) - 5002
                #end = int(value[1]) - 1
                #sequence = snpseq(seq, value, start, end, vcfdict)
                #sequence = replaceseq(sequence, replace)
                #o.write("##DOWNSTREAM##{down5k}".format(down5k=sequence))

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

