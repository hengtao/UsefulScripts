
import re,os,sys,time
import argparse
import logging

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def Baoliuxiaoshu(num, digit):
    digitnum = round(num, digit)
    return str(digitnum)

def StatSNPIndel(snpfile, snpoutstat, logger):
    typedict = dict()
    genotypedict = {
        '0/1': "hete",
        '1/1': "homo"
    }
    progencylist = list()
    progencydict = dict()
    groupnum = 0

    with open(snpfile) as f:
        for line in f:
            if line.startswith("#CHROM"):
                progencylist = line.strip().split("\t")[9:]
                groupnum = len(progencylist)
                for progency in progencylist:
                    progencydict[progency] = dict()
                    progencydict[progency]["intergenic"] = 0
                    progencydict[progency]["splicing"] = 0
                    progencydict[progency]["intron"] = 0
                    progencydict[progency]["upstream"] = 0
                    progencydict[progency]["downstream"] = 0
                    progencydict[progency]["updown"] = 0
                    progencydict[progency]["intergenic"] = 0
                    progencydict[progency]["synonymous"] = 0
                    progencydict[progency]["missense"] = 0
                    progencydict[progency]["stopgain"] = 0
                    progencydict[progency]["stoploss"] = 0
                    progencydict[progency]["UTR5"] = 0
                    progencydict[progency]["UTR3"] = 0
                    progencydict[progency]["frameshiftdeletion"] = 0
                    progencydict[progency]["frameshiftinsertion"] = 0
                    progencydict[progency]["nonframeshiftdeletion"] = 0
                    progencydict[progency]["nonframeshiftinsertion"] = 0
                    progencydict[progency]["Other"] = 0
            elif not line.startswith("#"):
                contentlist = line.strip().split("\t")
                progencygenolist = contentlist[9:]
                infolist = line.strip().split("\t")[7].split(";")
                FuncRefgene = infolist[17].replace("Func.refGene=", "")
                ExonicFuncRefGene = infolist[20].replace("ExonicFunc.refGene=", "")
                for i in range(groupnum):
                    genotype = progencygenolist[i].split(":")[0]
                    if genotype in genotypedict:
                        if FuncRefgene.endswith("splicing"):
                            progencydict[progencylist[i]]["splicing"] += 1
                        elif FuncRefgene.endswith("intronic"):
                            progencydict[progencylist[i]]["intron"] += 1
                        elif FuncRefgene.endswith("intergenic"):
                            progencydict[progencylist[i]]["intergenic"] += 1
                        elif FuncRefgene.startswith("upstream"):
                            if FuncRefgene.endswith("downstream"):
                                progencydict[progencylist[i]]["updown"] += 1
                            else:
                                progencydict[progencylist[i]]["upstream"] += 1
                        elif FuncRefgene.startswith("downstream"):
                            progencydict[progencylist[i]]["downstream"] += 1
                        elif FuncRefgene.startswith("UTR5"):
                            progencydict[progencylist[i]]["UTR5"] += 1
                        elif FuncRefgene.startswith("UTR3"):
                            progencydict[progencylist[i]]["UTR3"] += 1
                        elif FuncRefgene.startswith("exonic") and not FuncRefgene.endswith("splicing"):
                            if ExonicFuncRefGene.startswith("synonymous_SNV"):
                                progencydict[progencylist[i]]["synonymous"] += 1
                            elif ExonicFuncRefGene.startswith("nonsynonymous_SNV"):
                                progencydict[progencylist[i]]["missense"] += 1
                            elif ExonicFuncRefGene.startswith("stopgain"):
                                progencydict[progencylist[i]]["stopgain"] += 1
                            elif ExonicFuncRefGene.startswith("stoploss"):
                                progencydict[progencylist[i]]["stoploss"] += 1
                            elif ExonicFuncRefGene.startswith("frameshift_deletion"):
                                progencydict[progencylist[i]]["frameshiftdeletion"] += 1
                            elif ExonicFuncRefGene.startswith("frameshift_insertion"):
                                progencydict[progencylist[i]]["frameshiftinsertion"] += 1
                            elif ExonicFuncRefGene.startswith("nonframeshift_deletion"):
                                progencydict[progencylist[i]]["nonframeshiftdeletion"] += 1
                            elif ExonicFuncRefGene.startswith("nonframeshift_insertion"):
                                progencydict[progencylist[i]]["nonframeshiftinsertion"] += 1
                        else:
                            progencydict[progencylist[i]]["Other"] += 1
    logger.info("Starting to write stat result ...")
    with open(snpoutstat, 'w+') as o:
        o.write("SampleID\tintergenic\tintron\tupstream\tdownstream\tupdown\tsplicing\tUTR5\tUTR3\tmissense\tsynonymous\tstopgain\tstoploss\tframeshiftdeletion\tframeshiftinsertion\tnonframeshiftdeletion\tnonframeshiftinsertion\tOther\n")
        for progency in progencydict:
            o.write(progency + "\t" + str(progencydict[progency]["intergenic"]) + "\t" + str(progencydict[progency]["intron"]) + "\t" + str(progencydict[progency]["upstream"]) + "\t" + str(progencydict[progency]["downstream"])  + "\t" + str(progencydict[progency]["updown"]) + "\t" + str(progencydict[progency]["splicing"]) +
                    "\t" + str(progencydict[progency]["UTR5"]) + "\t" + str(progencydict[progency]["UTR3"]) + "\t" + str(progencydict[progency]["missense"]) + "\t" + str(progencydict[progency]["synonymous"]) + "\t" + str(progencydict[progency]["stopgain"]) + "\t" + str(progencydict[progency]["stoploss"]) + "\t" +
                    str(progencydict[progency]["frameshiftdeletion"]) + "\t" + str(progencydict[progency]["frameshiftinsertion"]) + "\t" + str(progencydict[progency]["nonframeshiftdeletion"]) + "\t" + str(progencydict[progency]["nonframeshiftinsertion"]) + "\t" + str(progencydict[progency]["Other"]) + "\n")

def main(args):
    snpvcf = args.snpvcf
    indelvcf = args.indelvcf
    outstat = args.statfile
    logfile = args.logfile

    ## Import logger, 获取logger实例
    logger    = logging.getLogger(__name__)
    ## 指定logger输出格式
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    #logger.setLevel(level = logging.INFO)
    handler = logging.FileHandler(logfile)
    handler.formatter = formatter
    ## 控制台日志
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.formatter = formatter
    ## 为logger添加日志处理器
    logger.addHandler(handler)
    logger.addHandler(console_handler)
    ## 指定日志的最低输出级别，默认为WARN
    logger.setLevel(level = logging.INFO)

    logger.info("Starting to process SNP VCF...")
    StatSNPIndel(snpvcf, outstat, logger)
    logger.info("Starting to process Indel VCF...")
    StatSNPIndel(indelvcf, outstat, logger)
    logger.info("Finishing stat...")

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/MergeIndelSNPVcf.py <snp.vcf> <indel.vcf> <Merged.vcf>
'''.format(scriptpath = scriptpath))

    parse.add_argument('-snpvcf', '--snpvcf', required = True, dest = "snpvcf", help = "filtered snp vcf contain only HQ SNPs", type = str, nargs = '?')
    parse.add_argument('-indelvcf', '--indelvcf', required = True, dest = "indelvcf", help = "filtered snp vcf contain only HQ SNPs", type = str, nargs = '?')
    parse.add_argument('-statfile', '--statfile', required=True, dest="statfile", help="file that stat result will be written in ", type=str, nargs='?')
    parse.add_argument('-logfile', '--logfile', required=True, dest="logfile",
                       help="Log file to record procedures of processing of this script",
                       metavar="Log file to record procedures of processing of this script", type=str, nargs='?')
    args = parse.parse_args()

    main(args)
