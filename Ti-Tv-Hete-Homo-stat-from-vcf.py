import re,os,sys,time
import argparse
import logging

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def Baoliuxiaoshu(num, digit):
    digitnum = round(num, digit)
    return str(digitnum)

def StatFromVcf(snpvcf, outstatfile, logger):
    purinelist = ['G', 'A']
    pyridinelist = ['C', 'T']
    genotypedict = {
        '0/1': "hete",
        '1/1': "homo"
    }
    transition = ""
    progencydict = dict()
    logger.info("Starting to process VCF...")
    with open(snpvcf, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                progencylist = line.strip().split("\t")[11:]
                for progency in progencylist:
                    progencydict[progency] = dict()
                    progencydict[progency]["TotalSNP"] = 0
                    progencydict[progency]["Tinum"] = 0
                    progencydict[progency]["Tvnum"] = 0
                    progencydict[progency]["TotalHete"] = 0
                    progencydict[progency]["TotalHomo"] = 0
            elif not line.startswith("#"):
                contentlist = line.strip().split("\t")
                ref = contentlist[3]
                alt = contentlist[4]
                if ref in purinelist and alt in pyridinelist:
                    transition = "Tv"
                elif ref in pyridinelist and alt in purinelist:
                    transition = "Tv"
                else:
                    transition = "Ti"
                progencygenolist = contentlist[11:]
                progencynum = len(progencygenolist)
                for i in range(progencynum):
                    genotype = progencylist[i].split(":")[0]
                    if genotype in genotypedict:
                        progencydict[progencylist[i]]["TotalSNP"] += 1
                        if genotype == "0/1":
                            progencydict[progencylist[i]]["TotalHete"] += 1
                        else:
                            progencydict[progencylist[i]]["TotalHomo"] += 1
                        if transition == "Tv":
                            progencydict[progencylist[i]]["Tvnum"] += 1
                        else:
                            progencydict[progencylist[i]]["Tinum"] += 1
    logger.info("Starting to write stat result ...")
    with open(outstatfile , 'r+') as o:
        for progency in progencydict:
            o.write(progency + "\t" + progencydict[progency]["TotalSNP"] + "\t" +  progencydict[progency]["Tinum"] + "\t" + progencydict[progency]["Tvnum"] + "\t" + Baoliuxiaoshu(progencydict[progency]["Tinum"]/progencydict[progency]["Tvnum"], 2) + progencydict[progency]["TotalHete"] + "\t" + progencydict[progency]["TotalHomo"] + "\n" )

def main(args):
    snpvcf = args.snpvcf
    outstatfile = args.statfile


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

    StatFromVcf(snpvcf, outstatfile, logger)
    logger.info("Finishing stat...")

if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:
python3 {scriptpath}/MergeIndelSNPVcf.py <snp.vcf> <indel.vcf> <Merged.vcf>
'''.format(scriptpath = scriptpath))

    parse.add_argument('-snpvcf', '--snpvcf', required = True, dest = "snpvcf", help = "filtered snp vcf contain only HQ SNPs", type = str, nargs = '?')
    parse.add_argument('-statfile', '--statfile', required=True, dest="statfile", help="file that stat result will be written in ", type=str, nargs='?')
    parse.add_argument('-logfile', '--logfile', required=True, dest="logfile",
                       help="Log file to record procedures of processing of this script",
                       metavar="Log file to record procedures of processing of this script", type=str, nargs='?')
    args = parse.parse_args()

    main(args)
