
#!/usr/bin/env python3
import re,os,sys,time
import argparse
import logging

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def Intervalfromvcf(args, invcf, outvcf, logger):
    chromlist  = args.chrom
    startslist = args.starts
    endslist = args.ends
    posnum  = len(startslist)
    o = open(outvcf, 'w+')
    pattern = re.compile(r"$")
    exonic = re.sub(pattern , ".exonic", outvcf)
    logger.info("Starting to read input vcf...")
    with open(invcf, 'r') as f:
        for line in f:
            if line.startswith("#CH"):
                o.write(line)
            else:
                contentlist = line.strip().split("\t")
                chrom = contentlist[0]
                for i in range(posnum):
                    if chrom == chromlist[i]:
                #if chrom in chromlist:
                        posit = int(contentlist[1])
                    #for i in range(posnum):
                        if posit >= startslist[i]  and posit <= endslist[i]:
                            o.write(line)
    o.close()
    logger.info("Intervals info are processed successfully...")
    os.system("grep -v 'Func.refGene=intergenic' {outvcf} > {exonic}".format(outvcf = outvcf, exonic = exonic))
def main(args):
    outvcf = args.outvcf
    invcf = args.invcf
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

    Intervalfromvcf(args, invcf, outvcf, logger)


if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:

python3 {scriptpath}/IntervalFromVcf.py <invcf> <outvcf> <startpos> <endpos> <logfile>

'''.format(scriptpath = scriptpath))

    parse.add_argument('-I', '--invcf', required = True, dest = "invcf", help = "strandard input vcf", metavar = "strandard input vcf", type = str, nargs = '?')
    parse.add_argument('-O', '--outvcf', required=True, dest="outvcf", help="out vcf with intervals",
                       metavar="out vcf with intervals", type=str, nargs='?')
    parse.add_argument('-C', '--chrom', required=True, dest="chrom", help="chrom number list",
                       metavar="chrom number list", type=str, nargs='+')
    parse.add_argument('-S', '--startpos', required=True, dest="starts", help="start pos list",
                       metavar="start pos list", type=int, nargs='+')
    parse.add_argument('-E', '--endpos', required=True, dest="ends", help="end pos list",
                       metavar="end pos list", type=int, nargs='+')
    parse.add_argument('-logfile', '--logfile', required=True, dest="logfile",
                       help="Log file to record procedures of processing of this script",
                       metavar="Log file to record procedures of processing of this script", type=str, nargs='?')
    args = parse.parse_args()

    main(args)
