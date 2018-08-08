import re,sub,sys,time
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def Merge(snp, indel, merge):
    chromlist = list()
    positdict = dict()
    o = open(merge, 'w+')
    with open(snp, 'r') as f:
        for line in f:
            if line.startswith("#"):
                o.write(line)
            else:
                contentlist = line.strip().split("\t")
                chrom = contentlist[0]
                posit = int(contentlist[1])
                if chrom not in chromlist:
                    chromlist.append(chrom)
                    positdict[chrom] = dict()
                    positdict[chrom][posit] = line
                else:
                    positdict[chrom][posit] = line
    with open(indel, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                contentlist = line.strip().split("\t")
                chrom = contentlist[0]
                posit = int(contentlist[1])
                if chrom not in chromlist:
                    chromlist.append(chrom)
                    positdict[chrom] = dict()
                    positdict[chrom][posit] = line
                else:
                    positdict[chrom][posit] = line
    for chrom in chromlist:
        for key in sorted(positdict[chrom].keys()):
            o.write(positdict[chrom][key])

def main(args):
    snpvcf = args.snpvcf
    indelvcf = args.indelvcf
    mergevcf = args.mergevcf

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

    Merge(snpvcf, indelvcf, mergevcf)


if __name__ == "__main__":
    scriptpath = os.path.split(os.path.realpath(__file__))[0]
    parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:

python3 {scriptpath}/MergeIndelSNPVcf.py <snp.vcf> <indel.vcf> <Merged.vcf>

'''.format(scriptpath = scriptpath))

    parse.add_argument('-snpvcf', '--snpvcf', required = True, dest = "snpvcf", help = "pos sorted snp vcf to be merged", metavar = "snp vcf to be merged", type = str, nargs = '?')
    parse.add_argument('-indelvcf', '--indelvcf', required=True, dest="indelvcf", help="pos sorted indel vcf to be merged",
                       metavar="indel vcf to be merged", type=str, nargs='?')
    parse.add_argument('-mergevcf', '--mergevcf', required=True, dest="mergevcf", help="merged vcf file,also sorted in pos",
                       metavar="merged vcf file,also sorted in pos", type=str, nargs='?')
    parse.add_argument('-logfile', '--logfile', required=True, dest="logfile",
                       help="Log file to record procedures of processing of this script",
                       metavar="Log file to record procedures of processing of this script", type=str, nargs='?')
    args = parse.parse_args()

    main(args)
