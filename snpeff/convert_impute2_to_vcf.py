#!/usr/bin/env python

'''
Convert output from IMPUTE2 to VCF file format.
VCF file format description: http://samtools.github.io/hts-specs/VCFv4.2.pdf
'''

import sys

def get_args(cline = sys.argv[1:]):
    '''
    '''
    import argparse
    parser = argparse.ArgumentParser(description = '''Convert output from IMPUTE2 to
                                     VCF file format.''')
    parser.add_argument('impute_file', help = 'haplotype output file from IMPUTE2')
    parser.add_argument('-o', '--output', help = 'output VCF file [default: stdout]',
                        metavar = "fname.vcf")
    args = parser.parse_args(cline)
    return args

def add_header(fname, num_samples):
    '''Writes the proper VCF header to fname.

    Input
      fname: file handle
      num_samples: the number of samples with genotype information (integer)
    '''
    
    descrip = '''CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'''
    list_num = []
    for num in range(int(num_samples)):
        list_num.append("NA%0.5d"%(num + 1))

    list_genot = str.join('\t', list_num)
    descrip = descrip + list_genot
    
    header = '''##fileformat=VCFv4.2\n##source=IMPUTE2\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n##INFO=<ID=GC,Number=G,Type=Integer,Description="Genotype Counts">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#%s\n'''%(descrip)
    
    fname.write(header)

    
    
def convert_line(line, chr):
    '''Convert line of IMPUTE2 to VCF format.

    Input:
      line: a string of text
      chr: the chromosome (string)
    
    Example:
      convert_line('--- rs006 60006 T G 0 1 1 0 0 1 1 0 0 1', 'chr20')
      returns 'chr20	60006	rs006	T	G	PASS	.	NS=5;AF=0.5;GC=0,5,0	GT	0|1	1|0	0|1	1|0	0|1'
    '''
    
    line_info = line.split(" ")
    
    haplo_ref = 0
    haplo_alt = 0
    gen_HH = 0
    gen_Hh = 0
    gen_hh = 0
    SNP_genot = []

    for hap in range(5, len(line_info),2):

        if hap == len(line_info) - 2 :
            line_info[hap+1] = line_info[hap+1].split('\n')[0] 
        
        if line_info[hap] == '0':
            haplo_ref += 1
        elif line_info[hap] == '1':
            haplo_alt += 1
            
        if line_info[hap+1] == '0':
            haplo_ref += 1
        elif line_info[hap+1] == '1':
            haplo_alt += 1
        
        genot = line_info[hap] + '|'+ line_info[hap+1]
        if genot == '0|0':
            gen_HH += 1
        elif genot == '0|1' or genot == '1|0':
            gen_Hh += 1
        elif genot == '1|1':
            gen_hh += 1

        if genot == '.|.':
            genot = '.'

        SNP_genot.append(genot)

    genotype_number = gen_HH + gen_Hh + gen_hh
    all_freq = round((haplo_alt/float(haplo_ref + haplo_alt)), 4)
    GT =  str.join('\t', SNP_genot) 
    count = str.join(',', [str(gen_HH), str(gen_Hh), str(gen_hh)])
    INFO = str.join(';', ['NS=' + str(genotype_number), 'AF=' + str(all_freq), 'GC=' + count])
    line_info[1], line_info[2] = line_info[2], line_info[1]
    
    Nline = []
    Nline.append(chr)
    for p1 in range(1,5):
        Nline.append(line_info[p1])
        if line_info[3] == '-':
            line_info[3] = '-'+ line_info[4]
        if line_info[4] == '-':
            line_info[4] = '-' + line_info[3]
    
    nline = str.join('\t', Nline)
    new_line = str.join('\t', [nline, "PASS", '.', INFO, 'GT', GT])
    
    return new_line
    
    

def convert_impute2_to_vcf(impute_file, chrom, vcf_file = sys.stdout):
    '''Converts IMPUTE2 haplotype file to VCF format.

    Input:
      impute_file: file handle to read IMPUTE file
      vcf_file: file handle to write VCF output (default: sys.stdout)
    '''

    first_line = impute_file.readline()
    first_new_line = convert_line(first_line, chrom)
    genot_number = (first_new_line.split(';')[0]).split('=')[1]
    
    add_header(vcf_file, genot_number)
    vcf_file.write(first_new_line + '\n')
    
    while True:
        line = impute_file.readline()
        if not line:
            break
        vcf_file.write(convert_line(line, chrom) + '\n')

def get_chrom (file_name):

    chrom = (file_name.split('/')[-1]).split('.')[0]
    return chrom

def get_gzip(impute_file):
    
    import gzip 
    file_name = (impute_file.split('/')[-1]).split('.')[-1]
    if file_name == 'gz':
        content = gzip.open(impute_file, 'rb') 
    else:
        content = open(impute_file, 'rb')

    return content
    
#####################################################################################################
            
if __name__ == '__main__':
    args = get_args()
    chr_name = get_chrom(args.impute_file)
    impute_file = get_gzip(args.impute_file)
    if args.output:
        vcf_file = open(args.vcf_file, 'w')
        convert_impute2_to_vcf(impute_file,  chr_name, vcf_file)
    else:
        convert_impute2_to_vcf(impute_file,  chr_name)
