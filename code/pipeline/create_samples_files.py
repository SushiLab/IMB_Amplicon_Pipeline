import sys
import glob

data_dir = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
samples = sys.argv[4].split(',')


def dofwrev(samples, outfile, ori, data_dir):
    of = open(outfile, 'w')
    sample_2_file = {}
    for trimmed_file in glob.glob(data_dir + '/2filterAndTrim/' + '*/*' + ori):
        sample = trimmed_file.split('/')[-2]
        sample_2_file[sample] = trimmed_file
    for sample in sorted(sample_2_file.keys()):
        of.write(f'{sample}\t{sample_2_file[sample]}\n')
    of.close()

dofwrev(samples, r1_file, 'R1.fastq.gz', data_dir)
dofwrev(samples, r2_file, 'R2.fastq.gz', data_dir)

