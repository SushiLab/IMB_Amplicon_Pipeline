import sys
import glob
import csv
import collections

data_dir = sys.argv[1]
if len(sys.argv) > 2:
    out_file = sys.argv[2]
else:
    out_file = False


def parse_bbduk_stats(logfile):
    with open(logfile) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('Input:'):
                splits = line.split()
                return int(splits[1]), int(splits[3])
        return 0, 0

def readStats(f):
    sample_2_inserts = {}
    with open(f) as handle:
        for line in handle:
            sample_2_inserts[line.strip().split()[0]] = line.strip().split()[1]
    return sample_2_inserts

def readabfile(f, samples):
    with open(f) as handle:
        sample_2_inserts = collections.Counter()
        csvFile = csv.DictReader(handle, delimiter='\t')
        for line in csvFile:
            for sample in samples:
                sample_2_inserts[sample] += int(line[sample])
        return sample_2_inserts


steps = set()
sample_2_step_2_inserts = {}
read_stats_files = glob.glob(data_dir + '/*/*/*readstats.log')
offset = -2
if len(read_stats_files) < 1:
    read_stats_files = glob.glob(data_dir + '/*/*readstats.log')
    offset = -1
for read_stats_file in read_stats_files:
    sample = read_stats_file.split('/')[offset].split('.')[0]
    step = read_stats_file.split('/')[offset - 1]
    reads, bases = parse_bbduk_stats(read_stats_file)
    inserts = int(reads / 2)
    steps.add(step)
    if sample not in sample_2_step_2_inserts:
        sample_2_step_2_inserts[sample] = {}
    sample_2_step_2_inserts[sample][step] = inserts

files = [('6bimeraRemoval', '6bimeraRemoval/seqtab.nobimera.stats'),
         ('4sampleInference_R1', '4sampleInference/R1.seqtab.stats'),
         ('4sampleInference_R2', '4sampleInference/R2.seqtab.stats'),
         ('5mergeReads', '5mergeReads/seqtab.mergereads.stats')]

for (step, f) in files:
    steps.add(step)
    sample_2_inserts = readStats(data_dir + '/' + f)
    for sample, inserts in sample_2_inserts.items():
        sample_2_step_2_inserts[sample][step] = inserts




otu_file = glob.glob(data_dir + '/*.otus.tsv')[0]
asv_file = glob.glob(data_dir + '/*.asvs.tsv')[0]
steps.add('8otu_counts')
steps.add('8asv_counts')

for sample, inserts in readabfile(otu_file, list(sample_2_step_2_inserts.keys())).items():
    sample_2_step_2_inserts[sample]['8otu_counts'] = inserts
for sample, inserts in readabfile(asv_file, list(sample_2_step_2_inserts.keys())).items():
    sample_2_step_2_inserts[sample]['8asv_counts'] = inserts

steps = sorted(list(steps))

tmp = '\t'.join(['sample'] + steps)

if out_file:
    sys.stdout = open(out_file, 'w')

sys.stdout.write(f'{tmp}\n')

for sample in sorted(sample_2_step_2_inserts.keys()):
    step_2_inserts = sample_2_step_2_inserts[sample]
    tmp = [sample]
    for step in steps:
        tmp.append(str(step_2_inserts[step]))
    tmp_str = '\t'.join(tmp)
    sys.stdout.write(f'{tmp_str}\n')
