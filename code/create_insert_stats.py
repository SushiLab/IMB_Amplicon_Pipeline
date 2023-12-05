import sys
import glob
import csv
import collections

# parse arguments
data_dir = sys.argv[1]
if len(sys.argv) > 2:
    out_file = sys.argv[2]
else:
    out_file = False

files = [('6bimeraRemoval', '6bimeraRemoval/seqtab.nobimera.stats'),
         ('4sampleInference_R1', '4sampleInference/R1.seqtab.stats'),
         ('4sampleInference_R2', '4sampleInference/R2.seqtab.stats'),
         ('5mergeReads', '5mergeReads/seqtab.mergereads.stats')]

steps = set()
sample_2_step_2_inserts = {}
otu_file = glob.glob(data_dir + '/*.otus.tsv')[0]
asv_file = glob.glob(data_dir + '/*.asvs.tsv')[0]

read_stats_files = glob.glob(data_dir + '/*/*/*readstats.log')
offset = -2
if len(read_stats_files) < 1:
    read_stats_files = glob.glob(data_dir + '/*/*readstats.log')
    offset = -1


def parse_bbduk_stats(logfile):
    """
    Parse BBDuk Statistics from Log File

    Parameters:
    - logfile (str): Path to the BBDuk log file containing statistics.

    Returns:
    - tuple: A tuple containing the input read count and retained read count.
             If unable to find statistics in the log file, it returns (0, 0).
    """
    with open(logfile) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('Input:'):
                splits = line.split()
                return int(splits[1]), int(splits[3])
        return 0, 0


def readStats(file_path):
    """
    Read Sample Insert Statistics

    Parameters:
    - file_path (str): Path to the file containing sample insert statistics.

    Returns:
    - dict: A dictionary containing sample IDs as keys and their corresponding insert counts.
    """
    sample_2_inserts = {}
    with open(file_path) as handle:
        for line in handle:
            sample_2_inserts[line.strip().split()[0]] = line.strip().split()[1]
    return sample_2_inserts


def readabfile(file_path, samples):
    """
    Read A/B File and Extract Sample Insert Counts

    Parameters:
    - file_path (str): Path to the file containing sample insert data.
    - samples (list): A list of sample names for which insert counts need to be extracted.

    Returns:
    - collections.Counter: A Counter object containing sample IDs as keys and their aggregated
                           insert counts as values.
    """
    with open(file_path) as handle:
        sample_2_inserts = collections.Counter()
        csvFile = csv.DictReader(handle, delimiter='\t')
        for line in csvFile:
            for sample in samples:
                sample_2_inserts[sample] += int(line[sample])
        return sample_2_inserts


# Iterating through read_stats_files list to extract sample information
for read_stats_file in read_stats_files:
    sample = read_stats_file.split('/')[offset].split('.')[0]
    step = read_stats_file.split('/')[offset - 1]
    reads, bases = parse_bbduk_stats(read_stats_file)
    inserts = int(reads / 2)
    steps.add(step)
    if sample not in sample_2_step_2_inserts:
        sample_2_step_2_inserts[sample] = {}
    sample_2_step_2_inserts[sample][step] = inserts

for (step, f) in files:
    steps.add(step)
    sample_2_inserts = readStats(data_dir + '/' + f)
    for sample, inserts in sample_2_inserts.items():
        sample_2_step_2_inserts[sample][step] = inserts

steps.add('8otu_counts')
steps.add('8asv_counts')

# Reading data from files (otu_file and asv_file) and updating inserts information
for sample, inserts in readabfile(otu_file, list(sample_2_step_2_inserts.keys())).items():
    sample_2_step_2_inserts[sample]['8otu_counts'] = inserts
for sample, inserts in readabfile(asv_file, list(sample_2_step_2_inserts.keys())).items():
    sample_2_step_2_inserts[sample]['8asv_counts'] = inserts

# Sorting steps and converting them to a list
steps = sorted(list(steps))

# Creating a header line for the output file
tmp = '\t'.join(['sample'] + steps)

# Write to output file
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
