import yaml
import sys
import random
import pathlib
import Bio.SeqIO.QualityIO as QualityIO
import gzip
import collections
import subprocess
import re
config_file = sys.argv[1]
threads = int(sys.argv[2])
def getPairsWritten(file):
    forward = file.split('/')[-1].replace('cutadapt_', '').replace('.log', '').split('---')[0] #cutadapt_515r_parada___v4---926r___v5.log
    reverse = file.split('/')[-1].replace('cutadapt_', '').replace('.log', '').split('---')[1]
    forward_name = forward.split('___')[0]
    forward_region = forward.split('___')[1]
    reverse_name = reverse.split('___')[0]
    reverse_region = reverse.split('___')[1]

    with open(file) as handle:
        for line in handle:
            if 'Pairs written (passing filters):' in line:
                return int(line.split()[-2].replace(',', ''))

def getSequences(count, file):
    data = []
    with gzip.open(file, 'rt') as handle:
        for cnt, (header, sequence, qual) in enumerate(QualityIO.FastqGeneralIterator(handle), 1):
            if cnt <= count:
                data.append((header, sequence, qual))
            else:
                break
    return data

 
yaml_data = {}
with open(config_file) as handle:
    yaml_data = yaml.safe_load(handle)
  

data_dir = yaml_data['data_dir']
samples = yaml_data['sample_file']
blocklist = yaml_data['blocklist']
primer_file = yaml_data['primers']


with open(primer_file) as handle:
    
    primers = {}
    primer_pairs = []
    for line in handle:
        splits = line.strip().split()
        if line.startswith('#'):
            forward_primer = primers[splits[1]]
            reverse_primer = primers[splits[-1]]
            primer_pairs.append((forward_primer, reverse_primer))
        else:
            primers[splits[2]] = splits[1:]

def whereAreTheNs(file, required_read_length):
    positions = collections.Counter()
    sequences_with_N = 0
    with open(file) as handle:
        for cnt, (header, sequence, qual) in enumerate(QualityIO.FastqGeneralIterator(handle), 1):
            sequence = sequence[:required_read_length]
            if 'N' in sequence:
                sequences_with_N += 1
                for position in [_.start() for _ in re.finditer('N', sequence)]:
                    positions[position] += 1
    total_N = sum(list(positions.values()))
    return total_N, sequences_with_N, positions

def getReadlength(file):
    reads = []
    with open(file) as handle:
        for cnt, (header, sequence, qual) in enumerate(QualityIO.FastqGeneralIterator(handle), 1):
            reads.append(len(sequence))
    return reads

def detreads_and_bases(f1, f2):
    r1_reads = []
    r2_reads = []
    with open(f1) as handle:
        for cnt, (header, sequence, qual) in enumerate(QualityIO.FastqGeneralIterator(handle), 1):
            r1_reads.append(len(sequence))
    with open(f2) as handle:
        for cnt, (header, sequence, qual) in enumerate(QualityIO.FastqGeneralIterator(handle), 1):
            r2_reads.append(len(sequence))
    reads = len(r1_reads) * 2
    inserts = len(r1_reads)
    bases = sum(r1_reads) + sum(r2_reads)
    return reads, inserts, bases

SAMPLENAMES = set()
BLOCKLISTNAMES = set()
with open(blocklist) as handle:
    for line in handle:
        BLOCKLISTNAMES.add(line.strip())

with open(samples) as handle:
    for line in handle:
        samplename = line.strip()
        if samplename not in BLOCKLISTNAMES:
            SAMPLENAMES.add(samplename)
            

temp_folder = data_dir + '/temp_folder/'

pathlib.Path(temp_folder).mkdir(parents=True, exist_ok=True)
raw_r1_file = temp_folder + 'raw_R1.fastq'
raw_r2_file = temp_folder + 'raw_R2.fastq'
total_reads = 0
if True:
    print('subsetting read files')
    dest_r1_file = open(raw_r1_file, 'w')
    dest_r2_file = open(raw_r2_file, 'w')
    for sample in SAMPLENAMES:
        source_r1_file = data_dir + '/0raw/' + sample + '/' + sample + '_R1.fastq.gz'
        source_r2_file = data_dir + '/0raw/' + sample + '/' + sample + '_R2.fastq.gz'
        maxcount = 5000
        r1_reads = getSequences(maxcount, source_r1_file)
        r2_reads = getSequences(maxcount, source_r2_file)
        for (header, sequence, qual) in r1_reads:
            dest_r1_file.write(f'@{header}\n{sequence}\n+\n{qual}\n')
            total_reads += 1
        for (header, sequence, qual) in r2_reads:
            dest_r2_file.write(f'@{header}\n{sequence}\n+\n{qual}\n')

    dest_r1_file.close()
    dest_r2_file.close()
print('FORWARD_PRIMER\tFORWARD_PRIMER_NAME\tREVERSE_PRIMER\tREVERSE_PRIMER_NAME\tRAW_INSERTS\tAVERAGE_READ_LENGTH_R1\tAVERAGE_READ_LENGTH_R2\tCUTADAPT_INSERTS\tPRIMER_HITS\tR1_READS_W_N\tR2_READS_W_N\tPARAM=TRUNCLENR1\tPARAM=TRUNCLENR2\tPARAM=QCMINLEN\tPARAM=MAXEE\tQC_INSERTS\tPERC_INSERTS\tPLANNED_OVERLAP\tESTIMATED_INSERT_SIZE')
raw_reads, raw_inserts, raw_bases = detreads_and_bases(raw_r1_file, raw_r2_file)
raw_readlength_r1 = getReadlength(raw_r1_file)
raw_readlength_r2 = getReadlength(raw_r2_file)
raw_readlength_r1_avg = int(sum(raw_readlength_r1)/len(raw_readlength_r1))
raw_readlength_r2_avg = int(sum(raw_readlength_r2)/len(raw_readlength_r2))
if True:
    for forward_primer, reverse_primer in primer_pairs:
        cutadapt_r1_file = temp_folder + f'cutadapt_{forward_primer[1]}___{forward_primer[2]}---{reverse_primer[1]}___{reverse_primer[2]}_R1.fastq'
        cutadapt_r2_file = temp_folder + f'cutadapt_{forward_primer[1]}___{forward_primer[2]}---{reverse_primer[1]}___{reverse_primer[2]}_R2.fastq'
        cutadapt_log_file = temp_folder + f'cutadapt_{forward_primer[1]}___{forward_primer[2]}---{reverse_primer[1]}___{reverse_primer[2]}.log'
        
        cutadapt_command_stub = f'cutadapt -O 12 --discard-untrimmed -g {forward_primer[0]} -G {reverse_primer[0]} -o {cutadapt_r1_file} -p {cutadapt_r2_file} {raw_r1_file} {raw_r2_file} -j {threads} --pair-adapters --minimum-length 75 &> {cutadapt_log_file}'
        if True:
            subprocess.check_call(cutadapt_command_stub, shell=True)
        pairs_written = getPairsWritten(cutadapt_log_file)
        perc_kept = int(100.0 * pairs_written/total_reads)
        cutadapt_reads, cutadapt_inserts, cutadapt_bases = detreads_and_bases(cutadapt_r1_file, cutadapt_r2_file)
        if perc_kept < 10.0:
            print(f'{forward_primer[0]}\t{forward_primer[1]}\t{reverse_primer[0]}\t{reverse_primer[1]}\t{raw_inserts}\t{raw_readlength_r1_avg}\t{raw_readlength_r2_avg}\t{cutadapt_inserts}\t{perc_kept}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA')            
            #print(f'------------------Discarding primer pair {forward_primer} {reverse_primer}. <10% ({perc_kept}%) reads made it through cutadapt')
            continue
        else:
            #print(f'Keeping primer pair {forward_primer} {reverse_primer}. >10% ({perc_kept}%) reads made it through cutadapt')
            readlength_r1 = getReadlength(cutadapt_r1_file)
            readlength_r2 = getReadlength(cutadapt_r2_file)
            readlength_r1_avg = int(sum(readlength_r1)/len(readlength_r1))
            readlength_r2_avg = int(sum(readlength_r2)/len(readlength_r2))
            estimated_insert_size = int(reverse_primer[1].split('r')[0]) - int(forward_primer[1].split('f')[0]) - len(forward_primer[0]) - len(reverse_primer[0]) + 1
            planned_overlap = 45
            required_read_length = int((2 * planned_overlap + estimated_insert_size)/2) 
            r1_required_read_length = required_read_length + 20
            r2_required_read_length = required_read_length - 20
            qcminlength = r2_required_read_length - 10
            maxees = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
            r2_total_N, r2_sequences_with_N, r2_N_positions = whereAreTheNs(cutadapt_r2_file, r2_required_read_length)
            r1_total_N, r1_sequences_with_N, r1_N_positions = whereAreTheNs(cutadapt_r1_file, r1_required_read_length)
            for maxee in maxees:
                r1_out_file = cutadapt_r1_file.replace('_R1.fastq', '_maxee-' + str(maxee) + '_R1.fastq')
                r2_out_file = cutadapt_r1_file.replace('_R1.fastq', '_maxee-' + str(maxee) + '_R2.fastq')
                command = f'Rscript /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAB_ANALYSIS_PAN/code/pipeline/dada2_filterandtrim.R {cutadapt_r1_file} {cutadapt_r2_file} {r1_out_file} {r2_out_file} {maxee} 2 0 FALSE {qcminlength} {r1_required_read_length} {r2_required_read_length} {threads} &> {r1_out_file}.log'
                    
                if True:
                    subprocess.check_call(command, shell=True)
                    pathlib.Path(r1_out_file).touch(exist_ok=True)
                    pathlib.Path(r2_out_file).touch(exist_ok=True)
                reads, inserts, bases = detreads_and_bases(r1_out_file, r2_out_file)
                perc_qc_inserts = int(100.0 * inserts/cutadapt_inserts)
                print(f'{forward_primer[0]}\t{forward_primer[1]}\t{reverse_primer[0]}\t{reverse_primer[1]}\t{raw_inserts}\t{raw_readlength_r1_avg}\t{raw_readlength_r2_avg}\t{cutadapt_inserts}\t{perc_kept}\t{r1_sequences_with_N}\t{r2_sequences_with_N}\t{r1_required_read_length}\t{r2_required_read_length}\t{qcminlength}\t{maxee}\t{inserts}\t{perc_qc_inserts}\t{planned_overlap}\t{estimated_insert_size}')

                #print(maxee, reads, estimated_insert_size,  bases, cutadapt_reads, cutadapt_bases)

