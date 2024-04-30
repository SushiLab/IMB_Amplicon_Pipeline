from pathlib import Path

'''
FILES + FOLDERS + PARAMS
'''
DATA_DIR = Path(config['data_dir'])
SAMPLE_FILE = Path(config['sample_file'])
BLOCKLIST_FILE = Path(config['blocklist'])
PRIMERS_FILE = Path(config['primers'])
SILVA_TRAINING_FILE = Path(config['silva_training'])
FORWARD_PRIMER_NAME = config['FORWARD_PRIMER_NAME']
REVERSE_PRIMER_NAME = config['REVERSE_PRIMER_NAME']
QC_MINLEN = config['QC_MINLEN']
QC_TRUNC_R1 = config['QC_TRUNC_R1']
QC_TRUNC_R2 = config['QC_TRUNC_R2']
QC_MAXEE = config['QC_MAXEE']
LE_NBASES = config['LE_NBASES']
USEARCH_DB = config['USEARCH_DB']
SCRIPTFOLDER = workflow.basedir + "/"

runCutadapt = config['runCutadapt']
allowUntrimmed = config['allowUntrimmed']
runQC = config['runQC']
runLearnErrors = config['runLearnErrors']
runInference = config['runInference']
runMergeReads = config['runMergeReads']
runRemoveBimeras = config['runRemoveBimeras']
runReadStats = config['runReadStats']
runASVTax = config['runASVTax']
runOTUTax = config['runOTUTax']
runUSEARCH = config['runUSEARCH']
runDefCom = config['runDefCom']


def get_primers(primer_file):
    """
    Reads a primer file and extracts primer sequences to create a dictionary of primer pairs.

    Args:
    - primer_file (str): The path to the file containing primer information.

    Returns:
    - primers (dict): A dictionary containing primer pairs. Keys are primer IDs, and values are lists of primer sequences.
                      The format is {primer_id: [forward_primer_sequence, reverse_primer_sequence]}.
    """
    primers = {}
    with open(primer_file) as handle:
        for line in handle:
            splits = line.strip().split()
            if line.startswith('#'):
                forward_primer = primers[splits[1]]
                reverse_primer = primers[splits[-1]]
            #primer_pairs.append((forward_primer, reverse_primer))
            else:
                primers[splits[2]] = splits[1:]
        return primers

# Get primer sequences
FORWARD_PRIMER_SEQUENCE = config.get('FORWARD_PRIMER_SEQUENCE')
if FORWARD_PRIMER_SEQUENCE is None:
    FORWARD_PRIMER_SEQUENCE = get_primers(PRIMERS_FILE)[FORWARD_PRIMER_NAME][0]
REVERSE_PRIMER_SEQUENCE = config.get('REVERSE_PRIMER_SEQUENCE')
if REVERSE_PRIMER_SEQUENCE is None:
    REVERSE_PRIMER_SEQUENCE = get_primers(PRIMERS_FILE)[REVERSE_PRIMER_NAME][0]

# Get reference sequence file and exit if none is provided and runDefCom should be run
REFERENCE_SEQUENCE_FILE = config.get('REFERENCE_SEQUENCE_FILE')
if (REFERENCE_SEQUENCE_FILE == None) and runDefCom:
    print("If you want to run runDefCom, please provide a reference sequence file.")
    exit(1)

# Define folder names
RAWFOLDER_NAME = '0raw'
CUTADAPTFOLDER_NAME = '1cutadapt'
FILTERANDTRIMFOLDER_NAME = '2filterAndTrim'
LEARNERRORFOLDER_NAME = '3learnerrors'
DADAFOLDER_NAME = '4sampleInference'
MERGEREADSFOLDER_NAME = '5mergeReads'
BIMERAREMOVALFOLDER_NAME = '6bimeraRemoval'
TAXONOMYFOLDER_NAME = '7taxonomy'
USEARCHFOLDER_NAME = "8uparsetax"

# Go through samples provided in the samples file and the blocklist.
SAMPLES = []
if SAMPLE_FILE.is_file():
    SAMPLES = set(SAMPLE_FILE.read_text().splitlines())
if not SAMPLE_FILE.is_file():
    print("Please provide a sample file. This error may be caused by running the snakemake command in the wrong directory.")
    exit(1)
if BLOCKLIST_FILE.is_file():
    BLOCKLIST = set(BLOCKLIST_FILE.read_text().splitlines())
    for blk in BLOCKLIST:
        SAMPLES.discard(blk)
if len(SAMPLES) == 0:
    print("No samples were found or all were discarded through the blocklist.")
    exit(1)

CUTADAPT_FILES = []
FILTERANDTRIM_FILES = []
LEARNERROR_FILES = []
DADA_FILES = []
MERGEREADS_FILES = []
BIMERA_FILES = []
TAXONOMY_FILES = []
OTU_ASV_FILES = []
STATS_FILES = []
STATS_FILES_TOT = []
DEFCOM_FILES = []

'''
END FILES + FOLDERS + PARAMS
'''

for sample in SAMPLES:
    '''
    SAMPLE SPECIFIC ANALYSIS FILES
    '''

    raw_stats = DATA_DIR.joinpath(RAWFOLDER_NAME,sample,sample + '.readstats.done')  #'{sample}.readstats.done'
    STATS_FILES.append(raw_stats)

    cutadapt_seqfile_1 = DATA_DIR.joinpath(
        CUTADAPTFOLDER_NAME,sample,sample + '_R1.fastq.gz')
    cutadapt_seqfile_2 = DATA_DIR.joinpath(
        CUTADAPTFOLDER_NAME,sample,sample + '_R2.fastq.gz')
    cutadapt_marker = DATA_DIR.joinpath(
        CUTADAPTFOLDER_NAME,sample,sample + '.cutadapt.done')
    cutadapt_stats = DATA_DIR.joinpath(CUTADAPTFOLDER_NAME,sample,sample + '.readstats.done')

    CUTADAPT_FILES.append(cutadapt_seqfile_1)
    CUTADAPT_FILES.append(cutadapt_seqfile_2)
    CUTADAPT_FILES.append(cutadapt_marker)
    STATS_FILES.append(cutadapt_stats)

    filterandtrim_seqfile_1 = DATA_DIR.joinpath(
        FILTERANDTRIMFOLDER_NAME,sample,sample + '_R1.fastq.gz')
    filterandtrim_seqfile_2 = DATA_DIR.joinpath(
        FILTERANDTRIMFOLDER_NAME,sample,sample + '_R2.fastq.gz')
    filterandtrim_marker = DATA_DIR.joinpath(
        FILTERANDTRIMFOLDER_NAME,sample,sample + '.filterandtrim.done')
    filterandtrim_stats = DATA_DIR.joinpath(FILTERANDTRIMFOLDER_NAME,sample,sample + '.readstats.done')

    FILTERANDTRIM_FILES.append(filterandtrim_seqfile_1)
    FILTERANDTRIM_FILES.append(filterandtrim_seqfile_2)
    FILTERANDTRIM_FILES.append(filterandtrim_marker)
    STATS_FILES.append(filterandtrim_stats)

    qc_stats_file = DATA_DIR.joinpath(
        FILTERANDTRIMFOLDER_NAME,sample,sample + '.qc.stats.pdf')
    STATS_FILES.append(qc_stats_file)

'''
PROJECT SPECIFIC ANALYSIS FILES
'''

# r1.samples & r2.samples
r1_samples_file = DATA_DIR.joinpath(FILTERANDTRIMFOLDER_NAME,'r1.samples')
r2_samples_file = DATA_DIR.joinpath(FILTERANDTRIMFOLDER_NAME,'r2.samples')
r1r2_marker_file = DATA_DIR.joinpath(
    FILTERANDTRIMFOLDER_NAME,'r1r2.samples.done')
LEARNERROR_FILES.append(r1_samples_file)
LEARNERROR_FILES.append(r2_samples_file)
LEARNERROR_FILES.append(r1r2_marker_file)

# Learn errors
r1_learnerrors_file = DATA_DIR.joinpath(
    LEARNERRORFOLDER_NAME,'R1.learnerrors.rds')
r1_learnerrors_marker = DATA_DIR.joinpath(
    LEARNERRORFOLDER_NAME,'R1.learnerrors.done')
r2_learnerrors_file = DATA_DIR.joinpath(
    LEARNERRORFOLDER_NAME,'R2.learnerrors.rds')
r2_learnerrors_marker = DATA_DIR.joinpath(
    LEARNERRORFOLDER_NAME,'R2.learnerrors.done')
LEARNERROR_FILES.append(r1_learnerrors_file)
LEARNERROR_FILES.append(r1_learnerrors_marker)
LEARNERROR_FILES.append(r2_learnerrors_file)
LEARNERROR_FILES.append(r2_learnerrors_marker)

# Inference
r1_dada_seqtab_file = DATA_DIR.joinpath(DADAFOLDER_NAME,'R1.seqtab.rds')
r1_dada_dd_file = DATA_DIR.joinpath(DADAFOLDER_NAME,'R1.dd.rds')
r1_dada_marker = DATA_DIR.joinpath(DADAFOLDER_NAME,'R1.dada.done')
r2_dada_seqtab_file = DATA_DIR.joinpath(DADAFOLDER_NAME,'R2.seqtab.rds')
r2_dada_dd_file = DATA_DIR.joinpath(DADAFOLDER_NAME,'R2.seqtab.rds')
r2_dada_marker = DATA_DIR.joinpath(DADAFOLDER_NAME,'R2.dada.done')
DADA_FILES.append(r1_dada_seqtab_file)
DADA_FILES.append(r1_dada_dd_file)
DADA_FILES.append(r1_dada_marker)
DADA_FILES.append(r2_dada_seqtab_file)
DADA_FILES.append(r2_dada_dd_file)
DADA_FILES.append(r2_dada_marker)

r1_seqtab_stats = DATA_DIR.joinpath(DADAFOLDER_NAME,'R1.seqtab.stats')
r1_seqtab_stats_marker = DATA_DIR.joinpath(DADAFOLDER_NAME,'R1.seqtab.stats.done')
r2_seqtab_stats = DATA_DIR.joinpath(DADAFOLDER_NAME,'R2.seqtab.stats')
r2_seqtab_stats_marker = DATA_DIR.joinpath(DADAFOLDER_NAME,'R2.seqtab.stats.done')
STATS_FILES.append(r1_seqtab_stats)
STATS_FILES.append(r1_seqtab_stats_marker)
STATS_FILES.append(r2_seqtab_stats)
STATS_FILES.append(r2_seqtab_stats_marker)

# Merge
merge_reads_seqtab = DATA_DIR.joinpath(
    MERGEREADSFOLDER_NAME,'seqtab.mergereads.rds')
merge_reads_dd = DATA_DIR.joinpath(MERGEREADSFOLDER_NAME,'dd.mergereads.rds')
merge_reads_marker = DATA_DIR.joinpath(
    MERGEREADSFOLDER_NAME,'mergereads.done')
MERGEREADS_FILES.append(merge_reads_seqtab)
MERGEREADS_FILES.append(merge_reads_dd)
MERGEREADS_FILES.append(merge_reads_marker)

merge_reads_seqtab_stats = DATA_DIR.joinpath(MERGEREADSFOLDER_NAME,'seqtab.mergereads.stats')
merge_reads_seqtab_stats_marker = DATA_DIR.joinpath(MERGEREADSFOLDER_NAME,'seqtab.mergereads.stats.done')
STATS_FILES.append(merge_reads_seqtab_stats)
STATS_FILES.append(merge_reads_seqtab_stats_marker)

# Bimera
nobimera_seqtab = DATA_DIR.joinpath(
    BIMERAREMOVALFOLDER_NAME,'seqtab.nobimera.rds')
nobimera_marker = DATA_DIR.joinpath(
    BIMERAREMOVALFOLDER_NAME,'seqtab.nobimera.done')

BIMERA_FILES.append(nobimera_marker)
BIMERA_FILES.append(nobimera_seqtab)
nobimera_seqtab_stats = DATA_DIR.joinpath(BIMERAREMOVALFOLDER_NAME,'seqtab.nobimera.stats')
nobimera_seqtab_stats_done = DATA_DIR.joinpath(BIMERAREMOVALFOLDER_NAME,'seqtab.nobimera.stats.done')
STATS_FILES.append(nobimera_seqtab_stats)
STATS_FILES.append(nobimera_seqtab_stats_done)

# Taxonomy
taxonomy_asv_file = DATA_DIR.joinpath(TAXONOMYFOLDER_NAME,'temp_asvs.tsv')
taxonomy_marker_file = DATA_DIR.joinpath(TAXONOMYFOLDER_NAME,'temp_asvs.done')

TAXONOMY_FILES.append(taxonomy_asv_file)
TAXONOMY_FILES.append(taxonomy_marker_file)

# OTUs
PROJECT_NAME = str(DATA_DIR).split('/')[-1]
asv_fasta_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.asvs.fasta')
asv_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.asvs.tsv')
otu_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.otus.tsv')
marker_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.done')
OTU_ASV_FILES.append(asv_fasta_file)
OTU_ASV_FILES.append(asv_file)
OTU_ASV_FILES.append(otu_file)
OTU_ASV_FILES.append(marker_file)

# Usearch
USEARCH_FILES = DATA_DIR.joinpath(USEARCHFOLDER_NAME,'uparse.done')
USEARCH_DIR = DATA_DIR.joinpath(USEARCHFOLDER_NAME)

# Final Insert stats
insert_counts_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.insert.counts')
STATS_FILES_TOT.append(insert_counts_file)

# Defined community
assignment_file = DATA_DIR.joinpath(f'{PROJECT_NAME}.refs.tsv')
DEFCOM_FILES.append(assignment_file)

if not runCutadapt:
    CUTADAPT_FILES = []
if not runQC:
    FILTERANDTRIM_FILES = []
if not runLearnErrors:
    LEARNERROR_FILES = []
if not runInference:
    DADA_FILES = []
if not runMergeReads:
    MERGEREADS_FILES = []
if not runRemoveBimeras:
    BIMERA_FILES = []
if not runReadStats:
    STATS_FILES = []
    STATS_FILES_TOT = []
if not runASVTax:
    TAXONOMY_FILES = []
if not runOTUTax:
    OTU_ASV_FILES = []
if not runUSEARCH:
    USEARCH_FILES = []
    USEARCH_DIR = []
if not runDefCom:
    DEFCOM_FILES = []

"""
In general the rules are run in the following order:
cutadapt, read_stats, dada2_filterAndTrim, qc_stats, r1r2_sample_files, dada2_learnErrors, dada2_inference,
seqtab_stats, dada2_mergeReads, bimera_removal, taxonomy, otu_asv, ref_assignment, insert_stats, uparse
"""

rule all:
    input:
        CUTADAPT_FILES,
        FILTERANDTRIM_FILES,
        STATS_FILES,
        LEARNERROR_FILES,
        DADA_FILES,
        MERGEREADS_FILES,
        BIMERA_FILES,
        TAXONOMY_FILES,
        OTU_ASV_FILES,
        STATS_FILES_TOT,
        USEARCH_FILES,
        DEFCOM_FILES


rule insert_stats:
    """
    Extract stats and write to file
    """
    input:
        stats_files=STATS_FILES,
        otu_asv_files=OTU_ASV_FILES,
    output:
        stats='{path}/' + PROJECT_NAME + '.insert.counts'
    shell:
        '''
        python {SCRIPTFOLDER}create_insert_stats.py {wildcards.path} {output.stats}
        '''


rule cutadapt:
    """
    Filter and trim fastq files, outputs trimmed reads which passed the filters
    """
    input:
        r1='{path}/' + RAWFOLDER_NAME + '/{sample}/{sample}_R1.fastq.gz',
        r2='{path}/' + RAWFOLDER_NAME + '/{sample}/{sample}_R2.fastq.gz',
    output:
        r1='{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}_R1.fastq.gz',
        r2='{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}_R2.fastq.gz',
        r1tmp=temp('{path}/' + CUTADAPTFOLDER_NAME +
                   '/{sample}/{sample}.temp1_R1.fastq'),
        r2tmp=temp('{path}/' + CUTADAPTFOLDER_NAME +
                   '/{sample}/{sample}.temp1_R2.fastq'),
        r1tmp2=temp('{path}/' + CUTADAPTFOLDER_NAME +
                    '/{sample}/{sample}.temp2_R1.fastq'),
        r2tmp2=temp('{path}/' + CUTADAPTFOLDER_NAME +
                    '/{sample}/{sample}.temp2_R2.fastq'),
        marker=touch('{path}/' + CUTADAPTFOLDER_NAME +
                     '/{sample}/{sample}.cutadapt.done')
    params:
        fwd_primer=FORWARD_PRIMER_SEQUENCE,
        rev_primer=REVERSE_PRIMER_SEQUENCE,
        minlength=QC_MINLEN,
        allow_untrimmed='--discard-untrimmed' if not allowUntrimmed else ''
    benchmark:
        '{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}.cutadapt.benchmark'
    threads:
        8
    log:
        command='{path}/' + CUTADAPTFOLDER_NAME + \
                '/{sample}/{sample}.cutadapt.command',
        log='{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}.cutadapt.log'
    shell:
        '''
        command="
        cutadapt -O 12 {params.allow_untrimmed} -g {params.fwd_primer} -G {params.rev_primer} -o {output.r1tmp} -p {output.r2tmp} {input.r1} {input.r2} -j {threads} --pair-adapters --minimum-length 75 &> {log.log}
        cutadapt -O 12 --times 5 -g {params.fwd_primer} -o {output.r1tmp2} -j {threads} {output.r1tmp} &>> {log.log}
        cutadapt -O 12 --times 5 -g {params.rev_primer} -o {output.r2tmp2} -j {threads} {output.r2tmp} &>> {log.log}
        cutadapt -o {output.r1} -p {output.r2} {output.r1tmp2} {output.r2tmp2} -j {threads} --minimum-length {params.minlength} &>> {log.log} 
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule dada2_filterAndTrim:
    """
    Truncate and discard reads using dada2 based on maxEE, quality score, number of N and length.
    """
    input:
        marker='{path}/' + CUTADAPTFOLDER_NAME + \
               '/{sample}/{sample}.cutadapt.done',
        fqgz1='{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}_R1.fastq.gz',
        fqgz2='{path}/' + CUTADAPTFOLDER_NAME + '/{sample}/{sample}_R2.fastq.gz',
    output:
        marker=touch('{path}/' + FILTERANDTRIMFOLDER_NAME +
                     '/{sample}/{sample}.filterandtrim.done'),
        fqgz1='{path}/' + FILTERANDTRIMFOLDER_NAME + \
              '/{sample}/{sample}_R1.fastq.gz',
        fqgz2='{path}/' + FILTERANDTRIMFOLDER_NAME + \
              '/{sample}/{sample}_R2.fastq.gz'
    params:
        maxee=QC_MAXEE,
        truncq='2',
        maxn='0',
        compress='TRUE',
        minlen=QC_MINLEN,
        trunc_R1=QC_TRUNC_R1,
        trunc_R2=QC_TRUNC_R2
    threads:
        4
    log:
        log='{path}/' + FILTERANDTRIMFOLDER_NAME + \
            '/{sample}/{sample}.filterandtrim.log',
        command='{path}/' + FILTERANDTRIMFOLDER_NAME + \
                '/{sample}/{sample}.filterandtrim.command'
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}dada2_filterandtrim.R {input.fqgz1} {input.fqgz2} {output.fqgz1} {output.fqgz2} {params.maxee} {params.truncq} {params.maxn} {params.compress} {params.minlen} {params.trunc_R1} {params.trunc_R2} {threads} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule qc_stats:
    """
    Plot a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file
    """
    input:
        r1_raw='{path}/' + RAWFOLDER_NAME + '/{sample}/{sample}_R1.fastq.gz',
        r2_raw='{path}/' + RAWFOLDER_NAME + '/{sample}/{sample}_R2.fastq.gz',
        r1_cut='{path}/' + CUTADAPTFOLDER_NAME + \
               '/{sample}/{sample}_R1.fastq.gz',
        r2_cut='{path}/' + CUTADAPTFOLDER_NAME + \
               '/{sample}/{sample}_R2.fastq.gz',
        marker_cut='{path}/' + CUTADAPTFOLDER_NAME + \
                   '/{sample}/{sample}.cutadapt.done',
        r1_fil='{path}/' + FILTERANDTRIMFOLDER_NAME + \
               '/{sample}/{sample}_R1.fastq.gz',
        r2_fil='{path}/' + FILTERANDTRIMFOLDER_NAME + \
               '/{sample}/{sample}_R2.fastq.gz',
        marker_fil='{path}/' + FILTERANDTRIMFOLDER_NAME + \
                   '/{sample}/{sample}.filterandtrim.done'
    output:
        pdf='{path}/' + FILTERANDTRIMFOLDER_NAME + \
            '/{sample}/{sample}.qc.stats.pdf'
    log:
        log='{path}/' + FILTERANDTRIMFOLDER_NAME + \
            '/{sample}/{sample}.qc.stats.log'
    shell:
        '''
        Rscript {SCRIPTFOLDER}qc_stats.R {input.r1_raw} {input.r2_raw} {input.r1_cut} {input.r2_cut} {input.r1_fil} {input.r2_fil} {output.pdf} &> {log.log}
        '''


rule r1r2_sample_files:
    """
    Create a file with sample names and the matching fastq file
    """
    input:
        FILTERANDTRIM_FILES
    output:
        r1_samples='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r1.samples',
        r2_samples='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r2.samples',
        marker=touch('{path}/' + FILTERANDTRIMFOLDER_NAME + '/r1r2.samples.done'),
    params:
        project_folder=DATA_DIR
    shell:
        '''
        python {SCRIPTFOLDER}create_samples_files.py {params.project_folder} {output.r1_samples} {output.r2_samples}
        '''

rule dada2_learnErrors:
    """
    Calculates error estimates for nucleotide transition probabilities
    """
    input:
        marker='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r1r2.samples.done',
        samples='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r{ori}.samples'
    output:
        marker=touch('{path}/' + LEARNERRORFOLDER_NAME +
                     '/R{ori}.learnerrors.done'),
        rds='{path}/' + LEARNERRORFOLDER_NAME + '/R{ori}.learnerrors.rds'
    params:
        nbases=LE_NBASES,
        folder=SCRIPTFOLDER
    threads:
        32
    log:
        log='{path}/' + LEARNERRORFOLDER_NAME + '/R{ori}.learnerrors.log',
        command='{path}/' + LEARNERRORFOLDER_NAME + \
                '/R{ori}.learnerrors.command'
    benchmark:
        '{path}/' + LEARNERRORFOLDER_NAME + '/R{ori}.learnerrors.benchmark'
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}dada2_learnErrors.R {input.samples} {output.rds} {threads} {params.nbases} {params.folder} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule dada2_inference:
    """
    Using the expectation of prior sequence variants created by dada2_learnErrors, we can use this information
    to increase sensitivity and detect singletons.
    """
    input:
        samples='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r{ori}.samples',
        marker_samples='{path}/' + \
                       FILTERANDTRIMFOLDER_NAME + '/r1r2.samples.done',
        marker='{path}/' + LEARNERRORFOLDER_NAME + '/R{ori}.learnerrors.done',
        errors='{path}/' + LEARNERRORFOLDER_NAME + '/R{ori}.learnerrors.rds',
    output:
        tab='{path}/' + DADAFOLDER_NAME + '/R{ori}.seqtab.rds',
        dd='{path}/' + DADAFOLDER_NAME + '/R{ori}.dd.rds',
        marker=touch('{path}/' + DADAFOLDER_NAME + '/R{ori}.dada.done')
    params:
        ref=REFERENCE_SEQUENCE_FILE,
        folder=SCRIPTFOLDER
    threads:
        32
    log:
        log='{path}/' + DADAFOLDER_NAME + '/R{ori}.dada2.log',
        command='{path}/' + DADAFOLDER_NAME + '/R{ori}.dada2.command'
    benchmark:
        '{path}/' + DADAFOLDER_NAME + '/R{ori}.dada2.benchmark'
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}dada2_inference.R {input.samples} {input.errors} {output.tab} {output.dd} {threads} {params.ref} {params.folder} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule dada2_mergeReads:
    """
    Merge the reads
    """
    input:
        marker_r1='{path}/' + DADAFOLDER_NAME + '/R1.dada.done',
        marker_r2='{path}/' + DADAFOLDER_NAME + '/R2.dada.done',
        samples_marker='{path}/' + \
                       FILTERANDTRIMFOLDER_NAME + '/r1r2.samples.done',
        samples_r1='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r1.samples',
        samples_r2='{path}/' + FILTERANDTRIMFOLDER_NAME + '/r2.samples',
        dd_r1='{path}/' + DADAFOLDER_NAME + '/R1.dd.rds',
        dd_r2='{path}/' + DADAFOLDER_NAME + '/R2.dd.rds'
    output:
        marker=touch('{path}/' + MERGEREADSFOLDER_NAME + '/mergereads.done'),
        merged_seqtab='{path}/' + MERGEREADSFOLDER_NAME + \
                      '/seqtab.mergereads.rds',
        merged_dd='{path}/' + MERGEREADSFOLDER_NAME + '/dd.mergereads.rds'
    threads:
        1
    log:
        log='{path}/' + MERGEREADSFOLDER_NAME + '/mergereads.log',
        command='{path}/' + MERGEREADSFOLDER_NAME + '/mergereads.command'
    benchmark:
        '{path}/' + MERGEREADSFOLDER_NAME + '/mergereads.benchmark'
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}dada2_mergeReads.R {input.samples_r1} {input.samples_r2} {input.dd_r1} {input.dd_r2} {output.merged_seqtab} {output.merged_dd} fw &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''


rule bimera_removal:
    """
    Remove bimeric sequences that are formed de novo during the sequencing process.
    """
    input:
        marker='{path}/' + MERGEREADSFOLDER_NAME + '/mergereads.done',
        merged_seqtab='{path}/' + MERGEREADSFOLDER_NAME + \
                      '/seqtab.mergereads.rds',
    output:
        tab_nobim='{path}/' + BIMERAREMOVALFOLDER_NAME + \
                  '/seqtab.nobimera.rds',
        marker=touch('{path}/' + BIMERAREMOVALFOLDER_NAME +
                     '/seqtab.nobimera.done')
    log:
        log='{path}/' + BIMERAREMOVALFOLDER_NAME + '/seqtab_nobimera.log',
        command='{path}/' + BIMERAREMOVALFOLDER_NAME + \
                '/seqtab_nobimera.command'
    threads:
        64
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}remove_bimera.R {input.merged_seqtab} {output.tab_nobim} {threads} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule taxonomy:
    """
    Load an RDS file containing a sequence table and assign the taxonomy of ASVs/OTUs using IDTAXA
    """
    input:
        tab_nobim='{path}/' + BIMERAREMOVALFOLDER_NAME + \
                  '/seqtab.nobimera.rds',
        marker='{path}/' + BIMERAREMOVALFOLDER_NAME + '/seqtab.nobimera.done'
    output:
        asv_table='{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.tsv',
        marker=touch('{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.done')
    log:
        log='{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.log',
        command='{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.command'
    benchmark:
        '{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.benchmark'
    threads:
        64
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}idtaxa_seqtab.R -i {input.tab_nobim} -s {SILVA_TRAINING_FILE} -t {threads} -o {output.asv_table} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''


rule otu_asv:
    """
    Loads an ASV table and uses the UPARSE algorithm to cluster ASVs into OTUs and produce an OTU table
    """
    input:
        asv_table='{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.tsv',
        marker='{path}/' + TAXONOMYFOLDER_NAME + '/temp_asvs.done'
    output:
        asv_fasta='{path}/' + PROJECT_NAME + '.asvs.fasta',
        asv_file='{path}/' + PROJECT_NAME + '.asvs.tsv',
        otu_file='{path}/' + PROJECT_NAME + '.otus.tsv',
        otu_fasta='{path}/' + PROJECT_NAME + '.otus.fasta',
        uparse='{path}/' + PROJECT_NAME + '.otus.uparse',
        marker=touch('{path}/' + PROJECT_NAME + '.done')
    log:
        log='{path}/' + PROJECT_NAME + '.log',
        command='{path}/' + PROJECT_NAME + '.command'
    benchmark:
        '{path}/' + PROJECT_NAME + '.benchmark'
    shell:
        '''
        command="
        Rscript {SCRIPTFOLDER}asvtab2otutab.R  -i {input.asv_table} -f {output.asv_fasta} -a {output.asv_file} -o {output.otu_file} -u {output.uparse} -x {output.otu_fasta} &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule read_stats:
    """
    Run bbmap: change quality encoding (qin) and run stats to check quality
    """
    input:
        fqz1='{sample}_R1.fastq.gz',
        fqz2='{sample}_R2.fastq.gz',
    output:
        marker=touch('{sample}.readstats.done'),
        bhist='{sample}.readstats.bhist',
        qhist='{sample}.readstats.qhist',
        qchist='{sample}.readstats.qchist',
        aqhist='{sample}.readstats.aqhist',
        bqhist='{sample}.readstats.bqhist',
        lhist='{sample}.readstats.lhist',
        gchist='{sample}.readstats.gchist'
    log:
        '{sample}.readstats.log'
    threads:
        4
    shell:
        '''
        reformat.sh -Xmx8G pigz=t bgzip=f threads={threads} qin=33 in1={input.fqz1} in2={input.fqz2} bhist={output.bhist} qhist={output.qhist} qchist={output.qchist} aqhist={output.aqhist} bqhist={output.bqhist} lhist={output.lhist} gchist={output.gchist} 2> {log}
       '''

rule seqtab_stats:
    """
    Create stats file from sequence table created by running dada2 inference.
    """
    input:
        rds_file='{sample}.rds'
    output:
        marker=touch('{sample}.stats.done'),
        stats='{sample}.stats'
    shell:
        '''
        Rscript {SCRIPTFOLDER}seqtab_stats.R {input.rds_file} {output.stats}
        '''


rule ref_assignment:
    """
    Perform sequence alignment between ASVs (Amplicon Sequence Variants) and a reference sequence database (assign ASVs
    to defined community members)
    """
    input:
        asvtab='{path}/' + PROJECT_NAME + '.asvs.tsv',
        asvs_fasta='{path}/' + PROJECT_NAME + '.asvs.fasta'
    output:
        asvs_ref='{path}/' + PROJECT_NAME + '.refs.tsv',
        marker=touch('{path}/' + PROJECT_NAME + '.refs.done')
    params:
        ref_fasta=REFERENCE_SEQUENCE_FILE,
        threads=32,
    log:
        log='{path}/' + PROJECT_NAME + '.refassign.log',
    threads:
        32
    shell:
        '''
        Rscript {SCRIPTFOLDER}assign_to_refs.R {input.asvtab} {input.asvs_fasta} {params.ref_fasta} 4 {output.asvs_ref} {params.threads}&> {log.log}
        '''

rule uparse:
    """
    Perform USEARCH sequence alignment against silva database for ASVs and OTUs and search for last common ancestor
    """
    input:
        input_otus_f=str(DATA_DIR / PROJECT_NAME) + '.otus.fasta',
        input_asvs_f=str(DATA_DIR / PROJECT_NAME) + '.asvs.fasta',
        database=USEARCH_DB
    output:
        done_file=touch(USEARCH_FILES),
        output_d=directory(USEARCH_DIR),
        output_lca_asvs=str(USEARCH_DIR) + PROJECT_NAME + '_asvs.lca',
        output_lca_otus=str(USEARCH_DIR) + PROJECT_NAME + '_otus.lca',
        output_tax_asvs=str(USEARCH_DIR) + PROJECT_NAME + '_asvs.tax',
        output_tax_otus=str(USEARCH_DIR) + PROJECT_NAME + '_otus.tax'
    shell:
        '''
        {SCRIPTFOLDER}uparse.sh -a {input.input_otus_f} -t {input.input_asvs_f} -o {output.output_d} -b {input.database} -n {PROJECT_NAME}
        '''
