
import glob


r1_files = glob.glob('0raw_original/*/*R1.fastq.gz')

forward_primer = "AGGAAGCCCCGGCTAACTC"
reverse_primer = "GACAGCCATGCAGCACCTC"

sample_2_r1_file = {}

marker_files = []
for r1_file in r1_files:
    marker = 'cutadapt_split/' + r1_file.split('/')[-1].replace('_R1.fastq.gz', '.cutadapt.done')
    sample = marker.split('/')[-1].split('.cutadapt')[0]
    sample_2_r1_file[sample] = r1_file
    marker_files.append(marker)



def getR1(wildcards):
    return sample_2_r1_file[wildcards.sample]

def getR2(wildcards):
    return getR1(wildcards).replace('_R1', '_R2')

rule all:
    input:
        marker_files
rule cutadapt:
    input:
        fqgzR1 = getR1,
        fqgzR2 = getR2
    output:
        marker = touch('cutadapt_split/{sample}.cutadapt.done'),
        fqgz1R1t = 'cutadapt_split/{sample}.fw_R1.fastq.gz',
        fqgz1R2t = 'cutadapt_split/{sample}.fw_R2.fastq.gz',
        fqgz2R1t = 'cutadapt_split/{sample}.rev_R1.fastq.gz',
        fqgz2R2t = 'cutadapt_split/{sample}.rev_R2.fastq.gz',
    log:
        log = 'cutadapt_split/{sample}.cutadapt.log',
        command = 'cutadapt_split/{sample}.cutadapt.command'
    threads:
        8
    shell:
        '''
        command="

        echo FW paired > {log.log};
        cutadapt -O 10 --action=retain --discard-untrimmed -g {forward_primer} -G {reverse_primer}  -o {output.fqgz1R1t} -p {output.fqgz1R2t} {input.fqgzR1} {input.fqgzR2} -j {threads} --pair-adapters --minimum-length 75 &>> {log.log};


        echo REV paired >> {log.log};
        cutadapt -O 10 --action=retain --discard-untrimmed -G {forward_primer} -g {reverse_primer}  -o {output.fqgz2R1t} -p {output.fqgz2R2t} {input.fqgzR1} {input.fqgzR2} -j {threads} --pair-adapters --minimum-length 75 &>> {log.log};


        ";
        echo "$command" > {log.command};
        eval "$command"
        '''
