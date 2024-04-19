import glob
import os

source_folder = '/path/to/your/data/' # where your raw data is located
dest_folder = '/path/to/your/analysis/0raw/' # 0raw folder within folder where you run your analysis

print('set -euxo pipefail')
r1_files = glob.glob(source_folder + '*/*raw_1.fq.gz') # change this to the correct file ending. Possibly add or remove "*/".
for r1_file in r1_files:
    r2_file = r1_file.replace('raw_1.fq.gz', 'raw_2.fq.gz') # change this to the correct file ending (same file ending as chosen for r1_files)
    if r1_file == r2_file:
        print(f'r1 and r2 file are identical, did you rename correctly?')
        exit(1)
    if not os.path.exists(r2_file):
        print(f'r2 file doesnt exist. Quitting.')
        exit(1)
    samplename = r1_file.split('/')[-1].replace('_raw_1.fq.gz', '') # change this to the correct file ending (same file ending as chosen for r1_files)

    new_sample_name = 'NAME24-1_' + samplename + '_METAB' # change NAME24-1
    out_folder = dest_folder + new_sample_name + '/'
    print(f'mkdir -p {out_folder}')
    dest_r1_file = out_folder + new_sample_name + '_R1.fastq.gz'
    dest_r2_file = out_folder + new_sample_name + '_R2.fastq.gz'
    print(f'ln -s {r1_file} {dest_r1_file}')
    print(f'ln -s {r2_file} {dest_r2_file}')
