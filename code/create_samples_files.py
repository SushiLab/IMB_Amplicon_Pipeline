import sys
import glob

data_dir = sys.argv[1]
r1_file = sys.argv[2]
if len(sys.argv) == 4:
    r2_file = sys.argv[3]


def create_file_with_samples_and_file_paths(out_file, orientation, data_directory):
    """
    Create a file with sample names and the matching fastq file. (needed by dada)

    Parameters:
    - out_file (str): name of the file that will be created
    - orientation (str): orientation of fastq file
    - data_directory (str): directory of the data
    """
    of = open(out_file, 'w')
    sample_2_file = {}
    offset = -2

    # Find filtered and trimmed files
    files = glob.glob(data_directory + '/2filterAndTrim/' + '*/*' + orientation)
    if len(files) < 1:
        files = glob.glob(data_directory + '/2filterAndTrim/' + '*' + orientation)
        offset = -1

    # Get needed information (sample name and file path)
    for trimmed_file in files:
        sample = trimmed_file.split('/')[offset]
        sample_2_file[sample] = trimmed_file

    # Write output files
    for sample in sorted(sample_2_file.keys()):
        of.write(f'{sample}\t{sample_2_file[sample]}\n')

    of.close()


create_file_with_samples_and_file_paths(r1_file, 'R1.fastq.gz', data_dir)
if len(sys.argv) == 4:
    create_file_with_samples_and_file_paths(r2_file, 'R2.fastq.gz', data_dir)
