import os
import sys
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from basecall_nanopore_methods import Methods
import pkg_resources
import shutil


__author__ = 'duceppemo'
__version__ = '0.1'


class Basecaller(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Guppy related
        self.device = args.device
        self.barcode_description = args.barcode_description
        self.barcode_kit = args.barcode_kit
        self.guppy_conf = args.guppy_conf
        self.recursive = args.recursive
        self.workflows = pkg_resources.resource_filename('data', 'workflows.tsv')

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        # Check input folder
        Methods.check_input(self.input)
        Methods.check_fast5(self.input)

        # Check if Guppy is installed
        Methods.check_guppy()

        # Check flowcell

        # Check library kit

        # Check

        print('\tAll good!')

        ##################
        #
        # Preparing outputs
        #
        ##################

        # Step completion report files
        done_basecalling = self.output_folder + '/done_basecalling'
        done_trimming = self.output_folder + '/done_trimming'
        done_filtering = self.output_folder + '/done_filtering'

        # Output folders to create
        basecalled_folder = self.output_folder + '/1_basecalled/'
        trimmed_folder = self.output_folder + '/2_trimmed/'
        filtered_folder = self.output_folder + '/3_filtered/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        ##################
        #
        # 1- Basecalling
        #
        ##################

        if not os.path.exists(done_basecalling):
            # Basecall fast5 to fastq
            Methods.run_guppy(self.input, basecalled_folder, self.guppy_conf, self.recursive,
                              self.device, self.barcode_kit)

            # Merge all fastq per barcode, if more than one file present
            Methods.merge_rename_fastq(basecalled_folder, self.barcode_kit)

            if self.barcode_description:
                sample_dict = Methods.parse_samples(self.barcode_description)
                Methods.rename_barcode(sample_dict, basecalled_folder)

            # Remove dump folder
            shutil.rmtree(basecalled_folder + 'guppy_basecaller-core-dump-db', ignore_errors=False, onerror=None)

            # Create "done" file for resuming purposes
            Methods.flag_done(done_basecalling)
        else:
            print('Skipping basecalling. Already done.')

        # Update sample_dict after extracting, only keep "pass" files
        self.sample_dict['basecalled'] = Methods.get_files(basecalled_folder, 'pass.fastq.gz')

        # Remove "unclassified" for next step if barcodes used
        if self.barcode_kit:
            self.sample_dict['basecalled'].pop('unclassified')

        ##################
        #
        # 2- Trim reads
        #
        ##################

        if not os.path.exists(done_trimming):
            print('Removing Nanopore adapters with Porechop...')
            Methods.run_porechop_parallel(self.sample_dict['basecalled'], trimmed_folder, self.cpu, self.parallel)
            Methods.flag_done(done_trimming)
        else:
            print('Skipping trimming. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['trimmed'] = Methods.get_files(trimmed_folder, '.fastq.gz')

        ##################
        #
        # 3- Filter reads
        #
        ##################

        # Get reference size
        if not os.path.exists(done_filtering):
            print('Filtering lower quality reads with Filtlong...')
            Methods.run_filtlong_parallel(self.sample_dict['trimmed'], filtered_folder, self.parallel)
            Methods.flag_done(done_filtering)
        else:
            print('Skipping filtering. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['filtered'] = Methods.get_files(filtered_folder, '.fastq.gz')

        print('DONE!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Basecall Nanopore raw data to fastq.')
    parser.add_argument('-i', '--input', metavar='/path/to/input_folder/',
                        required=True, type=str,
                        help='Folder that contains the fast5 files. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output_folder/',
                        required=True, type=str,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-f', '--flowcell', metavar='FLO-MIN106',
                        required=True, type=str,
                        help='Flowcell type used for sequencing. Mandatory')
    parser.add_argument('-l', '--library-kit', metavar='SQK-LSK109',
                        required=False, type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-k', '--barcode-kit', metavar='EXP-NBD104',
                        required=False, type=str,
                        help='Barcoding kit used. Use "unknown" if you know barcodes were used, but do not know '
                             'which kit. Not using this option will not perform barcode splitting. Optional')
    parser.add_argument('-b', '--barcode-description', metavar='/path/to/barcode_description.tsv',
                        required=False, type=str,
                        help='Tab-separated file with two columns with barcode assignments. '
                             'First column contains barcode names [barcode01, barcode02, etc.]. '
                             'Second column contains sample name. Avoid using special character. '
                             'Sample file in data folder. Optional.')
    parser.add_argument('-r', '--recursive',
                        action='store_true',
                        help='Look for fast5 recursively. Useful if fast5 are in multiple sub-folders. Optional')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False, type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-d', '--device', metavar='"cuda:0"',
                        required=True, type=str, default='"cuda:0"',
                        help='GPU device tp use. Typically "cuda:0" is one compatible graphics card is installed. '
                             'Use "cuda:0 cuda:1" (including the quotes) to use two graphics cards. '
                             'Default is "cuda:0". Mandatory.')
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False, type=int, default=2,
                        help='Number of samples to process in parallel for trimming and filtering. '
                             'Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False, type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({}). Optional.'.format(max_mem))
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Basecaller(arguments)
