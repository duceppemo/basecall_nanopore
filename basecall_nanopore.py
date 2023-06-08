import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from basecall_nanopore_methods import Methods
import pkg_resources
import shutil
from kits import Kits


__author__ = 'duceppemo'
__version__ = '0.1'

# TODO: add option to use multiple barcode kits
#       Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes.
# TODO: check that is the GPU version of guppy installed
#       Maybe by running "nvidia-smi" and checkking 1) no error (gpu driver working properly)
#       2) if "guppy_basecall_server" is running (printed in the ouptut of the command)


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
        self.gpu = args.gpu
        self.description = args.description
        self.barcode_kit = args.barcode_kit[0].split()
        self.sequencer = args.sequencer
        self.config = args.config
        self.flowcell = args.flowcell
        self.library_kit = args.library_kit
        self.recursive = args.recursive
        # self.accuracy = args.accuracy
        self.workflows = pkg_resources.resource_filename('data', 'workflows.tsv')

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):

        ##################
        #
        # Checks
        #
        ##################

        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        # Check input folder
        Methods.check_input(self.input)
        Methods.check_fast5(self.input)

        # Check if Guppy is installed
        Methods.check_guppy()

        # Check for config file
        Methods.check_config(self.config, self.flowcell, self.sequencer, self.library_kit)

        # Check barcodes
        Methods.check_barcode(self.barcode_kit, self.description)

        print('\tAll good!')

        ##################
        #
        # Preparing outputs
        #
        ##################

        # Step completion report files
        done_basecalling = self.output_folder + '/done_basecalling'
        done_qc = self.output_folder + '/done_QC'
        done_trimming = self.output_folder + '/done_trimming'
        done_filtering = self.output_folder + '/done_filtering'

        # Output folders to create
        basecalled_folder = self.output_folder + '/1_basecalled/'
        qc_folder = self.output_folder + '/2_qc/'
        trimmed_folder = self.output_folder + '/3_trimmed/'
        filtered_folder = self.output_folder + '/4_filtered/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        ##################
        #
        # 1- Basecalling
        #
        ##################

        if not os.path.exists(done_basecalling):
            # Retrieve proper configuration file
            if not self.config:
                guppy_conf = Methods.get_guppy_config(self.flowcell, self.library_kit, self.sequencer, self.workflows)
            else:
                guppy_conf = self.config

            # Basecall fast5 to
            Methods.run_guppy(self.input, basecalled_folder, guppy_conf, self.recursive,
                              self.gpu, self.barcode_kit)

            # Merge all fastq per barcode, if more than one file present
            Methods.merge_rename_fastq(basecalled_folder, self.barcode_kit)

            if self.description:
                sample_dict = Methods.parse_samples(self.description)
                Methods.rename_barcode(sample_dict, basecalled_folder)  # Also remove extra barcode folders

            # Remove dump folder
            dump_file = basecalled_folder + 'guppy_basecaller-core-dump-db'
            if os.path.exists(dump_file):
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
        # 2- QC
        #
        ##################

        if not os.path.exists(done_qc):
            print('Performing read QC with PycoQC...')
            Methods.run_pycoqc(basecalled_folder, qc_folder)
            Methods.flag_done(done_qc)
        else:
            print('Skipping QC. Already done.')

        ##################
        #
        # 3- Trim reads
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
        # 4- Filter reads
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

        ##################
        #
        # Done
        #
        ##################

        # Remove 'guppy_basecaller-core-dump-db' ?
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
    parser.add_argument('-s', '--sequencer',
                        required=False, type=str,
                        choices=['minion', 'promethion'],
                        help='Sequencer used. "minion" includes all sequencers except "promethion". '
                             'Optional. Default is "minion".')
    parser.add_argument('-c', '--config', metavar='dna_r9.4.1_450bps_sup.cfg',
                        required=False, type=str,
                        help='Guppy config file. Typically found in "/opt/ont/guppy/data" ending with ".cfg". '
                             'This is the prefered methods over choosing a "library kit/flowcell" combination. '
                             'Both methods are uncompatible. Optional.')
    parser.add_argument('-f', '--flowcell', metavar='FLO-MIN106',
                        required=False, type=str,
                        help='Flowcell type used for sequencing. Optional.')
    parser.add_argument('-l', '--library-kit', metavar='SQK-LSK109',
                        required=False, type=str,
                        help='Library kit used. Optional.')
    parser.add_argument('-b', '--barcode-kit', metavar='EXP-NBD104',
                        required=False, type=str, nargs='+',
                        help='Barcoding kit(s) used. Use "unknown" if you know barcodes were used, but do not know '
                             'which kit. Not using this option will not perform barcode splitting. For multiple '
                             'barcoding kits, use double quotes and space like this: "EXP-NBD104 EXP-NBA114". Optional')
    parser.add_argument('-d', '--description', metavar='/path/to/barcode_description.tsv',
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
    parser.add_argument('-g', '--gpu', metavar='"cuda:0"',
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
