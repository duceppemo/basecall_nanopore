import subprocess
import os
import sys
from concurrent import futures
import pathlib
from psutil import virtual_memory
from multiprocessing import cpu_count
import gzip
from glob import glob
import shutil
import pandas as pd
from kits import Kits


# mamba create -n nanopore -y -c bioconda \
# flye samtools parallel bbmap shasta porechop filtlong bandage minimap2 blast psutil pandas pycoqc pysam


class Methods(object):
    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def check_input(input_folder):
        if not os.path.exists(input_folder):
            raise Exception('Please select an existing folder as input.')

        # Check if folder
        if not os.path.isdir(input_folder):
            raise Exception('Please select a folder as input.')

    @staticmethod
    def check_fast5(input_folder):
        # Check if input folder contains fast5
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith('.fast5'):  # accept a tuple or string
                    return  # At least one fast5 file is present
        raise Exception('No fast5 files detected in provided input folder.')

    @staticmethod
    def check_guppy():
        cmd = ['guppy_basecaller', '--version']
        status = subprocess.getstatusoutput(' '.join(cmd))
        if status[0] != 0:
            raise Exception('The GPU version of Guppy must be installed first.')
        else:
            guppy_version = status[1]
            guppy_version = guppy_version.split(',')[1]  # Works with version 6.3.8
            guppy_version = '.'.join(guppy_version.split('.')[1:4])
            guppy_version = guppy_version.split('+')[0]
            print('Running Guppy{}'.format(guppy_version))

    @staticmethod
    def check_from_list(my_category, my_item, my_list):
        if my_item not in my_list:
            print('Please use of the following choice for {}: {}'.format(my_category, my_list))
            sys.exit()

    @staticmethod
    def check_version(log_file):
        # Not being used right now because versions are captured in the requirements.txt file
        with open(log_file, 'w') as f:
            # Python
            p = subprocess.Popen(['python', '--version'])
            stderr, stdout = p.communicate()

            # Porechop
            p = subprocess.Popen(['porechop', '--version'])
            stderr, stdout = p.communicate()

            # Filtlong
            p = subprocess.Popen(['filtlong', '--version'])
            stderr, stdout = p.communicate()

    @staticmethod
    def check_config(config, flowcell, sequencer, library_kit):
        if config and (flowcell or library_kit):
            raise Exception('Please chose a configuration file or a "library kit/flowcell/sequencer" combination, '
                            'not both.')
        if config:
            if config not in Kits.configuration_file_list:
                raise Exception('Please use one of the following supported configuration file: {}'
                                .format(Kits.configuration_file_list))
        if not config:
            if flowcell and library_kit and sequencer:
                # Check flowcell
                Methods.check_from_list('flowcell', flowcell, Kits.flowcell_list)
                # Check library kit
                Methods.check_from_list('library kit', library_kit, Kits.library_kit_list)
            else:
                raise Exception('Please make sure you are selection a library kit, a flowcell and a sequencer if '
                                'you are not using a configuration file')

    @staticmethod
    def check_barcode(barcode_kit, barcode_description):
        if barcode_kit:
            Methods.check_from_list('barcoding kit', barcode_kit, Kits.barcoding_kit_list)

    @staticmethod
    def get_guppy_config(flowcell, library, sequencer, workflows):
        # Parse workflow.tsv to pandas df
        df = pd.read_csv(workflows, sep='\t', header=0)

        # Return the config (4th column) matching the requested flowcell, library, accuracy and sequencer
        guppy_conf_list = df.loc[(df['flowcell'] == flowcell) & (df['kit'] == library), 'config_name'].values

        if sequencer == 'promethion':
            guppy_conf_list = [x for x in guppy_conf_list if 'prom' in x]  # Assume it will only return one...

        if len(guppy_conf_list) == 1:
            conf = guppy_conf_list[0]
        else:
            raise Exception('More than one config file matching '
                            '{}, {} and {}'.format(flowcell, library, sequencer))

        if not conf:
            raise Exception('Could not find a configuration for the requested flowcell, library kit, accuracy and '
                            'sequencer ({}, {}, {})'.format(flowcell, library, sequencer))
        else:
            return conf + '.cfg'

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if they do not exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder, ext):
        sample_dict = dict()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(ext):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]
                    sample_dict[sample] = file_path
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob(folder + '/*' + extension)

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass

    @staticmethod
    def gzipped_file_size(gzipped_file):
        with gzip.open(gzipped_file, 'rb') as f:
            return f.seek(0, whence=2)

    @staticmethod
    def merge_files(file_list, merged_file):
        with open(merged_file, 'wb') as wfd:
            for f in file_list:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

    @staticmethod
    def delete_unmerged(file_list):
        for f in file_list:
            os.remove(f)

    @staticmethod
    def merge_rename_fastq(fastq_folder, barcode_kit):
        # List all folders
        for i in ['pass', 'fail']:
            if not barcode_kit:
                fastq_list = glob(fastq_folder + i + '/fastq_runid_*.fastq.gz')
                merged_fastq = fastq_folder + i + '/' + i + '.fastq.gz'
                Methods.merge_files(fastq_list, merged_fastq)
                Methods.delete_unmerged(fastq_list)
            else:
                # List directory (each barcode)
                folder_list = glob(fastq_folder + i + '/*/')
                for barcode_folder in folder_list:
                    fastq_list = glob(barcode_folder + '/fastq_runid_*.fastq.gz')
                    barcode_name = barcode_folder.split('/')[-2]
                    merged_fastq = fastq_folder + i + '/' + barcode_name + '/' + barcode_name + '_' + i + '.fastq.gz'
                    Methods.merge_files(fastq_list, merged_fastq)
                    Methods.delete_unmerged(fastq_list)

    @staticmethod
    def parse_samples(barcode_desc):
        sample_dict = dict()

        with open(barcode_desc, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                try:
                    barcode, sample_name = line.split('\t')
                except IndexError:
                    raise Exception('The sample decription file must be a tab-separated file with two columns where '
                                    'the first one is the barcode (e.g. barcode01) name and the second one the sample '
                                    'name (e.g. my_sample')
                sample_dict[barcode] = sample_name

        return sample_dict

    @staticmethod
    def rename_barcode(sample_dict, basecalled_folder):
        for i in ['pass', 'fail']:
            # Rename folders
            folder_list = glob(basecalled_folder + i + '/*/')
            for barcode_folder in folder_list:
                barcode_name = barcode_folder.split('/')[-2]
                if barcode_name in sample_dict:
                    folder_new_name = '/'.join(barcode_folder.split('/')[:-2]) + '/' + sample_dict[barcode_name] + '/'
                    fastq_current_name = folder_new_name + '/' + barcode_name + '_' + i + '.fastq.gz'
                    fastq_new_name = folder_new_name + '/' + sample_dict[barcode_name] + '_' + i + '.fastq.gz'
                    os.rename(barcode_folder, folder_new_name)  # Rename folder
                    os.rename(fastq_current_name, fastq_new_name)  # Rename fastq
                elif barcode_name == 'unclassified':
                    continue
                else:  # Delete barcodes found but not present en description file. Not supposed to be there
                    shutil.rmtree(barcode_folder, ignore_errors=False, onerror=None)  # Delete non-empty folder

    @staticmethod
    def run_guppy(fast5_folder, basecalled_folder, guppy_conf, recursive, device, barcode_kit):
        Methods.make_folder(basecalled_folder)
        os.chdir(basecalled_folder)  # avoid "guppy_basecaller-core-dump-db" folder created in script location

        cmd = ['guppy_basecaller',
               '--config', guppy_conf,
               '--input_path', fast5_folder,
               '--save_path', basecalled_folder,
               '--calib_detect',
               '--records_per_fastq', str(0),
               '--compress_fastq',
               '--disable_pings',
               '--gpu_runners_per_device', str(2),
               '--chunk_size', str(1000),
               '--chunks_per_runner', str(128),
               '--device', device,
               '--detect_adapter',
               '--detect_primer',
               '--trim_adapters',
               '--trim_primers']
        if recursive:
            cmd += ['--recursive']
        if barcode_kit:
            cmd += ['--detect_barcodes']
            if barcode_kit != 'unknown':
                cmd += ['--barcode_kits', barcode_kit]

        subprocess.run(cmd)

    @staticmethod
    def rename_basecalled(basecalled_folder, sample_dict):
        """
        basecalled_folder
            |-pass
                |-barcode01
                    |-file1.fastq.gz
                    |-...
                |-...
                |-unclassified
            |-fail
                |-barcode01
                    |-file1.fastq.gz
                    |-...
                |-...
                |-unclassified
        """
        # Merge fastq so there is only one file per sample
        for root, directories, filenames in os.walk(basecalled_folder):
            for folder in directories:
                pass
            for filename in filenames:
                if filename.endswith('.fastq.gz'):  # accept a tuple or string
                    file_path = os.path.join(root, filename)

                    sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]
                    sample_dict[sample] = file_path
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        Methods.get_files(basecalled_folder)

    @staticmethod
    def run_pycoqc(basecalled_folder, report_folder):
        Methods.make_folder(report_folder)
        cmd = ['pycoQC',
               '-f', basecalled_folder + 'sequencing_summary.txt',
               '-o', report_folder + 'pycoQC_output.html']
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_porechop(sample, input_fastq, trimmed_folder, cpu):
        cmd = ['porechop',
               '-i', input_fastq,
               '-o', trimmed_folder + sample + '.fastq.gz',
               '--threads', str(cpu),
               '--check_reads', str(1000)]  # Only check adapter from 1,000 reads instead of 10,000

        print('\t{}'.format(sample))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_porechop_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_porechop(*x), args):
                pass

    @staticmethod
    def run_filtlong(sample, input_fastq, filtered_folder):
        print('\t{}'.format(sample))

        cmd = ['filtlong',
               '--keep_percent', str(95),  # Drop bottom 5% reads
               input_fastq]

        # Filtlong writes to stdout
        filtered_fastq = filtered_folder + sample + '.fastq.gz'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        with gzip.open(filtered_fastq, 'wb') as f:
            f.write(p.communicate()[0])

    @staticmethod
    def run_filtlong_parallel(sample_dict, output_folder, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_filtlong(*x), args):
                pass
