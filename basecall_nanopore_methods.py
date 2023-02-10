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


# mamba create -n nanopore -y -c bioconda \
# flye samtools parallel bbmap shasta porechop filtlong bandage minimap2 blast psutil


class Methods(object):
    barcoding_kit_list = ['EXP-NBD104', '', 'EXP-NBA114', 'SQK-RBK114.96', 'SQK-RBK114.24', 'SQK-NBD114.96',
                      'SQK-NBD114.24', 'SQK-PCB111.24', 'SQK-RPB004', 'SQK-PBK004', 'SQK-RBK110.96',
                      'SQK-RBK004', 'SQK-16S024', 'SQK-RAB204', 'EXP-PBC096', 'EXP-PBC001', 'EXP-NBD114',
                      'EXP-NBD104', 'EXP-NBD196', 'unknown']

    library_kit_list = ['SQK-16S024', 'SQK-CS9109', 'SQK-DCS108', 'SQK-DCS109', 'SQK-LRK001', 'SQK-LSK108',
                        'SQK-LSK109', 'SQK-LSK109-XL', 'SQK-LSK110', 'SQK-LSK110-XL', 'SQK-LSK111', 'SQK-LSK111-XL',
                        'SQK-LSK112', 'SQK-LSK112-XL', 'SQK-LSK114', 'SQK-LSK114-XL', 'SQK-LSK308', 'SQK-LSK309',
                        'SQK-LSK319', 'SQK-LWB001', 'SQK-LWP001', 'SQK-MLK111-96-XL', 'SQK-NBD111-24', 'SQK-NBD111-96',
                        'SQK-NBD112-24', 'SQK-NBD112-96', 'SQK-NBD114-24', 'SQK-NBD114-96', 'SQK-PBK004', 'SQK-PCB109',
                        'SQK-PCB110', 'SQK-PCB111-24', 'SQK-PCS108', 'SQK-PCS109', 'SQK-PCS111', 'SQK-PSK004',
                        'SQK-RAB201', 'SQK-RAB204', 'SQK-RAD002', 'SQK-RAD003', 'SQK-RAD004', 'SQK-RAD112',
                        'SQK-RAD114', 'SQK-RAS201', 'SQK-RBK001', 'SQK-RBK004', 'SQK-RBK110-96', 'SQK-RBK111-24',
                        'SQK-RBK111-96', 'SQK-RBK112-24', 'SQK-RBK112-96', 'SQK-RBK114-24', 'SQK-RBK114-96',
                        'SQK-RLB001', 'SQK-RLI001', 'SQK-RNA001', 'SQK-RNA002', 'SQK-RPB004', 'SQK-ULK001',
                        'SQK-ULK114', 'VSK-PTC001', 'VSK-VBK001', 'VSK-VMK001', 'VSK-VMK004', 'VSK-VPS001',
                        'VSK-VSK001', 'VSK-VSK003', 'VSK-VSK004']

    flowcell_list = ['FLO-FLG001', 'FLO-FLG111', 'FLO-FLG114', 'FLO-MIN106', 'FLO-MIN107', 'FLO-MIN110', 'FLO-MIN111',
                     'FLO-MIN112', 'FLO-MIN114', 'FLO-MINSP6', 'FLO-PRO001', 'FLO-PRO002', 'FLO-PRO002-ECO',
                     'FLO-PRO002M', 'FLO-PRO111', 'FLO-PRO112', 'FLO-PRO112M', 'FLO-PRO114', 'FLO-PRO114M']

    configuration_file_list = ['dna_r10.3_450bps_hac', 'dna_r10.3_450bps_hac_prom', 'dna_r10.4.1_e8.2_260bps_hac',
                               'dna_r10.4.1_e8.2_260bps_hac_prom', 'dna_r10.4.1_e8.2_400bps_hac',
                               'dna_r10.4.1_e8.2_400bps_hac_prom', 'dna_r10_450bps_hac', 'dna_r10.4_e8.1_hac',
                               'dna_r10.4_e8.1_hac_prom', 'dna_r9.4.1_450bps_hac', 'dna_r9.4.1_450bps_hac_prom',
                               'dna_r9.4.1_e8.1_hac', 'dna_r9.4.1_e8.1_hac_prom', 'dna_r9.5_450bps',
                               'rna_r9.4.1_70bps_hac', 'rna_r9.4.1_70bps_hac_prom']

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
                else:  # not supposed to be there
                    shutil.rmtree(barcode_folder, ignore_errors=False, onerror=None)  # Delete non-empty folder

    @staticmethod
    def remove_extra_barcodes(sample_dict, basecalled_folder):
        pass

    @staticmethod
    def run_guppy(fast5_folder, basecalled_folder, guppy_conf, recursive, device, barcode_kit):
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
