import subprocess
import gzip
from statistics import median
import csv
from Bio import SeqIO
from nanomath import ave_qual
import os



def statistics(file, output):
    '''
    Calculate the mean Phred quality score distribution, N50 and other read statistics

    :param file: Path to the input file
    :param output: Path to the output file
    :return: None. Saves the statistics to a TSV file
    '''

    qualities = [10,11,12,13,14,15,16,17,18]
    phred_count = {q: 0 for q in qualities}
    total_reads = 0
    total_length = 0
    lengths = []


    with open_fastq(file) as fastq:
        for record in SeqIO.parse(fastq,'fastq'):
            qualscores = record.letter_annotations['phred_quality']
            mean_qual = ave_qual(qualscores)

            total_reads += 1
            total_length += len(record.seq)
            lengths.append(len(record.seq))

            for q in qualities:
                if mean_qual >= q:
                    phred_count[q] += 1

    mean_length = round(total_length / total_reads, 2)
    median_length = median(lengths)
    N50 = calculate_n50(lengths)

    stats = {}

    stats['total_bases'] = total_length
    stats['total_reads'] = total_reads
    stats['mean_read_length'] = mean_length
    stats['median_read_length'] = median_length
    stats['N50'] = N50

    for q in qualities:
        stats[f'Q{q}'] = str(round((phred_count[q] / total_reads) * 100, 2)) + '%'


    with open(output / 'statistics.tsv','w', newline='') as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerows(stats.items())



def trimming(file, output):
    '''
    Trimms the sequence using Porechop_abi

    :param file: Path to the input file
    :param output: Path to the output file
    :return: None. Saves the trimmed reads to a compressed FASTQ file
    '''

    subprocess.run(f'porechop_abi -abi -i {file} -o {output / "sequence_trimmed.fastq.gz"}',
                    shell=True, check=True)



def filtering(file, output):
    '''
    Filter reads based on the selected minimum Phred quality of 10 using Chopper.

    :param file: Path to the input file
    :param output: Path to the output file
    :return: None. Saves the filtered reads to a compressed FASTQ file.
    '''

    if file.suffix == '.gz':
        subprocess.run(f'gunzip -c {file} | chopper -q 10 | gzip > {output / "sequence_filtered.fastq.gz"}',
                   shell=True, check=True)
    else:
        subprocess.run(f'cat -c {file} | chopper -q 10 | gzip > {output / "sequence_filtered.fastq.gz"}',
                   shell=True, check=True)



def open_fastq(file):
    '''
    Open a FASTQ file in text mode, automatically detecting gzip compression.

    :param file: Path to the input file
    :return: File handle for reading the FASTQ records.
    '''
    if file.suffix == '.gz':
        return gzip.open(file, 'rt')
    else:
        return open(file, 'rt')



def calculate_n50(lengths):
    '''
    Calculate the N50 statistic from a list of read lengths.

    :param lengths: List of read lengths
    :return: N50 value
    '''
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    half_length = sum(sorted_lengths) / 2
    cumulative_length = 0

    for length in sorted_lengths:
        cumulative_length += length

        if cumulative_length >= half_length:
            return length