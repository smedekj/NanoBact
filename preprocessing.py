import subprocess
import gzip
from Bio import SeqIO
from nanomath import ave_qual


def statistics(file):
    '''
    Calculate the mean Phred quality score distribution, N50 statistic and average read length.

    :param file: Path to the input file
    :return: Dictionary with percentages of reads with each Phred quality, N50 statistic and the average read length
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
    N50 = calculate_n50(lengths)

    phred_retention = {q: 0 for q in qualities}
    for q in qualities:
        phred_retention[q] = round((phred_count[q] / total_reads) * 100, 2)

    stats = phred_retention
    stats['total_reads'] = total_reads
    stats['mean_read_length'] = mean_length
    stats['total_bases'] = total_length
    stats['N50'] = N50

    return stats



def filtering(file, stats):
    '''
    Filter reads based on the selected minimum Phred quality threshold using Chopper.

    :param file: Path to the input file
    :param stats: Dictionary containing read quality statistics.
    :return: None. Saves the filtered reads to a compressed FASTQ file.
    '''
    benchmarks = filtering_benchmarks(stats)

    if file.suffix == '.gz':
        subprocess.run(f'gunzip -c {file} | chopper -q {benchmarks['phred']} | gzip > sequence_filtered.fastq.gz',
                   shell=True, check=True)
    else:
        subprocess.run(f'cat -c {file} | chopper -q {benchmarks['phred']} | gzip > sequence_filtered.fastq.gz',
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



def filtering_benchmarks(stats):
    '''
    Determine the minimum Phred quality threshold based on read retention.

    :param stats: Dictionary containing the percentage of reads meeting each Phred quality threshold.
    :return: Dictionary containing the selected minimum Phred threshold.
    '''
    benchmarks = {}
    for i in range(10,19):
        if stats[i] <= 92:
            benchmarks['phred'] = i
            break

    return benchmarks



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