import subprocess
import shutil


def assembly(file, assembler, output, threads):
    """
    Perform de novo genome assembly using the selected assembler.

    Supports Flye and Raven assemblers.

    :param file: Path to the input file containing sequencing reads.
    :param assembler: Name of the assembler ('flye' or 'raven').
    :param output: Path to the main output directory.
    :param threads: Number of CPU threads to use.
    :return: None. Saves the assembled genome as assembly/assembly.fasta.
    """

    output = output / 'assembly'
    output.mkdir(parents=True, exist_ok=True)

    if assembler == 'flye':
        flye_logs = output / 'flye_logs'
        subprocess.run(['flye', '--nano-hq', file, '--out-dir', flye_logs, '--threads', str(threads)],
            check=True
        )

        assembly_file = output / 'assembly.fasta'
        shutil.move(flye_logs / 'assembly.fasta', assembly_file)


    elif assembler == 'raven':
        with open(output / 'assembly.fasta','w') as fasta:
            subprocess.run(['raven', file, '--threads', str(threads)],
                stdout=fasta, check=True
            )



def polishing(file, assembly, output, threads):
    """
    Polish a genome assembly using Medaka.

    Uses the original reads and the assembly FASTA file.

    :param file: Path to the original input file.
    :param assembly: Path to the assembly FASTA file.
    :param output: Path to the main output directory.
    :param threads: Number of CPU threads to use.
    :return: None. Saves the polished assembly and intermediate files
             in the polishing directory.
    """

    output = output / 'polishing'

    subprocess.run(['medaka_consensus', '-i', file, '-d', assembly, '-o', output, '-t', str(threads), '--bacteria'],
        check=True
    )