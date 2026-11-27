from argument_processing import arguments_process
from preprocessing import statistics, trimming, filtering
from assembly_polishing import assembly, polishing
import sys
from pathlib import Path

def main():
    """
    Run the selected NanoBact processing steps.

    Processing steps are executed in a predefined order:
    statistics, trimming, filtering, assembly and polishing.

    Each processing step uses the output of the previous
    applicable step as its input.
    """

    output_dir = Path('output')
    output_dir.mkdir(parents=True, exist_ok=True)
    args = arguments_process(sys.argv)

    if args['whole_pipeline']:
        steps = {
            'statistics': False,
            'trimming': True,
            'filtering': True,
            'assembly': True,
            'polishing': True
        }
        print('Running the whole NanoBact pipeline')

    else:
        steps = {
            'statistics': args['statistics'],
            'trimming': args['trimming'],
            'filtering': args['filtering'],
            'assembly': args['assembly'],
            'polishing': args['polishing']
        }

    current_file = args['file']

    if steps['trimming']:
        print('Trimming reads')
        trimming(current_file, output_dir)

        current_file = output_dir / 'sequence_trimmed.fastq.gz'

    if steps['filtering']:
        print('Filtering low quality reads')
        filtering(current_file, output_dir)

        current_file = output_dir / 'sequence_filtered.fastq.gz'

    if steps['statistics']:
        print('Gathering read statistics')
        statistics(current_file, output_dir)

    assembly_file = None

    if steps['assembly']:
        assembler = args['assembler']

        print(f'Assembling using {assembler}')
        assembly(current_file, assembler, output_dir, args['threads'])

        assembly_file = output_dir / 'assembly' / 'assembly.fasta'

    if steps['polishing']:
        print('Polishing using Medaka')
        polishing(current_file, assembly_file, output_dir, args['threads'])


if __name__ == '__main__':
    main()