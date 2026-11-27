import gzip
import os
import shutil
import sys
from math import ceil
from pathlib import Path



def arguments_process(commands):
    '''
    Processes commands from command line and saves them.

    :param Commands: commands from command line (list)
    :return: Dictionary of arguments
    '''

    arguments = {
        'whole_pipeline': False,
        'file': None,
        'threads': None,
        'statistics': False,
        'trimming': False,
        'filtering': False,
        'assembly': False,
        'assembler': None,
        'polishing': False
    }

    i = 1

    while i < len(commands):
        arg = commands[i]
        if arg in ['-h', '--help']:
            help_text = """
            NanoBact - Tool for routine processing of nanopore data from bacterial isolates

            Usage:
                NanoBact -i <input> [options]

            Description:
                NanoBact is an automated pipeline for bacterial genome analysis.
                By default, all pipeline steps are executed.

            Required arguments:
                -i, --input         Input FASTQ file or directory containing FASTQ files.

            Optional arguments:
                -h, --help          Display this help message.

            Pipeline steps:
                -th, --threads      Number of threads used, Default: Available/2
                -s, --statistics    Calculate read statistics.
                -t, --trimming      Trim adapters.
                -f, --filtering     Filter low-quality reads.
                -a, --assembly      Perform genome assembly. Optionally specify assembler (flye or raven). Default: flye.
                -p, --polishing     Polishes assembled genome. Works only with assembly

            Examples:
                NanoBact -i sample.fastq.gz
                    Run the complete pipeline.

                NanoBact -i sample.fastq.gz -a
                    Run only genome assembly.

                NanoBact -i sample.fastq.gz -f -s
                    Run filtering and calculate statistics.

            For more information, visit the NanoBact documentation.
            Github: https://github.com/smedekj/NanoBact
            """
            print(help_text)
            sys.exit(0)



        elif arg in ['-i','--input']:

            if arguments['file'] is not None:
                raise ValueError('Input was specified more than once')

            if i + 1 >= len(commands) or commands[i+1].startswith('-'):
                raise ValueError(f"Missing value for {arg}")


            input_path = Path(commands[i+1])

            if input_path.is_file():
                if not is_fastq(input_path):
                    raise ValueError(f'Not a fastq file: {input_path}')
                output_path = input_path

            elif input_path.is_dir():

                fastq_files = []
                for file in sorted(input_path.iterdir()):
                    if file.is_file() and is_fastq(file):
                        if file.name == 'merged.fastq.gz':
                            continue
                        fastq_files.append(file)

                if not fastq_files:
                    raise ValueError(f'No fastq files found in {input_path}')

                merged_fastq = input_path / 'merged.fastq.gz'

                with gzip.open(merged_fastq,'wb') as output_file:
                    for fastq in fastq_files:
                        if fastq.name.endswith('.gz'):
                            with gzip.open(fastq, 'rb') as infile:
                                shutil.copyfileobj(infile, output_file)
                        else:
                            with open(fastq, 'rb') as infile:
                                shutil.copyfileobj(infile, output_file)

                output_path = merged_fastq

            else:
                raise ValueError(f'Input path does not exist: {input_path}')

            arguments['file'] = output_path
            i += 2
            print(output_path)

        elif arg in ['-th','--threads']:
            if arguments['threads'] is not None:
                raise ValueError('Threads was specified more than once')

            if i + 1 >= len(commands) or commands[i+1].startswith('-'):
                raise ValueError(f"Missing value for {arg}")

            try:
                threads = int(commands[i + 1])
            except ValueError:
                raise ValueError('Threads must be an integer')

            max_threads = (
                len(os.sched_getaffinity(0))
                if hasattr(os, 'sched_getaffinity')
                else (os.cpu_count() or 1)
            )

            if threads < 1:
                raise ValueError('Threads must be at least 1')

            if threads > max_threads:
                raise ValueError(f'Requested more threads than available ({max_threads})')

            arguments['threads'] = threads
            i += 2


        elif arg in ['-s','--statistics']:
            if arguments['statistics'] is True:
                raise ValueError('Statistics was specified more than once')

            arguments['statistics'] = True

            i += 1

        elif arg in ['-t','--trimming']:
            if arguments['trimming'] is True:
                raise ValueError('Trimming was specified more than once')

            arguments['trimming'] = True

            i += 1

        elif arg in ['-f','--filtering']:
            if arguments['filtering'] is True:
                raise ValueError('Filtering was specified more than once')

            arguments['filtering'] = True

            i += 1

        elif arg in ['-a','--assembly']:
            if arguments['assembly'] is True:
                raise ValueError('Assembly was specified more than once')

            arguments['assembly'] = True
            assemblers = ['flye', 'raven']

            if i + 1 >= len(commands) or commands[i+1].startswith('-'):
                arguments['assembler'] = 'flye'

                i += 1

            else:
                assembler = commands[i + 1]

                if assembler.lower() in assemblers:
                    arguments['assembler'] = assembler.lower()
                    i += 2
                else:
                    raise ValueError(f'Assembler "{assembler}" does not exist. Choose from: {assemblers}')

        elif arg in ['-p', '--polishing']:
            if arguments['polishing'] is True:
                raise ValueError('Polishing was specified more than once')

            arguments['polishing'] = True

            i += 1


        else:
            raise ValueError(f'Unknown argument: {arg}')

    if arguments['file'] is None:
        raise ValueError('Input file is required')

    if arguments['polishing'] and not arguments['assembly']:
        raise ValueError('Polishing works only with assembly')

    steps = ['statistics', 'trimming', 'filtering', 'assembly', 'polishing']

    if not any(arguments[step] for step in steps):
        for step in steps:
            arguments[step] = True
        arguments['whole_pipeline'] = True

    if arguments['threads'] is None:
        max_threads = (
            len(os.sched_getaffinity(0))
            if hasattr(os, 'sched_getaffinity')
            else (os.cpu_count() or 1)
        )
        arguments['threads'] = max(1, ceil(max_threads/2))


    return arguments



def is_fastq(file):
    '''
    Controls if file ends with a FASTQ suffix.

    :param file: file name (str)
    :return: bool
    '''
    suffixes = [suffix.lower() for suffix in file.suffixes]
    return (
        suffixes[-1:] == ['.fastq']
        or suffixes[-1:] == ['.fq']
        or suffixes[-2:] == ['.fastq', '.gz']
        or suffixes[-2:] == ['.fq', '.gz']
    )
