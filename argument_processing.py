import gzip
import shutil
import sys
from pathlib import Path



def arguments_process(commands):
    '''
    Processes commands from command line and saves them.

    :param Commands: commands from command line (list)
    :return: Dictionary of arguments
    '''

    arguments = {
        'file': None,
        'statistics': False,
        'trimming': False,
        'filtering': False,
        'assembly': False
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
                -s, --statistics    Calculate sequencing statistics.
                -t, --trimming      Trim adapters.
                -f, --filtering     Filter low-quality reads.
                -a, --assembly      Perform genome assembly.

            Examples:
                NanoBact -i sample.fastq.gz
                    Run the complete pipeline.

                NanoBact -i sample.fastq.gz -a
                    Run only genome assembly.

                NanoBact -i sample.fastq.gz -a -s
                    Run assembly and calculate statistics.

            For more information, visit the NanoBact documentation.
            Github: https://github.com/smedekj/NanoBact
            """
            print(help_text)
            sys.exit(0)



        elif arg in ['-i','--input']:

            if arguments['file'] is not None:
                raise ValueError('Input was specified more than once')

            if i + 1 >= len(commands):
                raise ValueError(f"Missing value for {arg}")


            input_path = Path(commands[i+1])

            if input_path.is_file():
                if not is_fastq(input_path):
                    raise ValueError(f'Not a fastq file: {input_path}')
                output_path = input_path

            elif input_path.is_dir():

                fastq_files = []
                for file in input_path.iterdir():
                    if file.is_file() and is_fastq(file):
                        if file.name == 'merged.fastq.gz':
                            continue
                        fastq_files.append(file)

                if not fastq_files:
                    raise ValueError(f'No fastq files found in {input_path}')

                merged_fastq = input_path / 'merged.fastq.gz'

                with gzip.open(merged_fastq,'wb') as output_file:
                    print('pred gzip')
                    for fastq in fastq_files:
                        if fastq.name.endswith('.gz'):
                            with gzip.open(fastq, 'rb') as infile:
                                shutil.copyfileobj(infile, output_file)
                        else:
                            with open(fastq, 'rb') as infile:
                                shutil.copyfileobj(infile, output_file)
                        print('cykl gzip')

                output_path = merged_fastq

            else:
                raise ValueError(f'Input path does not exist: {input_path}')

            arguments['file'] = output_path
            i += 2
            print(output_path)


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

            i += 1


        else:
            raise ValueError(f'Unknown argument: {arg}')

    if arguments['file'] is None:
        raise ValueError('Input file is required')

    steps = ['statistics', 'trimming', 'filtering', 'assembly']

    if not any(arguments[step] for step in steps):
        for step in steps:
            arguments[step] = True

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
