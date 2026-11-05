from pathlib import Path



def arguments_process(commands):
    '''
    Processes commands from command line and saves them.

    :param Commands: commands from command line (list)
    :return: Dictionary of arguments
    '''

    arguments = {
        'file': None
    }

    i = 1

    while i < len(commands):
        arg = commands[i]
        if arg in ['-i','-input']:
            if i + 1 >= len(commands):
                raise ValueError(f"Missing value for {arg}")

            path = Path(commands[i+1])

            if not path.is_file():
                raise ValueError(f'File does not exist: {path}')
            if not is_fastq(path):
                raise ValueError(f'Not a fastq file: {path}')

            arguments['file'] = path
            i += 2
        else:
            raise ValueError(f'Unknown argument: {arg}')

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
