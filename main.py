from argument_processing import arguments_process
from preprocessing import statistics, filtering
import sys

def main():
    args = arguments_process(sys.argv)

    print('Gathering filtering statistics')
    nanostatistics = statistics(args['file'])

    print('Filtering low quality reads')
    filtering(args['file'],nanostatistics)



if __name__ == '__main__':
    main()