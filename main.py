from argument_processing import arguments_process
from preprocessing import statistics, trimming, filtering
import sys
import os
from pathlib import Path

def main():
    output_dir = Path('output')
    output_dir.mkdir(parents=True, exist_ok=True)
    args = arguments_process(sys.argv)

    if args['statistics']:
        print('Gathering filtering statistics')
        statistics(args['file'], output_dir)

    if args['trimming']:
        print('Trimming reads')
        trimming(args['file'], output_dir)

    if args['filtering']:
        print('Filtering low quality reads')
        filtering(args['file'], output_dir)



if __name__ == '__main__':
    main()