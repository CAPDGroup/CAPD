import os
import argparse
import logging
from typing import List
from enum import Enum

trace = logging.getLogger(__name__)

def find_header_files(directory : str) -> List[str]:
    header_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(('.h', '.hpp')):
                header_files.append(os.path.join(root, file))
    return header_files


class State(Enum):
    PREAMBLE = 0
    INCLUSION_GUARD_IFNDEF = 1
    INCLUSION_GUARD_DEFINE = 2
    INCLUSION_GUARD_ENDIF = 3
    ERROR = -1


def process(file : str):
    trace.info(f'Processing {file}')

    file_tmp = file + '.tmp'
    with open(file, 'r') as ifs:
        with open(file_tmp, 'w') as ofs:

            try:
                state = State.PREAMBLE

                while True:
                    line = ifs.readline()

                    if line == '':
                        break

                    ofs.write(line)
                
            except:
                return None
            finally:
                pass
    
    os.remove(file)
    os.rename(file_tmp, file)
    os.chmod(file, 0o755)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Find .h and .hpp files in a directory")
    parser.add_argument("directory", type=str, help="Directory to search")
    args = parser.parse_args()

    if os.path.isdir(args.directory):
        header_files = find_header_files(args.directory)
        for file in header_files:
            process(file)
    else:
        print("Invalid directory")
        
