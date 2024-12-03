import os
import argparse
import logging
import re
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


class HeaderFile:
    def __init__(self, path):
        self.path = path
        self.header_guard = ''
        self.state = State.PREAMBLE

    def process_line(self, line : str) -> str:

        if self.state == State.PREAMBLE:
            result = re.match(r'#ifndef ([A-Z_]*)', line)
            if result:
                self.header_guard = result.group(1)
                trace.info(self.header_guard)
                self.state = State.INCLUSION_GUARD_IFNDEF
        
        elif self.state == State.INCLUSION_GUARD_IFNDEF:
            if line == f'#define {self.header_guard}':
                self.state = State.INCLUSION_GUARD_DEFINE
        
        elif self.state == State.INCLUSION_GUARD_DEFINE:
            return line
        
        elif self.state == State.INCLUSION_GUARD_ENDIF:
            return line
        
        else:
            self.state = State.ERROR

        return line


def process(path : str) -> HeaderFile:
    trace.info(f'Processing {path}')

    header_file = HeaderFile(path)

    file_tmp = path + '.tmp'
    with open(path, 'r') as ifs:
        with open(file_tmp, 'w') as ofs:

            try:
                while True:
                    line = ifs.readline()

                    # line = header_file.process_line(line)

                    if line == '':
                        break

                    ofs.write(line)
                
            except:
                return header_file
            finally:
                pass
    
    status = os.stat(path)
    permissions = status.st_mode & 0o777
    os.remove(path)
    os.rename(file_tmp, path)
    os.chmod(path, permissions)
    return header_file


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Find .h and .hpp files in a directory")
    parser.add_argument("directory", type=str, help="Directory to search")
    args = parser.parse_args()

    if os.path.isdir(args.directory):
        header_file_paths = find_header_files(args.directory)
        for path in header_file_paths:
            header_file = process(path)
            print(header_file.state)
    else:
        print("Invalid directory")
        
