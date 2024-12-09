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
            if self.__match_ifndef(line):
                self.state = State.INCLUSION_GUARD_IFNDEF
            elif self.__match_online_multiline_comment(line):
                pass
            elif self.__match_comment_or_empty_line(line):
                pass
            else:
                self.state = State.ERROR
        
        elif self.state == State.INCLUSION_GUARD_IFNDEF:
            if line == f'#define {self.header_guard}':
                self.state = State.INCLUSION_GUARD_DEFINE
        
        elif self.state == State.INCLUSION_GUARD_DEFINE:
            pass
        
        else:
            self.state = State.ERROR

        return line


    def __match_ifndef(self, line : str) -> bool:
        result = re.match(r'#ifndef\s+([A-Za-z0-9_]*)', line)
        if result:
            self.header_guard = result.group(1)
            trace.debug(self.header_guard)
            return True
        else:
            return False


    def __match_comment_or_empty_line(self, line : str) -> bool:
        return re.match(r'^\s*(//.*)?$', line)
    
    
    def __match_online_multiline_comment(self, line : str) -> bool:
        return re.match(r'^\s*/\*.*\*/\s*$', line)


def process(path : str) -> HeaderFile:
    trace.debug(f'Processing {path}')

    header_file = HeaderFile(path)

    file_tmp = path + '.tmp'
    with open(path, 'r', newline='') as ifs:
        with open(file_tmp, 'w', newline='') as ofs:

            try:
                while True:
                    line = ifs.readline()

                    line = header_file.process_line(line)

                    if line == '':
                        break

                    ofs.write(line)
                
            except Exception as e:
                trace.error(e)
                return header_file
            finally:
                pass
    
    status = os.stat(path)
    permissions = status.st_mode & 0o777
    os.remove(path)
    os.rename(file_tmp, path)
    os.chmod(path, permissions)
    return header_file


def check_unique_header_guard_names(header_files : List[HeaderFile]):

    header_guards = []
    for h in header_files:
        if h.header_guard in header_guards:
            trace.warning(f'Duplicate header guard label: {h.header_guard} in {h.path}')
        else:
            header_guards.append(h.header_guard)


if __name__ == "__main__":

    logging.basicConfig(level=logging.DEBUG)

    parser = argparse.ArgumentParser(description="Find .h and .hpp files in a directory")
    parser.add_argument('dirs', nargs='+', help="Directory/Directories to search")
    args = parser.parse_args()

    if all( [os.path.isdir(dir) for dir in args.dirs] ):

        header_file_paths = sum([find_header_files(dir) for dir in args.dirs], [])
        header_files = [process(path) for path in header_file_paths]

        trace.info(f'Files analyzed: {len(header_files)}')
        trace.info(f'Files with err: {sum([header_file.state == State.ERROR for header_file in header_files])}')

        for header_file in header_files:
            if header_file.state == State.ERROR:
                trace.warning(header_file.path)

        check_unique_header_guard_names(header_files)
        

    else:
        print("Invalid directory")
        
