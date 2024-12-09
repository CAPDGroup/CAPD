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
    def __init__(self, path : str):
        self.path = path
        self.header_guard = ''
        self.header_guard_new = ''
        self.state = State.PREAMBLE
        self.endif_counter = 0


    def process_line(self, line : str) -> str:

        if self.state == State.PREAMBLE:
            if self.__match_ifndef(line):
                self.state = State.INCLUSION_GUARD_IFNDEF
                self.endif_counter += 1
                # line = f'#ifndef {self.header_guard_new}'
            elif self.__match_online_multiline_comment(line):
                pass
            elif self.__match_comment_or_empty_line(line):
                pass
            else:
                self.state = State.ERROR
        
        elif self.state == State.INCLUSION_GUARD_IFNDEF:
            if self.__match_define(line):
                self.state = State.INCLUSION_GUARD_DEFINE
                # line = f'#define {self.header_guard_new}'
        
        elif self.state == State.INCLUSION_GUARD_DEFINE:
            if self.__match_if_or_ifdef(line):
                self.endif_counter += 1
            elif self.__match_endif(line):
                self.endif_counter -= 1

            if self.endif_counter == 0:
                self.state = State.INCLUSION_GUARD_ENDIF
                # line = f'#endif // {self.header_guard_new}'
        
        elif self.state == State.INCLUSION_GUARD_ENDIF:
            if self.__match_comment_or_empty_line(line):
                pass
            else:
                print(line)
                self.state = State.ERROR

        else:
            self.state = State.ERROR

        return line


    def __match_ifndef(self, line : str) -> bool:
        result = re.match(r'^#ifndef\s+([A-Za-z0-9_]*)', line)
        if result:
            self.header_guard = result.group(1)
            trace.debug(self.header_guard)

            self.header_guard_new = self.__process_header_guard(self.header_guard)
            trace.debug(self.header_guard_new)

            return True
        else:
            return False


    def __match_define(self, line : str) -> bool:
        return re.match(r'^#define\s+' + self.header_guard, line)
    

    def __match_if_or_ifdef(self, line : str) -> bool:
        return re.match(r'^\s*#if\(?\s+', line) or re.match(r'^\s*#ifdef\s+', line) or re.match(r'^\s*#ifndef\s+', line)
    

    def __match_endif(self, line : str) -> bool:
        return re.match(r'^\s*#endif\s*', line)


    def __match_comment_or_empty_line(self, line : str) -> bool:
        return re.match(r'^\s*(//.*)?$', line)
    
    
    def __match_online_multiline_comment(self, line : str) -> bool:
        return re.match(r'^\s*/\*.*\*/\s*$', line)
    

    def __process_header_guard(self, header_guard : str):

        # remove underscores
        header_guard = header_guard.strip('_')

        # set to uppercase
        header_guard = header_guard.upper()

        # ensure start with single CAPD_
        if header_guard.startswith('CAPD_CAPD_'):
            header_guard = header_guard[5:]
        elif not header_guard.startswith('CAPD_'):
            header_guard = 'CAPD_' + header_guard

        #ensure ending with _H or _HPP
        if self.path.endswith('.h'):
            if header_guard.endswith('_HPP'):
                header_guard = header_guard[:-2]
            elif not header_guard.endswith('_H'):
                header_guard = header_guard + '_H'
        elif self.path.endswith('.hpp'):
            if header_guard.endswith('_H'):
                header_guard = header_guard + 'PP'
            elif not header_guard.endswith('_HPP'):
                header_guard = header_guard + '_HPP'
        
        return header_guard


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
    
    trace.debug(header_file.state)

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
        trace.info(f'Files with err: {sum([header_file.state != State.INCLUSION_GUARD_ENDIF for header_file in header_files])}')

        for header_file in header_files:
            if header_file.state != State.INCLUSION_GUARD_ENDIF:
                trace.warning(header_file.path)

        check_unique_header_guard_names(header_files)
        

    else:
        print("Invalid directory")
        
