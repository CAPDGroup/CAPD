import os
import argparse

def find_header_files(directory):
    header_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(('.h', '.hpp')):
                header_files.append(os.path.join(root, file))
    return header_files

def main():
    parser = argparse.ArgumentParser(description="Find .h and .hpp files in a directory")
    parser.add_argument("directory", type=str, help="Directory to search")
    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print("Invalid directory")
        return

    header_files = find_header_files(args.directory)
    for file in header_files:
        print(file)

if __name__ == "__main__":
    main()
