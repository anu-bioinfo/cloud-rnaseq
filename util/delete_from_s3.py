#!/usr/bin/env python
import re
import argparse
import sys
import subprocess

def main():

    parser = argparse.ArgumentParser(
        description = "delete list of files listed in STDIN under the S3_PATH. \
            Use -f flag to force delete from aws. otherwise, this just prints out command")
    parser.add_argument('-s', action="store", dest='s3_path', default=False)
    parser.add_argument('-f', action="store_true", dest='force', default=False)
    results = parser.parse_args()
    if results.s3_path:
        count = 0
        for line in sys.stdin:
            count += 1
            line = line.rstrip()
            m = re.match(".*?([^ ]*)$", line)
            if m:
                command = "aws s3 rm %s/%s &" % (results.s3_path, m.group(1))
                if results.force:
                    subprocess.check_output(command, shell=True)
                else:
                    print command
            if count % 6 == 0:
                command = "sleep 5"
                if results.force:
                    output = subprocess.check_output(command, shell=True)
                else:
                    print command

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
