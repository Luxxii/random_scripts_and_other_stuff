import sys
import re

regex_pg_header = r"[A-Z0-9-_]+?\(.*?\)"

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file, "r") as in_fasta, open(out_file, "w") as out_fasta:
    write_out = False

    for line in in_fasta:
        if line.startswith(">"):
            desc = line.split("|", 2)[-1]

            line = line.replace(" ", "")
            matches_regex = re.finditer(regex_pg_header, line, re.MULTILINE)
            matches = [x.group() for x in matches_regex]

            if len(matches) == 1:
                write_out = True
            else:
                write_out = False

        if write_out:
            out_fasta.write(line)
