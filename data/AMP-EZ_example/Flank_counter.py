import sys
import datetime
import numpy as np
import pandas as pd
from numpy import savetxt

name_list = []
sequence_list = []

query_file = sys.argv[1] ##eg. "Bat_ACE2_aligned.fasta"

with open(query_file, 'r') as datafile:
    for line in datafile:
        if re.search(sys.argv[0], line):
            temp_sequence_variable = temp_sequence_variable + line.strip()
        if line[0] == ">":
            name_list.append(line.strip())