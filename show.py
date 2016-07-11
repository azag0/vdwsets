#!/usr/bin/env python3
import sys
import json
from tabulate import tabulate

print(tabulate(
    [{**row, 'idx': idx+1} for idx, row in enumerate(json.load(open(sys.argv[1])))],
    headers='keys'
))
