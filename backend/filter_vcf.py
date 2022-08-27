import sys
import json

MIN_MAX_AF = 0.001 # 0.1%

for line in sys.stdin:
    if line.startswith('#'):
        print(line)
        continue
    ele = line.strip().split()
    if ele[6] != "PASS":
        continue
    tmp = [x for x in ele[-1].split(";") if x.startswith('AF')]
    info = {x.split('=')[0]: float(x.split('=')[1]) for x in tmp}
    if max(info.values()) > MIN_MAX_AF:
        print(line)

