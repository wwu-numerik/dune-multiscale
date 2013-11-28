#!/usr/bin/env python
import os
import sys
import math
import csv

MIN_POW = 0
MAX_POW = 12

DATADIR = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()

HEADER_FN = os.path.join(DATADIR, 'header.csv')
IN_DELIMTER = ','
OUT_DELIMITER = ','

# READ header instead
HEADER = open(HEADER_FN, 'rb').readline().strip()
IN_FIELDS = ['procs'] + HEADER.split(IN_DELIMTER)
letters = [chr(c) for c in range(ord('A'), ord('Z')+1)] 
LETTERS = letters + ['A'+l  for l in letters] + ['B'+l  for l in letters]#itertools.prduct, letters, lterrs) ?

def result_file(exclude, out_fn):

    out_fields = [f for f in IN_FIELDS if exclude not in f and 'run' not in f]
    mean_rows = []
    with open(out_fn, 'wb') as outfile:
        out_row_count = 2
        out_csv = csv.DictWriter(outfile, out_fields, delimiter=OUT_DELIMITER)
        out_csv.writeheader()
        result_count = 0
        pows = []
        for procs in [int(math.pow(2, i)) for i in range(MIN_POW, MAX_POW + 1)]:

            for result_count, sub in enumerate([f for f in os.listdir(DATADIR) if os.path.isdir(os.path.join(DATADIR, f))
                                                                        and '-%d.' % (procs) in f]):
                with open(os.path.join(DATADIR, sub, 'profiler.csv')) as infile:
                    csv_reader = csv.DictReader(infile, delimiter=IN_DELIMTER)
                    in_row = csv_reader.next()
                    out_row = {key: in_row[key] for key in in_row.keys() if key in out_fields}
                    out_row['procs'] = procs
                    out_csv.writerow(out_row)
                    out_row_count += 1
            if result_count > 0:
                pows.append(procs)

            # avg row
            if result_count > 0:
                avg_row = { key: '=MITTELWERT({0}{1}:{0}{2})/1000'.format(LETTERS[i], out_row_count - result_count - 1,
                                                                           out_row_count - 1)
                           for i, key in enumerate(out_fields)}
                avg_row['procs'] = 'Mean'
                out_csv.writerow(avg_row)
                mean_rows.append(out_row_count)
                out_row_count += 1
                result_count = 0

        def w_row(base_idx, idx):
                fac_row = { key: '={0}{1}/{0}{2}'.format(LETTERS[i], mean_rows[base_idx], mean_rows[idx])
                           for i, key in enumerate(out_fields)}
                fac_row['procs'] = 'Factor {1}/{0}'.format(pows[base_idx], pows[idx])
                out_csv.writerow(fac_row)

        for l in range(1, len(mean_rows)):
            w_row(l - 1, l)
            out_row_count += 1
        for l in range(1, len(mean_rows)):
            w_row(0, l)
            out_row_count += 1

result_file('wall', out_fn=os.path.join(DATADIR, 'usr_result.csv'))
result_file('usr', out_fn=os.path.join(DATADIR, 'wall_result.csv'))
                
print('DONE')
