#!/usr/bin/env python
import os
import sys
import math
import csv

MIN_POW = 0
MAX_POW = 12

DATADIR = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()

out_fn = os.path.join(DATADIR, 'result.csv')
HEADER_FN = 'header.csv'
IN_DELIMTER = ','
OUT_DELIMITER = ','

# READ header instead
HEADER = 'run,fem.assembleMatrix_avg_usr,fem.assembleMatrix_max_usr,fem.assembleMatrix_avg_wall,fem.assembleMatrix_max_wall,fem.assembleRHS_avg_usr,fem.assembleRHS_max_usr,fem.assembleRHS_avg_wall,fem.assembleRHS_max_wall,fem.solve_avg_usr,fem.solve_max_usr,fem.solve_avg_wall,fem.solve_max_wall,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero_avg_usr,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero_max_usr,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero_avg_wall,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero_max_wall,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part_avg_usr,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part_max_usr,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part_avg_wall,msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part_max_wall,msfem.all_avg_usr,msfem.all_max_usr,msfem.all_avg_wall,msfem.all_max_wall,msfem.assembleMatrix_avg_usr,msfem.assembleMatrix_max_usr,msfem.assembleMatrix_avg_wall,msfem.assembleMatrix_max_wall,msfem.assembleRHS_avg_usr,msfem.assembleRHS_max_usr,msfem.assembleRHS_avg_wall,msfem.assembleRHS_max_wall,msfem.localproblemsolver.assemble_all_avg_usr,msfem.localproblemsolver.assemble_all_max_usr,msfem.localproblemsolver.assemble_all_avg_wall,msfem.localproblemsolver.assemble_all_max_wall,msfem.solveCoarse_avg_usr,msfem.solveCoarse_max_usr,msfem.solveCoarse_avg_wall,msfem.solveCoarse_max_wall,msfem.subgrid_list_avg_usr,msfem.subgrid_list_max_usr,msfem.subgrid_list_avg_wall,msfem.subgrid_list_max_wall,msfem.subgrid_list.create_avg_usr,msfem.subgrid_list.create_max_usr,msfem.subgrid_list.create_avg_wall,msfem.subgrid_list.create_max_wall,msfem.subgrid_list.create.finalize_avg_usr,msfem.subgrid_list.create.finalize_max_usr,msfem.subgrid_list.create.finalize_avg_wall,msfem.subgrid_list.create.finalize_max_wall,msfem.subgrid_list.identify_avg_usr,msfem.subgrid_list.identify_max_usr,msfem.subgrid_list.identify_avg_wall,msfem.subgrid_list.identify_max_wall,none.local_problem_solution_avg_usr,none.local_problem_solution_max_usr,none.local_problem_solution_avg_wall,none.local_problem_solution_max_wall'
IN_FIELDS = ['procs'] + HEADER.split(IN_DELIMTER)

OUT_FIELDS = [f for f in IN_FIELDS if 'usr' not in f and 'run' not in f]

letters = [chr(c) for c in range(ord('A'), ord('Z')+1)] 
LETTERS = letters + ['A'+l  for l in letters] + ['B'+l  for l in letters]#itertools.prduct, letters, lterrs) ?

mean_rows = []
with open(out_fn,'wb') as outfile:
    out_row_count = 2
    out_csv = csv.DictWriter(outfile, OUT_FIELDS, delimiter=OUT_DELIMITER )
    out_csv.writeheader()
    result_count = 0    
    pows = []
    for procs in [int(math.pow(2, i)) for i in range(MIN_POW, MAX_POW+1)]:
         
        for result_count,sub in enumerate([f for f in os.listdir(DATADIR) if os.path.isdir(os.path.join(DATADIR,f)) 
                                                                    and '-%d.'%(procs) in f]):
            with open( os.path.join(DATADIR,sub,'profiler.csv')) as infile:
                csv_reader = csv.DictReader(infile, delimiter=IN_DELIMTER)
                in_row = csv_reader.next()
                out_row = {key: in_row[key] for key in in_row.keys() if key in OUT_FIELDS} 
                out_row['procs'] = procs
                out_csv.writerow(out_row)
                out_row_count += 1
        if result_count > 0:
	    pows.append(procs)
                
        # avg row
        if result_count > 0:
            avg_row = { key: '=MITTELWERT({0}{1}:{0}{2})/1000'.format(LETTERS[i], out_row_count - result_count - 1, out_row_count - 1) 
                       for i,key in enumerate(OUT_FIELDS)}
            avg_row['procs'] = 'Mean'
            out_csv.writerow(avg_row)
            mean_rows.append(out_row_count)
            out_row_count += 1
            result_count = 0
    for idx in range(1,len(mean_rows)):
                fac_row = { key: '={0}{1}/{0}{2}'.format(LETTERS[i], mean_rows[idx-1],mean_rows[idx]) 
                           for i,key in enumerate(OUT_FIELDS)}
                fac_row['procs'] = 'Factor {1}/{0}'.format(pows[idx-1], pows[idx])
                out_csv.writerow(fac_row)
                out_row_count += 1
            
            
                
print('DONE')
