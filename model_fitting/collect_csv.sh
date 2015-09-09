#!/bin/bash

for p in  0.0 0.05 0.1 0.15 0.2 0.25 0.3
	do python collect_result_into_csv.py  -o   ~/work/data/lims2/modified_neurons_p/ephys_fit_result_r1.0_x1.0_y1.0_z1.0_p$p
done

for r in  0.5 0.7 0.9 1.0 1.1 1.3 1.5
	do python collect_result_into_csv.py  -o   ~/work/data/lims2/modified_neurons/ephys_fit_result_r$r\_x1.0_y1.0_z1.0_p0
done
