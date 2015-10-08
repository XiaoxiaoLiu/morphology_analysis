#!/bin/bash



#for p in  0.0 0.05 0.1 0.15 0.2 0.25 0.3
#   do python batch_fit.py  -l /home/xiaoxiaol/work/data/lims2/modified_neurons_p/r1.0_x1.0_y1.0_z1.0_p$p -o /home/xiaoxiaol/work/data/lims2/modified_neurons_p/ephys_fit_result_r1.0_x1.0_y1.0_z1.0_p$p &
#done


#for r in  0.5 0.7 0.9 1.0 1.1 1.3 1.5
for r in 0.6 0.8 1.2 1.4
	 do python batch_fit.py  -l /home/xiaoxiaol/work/data/lims2/modified_neurons/r$r\_x1.0_y1.0_z1.0_p0 -o /home/xiaoxiaol/work/data/lims2/modified_neurons/ephys_fit_result_r$r\_x1.0_y1.0_z1.0_p0 &
done





