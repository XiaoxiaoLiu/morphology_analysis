#!/bin/sh

/home/nathang/analysis/neuron_passive_fit.py $1 $2 $3
# /home/nathang/analysis/neuron_passive_fit2.py $1 $2 $3
/home/nathang/analysis/neuron_passive_fit_elec.py $1 $2 $3 $4 1
/home/nathang/analysis/neuron_passive_fit_elec.py $1 $2 $3 $4 10
/home/nathang/analysis/neuron_passive_fit_elec.py $1 $2 $3 $4 20


