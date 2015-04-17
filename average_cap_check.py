#!/usr/bin/env python

import argparse
import psycopg2
import h5py
import sys
import webbrowser
import matplotlib.pyplot as plt
import lims_orca_utils
import numpy as np

def get_sweep_from_orca(orca_path, sweep_num):
    exp_f = h5py.File(orca_path, 'r')
    exp_tag = 'epochs/Experiment_%d' % sweep_num
    if exp_tag not in exp_f:
        return (None, None, None)
    root = exp_f[exp_tag]
    v = root["response/timeseries/data"][0:]
    v *= 1e3 # to mV
    i = root["stimulus/timeseries/data"][0:]
    i *= 1e12 # to pA
    dt = 1.0 / root["stimulus/timeseries/sampling_rate"].value
    t = np.arange(0, len(v)) * dt
    return (v, i, t)    

def get_cap_check_indices(i):
    # Assumes that there is a test pulse followed by the stimulus pulses (downward first)
    di = np.diff(i)
    up_idx = np.flatnonzero(di > 0)
    down_idx = np.flatnonzero(di < 0)
    
    return up_idx[2::2], down_idx[1::2]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='analyze cap check sweep')
    parser.add_argument('specimen')
    parser.add_argument('--orca')
    args = parser.parse_args()
    
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()
    specimen_name, ephys_roi_result_id, specimen_id = lims_orca_utils.get_specimen_info_from_lims(args.specimen)

    if args.orca:
        orca_path = args.orca
    else: 
        orca_path = lims_orca_utils.get_orca_path_from_lims(ephys_roi_result_id)

    cur.execute("SELECT sw.sweep_number, sw.bridge_balance_mohm FROM ephys_sweeps sw JOIN ephys_stimuli stim \
         ON stim.id = sw.ephys_stimulus_id \
         WHERE sw.specimen_id = %s AND stim.description LIKE '%%C1SQCAPCHK%%'",
         (specimen_id,))
    
    sweeps_data = cur.fetchall()
    sweeps = [s[0] for s in sweeps_data]
    bridge_balances = [s[1] for s in sweeps_data]
    print "bridge avg {:.2f}".format(np.array(bridge_balances).mean())
    
    if len(sweeps) == 0:
        url = 'http://bouton/trace_viewer/?path=' + orca_path
        webbrowser.get('firefox').open(url)
        sweeps = map(int, raw_input("Enter the sweep numbers to analyze: ").split())

    cur.close()
    conn.close()

    up_avgs = {}
    down_avgs = {}
    
    plt.style.use('ggplot')
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    initialized = False
    for idx, s in enumerate(sweeps):
        v, i, t = get_sweep_from_orca(orca_path, s)
        if v is None:
            continue
        up_idxs, down_idxs = get_cap_check_indices(i)

        axes[0].plot(t, v)


        down_idx_interval = down_idxs[1] - down_idxs[0]
        skip_count = 0
        for j in range(len(up_idxs)):
            if j == 0:
                avg_up = v[(up_idxs[j] - 400):down_idxs[j + 1]]
                avg_down = v[(down_idxs[j] - 400):up_idxs[j]]
            elif j == len(up_idxs) - 1:
                avg_up = avg_up + v[(up_idxs[j] - 400):-2]
                avg_down = avg_down + v[(down_idxs[j] - 400):up_idxs[j]]
            else:
                avg_up = avg_up + v[(up_idxs[j] - 400):down_idxs[j + 1]]
                avg_down = avg_down + v[(down_idxs[j] - 400):up_idxs[j]]
        avg_up /= len(up_idxs) - skip_count
        avg_down /= len(up_idxs) - skip_count
        if not initialized:
            grand_up = avg_up - avg_up[0:400].mean()
            grand_down = avg_down - avg_down[0:400].mean()
            initialized = True
        else:
            grand_up = grand_up + (avg_up - avg_up[0:400].mean())
            grand_down = grand_down + (avg_down - avg_down[0:400].mean())

    # Show the test pulse
    axes[0].set_xlim(0, 0.04)
    
    grand_up /= len(sweeps)
    grand_down /= len(sweeps)
    
    t = 0.005 * np.arange(len(grand_up))
    with open("./results/data/" + args.specimen + '_upbase.dat', 'w') as f:
        np.savetxt(f, np.column_stack((t, grand_up)))
    with open("./results/data/" + args.specimen + '_downbase.dat', 'w') as f:
        np.savetxt(f, np.column_stack((t, grand_down)))

    # Plot the traces
    axes[1].plot(t, grand_up)
    axes[1].plot(t, -grand_down)
    axes[1].set_xlim(1, 4)
    plt.suptitle(args.specimen)
    plt.savefig("./results/test_pulse_plots/" + args.specimen + "_testpulseplot.png", bbox_inches="tight")
    plt.show()
    
    plt.figure()
    plt.plot(t, grand_up)
    plt.plot(t, -grand_down)
    plt.ylim(0.01, 10)
    plt.xlim(0, 100)
    plt.yscale('log')
    plt.show()
