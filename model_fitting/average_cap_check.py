#!/usr/bin/env python

from utilities.lims_orca_utils import *
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
    parser.add_argument('specimen_id')
    parser.add_argument('-o', dest='outputDir', default='./result')
    parser.add_argument('--orca')
    parser.add_argument('--noshow',dest="noshow",default=0)
    args = parser.parse_args()
    
    outputDir = args.outputDir
    output_png_path = outputDir + '/test_pulse_plots'
    if not os.path.isdir(output_png_path):
          os.system('mkdir -p '+ output_png_path)
    
    output_data_path = outputDir + '/data'
    if not os.path.isdir(output_data_path):
          os.system('mkdir -p '+ output_data_path)
    
    
    NOSHOW = args.noshow
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()
    specimen_id = args.specimen_id;
    specimen_name, ephys_roi_result_id = lims_orca_utils.get_ephys_id_from_lims(args.specimen_id)

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
    with open(output_data_path+"/" + args.specimen_id + '_upbase.dat', 'w') as f:
        np.savetxt(f, np.column_stack((t, grand_up)))
    with open(output_data_path +"/" + args.specimen_id + '_downbase.dat', 'w') as f:
        np.savetxt(f, np.column_stack((t, grand_down)))

    # Plot the traces
    axes[1].plot(t, grand_up,'r')
    axes[1].plot(t, -grand_down,'b')
    axes[1].set_xlim(1, 4)
    plt.suptitle(args.specimen_id)
    plt.savefig(output_png_path +"/" + args.specimen_id + "_testpulseplot.png", bbox_inches="tight")
    if not NOSHOW:
      plt.show()
    
    plt.figure()
    plt.plot(t, grand_up,'r')
    plt.plot(t, -grand_down,'b')
    plt.ylim(0.01, 10)
    plt.xlim(0, 100)
    plt.yscale('log')
    if not NOSHOW:
      plt.show()
    
