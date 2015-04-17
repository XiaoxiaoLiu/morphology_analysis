#!/usr/bin/env python

import psycopg2

def get_ephys_id_from_lims(specimen_id):
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()
    cur.execute("SELECT s.name, s.ephys_roi_result_id FROM specimens s WHERE s.id=%s", (specimen_id,))
    result = cur.fetchone()
    if not result:
        print "Could not find specimen ", specimen
        return (None, None, None)

    print "Specimen: " + result[0]
    print "EphysRoiResult: " + str(result[1])
    
    cur.close()
    conn.close()
    
    return result

def wet_specimen_info_from_lims(specimen):
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()
    cur.execute("SELECT s.name, s.ephys_roi_result_id, s.id FROM specimens s WHERE s.name LIKE %s", ('%' + specimen,))
    result = cur.fetchone()
    if not result:
        print "Could not find specimen ", specimen
        return (None, None, None)

    print "Specimen: " + result[0]
    print "EphysRoiResult: " + str(result[1])
    
    cur.close()
    conn.close()
    
    return result
    
def get_orca_path_from_lims(ephys_roi_result):
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()

    cur.execute("SELECT f.filename, f.storage_directory FROM well_known_files f \
             WHERE f.attachable_type = 'EphysRoiResult' AND f.attachable_id = %s AND f.filename LIKE '%%orca'", (ephys_roi_result,))
    result = cur.fetchone()
    
    if not result:
        print "Cannot find orca file"
        return ""
    
    orca_path = result[1] + result[0]
    print "Orca file: " + orca_path
    
    cur.close()
    conn.close()
    return orca_path
    
    
def get_swc_from_lims(specimen_id):
    conn = psycopg2.connect('host=limsdb2 dbname=lims2 user=limsreader password=limsro')
    cur = conn.cursor()

    SQL = "SELECT f.filename, f.storage_directory FROM \
     neuron_reconstructions n JOIN well_known_files f ON n.id = f.attachable_id \
     AND n.specimen_id = %s AND n.manual AND NOT n.superseded"
    cur.execute(SQL, (specimen_id,))
    result = cur.fetchone()
    swc_filename = result[0]
    swc_path = result[1] + result[0]
    print "SWC file: " + swc_path
    
    cur.close()
    conn.close()
    return swc_filename, swc_path
