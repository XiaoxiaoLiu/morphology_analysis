select sp.id specimen_id, sp.name as specimen_name, a3d.*, wkf.storage_directory || wkf.filename as orca_path from specimens sp
left join structures str on sp.structure_id = str.id
join alignment3ds a3d on sp.alignment3d_id = a3d.id
join ephys_roi_results err on err.id = sp.ephys_roi_result_id
join neuron_reconstructions nr on nr.specimen_id = sp.id
join well_known_files wkf on wkf.attachable_id = nr.id
where err.workflow_state = 'manual_passed'
and nr.superseded is false and nr.manual is true
order by sp.name;


#### new one suggestsed from  wayne

select DISTINCT sp.id specimen_id, sp.name as specimen_name, a3d.*, wkf.storage_directory || wkf.filename as orca_path
from specimens sp
left join structures str on sp.structure_id = str.id
LEFT join alignment3ds a3d on sp.alignment3d_id = a3d.id
join ephys_roi_results err on err.id = sp.ephys_roi_result_id
join neuron_reconstructions nr on nr.specimen_id = sp.id and nr.superseded is false and nr.manual is true
join well_known_files wkf on wkf.attachable_id = nr.id
where err.workflow_state = 'manual_passed'
order by sp.name;