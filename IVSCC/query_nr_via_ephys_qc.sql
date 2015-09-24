with dendrite_type as 
(
select sts.specimen_id, st.name 
from specimen_tags_specimens sts
join specimen_tags st on sts.specimen_tag_id = st.id
where st.name like 'dendrite type%'
)
select sp.id as specimen_id, sp.name as specimen_name, dt.name as dendrite_type, str.name as region_info, wkf.storage_directory || wkf.filename as orca_path from specimens sp
left join structures str on sp.structure_id = str.id
left join dendrite_type dt on dt.specimen_id = sp.id
join ephys_roi_results err on err.id = sp.ephys_roi_result_id
join neuron_reconstructions nr on nr.specimen_id = sp.id
join well_known_files wkf on wkf.attachable_id = nr.id
where err.workflow_state = 'manual_passed'
and nr.superseded is false and nr.manual is true
order by sp.name;




select sp.name from specimens sp
join ephys_roi_results err on err.id = sp.ephys_roi_result_id
join neuron_reconstructions nr on nr.specimen_id = sp.id
where nr.superseded IS FALSE AND nr.manual IS TRUE
and err.workflow_state = 'manual_passed';