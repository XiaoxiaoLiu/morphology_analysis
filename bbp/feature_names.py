
meta_names=['swc_file_name',
    'layer',
    'm-type',
    'shape_type',
    'dendrite_type',
    'detailed_type_name'
    ]


vaa3d_tree_features =[
    'average_bifurcation_angle_local',
    'average_bifurcation_angle_remote',
    'average_contraction',
    'average_diameter',
    'average_fragmentation',
    'average_parent_daughter_ratio',
    'max_branch_order',
    'max_euclidean_distance',
    'max_path_distance',
    'moment1',
    'moment10',
    'moment11',
    'moment12',
    'moment13',
    'moment14',
    'moment2',
    'moment3',
    'moment4',
    'moment5',
    'moment6',
    'moment7',
    'moment8',
    'moment9',
    'nodes_over_branches',
    'number_of_bifurcations',
    'number_of_branches',
    'number_of_nodes',
    'number_of_stems',
    'number_of_tips',
    'overall_depth',
    'overall_height',
    'overall_width',
    'soma_surface',
    'total_length',
    'total_surface',
    'total_volume'
    ]

non_gmi_features =[
    'average_bifurcation_angle_local',
    'average_bifurcation_angle_remote',
    'average_contraction',
    'average_diameter',
    'average_fragmentation',
    'average_parent_daughter_ratio',
    'max_branch_order',
    'max_euclidean_distance',
    'max_path_distance',
    'nodes_over_branches',
    'number_of_bifurcations',
    'number_of_branches',
    'number_of_nodes',
    'number_of_stems',
    'number_of_tips',
    'overall_depth',
    'overall_height',
    'overall_width',
    'soma_surface',
    'total_length',
    'total_surface',
    'total_volume'
    ]

apical_kg_features =[
    'kg_branch_centroid_distance_z_apical',
    'kg_branch_mean_from_centroid_z_apical',
    'kg_branch_stdev_from_centroid_z_apical',
    'kg_centroid_over_farthest_branch_apical',
    'kg_centroid_over_farthest_neurite_apical',
    'kg_centroid_over_radial_dist_apical', 'kg_mean_over_centroid',
    'kg_mean_over_farthest_branch_apical',
    'kg_mean_over_farthest_neurite_apical',
    'kg_mean_over_radial_dist_apical', 'kg_mean_over_stdev',
    'kg_num_branches_over_radial_dist_apical',
    'kg_num_outer_apical_branches', 'kg_outer_mean_from_center_z_apical',
    'kg_outer_mean_over_stdev', 'kg_outer_stdev_from_center_z_apical',
    'kg_peak_over_moment_z_apical', 'kg_radial_dist_over_moment_z_apical',
    'kg_soma_depth'
    ]

#71
BBP_features=['bb_first_moment_apical',
    'bb_first_moment_basal',
    'bb_first_moment_dendrite',
    'bb_first_moment_x_apical',
    'bb_first_moment_x_basal',
    'bb_first_moment_x_dendrite',
    'bb_first_moment_y_apical',
    'bb_first_moment_y_basal',
    'bb_first_moment_y_dendrite',
    'bb_first_moment_z_apical',
    'bb_first_moment_z_basal',
    'bb_first_moment_z_dendrite',
    'bb_max_branch_order_apical',
    'bb_max_branch_order_basal',
    'bb_max_branch_order_dendrite',
    'bb_max_path_length_apical',
    'bb_max_path_length_basal',
    'bb_max_path_length_dendrite',
    'bb_max_radial_distance_apical',
    'bb_max_radial_distance_basal',
    'bb_max_radial_distance_dendrite',
    'bb_mean_trunk_diameter_apical',
    'bb_mean_trunk_diameter_basal',
    'bb_mean_trunk_diameter_dendrite',
    'bb_first_moment_x_basal',
    'bb_first_moment_x_dendrite',
    'bb_first_moment_y_apical',
    'bb_first_moment_y_basal',
    'bb_first_moment_y_dendrite',
    'bb_first_moment_z_apical',
    'bb_first_moment_z_basal',
    'bb_first_moment_z_dendrite',
    'bb_max_branch_order_apical',
    'bb_max_branch_order_basal',
    'bb_max_branch_order_dendrite',
    'bb_max_path_length_apical',
    'bb_max_path_length_basal',
    'bb_max_path_length_dendrite',
    'bb_max_radial_distance_apical',
    'bb_max_radial_distance_basal',
    'bb_max_radial_distance_dendrite',
    'bb_mean_trunk_diameter_apical',
    'bb_mean_trunk_diameter_basal',
    'bb_mean_trunk_diameter_dendrite',
    'bb_number_branches_apical',
    'bb_number_branches_basal',
    'bb_number_branches_dendrite',
    'bb_number_neurites_apical',
    'bb_number_neurites_basal',
    'bb_number_neurites_dendrite',
    'bb_second_moment_apical',
    'bb_second_moment_basal',
    'bb_second_moment_dendrite',
    'bb_second_moment_x_apical',
    'bb_second_moment_x_basal',
    'bb_second_moment_x_dendrite',
    'bb_second_moment_y_apical',
    'bb_second_moment_y_basal',
    'bb_second_moment_y_dendrite',
    'bb_second_moment_z_apical',
    'bb_second_moment_z_basal',
    'bb_second_moment_z_dendrite',
    'bb_total_length_apical',
    'bb_total_length_basal','bb_total_length_dendrite',
    'bb_total_surface_area_apical',
    'bb_total_surface_area_basal',
    'bb_total_surface_area_dendrite',
    'bb_total_volume_apical',
    'bb_total_volume_basal',
    'bb_total_volume_dendrite']


#################################
#axon features: non_gmi
AXON_features=[]
for i in range (len(non_gmi_features)):
  AXON_features.append('axon_'+non_gmi_features[i])


# basl features:non_gmi
BASAL_features =[]
for i in range (len(non_gmi_features)):
  BASAL_features.append('basal_'+non_gmi_features[i])


#apical features: non_gmi  + kg
APICAL_features =[]
for i in range (len(non_gmi_features)):
  APICAL_features.append('apical_'+non_gmi_features[i])
APICAL_features.extend(apical_kg_features)


#all features reported: for feature specific analysis
ALL_FEATURES=[]
ALL_FEATURES.extend(AXON_features)
ALL_FEATURES.extend(BASAL_features)
ALL_FEATURES.extend(APICAL_features)
ALL_FEATURES.extend(BBP_features)


EXCITATORY_cols=[]
EXCITATORY_cols.extend(meta_names)
EXCITATORY_cols.extend(APICAL_features)
EXCITATORY_cols.extend(BASAL_features)


INHIBITORY_cols=[]
INHIBITORY_cols.extend(meta_names)
INHIBITORY_cols.extend(AXON_features)
INHIBITORY_cols.extend(BASAL_features)



