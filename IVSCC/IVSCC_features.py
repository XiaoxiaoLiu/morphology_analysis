__author__ = 'xiaoxiaoliu'
import numpy as np

basal_features = [
       #u'basal_dendrite_average_diameter',
       u'basal_dendrite_bifurcation_angle_local',
       u'basal_dendrite_bifurcation_angle_remote',
       u'basal_dendrite_bifurcation_centroid_over_distance_x',
       u'basal_dendrite_bifurcation_centroid_over_distance_y',
       u'basal_dendrite_bifurcation_centroid_over_distance_z',
       u'basal_dendrite_bifurcation_stdev_over_centroid',
       u'basal_dendrite_bifurcation_variance_over_distance_x',
       u'basal_dendrite_bifurcation_variance_over_distance_y',
       u'basal_dendrite_bifurcation_variance_over_distance_z',
       u'basal_dendrite_compartment_centroid_over_distance_x',
       u'basal_dendrite_compartment_centroid_over_distance_y',
       u'basal_dendrite_compartment_centroid_over_distance_z',
       u'basal_dendrite_compartment_stdev_over_centroid',
       u'basal_dendrite_compartment_variance_over_distance_x',
       u'basal_dendrite_compartment_variance_over_distance_y',
       u'basal_dendrite_compartment_variance_over_distance_z',
       u'basal_dendrite_contraction', u'basal_dendrite_depth',
       u'basal_dendrite_first_bifurcation_moment_x',
       u'basal_dendrite_first_bifurcation_moment_y',
       u'basal_dendrite_first_bifurcation_moment_z',
       u'basal_dendrite_first_compartment_moment_x',
       u'basal_dendrite_first_compartment_moment_y',
       u'basal_dendrite_first_compartment_moment_z', u'basal_dendrite_height',
       u'basal_dendrite_max_branch_order',
       u'basal_dendrite_max_euclidean_distance',
       u'basal_dendrite_max_path_distance',
       u'basal_dendrite_mean_fragmentation',
       #u'basal_dendrite_mean_parent_daughter_ratio',
       u'basal_dendrite_neurites_over_branches',
       u'basal_dendrite_num_bifurcations', u'basal_dendrite_num_branches',
       u'basal_dendrite_num_neurites',
       #u'basal_dendrite_num_nodes',
       u'basal_dendrite_num_stems', u'basal_dendrite_num_tips',
       #u'basal_dendrite_parent_daughter_ratio',
       u'basal_dendrite_second_bifurcation_moment_x',
       u'basal_dendrite_second_bifurcation_moment_y',
       u'basal_dendrite_second_bifurcation_moment_z',
       u'basal_dendrite_second_compartment_moment_x',
       u'basal_dendrite_second_compartment_moment_y',
       u'basal_dendrite_second_compartment_moment_z',
       #u'basal_dendrite_soma_surface',
       u'basal_dendrite_total_length',
       #u'basal_dendrite_total_surface',
       #u'basal_dendrite_total_volume',
       u'basal_dendrite_width']


apical_features = [
       #u'apical_dendrite_average_diameter',
       u'apical_dendrite_bifurcation_angle_local',
       u'apical_dendrite_bifurcation_angle_remote',
       u'apical_dendrite_bifurcation_centroid_over_distance_x',
       u'apical_dendrite_bifurcation_centroid_over_distance_y',
       u'apical_dendrite_bifurcation_centroid_over_distance_z',
       u'apical_dendrite_bifurcation_stdev_over_centroid',
       u'apical_dendrite_bifurcation_variance_over_distance_x',
       u'apical_dendrite_bifurcation_variance_over_distance_y',
       u'apical_dendrite_bifurcation_variance_over_distance_z',
       u'apical_dendrite_compartment_centroid_over_distance_x',
       u'apical_dendrite_compartment_centroid_over_distance_y',
       u'apical_dendrite_compartment_centroid_over_distance_z',
       u'apical_dendrite_compartment_stdev_over_centroid',
       u'apical_dendrite_compartment_variance_over_distance_x',
       u'apical_dendrite_compartment_variance_over_distance_y',
       u'apical_dendrite_compartment_variance_over_distance_z',
       u'apical_dendrite_contraction', u'apical_dendrite_depth',
       u'apical_dendrite_first_bifurcation_moment_x',
       u'apical_dendrite_first_bifurcation_moment_y',
       u'apical_dendrite_first_bifurcation_moment_z',
       u'apical_dendrite_first_compartment_moment_x',
       u'apical_dendrite_first_compartment_moment_y',
       u'apical_dendrite_first_compartment_moment_z',
       u'apical_dendrite_height', u'apical_dendrite_max_branch_order',
       u'apical_dendrite_max_euclidean_distance',
       u'apical_dendrite_max_path_distance',
       u'apical_dendrite_mean_fragmentation',
       #u'apical_dendrite_mean_parent_daughter_ratio',
       u'apical_dendrite_neurites_over_branches',
       u'apical_dendrite_num_bifurcations', u'apical_dendrite_num_branches',
       u'apical_dendrite_num_neurites',
       #u'apical_dendrite_num_nodes',
       #u'apical_dendrite_num_stems',
       u'apical_dendrite_num_tips',
       #u'apical_dendrite_parent_daughter_ratio',
       u'apical_dendrite_second_bifurcation_moment_x',
       u'apical_dendrite_second_bifurcation_moment_y',
       u'apical_dendrite_second_bifurcation_moment_z',
       u'apical_dendrite_second_compartment_moment_x',
       u'apical_dendrite_second_compartment_moment_y',
       u'apical_dendrite_second_compartment_moment_z',
       # u'apical_dendrite_soma_surface',
       u'apical_dendrite_total_length',
       #u'apical_dendrite_total_surface', u'apical_dendrite_total_volume',
       u'apical_dendrite_width']



axon_features = [
       #u'axon_average_diameter',
       u'axon_bifurcation_angle_local',
       u'axon_bifurcation_angle_remote',
       u'axon_bifurcation_centroid_over_distance_x',
       u'axon_bifurcation_centroid_over_distance_y',
       u'axon_bifurcation_centroid_over_distance_z',
       u'axon_bifurcation_stdev_over_centroid',
       u'axon_bifurcation_variance_over_distance_x',
       u'axon_bifurcation_variance_over_distance_y',
       u'axon_bifurcation_variance_over_distance_z',
       u'axon_compartment_centroid_over_distance_x',
       u'axon_compartment_centroid_over_distance_y',
       u'axon_compartment_centroid_over_distance_z',
       u'axon_compartment_stdev_over_centroid',
       u'axon_compartment_variance_over_distance_x',
       u'axon_compartment_variance_over_distance_y',
       u'axon_compartment_variance_over_distance_z', u'axon_contraction',
       u'axon_depth', u'axon_first_bifurcation_moment_x',
       u'axon_first_bifurcation_moment_y', u'axon_first_bifurcation_moment_z',
       u'axon_first_compartment_moment_x', u'axon_first_compartment_moment_y',
       u'axon_first_compartment_moment_z', u'axon_height',
       u'axon_max_branch_order', u'axon_max_euclidean_distance',
       u'axon_max_path_distance', u'axon_mean_fragmentation',
       #u'axon_mean_parent_daughter_ratio',
       u'axon_neurites_over_branches',
       u'axon_num_bifurcations', u'axon_num_branches', u'axon_num_neurites',
       #u'axon_num_nodes',
       #u'axon_num_stems',
       u'axon_num_tips',
       #u'axon_parent_daughter_ratio',
       u'axon_second_bifurcation_moment_x',
       u'axon_second_bifurcation_moment_y',
       u'axon_second_bifurcation_moment_z',
       u'axon_second_compartment_moment_x',
       u'axon_second_compartment_moment_y',
       u'axon_second_compartment_moment_z',
       #u'axon_soma_surface',
       u'axon_total_length',
       #u'axon_total_surface', u'axon_total_volume',
       u'axon_width']


def get_feature_names(type='basal'):
    if type =="all_dendrite":
        all_dendrite_features=[]
        all_dendrite_features.extend(basal_features)
        all_dendrite_features.extend(apical_features)
        all_dendrite_features.extend(axon_features)
        return  all_dendrite_features
    if type == "basal":
        return basal_features
    if type == "apical":
        print "apical features:", len(apical_features)
        return apical_features
    if type == "axon":
        return axon_features
    if type =="spiny_dendrite":
        # apical + basal
       spiny_dendrite_features=[]
       spiny_dendrite_features.extend(basal_features)
       spiny_dendrite_features.extend(apical_features)
       return  spiny_dendrite_features
    if type =="aspiny":
       # np.append(gl_feature_names, gmi_feature_names)
       return basal_features
    if type =="spiny_dendrite_no_z":
       spiny_dendrite_features=[]
       spiny_dendrite_features.extend(basal_features)
       spiny_dendrite_features.extend(apical_features)
       for feature in spiny_dendrite_features:
            if "_z" in feature:
                  spiny_dendrite_features.remove(feature)
                  print "remove ", feature
       return  spiny_dendrite_features
    if type == "apical_no_z":
        apical_no_z_features = apical_features[:] ### be careful! pass by copying
        for feature in apical_no_z_features:
            if "_z" in feature:
                  apical_no_z_features.remove(feature)
                  print "remove ", feature
        print "apical( no z) features:", len(apical_no_z_features)
        return apical_no_z_features