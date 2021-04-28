# hippo_nets
A repository for the data and and scripts used in the study: https://www.biorxiv.org/content/10.1101/2020.06.09.142166v1

# data
## matlab workspace
Within this workspace.mat file there are several arrays that were used when performing community detection (wholebrain= louv_comm_det_iteration_fmriprep_9p.m; DMN  = louv_comm_det_iteration_fmriprep_9p_DMN.m), extraction of hippocampal functional connectivity to all the large scale network (calculate_hippo_fc_subj_fmriprep_9p.m) and to the cortico-hippocampal networks (calc_hippo_fc_subj_DMN_MTL.m), and calculation of participation coefficient (calc_participation_strength_hipp_cortico_net.m).

Z_all is an ROIxROIxSubject array holding pairwise Fisher Z transformed functional connectivity. The names of the ROIs are held in the array names

Z_fmriprep_9p_cortical is an ROIxROIxSubject array holding pairwise Fisher Z transformed functional connectivity. The names of the ROIs are held in the array names_358.

Z_fmriprep_9p_avg is an ROIxROI array that is the Z_fmriprep_9p_cortical averaged along the subject dimension.The names of the ROIs are held in the array names_358.

The RdBl_map is the colormap for the imagesc figures. Colormap RGB coordinates used from Diverging Color Maps for Scientific Visualization
(Expanded). Kenneth Moreland.(https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf).

## Force-directed graphs
fmriprep_9p_pos_links_358.csv & fmriprep_9p_pos_nodes_358.csv
- the fmriprep_9p_pos_links_358.csv datafile contains the average functional connectivity matrix from the ROIs from the HCP-MMP1.0 atlas for the sample
- the fmriprep_9p_pos_nodes_358.csv contains the names of each ROI
- These csv files are a pair to be used together to plot the force-directed graphs using net_vis_fmriprep_9p.Rmd

## Path Length analysis
script for calculating pathlength after removal of each network is calc_pathlength_delete_nets.m and requires community labels from the Louvain community detection
plotting of this data done in path_length_analysis.ipynb using data from pathlength_data_all.xlsx and pathlength_data_all_hippo.xlsx which is obtained from running calc_pathlength_delete_nets.m 

## Participation coefficient data
fmriprep_9p_DMN_nodes_PC.csv
- This csv file contains the participation coefficient values for the DMN nodes to the whole brain

PC_DMN.xlsx
- This csv file contains the participation coefficient values for the DMN nodes for just within the DMN
These are plotted and stats are run using PC_figures.ipynb

## hippocampal functional connectivity

hippo_net_con.xlsx
- this spreadsheet contains the functional connectivity of the hippocampal regions of interest to the cortico-hippocampal networks

This data is plotted with hippo_nets_conn.ipynb

## Representional Profile similarity data

within_between_RDM_sim.xlsx
- this spreadsheet contains the within and between network representational profile similarity values for each subject
net2net_similarity.xlsx

net2net_similarity.xlsx
- this spreadsheet contains the averaged representational profile similarity values for every cortico-hippocampal network to every other cortico-hippocampal network for each subjects

These are plotted and stats are run using RDM_sim_plotting.ipynb

## neurosynth

nifti files used to decode terms: glasser_PMAT_net1.nii.gz, glasser_PMAT_net2.nii.gz, glasser_PMAT_net3.nii.gz, glasser_context_net.nii.gz

Terms are estimated using neurosynth_decoding.ipynb

curated terms for wordclouds: context_terms_358.xlsx, DMN1_terms_358.xlsx, DMN2_terms_358.xlsx, DMN3_terms_358.xlsx

Terms are plotted using wordcloud_fig.ipynb

# Figures

## Figure 1

Figure 1B
Functional connectivity matrix for the array Z_fmriprep_9p_avg_rearr in rearranged_matrices.mat
This data is from the array Z_fmriprep_9p_avg reorganized according to labels produced from the script louv_comm_det_iteration_fmriprep_9p.m
Data for this matrix was produced from CONN:Connectivity Toolbox (Whitfield-Gabrielli & Nieto-Castanon, 2012) version 18b that was used to estimate region-to-region functional connectivity from participant data that was preprocessed using fMRIPrep 1.4.1 Esteban et al (2019).

Figure 1C
The data used to produce this figure are in fmriprep_9p_pos_links_358.csv & fmriprep_9p_pos_nodes_358.csv
- the fmriprep_9p_pos_links_358.csv datafile contains the average functional connectivity matrix from Z_fmriprep_9p_avg
- the fmriprep_9p_pos_nodes_358.csv contains the label information of each ROI
These csv files are a pair to be used together to plot the force-directed graphs using net_vis_fmriprep_9p.Rmd. The R package igraph version 1.2.4.1 along with ForceAtlas2 version 0.1 were used to produce this plot.

Figure 2B
plotting of this data done in path_length_analysis.ipynb using data from pathlength_data_all.xlsx and pathlength_data_all_hippo.xlsx which is obtained from running calc_pathlength_delete_nets.m 

Figure 3
Functional connectivity matrix for the array Z_fmriprep_9p_PMAT_avg_rearr in rearranged_matrices.mat
This data is the functional connectivity values of the regions identified as the default mode network (DMN) in the louv_comm_det_iteration_fmriprep_9p.m script.
It has then undergone another round of Louvain Community detection and has been reorganized according to the partitioning solution that produced the highest modularity weighted z-Rand score from the script louv_comm_det_iteration_fmriprep_9p_DMN.m

Figure 4A
Force-directed plot created using net_vis_fmriprep_9p.Rmd. The R package igraph version 1.2.4.1 along with ForceAtlas2 version 0.1 were used to produce this plot. 
The DMN network labels were replaced with subnetwork labels

Figure 4B
Participation coefficients for the DMN subnetworks for the diversity of connections outside of the DMN and, separately, inside the DMN to other subnetworks
Data was calculated using participation_coef.m from the brain connectivity toolbox (BCT) along with the network labels attained from Louvain community detection (above) and was plotted using the notebook PC_figures.ipynb

Figure 5
Functional connectivity of the anterior and posterior hippocampus to the DMN subnetworks extracted from each participant using the labaels attained from Louvain community detection (above) and plotted using hippo_nets_conn.ipynb

Figure 6 LEFT
These data are the average spearman correlations of the representational dissimilarity profiles of of pairs of regions of interest that are in the same cortico-hippocampal network or pairs that are in different cortico-hippocampal networks. 
The data are found in within_between_RDM_sim.xlsx and they are plotted using the RDM_sim_plotting.ipynb notebook

Figure 6 Right
These data are the average spearman correlations of the representational dissimilarity profiles of pairs of regions in each cortico-hippocampal network, to every other cortico-hippocampal. 
The data are found in net2net_similarity.xlsx and plotted using the RDM_sim_plotting.ipynb notebook

Figure 7
The size of the words in the word clouds in this plot are related to how strongly each term relates to the spatial map of its corresponding cortico-hippocampal network
The term decoding was done using the nifti files : glasser_PMAT_net1.nii.gz, glasser_PMAT_net2.nii.gz, glasser_PMAT_net3.nii.gz, glasser_context_net.nii.gz and the terms are estimated using neurosynth_decoding.ipynb
This produced curated terms for wordclouds: context_terms_358.xlsx, DMN1_terms_358.xlsx, DMN2_terms_358.xlsx, DMN3_terms_358.xlsx

Terms are plotted using wordcloud_fig.ipynb

Supplemental Figure 2
Functional connectivity matrix for the array Z_fmriprep_9p_avg_rearr_WashU in washu_workspace.mat
This data is from the array Z_fmriprep_9p_avg reorganized according to labels produced from the script louv_comm_det_iteration_fmriprep_9p.m
Data for this matrix was produced from CONN:Connectivity Toolbox (Whitfield-Gabrielli & Nieto-Castanon, 2012) version 18b that was used to estimate region-to-region functional connectivity from participant data that was preprocessed using fMRIPrep 1.4.1 Esteban et al (2019).


Supplemental Figure 3A
Displayed average tsnr that were below 30 across the UCD and washu120 samples plotted on an inflated brain surface. Data for this are found in washu_workspace.mat.
The arrays that correspond to the tsnr are labelled WashU_tsnr and UCD_tsnr

Supplemental Figure 3B
t-values for the differences between the UCD and washu120 samples in tsnr for the low tsnr regions plotted on an inflated brain surface. 
t-values for all regions can be found in the arrays WashU_v_UCD_L and WashU_v_UCD_R in washu_workspace.mat. 
Positive t-values are those in which the UCD dataset has higher tnsr.

Supplemental Figure 4

Supplemental Figure 5

Supplemental Figure 6
These data are the average spearman correlations of the representational dissimilarity profiles of pairs of regions in each cortico-hippocampal network, to every other network in the brain. 
The data are found in net2net_similarity_358_all_March2021_wholecortex.xlsx and plotted using the RDM_sim_plotting-wholecortex.ipynb









