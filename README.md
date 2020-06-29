# hippo_nets
A repository for the data and and scripts used in the study: https://www.biorxiv.org/content/10.1101/2020.06.09.142166v1

# data
## matlab workspace
Within this workspace.mat file there are several arrays that will be used when performing community detection (wholebrain= louv_comm_det_iteration_fmriprep_9p.m; DMN  = louv_comm_det_iteration_fmriprep_9p_DMN.m), extraction of hippocampal functional connectivity to all the large scale network (calculate_hippo_fc_subj_fmriprep_9p.m) and to the cortico-hippocampal networks (calc_hippo_fc_subj_DMN_MTL.m), and calculation of participation coefficient (calc_participation_strength_hipp_cortico_net.m).

Z_all is an ROIxROIxSubject array holding pairwise Fisher Z transformed functional connectivity. The names of the ROIs are held in the array names

Z_fmriprep_9p_cortical is an ROIxROIxSubject array holding pairwise Fisher Z transformed functional connectivity. The names of the ROIs are held in the array names_358.

Z_fmriprep_9p_avg is an ROIxROI array that is the Z_fmriprep_9p_cortical averaged along the subject dimension.The names of the ROIs are held in the array names_358.

The RdBl_map is the colormap for the imagesc figures. Colormap RGB coordinates used from Diverging Color Maps for Scientific Visualization
(Expanded). Kenneth Moreland.(https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf).

## Force-directed graphs
fmriprep_9p_pos_links_358.csv & fmriprep_9p_pos_nodes_358.csv
- These csv files are a pair to be used together to plot the force-directed graphs using net_vis_fmriprep_9p.Rmd

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
