% Jan 20 2020
%calculate participation and strength for the hippo cortical networks
% run after louv_comm_det_iteration_fmriprep_9p_DMN.m

index_MTL_DMN_net = find(M1_rearr==11|M1_rearr==10); % index both the MTL and DMN nodes


thresholds = .05:.01:.2; %threshold fc matrix from 5% sparcity to 20%

%calc pc for all nodes all networks across threshold range
for thresh = 1:size(thresholds,2)
    PC_whole(:,thresh) = participation_coef(threshold_proportional(Z_fmriprep_9p_avg_pos,thresholds(thresh)),M1_rearr);
end

%calculate whole brain strength
strength_wholebrain = strengths_und(Z_fmriprep_9p_avg_pos)';
strength_whole_brain_PMAT = strength_wholebrain(index_MTL_DMN_net,1);

%% calculate the participation coefficient for the nodes in the DMN when they are broken down by subnetwork


%get a threshold range that is proportional to the whole network
thresholds = 80:1:95;

for thresh = 1:size(thresholds,2)
    PC_PMAT(:,thresh) = participation_coef(threshold_absolute(Z_fmriprep_9p_PMAT_avg_rearr,prctile(Z_fmriprep_9p_avg_pos(:),thresholds(thresh))),community_matrix_order_DMN);
end


for i = 1:size(PC_PMAT,1)
    PC_PMAT_avg(i,1) = mean(PC_PMAT(i,PC_PMAT(i,:)>0));
end

strength_PMAT= strengths_und(Z_fmriprep_9p_PMAT_avg_rearr);


