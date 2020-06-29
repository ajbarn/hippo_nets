 %% extract hippocampal connectivity to each network for each subject
 % run after louv_comm_det_iteration_fmriprep_9p.m

hippo_net_connectivity = zeros(4,max(community_matrix_order),size(Z_all,3)); % initiate matrix

%index the regions of interest for each network

trans_mat = [26:203,205:384]; % convert Z_fmriprep_9p_cortical to Z_all
subcort_mat = [1:4]; % call on each hippocampal ROI in Z_all
for sub = 1:size(Z_all,3)
    for i = 1:max(community_matrix_order)
        for j = 1:size(subcort_mat,2)
            index = find(M1==i);
            temp_net_corr = Z_all(subcort_mat(j),trans_mat(index),sub);
            hippo_net_connectivity(j,i,sub) = mean(temp_net_corr);
        end
    end
end

figure
imagesc(mean(hippo_net_connectivity,3)./(sqrt(var(hippo_net_connectivity,[],3))/sqrt(40)))

mean(hippo_net_connectivity,3) % mean hippocampal activity for left/right, ant/post hippocampus to every network
sqrt(var(hippo_net_connectivity,[],3))./sqrt(40) %sd hippocampal activity for left/right, ant/post hippocampus to every network