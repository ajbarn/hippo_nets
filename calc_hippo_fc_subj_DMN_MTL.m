%% calculate hippocampal and MTL connectivity to DMN networks
% run after louv_comm_det_iteration_fmriprep_9p_DMN.m

Z_trans_mat = [26:203,384,205:383];

% index the labels for the MTN and DMN
% these labels are arbitrary as far as louvain community detection is
% concerned and will change run to run; enter accordingly

index_MTL_net = M1_rearr==11; % find the MTL nodes which were given the network label 11

index_DMN_net = find(M1_rearr==10); %extracting out JUST the DMN

index_MTL_DMN_net = find(M1_rearr==11|M1_rearr==10); % index both the MTL and DMN nodes

index_3net_DMN = zeros(size(Z_fmriprep_9p_cortical,1),1);
index_3net_DMN(Z_trans_mat(index_DMN_net(:)),1) = M1_PMAT;

labels = zeros(size(Z_all,1),1); %first make a label of the communities on the whole network
%labels(1:4) = 5; % post hippocampal ROIs labelled as 5
labels(1:2) = 5; % post hippocampal ROIs labelled as 5
labels(3:4) = 7; % ant hippocampal ROIs labelled as 7
labels(Z_trans_mat(index_MTL_net)) = 6; %MTL ROIs as 6
labels(index_3net_DMN>0,1)=index_3net_DMN(index_3net_DMN>0,1);




Z_hip_MTL_DMN = Z_all([1:4,Z_trans_mat(index_MTL_DMN_net(:))],[1:4,Z_trans_mat(index_MTL_DMN_net(:))],:); % get an array for hippocampus ROIs, MTL, and DMN ROIs for all subjects
names_hip_MTL_DMN_rearr = names(1,[1:4,Z_trans_mat(index_MTL_DMN_net(:))]); % get the names of these ROIs to cross reference
labels_rearr = labels([1:4,Z_trans_mat(index_MTL_DMN_net(:))],1); % get the community labels of the ROIs in Z_hip_MTL_DMN
Z_hip_MTL_DMN(isnan(Z_hip_MTL_DMN)) = 1; % turn diagonan NaNs into 1s
Z_hip_MTL_DMN_avg = mean(Z_hip_MTL_DMN,3); % average across subjects in dimension 3


%% change some of the labeling so that it can be plotted with hippo and MTL first
labels_rearr(labels_rearr==1) = 9;
labels_rearr(labels_rearr==2) = 10;
labels_rearr(labels_rearr==3) = 8;
labels_rearr(labels_rearr==6) = 3;
labels_rearr(labels_rearr==5) = 2; %post hipp
labels_rearr(labels_rearr==7) = 1; %ant hipp


%% plot the communities reordered
clear INDSORT
[X,Y,INDSORT] = grid_communities(labels_rearr);
figure
imagesc(Z_hip_MTL_DMN_avg(INDSORT,INDSORT))
hold on;
plot (X,Y,'b','linewidth',2);


names_hip_MTL_DMN_ordered = names_hip_MTL_DMN_rearr(1,INDSORT);

for x = 1:size(labels_rearr,1)
community_matrix_order_hip_MTL_DMN = labels_rearr(INDSORT,1);
end



%% calc hippocampal connectivity to each subnetwork
Z_hip_MTL_DMN_avg_rearr=Z_hip_MTL_DMN_avg(INDSORT,INDSORT);
Z_hip_MTL_DMN_rearr = Z_hip_MTL_DMN(INDSORT,INDSORT,:);



% individual subject matrix
community_list = unique(community_matrix_order_hip_MTL_DMN);

for sub = 1:size(Z_hip_MTL_DMN_rearr,3)
    for i = 1:size(community_list,1) % loop through the MTN and DMN subnetworks
        for j = 1:4 % loop through hippocampal ROIs (left/right; ant/post)
            index = find(community_matrix_order_hip_MTL_DMN==community_list(i));
            temp_net_corr = Z_hip_MTL_DMN_rearr(j,index,sub);
            hippo_sub_conn(j,i,sub) = mean(temp_net_corr);
        end
    end
end

hippo_sub_conn_reshape = reshape(hippo_sub_conn,size(hippo_sub_conn,1)*size(hippo_sub_conn,2),size(hippo_sub_conn,3),1);
hippo_sub_conn_reshape = hippo_sub_conn_reshape';  %% this data is subsequently analyzed in hippo_nets_conn.ipynb