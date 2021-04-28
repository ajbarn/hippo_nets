%% this script will run the path length analysis for the visual network to DMN and visual network to hippocampus

% this script uses the function distance_wei_floyd.m from the BCT
% (https://sites.google.com/site/bctnet/)
%% remove networks and then calculate pathlength between DMN and visual network

[X,Y,INDSORT] = grid_communities(M1_rearr); %M1_rearr is the labels produced from Louvain Community detection
Z_fmriprep_9p_new_rearr = Z_fmriprep_9p_cortical(INDSORT,INDSORT);

%loop through each network that isn't vis or DMN
delete_nets = [2,3,4,5,6,7,8,9,11];


SPL_vis_DMN = zeros(size(Z_fmriprep_9p_cortical,3),size(delete_nets,2));
for i=1:size(delete_nets,2)
    clear SPL_deleted
    for sub=1:size(Z_fmriprep_9p_new_rearr,3)
        index = M1_rearr(INDSORT)~=delete_nets(i);
        community_matrix_order_paper_deleted = community_matrix_order_paper(index,1);
        SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_new_rearr(index,index,sub),0),'inv');
        SPL_vis_DMN(sub,i) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==10,sub),...
            [1,numel(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==10,sub))]));
    end
end

%% remove networks and then calculate path length between hippocampus and visual network


Z_fmriprep_9p_hipp =Z_all([1:4,26:203,205:384],[1:4,26:203,205:384],:); % select only the HCP_MMP regions and remove the ?? regions
Z_fmriprep_9p_hipp=Z_fmriprep_9p_hipp([1:182,362,183:361],[1:182,362,183:361],:); % move p24 to proper place

M1_rearr_hippocampus(5:362) = M1_rearr;
M1_rearr_hippocampus(1:4) = 0;
M1_rearr_hippocampus = M1_rearr_hippocampus';

[X,Y,INDSORT] = grid_communities(M1_rearr_hippocampus);

community_matrix_order_paper_hipp(5:362,1) = community_matrix_order_paper;
community_matrix_order_paper_hipp(1:4,1) = 0;

%loop through each network that isn't vis or hipp
delete_nets = 2:11;
Z_fmriprep_9p_hipp_new_rearr = Z_fmriprep_9p_hipp(INDSORT,INDSORT,:);

SPL_vis_hipp = zeros(size(Z_fmriprep_9p_hipp,3),size(delete_nets,2));
for i=1:size(delete_nets,2)
    clear SPL_deleted
    for sub=1:size(Z_fmriprep_9p_hipp,3)
        index = M1_rearr_hippocampus(INDSORT)~=delete_nets(i);
        community_matrix_order_paper_deleted = community_matrix_order_paper_hipp(index,1);
        SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_hipp_new_rearr(index,index,sub),0),'inv');
        SPL_vis_hipp(sub,i) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==0,sub),...
            [1,numel(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==0,sub))]));
    end
end

%% specificity of the MTN compared to removal of 30 random nodes

[X,Y,INDSORT] = grid_communities(M1_rearr_hippocampus);

index_list = (1:362)';
index_to_remove = index_list(M1_rearr_hippocampus(INDSORT)>1)';
clear SPL_deleted


for i = 1:1000
    for sub=1:size(Z_fmriprep_9p_hipp,3)
        random_index = index_to_remove(randperm(length(index_to_remove),30))';
        index= ones(362,1);
        index(random_index)=0;
        community_matrix_order_paper_deleted = community_matrix_order_paper_hipp(index_list(index==1),1);
        SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_hipp_new_rearr(index_list(index==1),index_list(index==1),sub),0),'inv');
        SPL_vis_hipp_perm(sub,i) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==0,sub),...
            [1,numel(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==0,sub))]));
    end
end

%compare each subjects path length after DMN removal to the distribution
%following 76 random nodes

for sub = 1:size(Z_fmriprep_9p_hipp,3)
    MTN_PL_effect_hipp(sub,1) = (SPL_vis_hipp(sub,9)-mean(SPL_vis_hipp_perm(sub,:)))/std(SPL_vis_hipp_perm(sub,:));
end

[h,p,ci,stats] = ttest(MTN_PL_effect_hipp)

%% specificity: does removal of the MTN result in disproportionately high path length to all networks or is its mediation of fc specific to visual network?
[X,Y,INDSORT] = grid_communities(M1_rearr);
 

%loop through each network that isn't vis or DMN
delete_nets = [1,2,3,4,5,6,7,8,9,11];


SPL_DMN = zeros(size(Z_fmriprep_9p_new_rearr,3),size(delete_nets,2),size(delete_nets,2));
for target = 1:size(delete_nets,2)
    for i=1:size(delete_nets,2)
        clear SPL_deleted
        for sub=1:size(Z_fmriprep_9p_new_rearr,3)
            index = M1_rearr(INDSORT)~=delete_nets(i);
            community_matrix_order_paper_deleted = community_matrix_order_paper(index,1);
            SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_new_rearr(index,index,sub),0),'inv');
            SPL_DMN(sub,i,target) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==delete_nets(target),community_matrix_order_paper_deleted==10,sub),...
                [1,numel(SPL_deleted(community_matrix_order_paper_deleted==delete_nets(target),community_matrix_order_paper_deleted==10,sub))]));
        end
    end
end

%% specificity: does removal of the MTN result in disproportionately high path length compared to removal of 30 random nodes (number of nodes in MTN)

[X,Y,INDSORT] = grid_communities(M1_rearr);

index_list = (1:358)';
index_to_remove = index_list(M1_rearr(INDSORT)>1)';
clear SPL_deleted


for i = 1:1000
    for sub=1:size(Z_fmriprep_9p_new,3)
        random_index = index_to_remove(randperm(length(index_to_remove),30))';
        index= ones(358,1);
        index(random_index)=0;
        community_matrix_order_paper_deleted = community_matrix_order_paper(index_list(index==1),1);
        SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_new_rearr(index_list(index==1),index_list(index==1),sub),0),'inv');
        SPL_vis_DMN_perm(sub,i) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==10,sub),...
            [1,numel(SPL_deleted(community_matrix_order_paper_deleted==1,community_matrix_order_paper_deleted==10,sub))]));
    end
end

%compare each subjects path length after DMN removal to the distribution
%following 30 random nodes

for sub = 1:size(Z_fmriprep_9p_new,3)
    MTN_PL_effect(sub,1) = (SPL_vis_DMN(sub,9)-mean(SPL_vis_DMN_perm(sub,:)))/std(SPL_vis_DMN_perm(sub,:));
end

[h,p,ci,stats] = ttest(MTN_PL_effect);

%% specificity of hippocampus effects

M1_rearr_hippocampus(5:362) = M1_rearr;
M1_rearr_hippocampus(1:4) = 0;
M1_rearr_hippocampus = M1_rearr_hippocampus';

[X,Y,INDSORT] = grid_communities(M1_rearr_hippocampus);

community_matrix_order_paper_hipp(5:362,1) = community_matrix_order_paper;
community_matrix_order_paper_hipp(1:4,1) = 0;

%loop through each network that isn't vis or hipp
delete_nets = 1:11;
Z_fmriprep_9p_hipp_new_rearr = Z_fmriprep_9p_hipp(INDSORT,INDSORT,:);

SPL_hipp = zeros(size(Z_fmriprep_9p_hipp,3),size(delete_nets,2),size(delete_nets,2));
for target = 1:size(delete_nets,2)
    for i=1:size(delete_nets,2)
        clear SPL_deleted
        for sub=1:size(Z_fmriprep_9p_hipp,3)
            index = M1_rearr_hippocampus(INDSORT)~=delete_nets(i);
            community_matrix_order_paper_deleted = community_matrix_order_paper_hipp(index,1);
            SPL_deleted(:,:,sub) = distance_wei_floyd(threshold_absolute(Z_fmriprep_9p_hipp_new_rearr(index,index,sub),0),'inv');
            SPL_hipp(sub,i,target) = mean(reshape(SPL_deleted(community_matrix_order_paper_deleted==delete_nets(target),community_matrix_order_paper_deleted==0,sub),...
                [1,numel(SPL_deleted(community_matrix_order_paper_deleted==delete_nets(target),community_matrix_order_paper_deleted==0,sub))]));
        end
    end
end

delete_nets = [1:10];

for target=1
    fprintf('target is')
    target
    for i = [2,3,4,5,6,7,8,9,10]
        i
        [h,p,ci,stats] = ttest(SPL_hipp(:,11,target),SPL_hipp(:,delete_nets(i),target));
        stats
    end
end





%% boxplot 
figure
for i = 1:9
subplot(5,2,i)
boxplot(SPL_DMN(:,:,i))
end

figure
for i = 1:10
subplot(5,2,i)
boxplot(SPL_hipp(:,:,i))
end



