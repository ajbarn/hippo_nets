% extract out ROIs from the DMN and perform community detection across
% range of gamma resolutions

%uses community_louvain from the brain connectivity toolbox
%(https://sites.google.com/site/bctnet/)
%uses zrand from http://commdetect.weebly.com/


%% extract out just the DMN nodes
% the exact network number assigned to the DMN is random and will change
% every time
trans_mat = [26:203,205:384];
index_hippo_nets = find(M1_rearr==10); %extracting out JUST the DMN (the exact number assigned to the DMN is arbitrary and will change from solution to solution when running louvain community detection)

%select out the fc for just the DMN
Z_fmriprep_9p_PMAT_avg = Z_fmriprep_9p_avg(index_hippo_nets(:),index_hippo_nets(:)); Z_all([trans_mat(index_hippo_nets(:))],[trans_mat(index_hippo_nets(:))],:);
%change diagonal to 1's
Z_fmriprep_9p_PMAT_avg(isnan(Z_fmriprep_9p_PMAT_avg)) = 1;


%% perform community detection across range of gamma resolutions and track ZR
%and Q
gamma_range_PMAT = .75:.01:1.1;
conn_matrix = Z_fmriprep_9p_PMAT_avg;
iterations = 1000;
ZR_PMAT = zeros((iterations^2-iterations)/2,size(gamma_range_PMAT,2)); 
M_PMAT = zeros(size(conn_matrix,1),iterations,size(gamma_range_PMAT,2));
Q1_PMAT = zeros(iterations,size(gamma_range_PMAT,2));

gamma_number =1;
for gamma = gamma_range_PMAT 

    iter = 1;
    while iter < iterations+1
        % Iterative community finetuning.
        % W is the input connection matrix.
        n  = size(conn_matrix,1);             % number of nodes
        M_PMAT(:,iter,gamma_number)  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1_PMAT(iter,gamma_number) = 0;            % initialize modularity values
        while Q1_PMAT(iter,gamma_number)-Q0>1e-5           % while modularity increases
            Q0 = Q1_PMAT(iter,gamma_number);                % perform community detection
            [M_PMAT(:,iter,gamma_number), Q1_PMAT(iter,gamma_number)] = community_louvain(conn_matrix, gamma, M_PMAT(:,iter,gamma_number),'negative_asym'); %in the community_louvain function input the array of interest, resolution parameter, and initial community affiliations that will get tuned (Here I have M(:,iter) because there will be 1000 sets of community affiliations)
        end
        iter=iter+1;
    end
    
    idx=1;
    k=1;
    for i = 1:size(M_PMAT,2)
        for j = k:size(M_PMAT,2)
            if i~=j
                ZR_PMAT(idx,gamma_number) = zrand(M_PMAT(:,i,gamma_number),M_PMAT(:,j,gamma_number));
                idx = idx+1;
            end
        end
        k=k+1;
    end
    gamma_number = gamma_number +1;
end


%% calculate ZR max
ZR_gamma_PMAT = zeros(iterations,iterations,size(gamma_range_PMAT,2));

for gamma_number = 1:size(gamma_range_PMAT,2)
    k=1;
    idx = 1;
    for j = 1:iterations
        for i=k:iterations
            
            if i~=j
                
                ZR_gamma_PMAT(i,j,gamma_number) = ZR_PMAT(idx,gamma_number);
                idx=idx+1;
            end
        end
        k=k+1;
    end
end


for gamma_number = 1:size(gamma_range_PMAT,2)
    k=1;
    idx = 1;
    for j = 1:iterations
        for i=k:iterations
            
            if i~=j
                
                ZR_gamma_PMAT(j,i,gamma_number) = ZR_gamma_PMAT(i,j,gamma_number);
                idx=idx+1;
            end
        end
        k=k+1;
    end
end

ZR_PMAT_gamma_mean = nanmean(ZR_gamma_PMAT,2);
ZR_PMAT_gamma_mean = reshape(ZR_PMAT_gamma_mean,iterations,size(gamma_range_PMAT,2),1);
ZR_PMAT_max = max(ZR_PMAT_gamma_mean,[],1);

for i = 1:size(gamma_range_PMAT,2)
    max_ZR_PMAT_index(1,i) = find(ZR_PMAT_gamma_mean(:,i)==max(ZR_PMAT_gamma_mean(:,i),[],1),1,'first');
end

for i = 1:size(gamma_range_PMAT,2)
ZR_PMAT_max_mul_Q(1,i) = ZR_PMAT_max(1,i)*Q1_PMAT(max_ZR_PMAT_index(i),i);
end

figure
plot(gamma_range_PMAT, ZR_PMAT_max_mul_Q)

%% select the solution with highest ZR*Q for best_gamma_number
best_gamma_number = 11; 
M1_PMAT=M_PMAT(:,max_ZR_PMAT_index(1,best_gamma_number),best_gamma_number);

clear INDSORT
[X,Y,INDSORT] = grid_communities(M1_PMAT);
figure
imagesc(Z_fmriprep_9p_PMAT_avg(INDSORT,INDSORT))
hold on;
plot (X,Y,'b','linewidth',2);


% rearrange for publication
M1_PMAT_rear = zeros(size(M1_PMAT,1),1);
M1_PMAT_rear(M1_PMAT==3)=2;
M1_PMAT_rear(M1_PMAT==2)=3;
M1_PMAT_rear(M1_PMAT==1)=1;


clear INDSORT
[X,Y,INDSORT] = grid_communities(M1_PMAT_rear);
figure
imagesc(Z_fmriprep_9p_PMAT_avg(INDSORT,INDSORT)) % first ATN then PMN then MPN
hold on;
plot (X,Y,'b','linewidth',2);


Z_fmriprep_9p_PMAT_avg_rearr=Z_fmriprep_9p_PMAT_avg(INDSORT,INDSORT);

for x = 1:size(M1_PMAT,1)
community_matrix_order_DMN = M1_PMAT_rear(INDSORT,1);
end