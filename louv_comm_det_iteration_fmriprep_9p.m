%%run community detection on a range of
%%resolution values and track modularity (Q) and assignment consistency
%%(Z-rand)

%uses community_louvain from the brain connectivity toolbox
%(https://sites.google.com/site/bctnet/)
%uses zrand from http://commdetect.weebly.com/

Z_fmriprep_9p_avg_pos= threshold_absolute(Z_fmriprep_9p_avg,0); % remove negative weights
%% run community detection on the Z_fmriprep_9p_avg matrix, and track modularity (Q) and consistency, (AR/ZR)
gamma_range = 2:.005:2.5; % gamma is starting at 2 because prior solutions don't break apart somatomotor and auditory network
conn_matrix = Z_fmriprep_9p_avg_pos;
iterations = 1000;
SAR = zeros((iterations^2-iterations)/2,size(gamma_range,2)); 
ZR = zeros((iterations^2-iterations)/2,size(gamma_range,2)); 
M = zeros(size(conn_matrix,1),iterations,size(gamma_range,2));
Q1 = zeros(iterations,size(gamma_range,2));

h = waitbar(0);
pause(.01)
conn_matrix(conn_matrix<0)= 0; % remove negative weights
gamma_number =1;
for gamma = gamma_range 
    iter = 1;
    while iter < iterations+1
        % Iterative community finetuning.
        % W is the input connection matrix.
        n  = size(conn_matrix,1);             % number of nodes
        M(:,iter,gamma_number)  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1(iter,gamma_number) = 0;            % initialize modularity values
        while Q1(iter,gamma_number)-Q0>1e-5           % while modularity increases
            Q0 = Q1(iter,gamma_number);                % perform community detection
            [M(:,iter,gamma_number), Q1(iter,gamma_number)] = community_louvain(conn_matrix, gamma, M(:,iter,gamma_number),'negative_asym'); %in the community_louvain function input the array of interest, resolution parameter, and initial community affiliations that will get tuned (Here I have M(:,iter) because there will be 1000 sets of community affiliations)
        end
        iter=iter+1;
        
    end
    idx=1;
    k=1;
    for i = 1:size(M,2)
        for j = k:size(M,2)
            if i~=j
                [ZR(idx,gamma_number),SR(idx,gamma_number), SAR(idx,gamma_number)] = zrand(M(:,i,gamma_number),M(:,j,gamma_number));
                idx = idx+1;
            end
        end
        k=k+1;
    end
    gamma_number = gamma_number +1; 
    waitbar(gamma_number/size(gamma_range,2),h);
    
end

%% calculate highest mean ZRand solutions and 
%%and plot to select those the ZR*Q metric

ZR_gamma = zeros(iterations,iterations,size(gamma_range,2));

for gamma_number = 1:size(gamma_range,2)
    k=1;
    idx = 1;
    for j = 1:iterations
        for i=k:iterations
            
            if i~=j
                
                ZR_gamma(i,j,gamma_number) = ZR(idx,gamma_number);
                idx=idx+1;
            end
        end
        k=k+1;
    end
end


for gamma_number = 1:size(gamma_range,2)
    k=1;
    idx = 1;
    for j = 1:iterations
        for i=k:iterations
            
            if i~=j
                
                ZR_gamma(j,i,gamma_number) = ZR_gamma(i,j,gamma_number);
                idx=idx+1;
            end
        end
        k=k+1;
    end
end

ZR_gamma_mean = mean(ZR_gamma,2);
ZR_gamma_mean = reshape(ZR_gamma_mean,iterations,size(gamma_range,2),1);
ZR_max = max(ZR_gamma_mean,[],1);

% find the top solutions
for i = 1:size(gamma_range,2)
    max_ZR_index(1,i) = find(ZR_gamma_mean(:,i)==max(ZR_gamma_mean(:,i),[],1),1,'first');
end

%multiply the modularity for the top ZR solution with the ZR across all
%resolutions
for i = 1:size(gamma_range,2)
ZR_max_mul_Q(1,i) = ZR_max(1,i)*Q1(max_ZR_index(i),i);
end

%% plot ZR*Q across the range of gamma resolutions and select solution with the peak ZR*Q

figure 
plot (gamma_range(2:size(gamma_range,2)),ZR_max_mul_Q(1,2:size(gamma_range,2))) 

figure 
plot (gamma_range(2:size(gamma_range,2)),ZR_max(1,2:size(gamma_range,2))) 

i = 1;
%calc number of communities for the top network assignments
for best_gamma_number = 2:size(gamma_range,2)
    M1=M(:,max_ZR_index(1,best_gamma_number),best_gamma_number);
    net_num(i,1) = max(M1);
    i=i+1;
end

figure
plot(gamma_range(2:size(gamma_range,2)),net_num)

%% select highest ZR*Q
best_gamma_number =2; 
M1=M(:,max_ZR_index(1,best_gamma_number),best_gamma_number);


%plot fc matrix reordered according to community membership
[X,Y,INDSORT] = grid_communities(M1);
figure
imagesc(Z_fmriprep_9p_avg(INDSORT,INDSORT))
hold on;
plot (X,Y,'b','linewidth',2);

for x = 1:size(M1,1)
community_matrix_order = M1(INDSORT,1);
end

%% re-arrange M1 for publication

M1_rearr = zeros (size(M1,1),1);
M1_rearr(M1==2) = 4;
M1_rearr(M1==9) = 7;
M1_rearr(M1==10)=1;
M1_rearr(M1==5)= 2;
M1_rearr(M1==1) = 11;
M1_rearr(M1==6) = 5;
M1_rearr(M1==11) = 9;
M1_rearr(M1==8) = 3;
M1_rearr(M1==7) = 6;
M1_rearr(M1==4) = 8;
M1_rearr(M1==3) = 10;

[X,Y,INDSORT] = grid_communities(M1_rearr);
figure
imagesc(Z_fmriprep_9p_avg(INDSORT,INDSORT))
hold on;
%plot (X,Y,'b','linewidth',2);
names_reordered = names_358;(INDSORT);
Z_fmriprep_9p_avg_rearr = Z_fmriprep_9p_avg(INDSORT,INDSORT);
for i = 1:357
    Z_fmriprep_9p_avg_rearr(i,i) = 0;
end


Z_fmriprep_9p_avg_rearr = Z_fmriprep_9p_avg(INDSORT,INDSORT); % rearrange matrix to group regions according to community

