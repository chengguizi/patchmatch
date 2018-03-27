%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling

clear; close all;

A = double(imread('imageA.jpg')); % source
%A = imresize(A,0.5);
B = double(imread('imageB.jpg')); % target / to be reconstructed
%B = imresize(B,0.5);

p_size=5; % must be odd, patch size
hp_size=(p_size-1)/2;
num_pass = 4;

[row_A,column_A,channel_A] = size(A);
[row_B,column_B,channel_B] = size(B);

radius0 = max(row_A,column_A)/2;
alpha = 0.5;

% regeneration patch
regen_p_size = 5; % should be smaller than p_size
regen_hp_size = (regen_p_size-1)/2;
patch_construct = false;

%% memory declaration
tic
% create a padded version of B (target)
% non-active region is marked as NAN, hence excluded from SSD calculation
% later on
padded_B = nan(row_B+p_size-1,column_B+p_size-1,channel_B);
i_start = 1+hp_size;
i_end = row_B+hp_size;
j_start = 1+hp_size;
j_end = column_B+hp_size;

padded_B(i_start:i_end , j_start:j_end,:) = B;
%imshow(uint8(padded_B));

% create near-neighbour field matrix
NNF = cat(3, ...
    randi([1+hp_size,row_A-hp_size],[row_B,column_B]), ...
    randi([1+hp_size,column_A-hp_size],[row_B,column_B]));

%% one-time initialisation
NNF_ssd = inf(row_B,column_B);
for i= 1:row_B
    for j= 1:column_B
        NNF_ssd(i,j) = calc_ssd(NNF(i,j,1),NNF(i,j,2),A,padded_B,i,j,hp_size);
    end
end

%% iteration steps
visualise_correspondance_vectors(NNF,B,A,0,p_size);
for pass = 1:num_pass
    %% even/odd logic
    if mod(pass,2) ~= 0 % is odd
        fprintf('[%d]odd propagation...\n',pass);
        i_seq = 1:row_B;
        j_seq = 1:column_B;
    else
        fprintf('[%d]even propagation...\n',pass);
        i_seq = row_B:(-1):1;
        j_seq = column_B:(-1):1;
    end
    
    %% propagation
        
    for i=i_seq
        for j=j_seq

            if mod(pass,2) ~= 0 % is odd
                
                top_ssd = NNF_ssd(max(i-1,1),j);
                left_ssd = NNF_ssd(i,max(j-1,1));
                current_ssd = NNF_ssd(i,j);

                valid_top = (NNF(max(i-1,1),j,1) + 1 <= row_A - hp_size); % assert shifting 1px will not out-of-bound
                valid_left = (NNF(i,max(j-1,1),2) + 1 <= column_A - hp_size); % same as above

                if (valid_top && top_ssd<current_ssd && top_ssd < left_ssd) % propagate from top
                    NNF(i,j,1) = NNF(i-1,j,1)+1;
                    NNF(i,j,2) = NNF(i-1,j,2);
                    NNF_ssd(i,j) = calc_ssd(NNF(i,j,1),NNF(i,j,2),A,padded_B,i,j,hp_size);
                elseif (valid_left && left_ssd<top_ssd && left_ssd < current_ssd) % propagate from left
                    NNF(i,j,1) = NNF(i,j-1,1);
                    NNF(i,j,2) = NNF(i,j-1,2)+1;
                    NNF_ssd(i,j) = calc_ssd(NNF(i,j,1),NNF(i,j,2),A,padded_B,i,j,hp_size);
                end

            else % is even           

                bottom_ssd = NNF_ssd(min(i+1,row_B),j);
                right_ssd = NNF_ssd(i,min(j+1,column_B));
                current_ssd = NNF_ssd(i,j);

                valid_bottom = (NNF(min(i+1,row_B),j,1) - 1 > hp_size); % assert shifting 1px will not out-of-bound
                valid_right = (NNF(i,min(j+1,column_B),2) - 1 > hp_size); % same as above

                if (valid_bottom && bottom_ssd<current_ssd && bottom_ssd < right_ssd) % propagate from top
                    NNF(i,j,1) = NNF(i+1,j,1)-1;
                    NNF(i,j,2) = NNF(i+1,j,2);
                    NNF_ssd(i,j) = calc_ssd(NNF(i,j,1),NNF(i,j,2),A,padded_B,i,j,hp_size);
                elseif (valid_right && right_ssd<bottom_ssd && right_ssd < current_ssd) % propagate from left
                    NNF(i,j,1) = NNF(i,j+1,1);
                    NNF(i,j,2) = NNF(i,j+1,2)-1;
                    NNF_ssd(i,j) = calc_ssd(NNF(i,j,1),NNF(i,j,2),A,padded_B,i,j,hp_size);
                end
            end % end if
            
            %% random search
            r_pass = 0;
            ci = NNF(i,j,1);
            cj = NNF(i,j,2);
            while true
                r = floor(alpha^r_pass*radius0);
                if r < 1
                    break;
                end

                % setting up the search boundary

                i_min = max(ci-r,1+hp_size);
                i_max = min(ci+r,row_A-hp_size);
                j_min = max(cj-r,1+hp_size);
                j_max = min(cj+r,column_A-hp_size);

                newi = randi([i_min,i_max]);
                newj = randi([j_min,j_max]);

                new_ssd = calc_ssd(newi,newj,A,padded_B,i,j,hp_size);

                if (new_ssd<NNF_ssd(i,j)) % new minima found!
                    NNF(i,j,1) = newi;
                    NNF(i,j,2) = newj;
                    NNF_ssd(i,j) = new_ssd;
                end

                r_pass = r_pass + 1;
            end
        end % end j
    end % end i
        
    % reconstruction debug
    reconstruct_B = zeros(size(B));
    
    for i=1:row_B
        for j=1:column_B
            reconstruct_B(i,j,:) = A(NNF(i,j,1),NNF(i,j,2),:);
        end
    end
    figure(1);
    imshow(uint8(reconstruct_B));
    
    % debug
    visualise_correspondance_vectors(NNF,B,A,pass,p_size);
    
    drawnow
    
    
    
end % end for pass

%% reconstruction

td = toc;

fprintf('time elapsed: %.2f s\n',td);

reconstruct_B = zeros(size(B));

if patch_construct
    for i=1+regen_hp_size:regen_p_size:row_B
        for j=1+regen_hp_size:regen_p_size:column_B
            trim_i = max (i+regen_hp_size - row_B,0);
            trim_j = max (j+regen_hp_size - column_B,0);
            reconstruct_B( i-regen_hp_size:i+regen_hp_size-trim_i ,j-regen_hp_size: j+regen_hp_size-trim_j,:) ...
                = A( NNF(i,j,1)-regen_hp_size:NNF(i,j,1)+regen_hp_size-trim_i , NNF(i,j,2)-regen_hp_size:NNF(i,j,2)+regen_hp_size-trim_j,:);
        end
    end
else
    for i=1:row_B
        for j=1:column_B
            reconstruct_B(i,j,:) = A(NNF(i,j,1),NNF(i,j,2),:);
        end
    end
end
reconstruct_B_uint8 = uint8(reconstruct_B);
figure(1);
imshow(reconstruct_B_uint8);

imwrite(reconstruct_B_uint8,'reconstructed.png');
