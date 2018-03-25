%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling


function NNF_ssd = calc_ssd(ci,cj,A,padded_B,i,j,hp_size)
    % ci,cj patch centre in A

    diff = A(ci-hp_size:ci+hp_size,cj-hp_size:cj+hp_size,:) ...
         - padded_B(i+hp_size-hp_size:i+2*hp_size,j+hp_size-hp_size:j+2*hp_size,:);
    valid_diff = diff(~isnan(diff(:)));
    NNF_ssd = sum(valid_diff.^2)/length(valid_diff); % normalise according to valid pixel area
end