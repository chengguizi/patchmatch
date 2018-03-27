%%% Cheng Huimin
%%% A0138497M
%%% EE4212 Assignment: Non-Parametric Sampling

% NNF should has the same size as B
function visualise_correspondance_vectors(NNF,B,A,pass,p_size)

max_height = max(size(B,1),size(A,1));
max_width = max(size(B,2),size(A,2));
max_length = norm([max_height,max_width]);

hsv_map = ones(size(B));
% :1 = hue
%hsv_map(:,:,2) = hsv_map(:,:,2)*0; % :2 = saturation
hsv_map(:,:,3) = hsv_map(:,:,3)*0.7; % :3 = brightness

[row,column,~] = size(NNF);

for i=1:row
    for j=1:column
        offset = squeeze(NNF(i,j,:)) - [i;j];
        mag = norm(offset);
        dir = atan2(offset(2) , offset(1))+pi;
        hsv_map(i,j,2) = mag / max_length;
        hsv_map(i,j,1) = dir / (2*pi);
    end
end

rgb_map = hsv2rgb(hsv_map)*255;

figure(2)
imshow( uint8(rgb_map) );

filename = ['covec_map_pass',num2str(pass),'_psize',num2str(p_size),'.png'];
imwrite(uint8(rgb_map),filename);

end