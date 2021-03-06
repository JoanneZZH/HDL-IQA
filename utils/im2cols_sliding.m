function [ C ] = im2cols_sliding(I, sz_patch, sz_step)
%% im2cols Rearrange image blocks into columns.
%  [ C ] = im2cols(I, sz_patch) rearranges the image patches into columns.
%   The image patches are sampled using a sliding window. Only square
%   patch is supported.
%   Input
%       - I: 2D image matrix
%       - sz_patch: size of image patch
%   Output
%       - C: columnss
%   See also COLS2IM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use of this code is free for research purposes only.
%
%Author:  Yuhui Quan
%
%Last Revision: 22-Jun-2014
 %
 %
    if numel(sz_patch)>1, sz_patch = min(sz_patch); end
%     sz_step = 1;        %%%SLIDING
    [n_row, n_col] = size(I);
    C  = zeros( sz_patch^2,(ceil((n_row - sz_patch )/ sz_step) + 1) * ( ceil( (n_col - sz_patch) / sz_step) + 1) );
    cnt = 1; 
    
    for j = 1: sz_step: (ceil((n_col - sz_patch) / sz_step)*sz_step + 1)
        for i = 1: sz_step: (ceil((n_row - sz_patch) / sz_step)*sz_step + 1)
            P = I(i: i + sz_patch - 1, j: j + sz_patch - 1);
            C (:,cnt) = P(:); 
            cnt = cnt+1;
        end;
    end;
    
   C = C - repmat(mean(C,1),[size(C,1) 1]);

end