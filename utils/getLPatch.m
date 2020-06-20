function patch= getLPatch(im, patch_size, location)
% Sample gray-image patches

if size(im, 3) == 3,
    Img = rgb2gray(im);
else
    Img = im;
end

 xrow = location(:,1);
 ycol = location(:,2);

patch_num = length(xrow);

Img = double(Img);

patch = zeros(patch_size^2,  length(xrow));
 
for ii = 1:patch_num,    
    row = xrow(ii);
    col = ycol(ii);
    
    Hpatch = Img(row:row+patch_size-1,col:col+patch_size-1);
       
    patch(:,ii) = Hpatch(:)-mean(Hpatch(:));
end

