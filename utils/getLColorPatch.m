function patch= getLColorPatch(im, patch_size, location)
% Sample color-image patches

if size(im, 3) ~= 3,
	disp('Using inproper color function: getLColorPatch');
end

 xrow = location(:,1);
 ycol = location(:,2);

patch_num = length(xrow);

im = double(im);

patch = zeros((patch_size^2)*3,  length(xrow));
 
for ii = 1:patch_num,    
    row = xrow(ii);
    col = ycol(ii);
    
    I1 = im(row:row+patch_size-1,col:col+patch_size-1,1);
	I1 = I1(:)';
	I2 = im(row:row+patch_size-1,col:col+patch_size-1,2);
	I2 = I2(:)';
	I3 = im(row:row+patch_size-1,col:col+patch_size-1,3);
	I3 = I3(:)';
    Hpatch = [I1 I2 I3]';
       
    patch(:,ii) = Hpatch(:)-mean(Hpatch(:));
end
