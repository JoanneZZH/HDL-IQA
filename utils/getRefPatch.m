function [Patch, location] = getRefPatch(im, patch_size, patch_num,slt_choice,Tv)
% Sample gray-image patches

if size(im, 3) == 3,
    Img = rgb2gray(im);
else
    Img = im;
end

[nrow, ncol] = size(Img);

Img = double(Img);
Y = im2col(Img,[patch_size patch_size], 'distinct');
% Y = im2col(Img,[patch_size.^2 nrow*ncol ], 'distinct');
% pvars = var(Y,0,1);
pvars = var(Y);
meanVar = median(pvars);

Patch = [];
x1 = []; x2 = [];

num = (nrow+patch_size-1)*(ncol+patch_size-1);
xrow = floor(rand(num,1)*(nrow - patch_size))+1;
ycol = floor(rand(num,1)*(ncol - patch_size))+1;

idx = 1;
num = 0;

if (strcmp('rnd',slt_choice))
    for ii = 1:patch_num
        row = xrow(ii);
        col = ycol(ii);
        Hpatch = Img(row:row+patch_size-1,col:col+patch_size-1);
        patch = Hpatch(:) - mean(Hpatch(:));
        Patch = [Patch patch];
        x1 = [x1 row];
        x2 = [x2 col];
    end
    
else
    while 1
        row = xrow(idx);
        col = ycol(idx);
        Hpatch = Img(row:row+patch_size-1,col:col+patch_size-1);
        if (strcmp('etp',slt_choice))
            flag = (localEntropy(Hpatch)>0);
        else
            flag = (var(Hpatch,0,1)>(meanVar*Tv));
        end
        if(flag)
            patch = Hpatch(:)-mean(Hpatch(:));
            Patch = [Patch patch];
            x1 = [x1 row];
            x2 = [x2 col];
            num = num + 1;
        end
        
        if(num==patch_num)
            break;
        end
        if(idx==numel(xrow))
            pnum = size(Patch,2);
            if (pnum<patch_num)
                for kk = pnum+1:patch_num
                    row = xrow(kk);
                    col = ycol(kk);
                    Hpatch = Img(row:row+patch_size-1,col:col+patch_size-1);
                    patch = Hpatch(:) - mean(Hpatch(:));
                    Patch = [Patch patch];
                    x1 = [x1 row];
                    x2 = [x2 col];
                end
            end
            break;
        end      
        idx = idx + 1;
    end
end

location(:,1) = x1(:);
location(:,2) = x2(:);
