function [Patch, location] = getRefColorPatch(im, patch_size, patch_num,slt_choice,Tv)
% Sample color-image patches

if (ndims(im)==3)
Img = rgb2gray(im);
[nrow, ncol] = size(Img);
else
	disp('Using inproper color function: getRefColorPatch');
    exit();
end
im = double(im);
Y1 = im2col(im(:,:,1),[patch_size.^2 nrow*ncol ], 'distinct');
Y2 = im2col(im(:,:,2),[patch_size.^2 nrow*ncol ], 'distinct');
Y3 = im2col(im(:,:,3),[patch_size.^2 nrow*ncol ], 'distinct');
Y = [Y1' Y2' Y3']';
pvars = var(Y,0,1);
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
		I1 = im(row:row+patch_size-1,col:col+patch_size-1,1);
		I1 = I1(:)';
		I2 = im(row:row+patch_size-1,col:col+patch_size-1,2);
		I2 = I2(:)';
		I3 = im(row:row+patch_size-1,col:col+patch_size-1,3);
		I3 = I3(:)';
        Hpatch = [I1 I2 I3]';
        patch = Hpatch(:) - mean(Hpatch(:));
        Patch = [Patch patch];
        x1 = [x1 row];
        x2 = [x2 col];
    end
    
else
    while 1
        row = xrow(idx);
        col = ycol(idx);
		I1 = im(row:row+patch_size-1,col:col+patch_size-1,1);
		I1 = I1(:)';
		I2 = im(row:row+patch_size-1,col:col+patch_size-1,2);
		I2 = I2(:)';
		I3 = im(row:row+patch_size-1,col:col+patch_size-1,3);
		I3 = I3(:)';
        Hpatch = [I1 I2 I3]';
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
                   	I1 = im(row:row+patch_size-1,col:col+patch_size-1,1);
					I1 = I1(:)';
					I2 = im(row:row+patch_size-1,col:col+patch_size-1,2);
					I2 = I2(:)';
					I3 = im(row:row+patch_size-1,col:col+patch_size-1,3);
					I3 = I3(:)';
					Hpatch = [I1 I2 I3]';
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

