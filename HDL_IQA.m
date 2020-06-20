function [ score ] = HDL_IQA( I_ref_ori,I_dist_ori )
% This function includes four steps for the color image and three steps for
% the gray image;
N_ref = ndims(I_ref_ori);
N_dist = ndims(I_dist_ori);
if(N_ref~=N_dist)
    error('The dimension of two images are not matched.');
end

params = [];
if(N_ref == 2)   % if they are gray images
   [ref_dic, dist_dic] = getDic(I_ref_ori,I_dist_ori);
   [ref_grad, dist_grad] = getGrad(I_ref_ori,I_dist_ori);
   [ref_lu, dist_lu] = getLu(I_ref_ori,I_dist_ori);
else
  I_ref = rgb2gray(I_ref_ori);
  I_dist = rgb2gray(I_dist_ori);
  [ref_dic, dist_dic] = getDic(I_ref,I_dist);
  [ref_grad, dist_grad] = getGrad(I_ref,I_dist); 
  [ref_lu, dist_lu] = getLu(I_ref_ori,I_dist_ori);
  [ref_cr, dist_cr, ref_cb, dist_cb] = getCrCb(I_ref_ori,I_dist_ori);
  params.ref_cr = ref_cr;
  params.dist_cr = dist_cr;
  params.ref_cb = ref_cb;
  params.dist_cb = dist_cb;
end
   params.ref_dic = ref_dic;
   params.dist_dic = dist_dic;
   params.ref_grad = ref_grad;
   params.dist_grad = dist_grad;
   params.ref_lu = ref_lu;
   params.dist_lu = dist_lu;
   
   score = getScore(params);
end


%% This version for tarning the dist dictionary when fixing the coef of the ref
%% and guaratee the corressponding relation between atoms of ref_dic and dist_dic
function [ref_dic, dist_dic] = getDic( I_ref,I_dist )
% input: gray images

%   setting used parameters
    patch_size = 8;
    dic_size = 20;
    slide_step = 4;
    ratio = 1;  % default = 0.5(this is not true)
    errT = 0.1;

%    Initialization params for func KSVD
    numIteration = 10;   % default 10
    initialDictionary = InitialDicRealDFT( dic_size ,patch_size);
    paramksvd = [];
    paramksvd.initdict = initialDictionary;
    paramksvd.Edata = errT;
    paramksvd.dictsize = dic_size;
    paramksvd.iternum = numIteration;
    paramksvd.memusage = 'high';

    I_ref_grad  =  imresize(I_ref,ratio);  % just rescaling the size of images
    Y_ref = im2cols_sliding(I_ref_grad, patch_size,slide_step);
    paramksvd.data = Y_ref;
    [ref_dic,ref_coef] = ksvd(paramksvd);
    
    I_dist_grad  = imresize(I_dist,ratio); 
    Y_dist = im2cols_sliding(I_dist_grad, patch_size,slide_step);
    dist_dic = (Y_dist*ref_coef')/(ref_coef*ref_coef'); 

end

function [ ref_grad, dist_grad ] = getGrad( I_ref,I_dist )
     %% get gradient feature
    dx = [-1 0 1];
    dy = [-1 ; 0 ; 1];
    
    I_ref = im2double(I_ref);
    [M, N] = size(I_ref);
    f = max(1,round(min(M,N)/256));
    I_ref = imresize(I_ref, 1/f );
    g_x_r = conv2(I_ref, dx, 'same');
    g_y_r = conv2(I_ref, dy, 'same');
    ref_grad = sqrt(g_x_r.^2 + g_y_r.^2);

    I_dist = im2double(I_dist);
    I_dist = imresize(I_dist, 1/f );
    g_x = conv2(I_dist, dx, 'same');
    g_y = conv2(I_dist, dy, 'same');
    dist_grad = sqrt(g_x.^2 + g_y.^2);

end

function [ref_lu, dist_lu] = getLu(I_ref, I_dist)

  %   setting used parameters
    patch_size = 8;  Tv = 0.7; slt_choice = 'vrc'; patch_num = 1000;
    
    M = size(I_ref,1);
    N = size(I_ref,2);
    f = max(1,round(min(M,N)/256));
    
    I_ref = imresize(I_ref, 1/f );
    I_dist = imresize(I_dist, 1/f );
    n_dim = ndims(I_ref);
    % ------ SELECTION OF MEAN VALUE PAIRS ------


    if (n_dim == 2)
        [Xr, location] = getRefPatch(I_ref,patch_size,patch_num,slt_choice,Tv);
        Xd = getLPatch(I_dist,patch_size,location);
    else
        [Xr, location] = getRefColorPatch(I_dist,patch_size,patch_num,slt_choice,Tv);
        Xd = getLColorPatch(I_dist,patch_size,location);
    end
    Xr = im2double(Xr);        Xd = im2double(Xd);
    mXr = mean(abs(Xr));       mXd = mean(abs(Xd));     % mean value

    ref_lu = mXr-mean(mXr);
    dist_lu = mXd-mean(mXd);

end

function [ ref_cr, dist_cr, ref_cb, dist_cb ] = getCrCb(I_ref,I_dist)
%% get color hue and saturation feature  
  %   setting used parameters
    patch_size = 8;  Tv = 0.7; slt_choice = 'vrc'; patch_num = 1000;
    
    M = size(I_ref,1);
    N = size(I_ref,2);
    f = max(1,round(min(M,N)/256));
    
    I_ref = im2double(I_ref);
    I_ref = rgb2ycbcr(I_ref);
    I_ref = imresize(I_ref, 1/f );
    cb_ref = I_ref(:,:,2);
    cr_ref = I_ref(:,:,3);
    [Y_ref_cb, location] = getRefPatch(cb_ref,patch_size,patch_num,slt_choice,Tv);
    Y_ref_cr = getLPatch(cr_ref,patch_size,location);
    
    I_dist = im2double(I_dist);
    I_dist = rgb2ycbcr(I_dist);
    I_dist = imresize(I_dist, 1/f );
    cb = I_dist(:,:,2);
    cr = I_dist(:,:,3);
    Y_dist_cb = getLPatch(cb,patch_size,location);
    Y_dist_cr = getLPatch(cr,patch_size,location);

    Xr_cb = im2double(Y_ref_cb);     Xd_cb = im2double(Y_dist_cb);
    mXr_cb = mean(abs(Xr_cb));       mXd_cb = mean(abs(Xd_cb));     % mean value
    
    Xr_cr = im2double(Y_ref_cr);     Xd_cr = im2double(Y_dist_cr);
    mXr_cr = mean(abs(Xr_cr));       mXd_cr = mean(abs(Xd_cr));     % mean value
    
    ref_cb = (mXr_cb - mean(mXr_cb));
    dist_cb = (mXd_cb - mean(mXd_cb));
    ref_cr = (mXr_cr - mean(mXr_cr));
    dist_cr = (mXd_cr - mean(mXd_cr));

end

function [ score ] = getScore( params )
%GETSCORE Calculate the final score of a degraded image
ref_dic = params.ref_dic;
dist_dic = params.dist_dic;
ref_grad = params.ref_grad;
dist_grad = params.dist_grad;
ref_lu = params.ref_lu;
dist_lu = params.dist_lu;

%% calculate 6 sub-scores
% the score about dictionaries
score1 = mean(mean((2*ref_dic.*dist_dic+0.5) ./(ref_dic.^2 + dist_dic.^2+0.5)));
score5 = corr(ref_dic(:),dist_dic(:));
diff = ref_dic - dist_dic;
score7 = -log10(sqrt(sum(sum(diff.^2))));

% the score about gradient amplitude
C = 0.1;
c1 = 2*ref_grad.*dist_grad + C;
c2 = ref_grad.^2 + dist_grad.^2 + C;
grad_2 = mean(mean(c1./c2));

% the score about luminance information
Cm = 0.01;
l_score = (sum(ref_lu.*dist_lu)+Cm) / (sqrt(sum(ref_lu.^2)*sum(dist_lu.^2))+Cm);

num_feas = length(params);
if( num_feas==6 )
    score = 0.01*score7 + 0.2*score5 + 0.1*score1 + 0.8*grad_2 + 0.2*l_score;
else
   % the score about color saturation and hue for color images
    mXr_cb = params.ref_cb; mXd_cb = params.dist_cb;
    mXr_cr = params.ref_cr; mXd_cr = params.dist_cr;
    c = 0.01;
    c1 = 2*mXr_cb.*mXd_cb + c;
    c2 = mXr_cb.^2 + mXd_cb.^2 + c;
    c3 = c1./c2;
    c4 = 2*mXr_cr.*mXd_cr + c;
    c5 = mXr_cr.^2 + mXd_cr.^2 + c;
    c6 = c4./c5;
    crcb_score = mean(c3+c6);

    % calculate final score using one group of hyper
    score = 0.01*score7 + 0.2*score5 + 0.1*score1 + 0.8*grad_2 + 0.2*l_score+ 0.5*crcb_score;
    score = score/3;
end

end