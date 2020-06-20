function [ initialDictionary ] = InitialDicRealDFT( dic_size, patch_size)
% Initialize a dictionary

    Pn=ceil(sqrt(dic_size));
    DFT = dftmtx(patch_size);
    DFT = kron(DFT',DFT);
    rDFT = real(DFT);
    initialDictionary= rDFT(:,2:dic_size+1);
    % initialDictionary= DCT(:,1:dic_size);

end

