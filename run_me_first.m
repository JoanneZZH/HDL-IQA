basedir = pwd; 

str = computer ; 

if strcmpi(str(1:2),'pc')
    slash = '\' ; 
else
    slash = '/' ; 
end


addpath([basedir slash 'KSVD' slash 'ksvdbox13' ])   ;
addpath([basedir slash 'KSVD' slash 'ompbox10' ])   ;
addpath([basedir slash 'KSVD' slash 'ompbox10' ])   ;


cd([basedir slash 'KSVD' slash 'ksvdbox13' slash 'private']);
make; 

cd([basedir slash 'KSVD' slash 'ompbox10' slash 'private']);
make

cd(basedir) ; 

display('done!');







