function [dicom, dicomvid, inf, dimx, dimy, dimz, no_dyn] = readDicom(location,number)
%% Read DICOM files 
% Version 25092015 

%% Initiate
if nargin<2
    number=0;
end

str_location=[location,'/*.dcm'];
s=dir(str_location);
inf=dicominfo([location,'/',s(1).name]);
dimx=inf.Width;
dimy=inf.Height;
dimz=1;
avgsize=mean([s.bytes]);
ss=s;

for j=1:size(s,1)
    if s(j).bytes > 2*avgsize
        ss(j)=[];
    end
end

s=ss;
no_dyn=size(s,1);

if number==0
    dicom=zeros(dimx,dimy,1,size(s,1)); % preallocate memory
    
    for j=1:size(s,1)
        dicom(:,:,1,j)=dicomread([location,'/',s(j).name]);
    end
else
    dicom=dicomread([location,'/',s(number).name]); % Just read one dicom file
end


%% Video mode
dicomvid=dicom/max(max(max(dicom)));

%% END of function.
end
