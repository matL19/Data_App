function [Spectra, Wavenumbers, SpectraTitles, Filenames, ...
    SpectraComments] = LoadSpectra (varargin)

%
% LoadSpectra.m
%  
% Imports the absorbance data in .SPA spectrum files into a set of arrays 
% with data from the selected files stored in columns.
%
% Kurt Oldenburg - 06/28/16
% 
if nargin == 0
[Filenames,pathname]=uigetfile({'*.spa','Thermo Spectrum (*.spa)'}, ...
       'MultiSelect','on','Select Spectra Files...');
elseif nargin == 2
    pathname = varargin{1};
    Filenames = varargin{2};
elseif nargin == 3
    pathname = varargin{1};
    fileroot = varargin{2};
    nums = varargin{3};
    
    cd(pathname);
    
    files = dir(fileroot + "*");
    if isnumeric(nums)
        Filenames = {files(nums).name};
%     elseif isstring(nums)
%         for ii = 1:numel(files)
%             if contains(files(ii).name, fileroot + nums(ii))
%                Filenames{ii} = files(ii).name 
%             end
%         end
    end
else
    error("Error: Invalid set of arguments. ...If using arguments, enter (filepath,fileroot,nums)")
end

cd (pathname);  % Change to directory where the spectrum files are.

if ischar(Filenames)== 1           % If only 1 file is selected, Filenames
    NumSpectra = 1;                % is a char instead of a cell of chars,
else                               % which messes up fopen.
    NumSpectra =length(Filenames);
end

for i = 1:NumSpectra 

    DataStart=0;
    CommentStart=0;

    if NumSpectra == 1       
        fid=fopen(Filenames,'r');
    else
        fid=fopen(Filenames{i},'r');
    end;
    
    fseek(fid,30,'bof');
    SpectraTitles(i)={char(nonzeros(fread(fid,255,'uint8'))')};

    fseek(fid,564,'bof');
    Spectrum_Pts=fread(fid,1,'int32');

    fseek(fid,576,'bof');
    Max_Wavenum=fread(fid,1,'single');
    Min_Wavenum=fread(fid,1,'single');

    % The Wavenumber values are assumed to be linearly spaced between
    % between the Min and Max values. The array needs to be flipped 
    % around to get the order lined up with the absorbance data.
    
    Wavenumbers(:,i)=flipud(linspace(Min_Wavenum,...
        Max_Wavenum,Spectrum_Pts).')';

    
    % The starting byte location of the absorbance data is stored in the
    % header. It immediately follows a flag value of 3:
    
    Flag=0; 
    
    fseek(fid,288,'bof');
    
    while Flag ~= 3         
        Flag = fread(fid,1,'uint16');   
    end;
    
    DataPosition=fread(fid,1,'uint16')';
    fseek(fid,DataPosition,'bof');
    
    Spectra(:,i)=fread(fid,Spectrum_Pts,'single');

    % Same story goes for the Comments section with a flag of 4.
    % The size of the section is the difference between the two.
    
%     Flag=0;
%     
%     fseek(fid,288,'bof');
%     
%     while Flag ~= 4         
%         Flag = fread(fid,1,'uint16');   
%     end 
%     
%     CommentPosition=fread(fid,1,'uint16')';
%     SpectraComments(i)={char(nonzeros(fread(fid, ...
%         (DataPosition-DataPosition), 'uint8'))')};
    
    fclose(fid);
    
end;


 
