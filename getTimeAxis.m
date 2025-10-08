function t = getTimeAxis(varargin)
if nargin == 1
    filenames = uigetfile({'*.spa','Thermo Spectrum (*.spa)'}, ...
        'MultiSelect','on','Select Spectra Files...');
elseif nargin == 3
    pathname = varargin{1};
    fileroot = varargin{2};
    nums = varargin{3};
    
    cd(pathname);
    
    files = dir(fileroot + "*");
    filenames = {files(nums).name};
else
    error("Error: Invalid set of arguments. ...If using arguments, enter (filepath,fileroot,nums)")
end
times = [];
for ii = 1:numel(filenames)
    g = dir(filenames{ii});
    times = [times datetime(g.date)];
end
sec = [];
for ii = 1:numel(times)
    sec = [sec seconds(times(ii) - times(1))];
end

t = sec;

end