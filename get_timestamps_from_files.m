function time_points = get_timestamps_from_files(path,varargin)
cd(path)
if nargin == 1
    filenames = uigetfile({'*.spa','Thermo Spectrum (*.spa)'}, ...
        'MultiSelect','on','Select Spectra Files...');
else
    while numel(varargin) >= 2
        var = varargin{1};
        val = varargin{2};
        switch var
            case "relative"
                rel_time = val;
            case "fileroot"
                fileroot = val;
            case "filenums"
                filenums = val;
            otherwise
                error(var + " is not a valid name/value pair argument.")
        end
        varargin = varargin(3:end);
    end
    
    files = dir(fileroot + "*");
    filenames = {files(filenums).name};
end

if ~exist('rel_time','var')
    rel_time = false;
end

times = [];
for ii = 1:numel(filenames)
    g = dir(filenames{ii});
    times = [times datetime(g.date)];
end

if rel_time
    sec = [];
    for ii = 1:numel(times)
        sec = [sec seconds(times(ii) - times(1))];
    end
    time_points = sec;
else
    time_points = times;
end
end