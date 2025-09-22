function baseline_corrected_data = baselineCorrect(varargin)
% baseline corrects the data
% if 1 input, baselineCorrect(data)
%       corrects to the minimum point in the data
% if 2 inputs, baselineCorrect(x,data)
%       corrects to the minimum point in the data
% if >= 3 inputs, baselineCorrect(x,data,"name",value)
%       name/value pairs:
%          "value" - baseline corrects to the specific value of the x-axis

if nargin <= 0
    error("Function baselineCorrect requires at least 1 input.")
elseif nargin == 1
    y = varargin{1};
    method = "min";
elseif nargin == 2
    x = varargin{1};
    y = varargin{2};
    method = "min";
else
    x = varargin{1};
    y = varargin{2};
    varargin = varargin(3:end);
    
    method = "min";
    correction_pt = 0;
    while numel(varargin) >= 2
        var = varargin{1};
        val = varargin{2};
        switch var
            case "value"
                correction_pt = val;
                method = "point val";
            otherwise
                error("Invalid name/value pairs!")
        end
        varargin = varargin(3:end);
    end
end

if method == "min"
    baseline_corrected_data = y - min(y);
elseif method == "point val" && exist('x','var')
    dx = abs(x(2) - x(1));
    l = x > correction_pt - dx/2 & x < correction_pt + dx/2;
    y0 = y(l);
    if numel(y0) ~= 1
        error("Correction point value yielded more than 1 index.")
    end
    baseline_corrected_data = y - y0;
elseif method == "point val" && ~exist('x','var')
    error("An x-axis is needed to do a point-correction.")
end
end