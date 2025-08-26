function baseline_corrected_data = baselineCorrect(varargin)
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
    baseline_corrected_data = data - y0;
elseif method == "point val" && ~exist('x','var')
    error("An x-axis is needed to do a point-correction.")
end
end