function [c,ww] = smooth(varargin)
% Pared down version of Matlab curviefit version
% SMOOTH  Smooth data.
%   Z = SMOOTH(Y) smooths data Y using a 5-point moving average.
%
%   Z = SMOOTH(Y,SPAN) smooths data Y using SPAN as the number of points used
%   to compute each element of Z.

if nargin < 1
    error(message('curvefit:smooth:needMoreArgs'));
end

[varargin{:}] = convertStringsToChars(varargin{:});

if nargout > 1 % Called from the GUI cftool
    ws = warning('off', 'all'); % turn warning off and record the previous warning state.
    [lw,lwid] = lastwarn;
    lastwarn('');
else
    ws = warning('query','all'); % Leave warning state alone but save it so resets are no-ops.
end

% is x given as the first argument?
if nargin==1 || ( nargin > 1 && (length(varargin{2})==1 || ischar(varargin{2})) )
    % smooth(Y) | smooth(Y,span,...) | smooth(Y,method,...)
    is_x = 0; % x is not given
    y = varargin{1};
    y = y(:);
    x = (1:length(y))';
else % smooth(X,Y,...)
    is_x = 1;
    y = varargin{2};
    x = varargin{1};
    y = y(:);
    x = x(:);
end

% is span given?
span = [];
if nargin == 1+is_x || ischar(varargin{2+is_x})
    % smooth(Y), smooth(X,Y) || smooth(X,Y,method,..), smooth(Y,method)
    is_span = 0;
else
    % smooth(...,SPAN,...)
    is_span = 1;
    span = varargin{2+is_x};
end

t = length(y);
if t == 0
    c = y;
    ww = '';
    if nargout > 1
        ww = lastwarn;
        lastwarn(lw,lwid);
        warning(ws);  % turn warning back to the previous state.
    end
    return
elseif length(x) ~= t
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:XYmustBeSameLength'));
end

% realize span
if span <= 0
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:spanMustBePositive'));
end
if span < 1, span = ceil(span*t); end % percent convention
if isempty(span), span = 5; end % smooth(Y,[],method)

idx = 1:t;

sortx = any(diff(isnan(x))<0);   % if NaNs not all at end
if sortx || any(diff(x)<0) % sort x
    [x,idx] = sort(x);
    y = y(idx);
end

if islogical(y)
    y = double(y);
end

c = NaN(size(y),'like',y);

ok = ~isnan(x);
c(ok) = moving(x(ok),y(ok),span);
 
c(idx) = c;

if nargout > 1
    ww = lastwarn;
    lastwarn(lw,lwid);
    warning(ws);  % turn warning back to the previous state.
end

%--------------------------------------------------------------------
function c = moving(x,y, span)
% moving average of the data.

ynan = isnan(y);
span = floor(span);
n = length(y);
span = min(span,n);
width = span-1+mod(span,2); % force it to be odd
xreps = any(diff(x)==0);
if width==1 && ~xreps && ~any(ynan), c = y; return; end
if ~xreps && ~any(ynan)
    % simplest method for most common case
    c = filter(ones(width,1)/width,1,y);
    cbegin = cumsum(y(1:width-2));
    cbegin = cbegin(1:2:end)./(1:2:(width-2))';
    cend = cumsum(y(n:-1:n-width+3));
    cend = cend(end:-2:1)./(width-2:-2:1)';
    c = [cbegin;c(width:end);cend];
elseif ~xreps
    % with no x repeats, can take ratio of two smoothed sequences
    yy = y;
    yy(ynan) = 0;
    nn = double(~ynan);
    ynum = moving(x,yy,span);
    yden = moving(x,nn,span);
    c = ynum ./ yden;
else
    % with some x repeats, loop
    notnan = ~ynan;
    yy = y;
    yy(ynan) = 0;
    c = zeros(n,1,'like',y);
    for i=1:n
        if i>1 && x(i)==x(i-1)
            c(i) = c(i-1);
            continue;
        end
        R = i;                                 % find rightmost value with same x
        while(R<n && x(R+1)==x(R))
            R = R+1;
        end
        hf = ceil(max(0,(span - (R-i+1))/2));  % need this many more on each side
        hf = min(min(hf,(i-1)), (n-R));
        L = i-hf;                              % find leftmost point needed
        while(L>1 && x(L)==x(L-1))
            L = L-1;
        end
        R = R+hf;                              % find rightmost point needed
        while(R<n && x(R)==x(R+1))
            R = R+1;
        end
        c(i) = sum(yy(L:R)) / sum(notnan(L:R));
    end
end