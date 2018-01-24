function H = dfa(y,x,detrend_order,xstart_xend,box_sizes,plotting)
% Program DFluC : Detrended Fluctuation Analysis
%
%     H = DFA(y, x, detrend_order, [xstart xend], box_sizes, plotting)
%
% estimates the Hurst exponent H of a random walk-like process y (or
% part of it as specified by [xstart xend]).
%
% Input arguments (only required input is y):
%     y is an array of observed values, or may contain some NaN entries
%       (which are interpreted as missing values).
%
%     x is an array indicating the sample points (where the values of y
%       are sampled). If x is not provided, then y is assumed to be
%       regularly sampled.
%
%     detrend_order specifies the degree of polynomials to be used for
%       local detrending. Default is 2.
%
%     [xstart xend] defines the portion of data for which the Hurst
%       exponent is to be estimated. xstart and xend are assumed to have
%       the same unit as x. If [xstart xend] is not provided, then the
%       Hurst exponent will be estimated for the entire range of x.
%
%     box_sizes is an array indicating several box sizes to be used for
%       local polynomial detrending. Each entry is assumed to have the
%       same unit as x. If box_sizes is not provided, then 50 automatically
%       determined box sizes will be used.
%
%     plotting is true if a plot of log(root-mean-square fluctuation)
%       versus log(box size) is to be shown, false otherwise. The slope
%       of the least-squares fit line is the estimated Hurst exponent H.
%       Default is false.
%
%
% Author: Yishul Wei. All rights reserved.
%
% Program DFluC is intended to be an academic software program. Permission
% to use, copy, modify, and distribute the software and its documentation
% for not-for-profit purposes is granted to any person obtaining a copy of
% the source code, provided that this permission notice appear in all
% copies. For other uses, please contact the author (Y. Wei).
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
%
% If you publish results obtained with this program, you may include the
% following citation:
%
%     Colombo, Wei, Ramautar, Linkenkaer-Hansen, Tagliazucchi, Van Someren.
%       More severe insomnia complaints in people with stronger long-range
%       temporal correlations in wake resting-state EEG.
%       Frontiers in Physiology 7 (2016): 576. doi:10.3389/fphys.2016.00576
try
    assert(isnumeric(y));
    
    if (nargin < 2) || isempty(x)
        x = 1:numel(y);
    elseif (numel(x)~=numel(y))
        warning('Sample points do not match data lentgh.');
        H = nan; return
    else
        [x,ord] = sort(x(:),'ascend'); y = y(ord);
    end
    x = x(isfinite(y)); y = y(isfinite(y));
    x = x(:).'; y = y(:).';
    dx = diff(x); dx = min(dx(dx > eps));
    
    if (nargin < 3) || isempty(detrend_order)
        detrend_order = 2;
    else
        detrend_order = unique(detrend_order(~isnan(detrend_order)));
        if numel(detrend_order) ~= 1
            if isempty(detrend_order)
                detrend_order = 2;
            else
                detrend_order = max(detrend_order);
            end
            warning('Ambiguous detrending order specification. Detrending order %d will be used.', detrend_order);
        end
    end
    
    if (nargin < 4) || isempty(xstart_xend)
        xend = max(x); xstart = min(x) - dx;
    else
        xend = max(xstart_xend(:)); xstart = min(xstart_xend(:));
        if sum(~isnan(unique(xstart_xend))) > 2
            warning('Ambiguous data range specification. Range [%g %g] will be used.', xstart, xend);
        end
    end
    
    if all(abs(x-xstart) > eps)
        xbeg = xstart;
    else
        xbeg = xstart - dx;
    end
    msk = (x > xbeg) & (x <= xend);
    if sum(msk) < (3*(detrend_order + 2))
        if (nargin < 4) || isempty(xstart_xend)
            warning('Too few data points. Detrending order %g requires at least %g data points.', ...
                detrend_order, 3*(detrend_order+2));
        else
            warning('The range [%g %g] contains too few data points. Detrending order %g requires at least %g data points.', ...
                xstart, xend, detrend_order, 3*(detrend_order+2));
        end
        H = nan; return
    else
        x = x(msk); y = y(msk);
    end
    total_span = xend - xbeg;
    
    if (nargin < 5) || isempty(box_sizes)
        box_sizes = logspace(log10(iqr(x))-1.001, log10(total_span)-0.5, min(50,sum(msk)));
    else
        box_sizes = unique(box_sizes(isfinite(box_sizes) & (box_sizes > 0)));
        if numel(box_sizes) < 2
            warning('At least 2 valid distinct box sizes are needed.');
            H = nan; return
        end
    end
    box_sizes = box_sizes(:);
    
    if (nargin < 6) || isempty(plotting)
        plotting = false;
    else
        plotting = logical(plotting);
    end
catch
    error('Non-numeric input...');
end

nbox = floor(total_span ./ box_sizes);
logfluc = nan(numel(box_sizes),1);

for k=1:numel(box_sizes)
    box_bound = xbeg + ((total_span - (box_sizes(k)*nbox(k)))/2) + (box_sizes(k)*(0:nbox(k)));
    ibeg = arrayfun(@(bn)min([inf find(x>bn,1,'first')]), box_bound(1:(end-1)));
    iend = arrayfun(@(bn)max([-inf find(x<=bn,1,'last')]), box_bound(2:end));
    
    logfluc(k) = ...
        log10(mean(cell2mat(arrayfun(@(bi,ei)polyfit_residual_sq(x(bi:ei), y(bi:ei), detrend_order), ...
        ibeg(ibeg<iend), iend(ibeg<iend), 'UniformOutput',false)))) / 2;
end

fin = isfinite(logfluc);
if sum(fin) < 2
    warning('Fluctuation function cannot be estimated. Box sizes are probably unfit for the data range.');
    H = nan;
else
    log_box_sizes = log10(box_sizes);
    H = [log_box_sizes(fin) ones(sum(fin),1)] \ logfluc(fin);
    if plotting
        figure; plot(log_box_sizes,logfluc,'*'); pbaspect([1 1 1]); daspect([1 1 1]);
        xlabel 'Log_{10}(Box Size)'; ylabel 'Log_{10}(RMS Fluctuation)';
        mid = (min(log_box_sizes) + max(log_box_sizes)) / 2;
        lsline; text(mid, [mid 1]*H, sprintf('H=%.4f',H(1)), ...
            'HorizontalAlignment','center','Rotation',atan(H(1))*180/pi,'BackgroundColor','w');
    end
    H = H(1);
end
end

function resq = polyfit_residual_sq(x,y,detrend_order)
if numel(x) <= (detrend_order + 1)
    resq = [];
else
    [p,~,mu] = polyfit(x,y,detrend_order);
    resq = abs(y - polyval(p,x,[],mu)) .^ 2;
end
end
