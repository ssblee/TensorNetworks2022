function yr = KKi2r_Ex (xi, yi, varargin)
% < Description >
%
% yr = KKi2r_Ex (xi, yi [,option])
%
% Convert the imaginary part of a causal function into the real part, by
% using the Kramers-Kronig relation. If the input is the real part, then 
% yi = -KKi2r(xr,yr) yields the imaginary part.
% This function computes the Cauchy principal value (also called principal
% value integral)
%    yr (xr) = (1/pi) * P. \int_{-\infty}^{\infty} dx yi(x) / (x-xr) ,
% where P. means the principal value and yi(x), yr(xr) are the functions
% corresponding to the discrete data (xi, yi), (xi, yr), respectively. The
% integration is done numerically, assuming that the input function is the
% linear interpolation among the points (xi, yi) on the interval [xi(1),
% xi(end)], and have 1/x tails on (-\infty, xi(1)) and (xi(end), \infty)
% extrapolated from (xi(1), yi(1) and (xi(end), yi(end)), respectively.
%
% < Input >
% xi : [vector] Grid points of x values on which the imaginary part(s) of
%       the causal function(s) are defined.
% yi : [vector, matrix, or multi-dimensional array] The imaginary part(s)
%       of the causal function(s) on the grid given by "xi". Each column of
%       "yi" corresponds to a function, representing the imaginary part
%       of a function as a piecewise linear function connecting the (x,y)
%       pairs specified xi and the column elements. numel(xi) == size(yi,1)
%       should hold.
%
% < Output >
% yr : [vector, matrix, or multi-dimensional array] The real part(s) of the
%       causal function(s) on the grid specifed by xi. It has the same
%       structure as the input "yi"; each column of "yr" corresponds to the
%       column of yi with the same column index.
%
% Rewritten by S.Lee (Oct.28,2022): Revised for the TN lecture course at
%       SNU.

% % % sanity check
if isvector(xi)
    xi = xi(:);
else
    error('ERR: input frequency ''xi'' is not a vector.');
end
if isvector(yi)
    yi = yi(:);
end
if numel(xi) ~= size(yi,1)
    error('ERR: numel(xi) ~= size(yi,1)');
end

sz = size(yi); % save the size for future reshaping
yi = yi(:,:); % make yi to matrix; reshape later

if any(diff(xi) < 0)
    disptime('WRN: input ''xi'' is not in ascending order.');
    [xi,ids] = sort(xi,'ascend');
    yi = yi(ids,:);
end
if any(diff(xi) == 0)
    error('ERR: Input ''xi'' is not unique.');
end
% % % % 

dxi = diff(xi);
a = (yi(2:end,:) - yi(1:end-1,:))./dxi; % piecewise slope
b = (-xi(1:end-1).*yi(2:end,:)+xi(2:end).*yi(1:end-1,:))./dxi; % piecewise y intercept

ids = find(xi ~= 0); % to be used for adding the tail contributions

yr = zeros(size(yi)); % result array

% % % % TODO (start) % % % %

% perform the numerical integration, assuming that the input function is the
% linear interpolation among the points (xi, yi) on the interval [xi(1),
% xi(end)], and have 1/x tails on (-\infty, xi(1)) and (xi(end), \infty)
% extrapolated from (xi(1), yi(1) and (xi(end), yi(end)), respectively.

% % % % TODO (end) % % % %

if ~isequal(size(yr),sz)
    yr = reshape(yr,sz);
end

% catch e
%     disp2(getReport(e))
%     disp2('Let''s Debug!');
%     keyboard;
% end

end