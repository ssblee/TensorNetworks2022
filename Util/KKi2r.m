function yr = KKi2r (xi, yi, varargin)
% < Description >
%
% yr = KKi2r (xi, yi [,option])
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

% Contribution of a piecewise linear segment y = (a*x+b) over [x1,x2] to
% the point at x0 is given by the integral of (a*x+b)/(x-x0) over [x1,x2]:
%   a*(x2-x1) + (a*x0+b)*log(abs((x2-x0)/(x1-x0)))  --- (1)
% Since a = (y2-y1)/(x2-x1), y1 = a*x1+b, and y2 = a*x2+b, the sum of the
% first terms will be yi(end)-yi(1).

% contribution from the second term of Eq. (1)
for it = (1:numel(xi))
    yr(it,:) = sum( (a(1:it-2,:)*xi(it)+b(1:it-2,:)).* ...
                    log((xi(it)-xi(2:it-1,:))./(xi(it)-xi(1:it-2,:))) , 1) + ...
               sum( (a(it+1:end,:)*xi(it)+b(it+1:end,:)).* ...
                    log((xi(it+2:end,:)-xi(it))./(xi(it+1:end-1,:)-xi(it))) , 1);
end
yr(2:end-1,:) = yr(2:end-1,:) + yi(2:end-1,:).*log(abs(dxi(2:end)./dxi(1:end-1)));

% contribution from the first term of Eq. (1)
yr = yr + (yi(end,:)-yi(1,:)); 

% contribution from a 1/x tail (= y1*x1/x stretching from the point (x1,y1) at the edge) to point at x0:
% \int_{x1}^{inf} dx (y1*x1/x) * (1/(x-x0)) = (y1*x1/x0)*log(abs(x1/(x1-x0))) --- (2)
        
% from the right tail
yr(ids(1:end-1),:) = yr(ids(1:end-1),:) + ...
    xi(end)*log(xi(end)./(xi(end)-xi(ids(1:end-1)))).*(yi(end,:)./xi(ids(1:end-1)));
        
% from the left tail: end <-> 1 / takes opposite sign (-1) to the
% contribution to yr(1:end-1) due to opposite integration interval [-inf, x1]
yr(ids(2:end),:) = yr(ids(2:end),:) - ...
    xi(1)*log(xi(1)./(xi(1)-xi(ids(2:end)))).*(yi(1,:)./xi(ids(2:end)));
        
% at zero frequency
yr(xi == 0,:) = yr(xi == 0,:) + (yi(end,:)-yi(1,:));
        
% at the edges: the sum of the second term in Eq.(1) and the term in
% Eq.(2). The divergent terms are cancelled out.
yr(end,:) = yr(end,:) + yi(end,:)*log(xi(end)/dxi(end));
yr(1,:) = yr(1,:) - yi(1,:)*log(abs(xi(1)/dxi(1))); % opposite sign simliarly as for yr(2:end)

yr = yr/pi; % factor 1/pi due to the definition of KK relation

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