% ZeroX has been modified so as to avoid the following error "Error using griddedInterpolant
% The grid vectors are not strictly monotonic increasing."
%
% The following modifications are based off of this MatLab forum answer:
% https://www.mathworks.com/matlabcentral/answers/283718-grid-vectors-not-strictly-monotonic-increasing
function ZC = ZeroX(x,y)
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
for k1 = 1:numel(zxidx)
    idxrng = max([1 zxidx(k1)-1]):min([zxidx(k1)+1 numel(y)]);
    xrng = x(idxrng);
    yrng = y(idxrng);
    
    % Beginning of ZeroX2 modifications. The naming conventions follow
    % those in the referenced MatLab forum, except that "X" is "yrng" and
    % "Y" is "xrng".
    [yrng2, ~, jyrng] = unique(yrng); %yrng is a new array containing the unique values of yrng. jyrng contains the indices in yrng that correspond to the original vector. yrng = yrng2(jyrng)
    xrng2 = accumarray(jyrng, xrng, [], @mean); %This function creates a new array "xrng2" by applying the function "@mean" to all elements in "xrng" that have identical indices in "jyrng". Any elements with identical X values will have identical indices in jyrng. Thus, this function creates a new array by averaging values with identical X values in the original array.
        
    ZC(k1) = interp1( yrng2(:), xrng2(:), 0, 'linear', 'extrap' );
    
end
end
