function Vvar = vee(w)
%VEE Summary of this function goes here
%   Detailed explanation goes here
    tmp = [];
    if isa(w, 'symfun') && size(w,1) == 1 && size(w,2) == 1 % so(3)
        syms t
        tmp = w(t);
    elseif isa(w, 'sym') &&  size(w,1) == 3 && size(w,2) == 3 % so(3)
        tmp = w;
    end
    Vvar = [tmp(3,2); tmp(1,3); tmp(2,1)];
end