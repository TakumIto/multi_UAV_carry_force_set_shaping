function Wvar = wedge(v)
%WEDGE Summary of this function goes here
%   Detailed explanation goes here
    tmp = [];
    
%     if isa(v, 'sym') && xor(size(v,1) == 3, size(v,2) == 3) % so(3)
    if xor(size(v,1) == 3, size(v,2) == 3) % so(3)
        tmp = v;      
    end
    tmp = v; 

    Wvar = [0 -tmp(3) tmp(2);
        tmp(3) 0 -tmp(1);
        -tmp(2) tmp(1) 0];
end

