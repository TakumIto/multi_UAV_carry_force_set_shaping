function Ad = adtrans(R, p)
%ADTRANS Summary of this function goes here
%   Detailed explanation goes here
    Ad = [R wedge(p)*R;
        zeros(3) R];
end