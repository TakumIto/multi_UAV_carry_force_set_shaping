function dexpr = tderiv(expr, t, varargin)
%TDERIV Summary of this function goes here
%   Detailed explanation goes here

    vars = symvar(expr);
    names = arrayfun(@(x)char(x), vars, 'UniformOutput', false);
    
    prefix = 'd';
    order = 1;
    
    if nargin > 2
        vari = 1;
        
        if isnumeric(varargin{vari})
            order = uint32(varargin{vari});
            vari = vari + 1;
        end
        
        if nargin >= vari && strcmp(varargin{vari}, 'Prefix')
            prefix = varargin{vari+1};
            vari = vari + 2;
        end

        if nargin >= vari && strcmp(varargin{vari}, 'Exclude')
            for exvar=[varargin{vari+1}]
                exi = strcmp(names, char(exvar));
                vars(exi)=[];
                names(exi)=[];
            end
        end
    end
    
    dexpr = zeros(size(expr));
    for k=1:size(names, 2)
        var = vars(k);
        name = [prefix names{k}];
        evalin('caller', ['syms ' name]);
        dexpr = dexpr + diff(expr, var) .* evalin('caller', name);
    end
    
    if order > 1
        dexpr = tderiv(dexpr, t, order-1, varargin{2:end});
    end
end

