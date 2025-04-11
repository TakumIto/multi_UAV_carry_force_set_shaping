function HFS = compute_HFS(M_ftau_fun,f_max,opt_angles)
%COMPUTE_HFS Summary of this function goes here
%   Detailed explanation goes here

Mp_val = double(M_ftau_fun(opt_angles));
M_null = null(Mp_val(4:6,:));

facet_axes = nchoosek(1:16,13);
ells = table2array( ...
    combinations([0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]) ...
    );
HFS_vertex = NaN(2^16,16);
n_vrt = 0;

for i=1:length(facet_axes)
    facet_index = facet_axes(i,:);
    M_facet = M_null(facet_index, :);
    remain_index = 1:16;
    remain_index(facet_index) = [];
    M_remain = M_null(remain_index, :);
    if rank(M_facet) == 13
        % facet_index
        for ell = ells'
            a = M_facet \ ell;
            ell_remain = M_remain * a;
            if 0 <= min(ell_remain) && max(ell_remain) <= 1
                n_vrt = n_vrt + 1;
                HFS_vertex(n_vrt, facet_index) = ell;
                HFS_vertex(n_vrt, remain_index) = ell_remain;
            end
        end
    end
    % if rem(i, 50) == 0
    %     disp(i/length(facet_axes))
    % end
end
HFS_vertex = rmmissing(HFS_vertex);
HFS = (f_max*Mp_val(1:3,:)*HFS_vertex')';

end

