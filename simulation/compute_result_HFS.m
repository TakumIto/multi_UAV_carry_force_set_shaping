addpath .. function
ts = 0.1;
reduced_idx = 1:10:length(time);

Mp_tmp = @(gm) Mp_fun(gm);
HFSs = cell(1,length(reduced_idx));
parfor j=1:length(reduced_idx)
    i=reduced_idx(j);
    HFS = compute_HFS(Mp_tmp,f_max,gamma_(i,:)');
    HFSs{j} = HFS;
    disp(j/length(reduced_idx))
end

save('simlog/arc_HFS.mat',"HFSs")

%% 
Verts = cell(1,length(reduced_idx));
parfor j=1:length(reduced_idx)
    HFS = HFSs{j};
    try
        [k,~] = convhull(HFS(:,1),HFS(:,2),HFS(:,3),'Simplify',true);
        Verts{j} = k;
    catch
        Verts{j} = [];
    end
    disp(j/length(reduced_idx))
end

save('simlog/arc_HFS.mat',"HFSs","Verts")