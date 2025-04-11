%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ad_f^i g を計算します
% g としては vector か distribution (基底を並べた cell 配列) のいずれかを指定できます
% g として distiribution を指定した場合は, 
% g のすべての基底 g_j に対して ad_f^i g_j を計算し, 
% ad_f^i g_j を新たな基底とする distribution (cell 配列) を返します
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = ad(i, x, f, g)

if isa(g, 'double') || isa(g, 'sym')
	out = ad_vector(i, x, f, g);
elseif isa(g, 'cell')
	out = ad_distribution(i, x, f, g);
else
	out = nan;
end

end

function out = ad_vector(i, x, f, g)

if i > 0
	out = LieBracket(x, f, ad_vector(i-1, x, f, g));
else
	out = g;
end

end

% (ゼロベクトルを基底としてしまう可能性あり)
function out = ad_distribution(i, x, f, G)
dimG = size(G, 2);	% dimension of G

out = cell([1, dimG]);
for j = 1:dimG
	out{1, j} = ad_vector(i, x, f, G{1, j});
end

end
