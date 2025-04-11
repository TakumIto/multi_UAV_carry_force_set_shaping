%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% span(g1, g2, ...) を計算します
% 戻り値は distribution の基底を並べた cell 配列です
% 引数は次元の等しい vector か distribution (cell配列) とします
% 引数の個数に制限はなく, vector と distribution が混在していても問題ありませんが
% distribution を指定する場合は, その基底がすべて独立となるように選んでください(問題なくしたい)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spanD, exclusion] = span(varargin)

num = length(varargin);

out_tmp = varargin{1};
exclusion = nan;
if num >= 2
	for i=2:num
		[out_tmp, exc_tmp] = span2(out_tmp, varargin{i});
		if ~isnan(exc_tmp)
			exclusion = exc_tmp;
		end
	end
end
if isa(out_tmp, 'cell')
	spanD = out_tmp;
elseif isa(out_tmp, 'double') || isa(out_tmp, 'sym')
	spanD = {out_tmp};
else
	spanD = nan;
end

end

%%
function [spanD, exclusion] = span2(D1, D2)

if isa(D1, 'cell') && isa(D2, 'cell')
	% D1, D2: distribution
	[spanD, exclusion] = span_dd(D1, D2);
elseif isa(D1, 'cell') && (isa(D2, 'double') || isa(D2, 'sym'))
	% D1: distribution, D2: vector
	[spanD, exclusion] = span_dv(D1, D2);
elseif isa(D2, 'cell') && (isa(D1, 'double') || isa(D1, 'sym'))
	% D1: vector, D2: distribution
	[spanD, exclusion] = span_dv(D2, D1);
elseif (isa(D1, 'double') || isa(D1, 'sym')) ...
		&& (isa(D2, 'double') || isa(D2, 'sym'))
	% D1, D2: vector
	if rank([D1, D2]) == 2
		spanD = [{D1}, {D2}];
		exclusion = 1;
	else
		spanD = {D1};
		exclusion = nan;
	end
else
	spanD = nan;
	exclusion = nan;
end

end

%% D: distribution, f: vector
function [spanD, exclusion] = span_dv(D, f)
dimD = size(D, 2);	% dimension of D1
n = length(f);		% number of state variables
Dmat = vpa(zeros(n, dimD));
for i = 1:dimD
	Dmat(:, i) = D{1, i};
end
if rank([Dmat, f]) > dimD
	spanD = [D, {f}];
	exclusion = nan;
else
	spanD = D;
	exclusion = 1;
end

end

%% D1, D2: distribution
function [spanD, exclusion] = span_dd(D1, D2)
Dout = D1;
dimD1 = size(D1, 2);	% dimension of D1
dimD2 = size(D2, 2);	% dimension of D2
n = length(D1{1, 1});	% number of state variables

Dmat = vpa(zeros(n, dimD1));
for i = 1:dimD1
	Dmat(:, i) = D1{1, i};
end

count     = 0;
exc_tmp = zeros(dimD2);
for i=1:dimD2	% もう少し良いアルゴリズムを...
	f = D2{1, i};
	if rank([Dmat, f]) > dimD1+count
		count = count + 1;
		Dout{1, dimD1+count} = f;
		Dmat(:, dimD1+count) = f;
%		D1out = [D1out, {f}];	% この辺がキモい
%		D1mat = [D1mat, f];
	else
%		fprintf('i=%d\n', i);
		exc_tmp(i-count) = i;
	end
end

if count ~= dimD2
	exclusion = exc_tmp(1:dimD2-count);
else
	exclusion = nan;
end
spanD = Dout;
end

