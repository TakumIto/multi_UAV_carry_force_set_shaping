function out = invClosure(D, x)
% involutive closure を計算します
% out = invClosure(D, x);
% x: 座標ベクトル
% D: cell 配列

flag = isa(D, 'cell');

% もう少しアルゴリズムを改善したい
if flag
	Dout=D;
	c = size(Dout, 2);	% 要素の個数
	if c >= length(x)
		out = Dout;
		return;
	end
	for i = 1:c		% cell2mat が使えなかったので
		if i == 1
			if isa(Dout{1, i}, 'double')
				Dmat = vpa(Dout{1, i});
			else
				Dmat = Dout{1, i};
			end
		else
			Dmat(:, i) = Dout{1, i};
			%		Dmat = [Dmat, Dout{1, i}];
		end
	end
else
	Dmat=D;
	c = length(Dmat(1, :));	% 要素の個数
	if c >= length(x)
		out = Dmat;
		return;
	end
	for i = 1:c
		Dout{1, i} = Dmat(:, i);
	end
end
c_old = 0;
while c_old ~= c
	c_old = c;
	for i = 1:c-1
		for j = i+1:c
			X = Dout{1, i};
			Y = Dout{1, j};
			liebracket = LieBracket(x, X, Y);
			if rank([Dmat, liebracket]) > c		% not involutive
				c = c+1;
				fprintf('\t要素%dと要素%dの Lie Bracket を要素%dに加えます\n',i,j,c);
				Dout{1, c} = liebracket;
				Dmat(:, c) = liebracket;
%				Dout = [Dout, {liebracket}];
%				Dmat = [Dmat, liebracket];
				break;
			end
		end
		if c_old ~= c
			break;
		end
	end
end
if flag
	out = Dout;
else
	out = Dmat;
end
end
