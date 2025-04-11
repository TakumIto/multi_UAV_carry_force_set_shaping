function out = LieBracket(x, f, g)
%%LIEBRACKET - Lie bracket の計算
%	Lie Bracket [f, g] を計算します (= L_g f - L_f g)
%	
%	out = LieBracket(x, f, g)

out = LieD(g, x, f) - LieD(f, x, g);
end


