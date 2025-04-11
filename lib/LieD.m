function out = LieD(phi, x, f)
% Lie微分するだけです
% L_f phi を計算します
% xとfは, 次元が等しい縦ベクトルとしてください
% LieD(phi, x, f)

[m, n] = size(phi);
if m == 1 && n == 1
	out = transpose(gradient(phi, x))*f;
else
	out = jacobian(phi, x)*f;
end
end
