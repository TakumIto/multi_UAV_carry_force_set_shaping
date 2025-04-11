addpath('lib/');

% Ag=[   0 9.8 0;
%     -9.8   0 0;
%        0   0 0];
% AL=[zeros(3) eye(3) zeros(3);
%     zeros(3) zeros(3) zeros(3);
%     zeros(3) zeros(3) zeros(3)];
AL=[zeros(3) eye(3) zeros(3) zeros(3);
    zeros(3) zeros(3) zeros(3) zeros(3);
    zeros(3) zeros(3) zeros(3) eye(3);
    zeros(3) zeros(3) zeros(3) zeros(3)];

BL=[zeros(3) zeros(3);
    eye(3) zeros(3);
    zeros(3) zeros(3);
    zeros(3) eye(3)];

CL=[eye(3) zeros(3) zeros(3) zeros(3);
    zeros(3) zeros(3) eye(3) zeros(3)];


sys=ss(AL,BL,CL,zeros(6));
GL = zeros(6,9);
[KL,S,e]=lqi(sys,Q,R);
KL
rankABC = rank([AL,BL;CL,zeros(6)])

