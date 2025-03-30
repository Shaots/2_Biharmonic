% Warning:
%           X, Y, LAM are given in global order
%
% function is able to also return values of P, Px and Py
function [DelP]= HCT_PVALS(X,Y, LCOEFi, i, LAM)

order= [1,2,3,1,2];
x= [X(order(i:i+2))];
y= [Y(order(i:i+2))];
L= [LAM(order(i:i+2))];
L2= L.^2;
L3= L.^3;

% LV vector of local basis functions in barycentric coordinates (10 functions)

% LV= [ L3(1),L3(2),L3(3),...
%       L2(1)*L(3),L2(1)*L(2), L2(2)*L(1),L2(2)*L(3), L2(3)*L(2), L2(3)*L(1),...
%       L(1)*L(2)*L(3) ];

% DLV1= [ 3*L2(1),0,0, 2*L(1)*L(3),2*L(1)*L(2), L2(2),0,0, L2(3), L(2)*L(3) ];
% DLV2= [ 0,3*L2(2),0, 0, L2(1), 2*L(2)*L(1),2*L(2)*L(3), L2(3), 0, L(1)*L(3) ];
% DLV3= [ 0,0,3*L2(3), L2(1),0, 0,L2(2), 2*L(3)*L(2), 2*L(3)*L(1), L(1)*L(2) ];

% DX= (y(2)-y(3))*DLV1 + (y(3)-y(1))*DLV2 + (y(1)-y(2))*DLV3;
% DY= (x(3)-x(2))*DLV1 + (x(1)-x(3))*DLV2 + (x(2)-x(1))*DLV3;

DLV11= [ 6*L(1),0,0, 2*L(3),2*L(2), 0,0,0,0,0];
DLV22= [ 0,6*L(2),0, 0,0, 2*L(1),2*L(3), 0,0,0];
DLV33= [ 0,0,6*L(3), 0,0,0,0, 2*L(2),2*L(1), 0];

DLV12= [ 0,0,0,0, 2*L(1),2*L(2), 0,0,0, L(3)];
DLV23= [ 0,0,0,0, 0,0,2*L(2),2*L(3), 0, L(1)];
DLV31= [ 0,0,0, 2*L(1),0,0,0,0, 2*L(3), L(2)];

D1= (y(2)-y(3))^2*DLV11 + (y(3)-y(1))^2*DLV22 + (y(1)-y(2))^2*DLV33;
D3= (x(3)-x(2))^2*DLV11 + (x(1)-x(3))^2*DLV22 + (x(2)-x(1))^2*DLV33;

D2= 2*(y(2)-y(3))*(y(3)-y(1))*DLV12 + ...
    2*(y(3)-y(1))*(y(1)-y(2))*DLV23 + ...
    2*(y(1)-y(2))*(y(2)-y(3))*DLV31;

D4= 2*(x(3)-x(2))*(x(1)-x(3))*DLV12 + ...
    2*(x(1)-x(3))*(x(2)-x(1))*DLV23 + ...
    2*(x(2)-x(1))*(x(3)-x(2))*DLV31;

DET= x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2));
DelP= (LCOEFi*(D1+D2+D3+D4)')/(DET^2);

% PVAL(1)= (LCOEFi*LV');
% PVAL(2)= (LCOEFi*DX')/DET;
% PVAL(3)= (LCOEFi*DY')/DET;