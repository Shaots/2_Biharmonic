% EXPLICITLY VERIFIED by running
%
% HCT_LAMin2out([0,.075,0],[0,0.25,1],LAMin,1)
% LAMin= [1,0,0] -> LAMout= [1/3,1/3,1/3]
% LAMin= [0,1,0] -> LAMout= [0,1,0]
% LAMin= [0,0,1] -> LAMout= [0,0,1]
%
% HCT_LAMin2out([0,.075,0],[0,0.25,1],LAMin,2)
% LAMin= [1,0,0] -> LAMout= [1,0,0]
% LAMin= [0,1,0] -> LAMout= [1/3,1/3,1/3]
% LAMin= [0,0,1] -> LAMout= [0,0,1]
%
% and also by 'if' from below
function [LAMout,Xout,Yout]= HCT_LAMin2out(X,Y,LAMin,i)

order= [1,2,3,1,2];
ind= order(i:i+2);
x= [X(ind(1:3))]; 
y= [Y(ind(1:3))];
LAMBDA=[LAMin(ind(1:3))];

Xin= [mean(x),x(2),x(3)];
Yin= [mean(y),y(2),y(3)];

Xout= LAMBDA*Xin';
Yout= LAMBDA*Yin';

DET= x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2));
LAMout(ind(1))= ( (Xout-x(2))*(y(2)-y(3)) + (x(3)-x(2))*(Yout-y(2)) )/DET; 
LAMout(ind(2))= ( (Xout-x(3))*(y(3)-y(1)) + (x(1)-x(3))*(Yout-y(3)) )/DET;
LAMout(ind(3))= ( (Xout-x(1))*(y(1)-y(2)) + (x(2)-x(1))*(Yout-y(1)) )/DET;

if ( or( abs(Xout - LAMout*X') > 1e-15 , abs(Yout - LAMout*Y') > 1e-15 ) ) 
    disp('Warning in recalculation of lambda-vector from local to global')
end    