% U is exact solution with laplacian DelU
% Realization of non-conforming majorant
function [Et,Nt,Dt,Rt1,Rt2]= HCT_INTschP2(p,t,f,DelU,DelV,KS,KSx,KSy,YS1,YS2);
         
[dYS1x,tmp]= pdegrad(p,t,YS1);
[tmp,dYS2y]= pdegrad(p,t,YS2);

               % we do the same trick as global-to-local reindexation
               % indexation [i,i+1,i+2] can be used for i=1,2,3         
LAMBin= [0.5, 0.5, 0, 0.5, 0.5]; 

nt= size(t,2);
Et= zeros(1,nt);
Nt= zeros(1,nt);
Dt= zeros(1,nt);
Rt1= zeros(1,nt);
Rt2= zeros(1,nt);

for T=1:nt
    NODE= t(1:3,T);
    X=  p(1,NODE(1:3));
    Y=  p(2,NODE(1:3));
    S= 0.5 * abs( (X(1)-X(3))*(Y(2)-Y(3)) + (X(2)-X(3))*(Y(3)-Y(1)) );
    weight= S/3;
        
    KST= KS(NODE(1:3));
    YS1T= YS1(NODE(1:3));
    YS2T= YS2(NODE(1:3));    
                              % i is index of HCT subtriangles
    for i=1:3
        for j=1:3
               LAMin= [LAMBin(j:j+2)];
               [LAMout,Xout,Yout]= HCT_LAMin2out(X,Y,LAMin,i);
               KSVAL= LAMout*KST;               
                                               % SEE HCT_DELV 
               Dt(T)= Dt(T) + weight* (DelV(i,j,T) - KSVAL)^2;               
               Et(T)= Et(T) + weight* (DelV(i,j,T) - DelU(i,j,T))^2;               
               Nt(T)= Nt(T) + weight* DelU(i,j,T)^2;
        end
    end  
    
    YST_12= [ YS1T(1)+YS1T(2), YS2T(1)+YS2T(2) ]/2;
    YST_23= [ YS1T(2)+YS1T(3), YS2T(2)+YS2T(3) ]/2;
    YST_31= [ YS1T(3)+YS1T(1), YS2T(3)+YS2T(1) ]/2;        
    Rt12= (KSx(T)-YST_12(1))^2 + (KSy(T)-YST_12(2))^2 ;    
    Rt23= (KSx(T)-YST_23(1))^2 + (KSy(T)-YST_23(2))^2 ;        
    Rt31= (KSx(T)-YST_31(1))^2 + (KSy(T)-YST_31(2))^2 ;            
    Rt1(T)= weight* ( Rt12 + Rt23 + Rt31 );
    Rt2(T)= S*(dYS1x(T) + dYS2y(T)-f(T))^2;
end    