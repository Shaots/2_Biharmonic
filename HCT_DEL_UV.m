% THE MOST EXPENSIVE SUBROUTINE
%
function [DelU,DelV,intDelVPhi]= HCT_DEL_UV(p,t,LaplU,A1,A2,A3,V,GVx,GVy,check);

% we do the same trick as global-to-local reindexation
% indexation [i,i+1,i+2] can be used for i=1,2,3         
LAMBin= [0.5, 0.5, 0, 0.5, 0.5]; 

nt= size(t,2);
rand('state',sum(100*clock));
Tver= int32(1+rand(1).*(nt-1));

DelV= zeros(3,3,nt);
intDelVPhi= zeros(3,nt);
LCOEFV= zeros(3,10);

for T=1:nt
    NODE= t(1:3,T);
    X=  p(1,NODE(1:3));
    Y=  p(2,NODE(1:3)); 
    S= 0.5 * abs( (X(1)-X(3))*(Y(2)-Y(3)) + (X(2)-X(3))*(Y(3)-Y(1)) );
    Si= S/3;
    VT= V(NODE(1:3));
    VxT= GVx(NODE(1:3));
    VyT= GVy(NODE(1:3));       
    
% i is index of HCT subtriangle area of every subtriangle is S/3
% A1,A2,A3 are global and are defined for concrete example
% ORDER OF INDEXATION OF DelV(i,j,T) IS IMPORTANT IN HCT_INTschP2 - DO NOT CHANGE IT

    for i=1:3
        
% the next step is CHEEP
        [LCOEFV(i,1:10)]=HCT_LCOEFi(X,Y,VT,VxT,VyT,i);
        for j=1:3
            LAMin= [LAMBin(j:j+2)];
            [LAMout,Xout,Yout]= HCT_LAMin2out(X,Y,LAMin,i);
            
% the next step is most EXPENSIVE: verified by substituting 
% DelU(i,j,T)= DelV(i,j,T) + 0.01, instead.
% Reason - large amount of calculations for inline function
            if ( j <= 4-i )
                DelU(i,j,T)= LaplU(Xout,Yout,A1,A2,A3);
            else    
% special trick to reduce costs: 6 values instead of 9.                 
                DelU(i,j,T)= DelU(4-j,4-i,T);
            end    

            [DelV(i,j,T)]= HCT_PVALS(X,Y, LCOEFV(i,1:10), i, LAMout); 
            intDelVPhi(1:3,T)= intDelVPhi(1:3,T) + DelV(i,j,T)* Si/3 * LAMout(1:3)';
        end  
    end  
 
%=============================================================================        
% comparatively, 'if' costs nothing 
%=============================================================================
    if ( isequal(T,Tver) & check==1 )
        clc;
        fprintf('Randomly selected element for verification\n'); 
        fprintf('\n\nVerify continuity (the center of mass)\n');        
        LAM= [1/3,1/3,1/3]
        PVAL= [0,0,0];
        for i=1:3
            [PVAL,DelP]= HCT_PVALi(X,Y, LCOEFV(i,1:10), i, LAM);
            disp(sprintf('\nsubtriangle %d: P= %g, Px= %g, Py= %g',i,PVAL)); 
        end        

        fprintf('\n\nVerify that laplacian P is an element of P^1\n');
        for i=1:3
            fprintf('\nsubtriangle %d:\n',i);        
            [LAMout1,Xout,Yout]= HCT_LAMin2out(X,Y,[1,0,0],i);
            [LAMout2,Xout,Yout]= HCT_LAMin2out(X,Y,[0,1,0],i);
            [LAMout3,Xout,Yout]= HCT_LAMin2out(X,Y,[0,0,1],i);
            [LAMoutC,Xout,Yout]= HCT_LAMin2out(X,Y,[1/3,1/3,1/3],i);        
            [PVAL,DelP1]= HCT_PVALi(X,Y, LCOEFV(i,1:10), i, LAMout1);
            [PVAL,DelP2]= HCT_PVALi(X,Y, LCOEFV(i,1:10), i, LAMout2);
            [PVAL,DelP3]= HCT_PVALi(X,Y, LCOEFV(i,1:10), i, LAMout3);    
            [PVAL,DelPC]= HCT_PVALi(X,Y, LCOEFV(i,1:10), i, LAMoutC);        
            MeanDelP= (DelP1 + DelP2 + DelP3)/3;
            fprintf('\n%g - mean of values in the vertices\n',MeanDelP);             
            fprintf('\n%g - value in the center\n',DelPC);                         
        end 

        pause
        clc;    
        order= [1,2,3,1,2];
        PVAL1= [0,0,0];
        PVAL2= [0,0,0];
        fprintf('\n\nVerify continuity:\n');        
        for i=1:3            
            fprintf('\n\nLocal vertex i= %d\n',i);
            LAM([order(i),order(i+1),order(i+2)])= [1,0,0]
            [PVAL1,DelP1]= HCT_PVALi(X,Y, LCOEFV(order(i+1),1:10), order(i+1), LAM);
            [PVAL2,DelP2]= HCT_PVALi(X,Y, LCOEFV(order(i+2),1:10), order(i+2), LAM);
            disp(sprintf('\nP(i+1)= %g, P(i+1)x= %g, P(i+1)y= %g',PVAL1));             
            disp(sprintf('\nP(i+2)= %g, P(i+2)x= %g, P(i+2)y= %g',PVAL2));            
            fprintf('\n     V= %g,      Vx= %g,      Vy= %g\n',VT(i),VxT(i),VyT(i));
           
            fprintf('\n\nPoint on common edge of %d and %d subtriangles\n',order(i+1),order(i+2));
            LAM([order(i),order(i+1),order(i+2)])= [0.83,0.075,0.075]
            [PVAL1,DelP1]= HCT_PVALi(X,Y, LCOEFV(order(i+1),1:10), order(i+1), LAM);
            [PVAL2,DelP2]= HCT_PVALi(X,Y, LCOEFV(order(i+2),1:10), order(i+2), LAM);           
            disp(sprintf('\nP(i+1)= %g, P(i+1)x= %g, P(i+1)y= %g',PVAL1));             
            disp(sprintf('\nP(i+2)= %g, P(i+2)x= %g, P(i+2)y= %g',PVAL2));
            pause
            clc;
        end
        fprintf('verification is done: continue\n');
    end       
%=============================================================================
end