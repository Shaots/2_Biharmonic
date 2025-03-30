 function prnf2D( p,e,t, f2D, ftitle, fcolormap, fcolorbar )
 
 if ( length(f2D)==1 )                                                   
     f2D= f2D*ones(1,size(t,2));
 end    

 pdeplot(p,e,t,'xydata',f2D,'xystyle','flat','mesh','on','colorbar',fcolorbar);
 title(ftitle,'FontSize',21,'FontWeight','bold');
 colormap(fcolormap);
 axis equal;
 brighten(0.3);