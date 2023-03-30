clear all

NGauss=10;
Ndim = 8;

[x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,zeros(Ndim,1),ones(Ndim,1));

tic
myint= 0;
for i0=1:NGauss
for i1=1:NGauss
  for i2=1:NGauss
    for i3=1:NGauss
      for i4=1:NGauss
        for i5=1:NGauss
        for i6=1:NGauss
        for i7=1:NGauss
%        for i8=1:NGauss
%        for i9=1:NGauss
          %myint= myint + wn(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9);
          myint= myint + wn(i0,i1,i2,i3,i4,i5,i6,i7);
        end
      end
    end
  end
end
end
end
end
toc
%end
%end
myint 


syms x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
tic
g7 = w' * (x *x0*x1*x2*x3*x4*x5*x6)
g6 = w' *subs(g7 , {x6}, x )
g5 = w' *subs(g6 , {x5}, x )
g4 = w' *subs(g5 , {x4}, x )
g3 = w' *subs(g4 , {x3}, x )
g2 = w' *subs(g3 , {x2}, x )
g1 = w' *subs(g2 , {x1}, x )
g0 = w' *subs(g1 , {x0}, x )
toc
