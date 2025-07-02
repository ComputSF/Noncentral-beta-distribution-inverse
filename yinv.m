function [y,ier]=yinv(z, x, p, q)
%--------------------------------------------------
%Inversion of B_(x,y)(p,q)=z.
%Computes y when z, x, p, q are given.
%Newton's method is used for the main iterations
%--------------------------------------------------
% Outputs:
%     ier: error flag.
%          ier=0,computation succesful
%          ier=1,the maximum number of iterations 
%                has been reached. Possible loss of
%                accuracy in the y-value.
%     y:  approximation found. 
%---------------------------------------------------
eps=1.0e-14;
kmax=10;
ier=0;
r= p+q;
if p>9 && q>9 
  zeta= sqrt(2/r)*erfcinv(2*z);
  y= yzeta(zeta, x, p, q, eps);
else
  y0=(x+2*p)/(x+2*p+2*q);
  y=y0;
  [ib,ierr]=Bpqxy(x,y,p,q);
  if ib<z && ierr==0
    y=illinois(1,ib,y0,z,x,p,q);
  else
    y=illinois(0,ib,y0,z,x,p,q);
  end
end  
d= 1; k= 1;
while abs(d) > eps &&k < kmax
  ib=Bpqxy(x,y,p,q);
  dbdy= Bpqxydery(x,y,p,q);
  if dbdy==0 || isinf(dbdy)>0
    k=kmax;
  else  
    y1= y-(ib-z)/dbdy; 
    d= y/y1-1;
    k= k+1;
    if y1>1 || y1<0
      k=kmax;
    else 
      y= y1;
    end  
  end  
end
if k==kmax
  ier=1;
end  
end
function y=illinois(icho,ib,y,z,x,p,q)
if icho==0
  a=0;
  b=y;
  fa=-z;
  fb=ib-z;
else  
  a=y;
  b=1;
  fa=ib-z;
  fb=1-z;
end  
eps=1.e-5;
erro=1;
n=1;
while erro>eps && n<5
  n=n+1;
  c=b-fb*(b-a)/(fb-fa);
  erro=abs(1-b/c);
  [B3,ier3]=Bpqxy(x,c,p,q);
  if ier3==0
    fc=B3-z;
    if fb*fc<0
      a=b;fa=fb;
    else
      fa=fa/2;  
    end
    b=c;fb=fc;
  end
end
y=c;
end      
function yz=yzeta(zeta, x, p, q, eps)
r= p+q; s2= q/r; c2= p/r; w= x+2*r; 
r2= r*r; r3= r*r2;
w2= w*w; w3= w*w2; w4= w*w3; 
w5= w*w4; w6= w*w5; w7= w*w6; w8=  w*w7;
w9= w*w8; w10=w*w9; w12= w2*w10; w14= w7*w7; w21= w14*w7;
s22= s2*s2; s23=s22*s2; s24= s2*s23; s25= s2*s24; s26= s25*s2;
y0= (x+2*p)/w; 
yk(1)= -2*sqrt(-(4*s2-w2)*s2)/w2; 
if abs(zeta)<0.01
  y12= yk(1)*yk(1); y15= y12*y12*yk(1); y18= y12*y15*yk(1);
  yk(2)= -(8/3)*(48*s22*r-12*w2*s2*r-8*w*s22+w4*r)*s2/(y12*w7);
  yk(3)= (16/9)*(w8*r2+9216*s24*r2-5760*s23*w2*r2+64*s24*w2+...
		 1920*s23*w3*r+960*s22*w4*r2-144*s23*w4-160*s22*w5*r...
         -60*s2*w6*r2-3072*s24*w*r)*s22/(y15*w14);
  yk(4)= -(256/135)*(1105920*s26*r3-w12*r3-218880*s24*w5*r2-...
		 46080*s25*w4*r+414720*s24*w4*r3-552960*s26*w*r2...
         +552960*s25*w3*r2+512*s26*w3-3456*s25*w5...
         -72*s2*w10*r3+46080*s26*w2*r+24480*s23*w7*r2...
         -840*s22*w9*r2-60480*s23*w6*r3+32640*s24*w6*r-...
		 1105920*s25*w2*r3-1728*s24*w7+3744*s22*w8*r3...
         -2160*s23*w8*r)*s23/(y18*w21);
  y= y0; zetak= 1; k= 0; d= 1;
  while abs(d) > eps && k < 4
	k= k+1;  zetak= zeta*zetak; 
    term= yk(k)*zetak; 
	  y= y+term; 
	  d= term/y;
  end
elseif zeta < -0.5 
  xi= x/(2*r);
  w= sqrt((xi-c2)*(xi-c2)+4*xi);
  if xi > c2 
    t0= (xi-c2+w)/(2*xi);
  else
    t0= 2/(c2-xi+w);
  end  
  y= 1-exp(-(zeta*zeta/2+ft(t0,s2,xi))/s2);
elseif zeta > 0.5
  y= exp(-5*zeta*zeta);
else
  w= y0*(y0-1)/yk(1);     
  if w+zeta>0
	  y= 2*w*y0/(w+zeta+sqrt((w+zeta)*(w+zeta)-4*zeta*y0*w));
  else
	  y= (w+zeta-sqrt((w+zeta)*(w+zeta)-4*zeta*y0*w))/(2*zeta);
  end
end  
y1= y; 
d= 1; k= 0;
while abs(d) > eps && k < 10
  k= k+1; [zy,zetak]= zetaxypq(x, y1, p, q);
  y= y1-(zetak-zeta)/zy;
  d= y/y1-1; 
  y1= y;
end
yz= y;
end
function [dzetady, zetaxy]=zetaxypq(x,y,p,q)
r= p+q; xi= x*y/(2*r); c2= p/r; s2= q/r;
w= sqrt((xi-c2)*(xi-c2)+4*xi);
tp= 1/y; 
if xi > c2 
  t0= (xi-c2+w)/(2*xi);
else
  t0= 2/(c2-xi+w);
end  
if abs(tp-t0) < 1.0e-8
  zeta1= sqrt(w*(c2*(w+c2-xi)+2*xi)/(2*s2));
  zeta= zeta1*(tp-t0);
  y1= -2*sqrt(r*q*(4*r*(p+x)+x*x))/((2*r+x)*(2*r+x));
  dzetady= (x*y*(1-y)+(2*r+x)*y1*zeta1)/(2*y*r*(1-y)*zeta1);
else
  ftp= ft(tp, s2, xi); 
  ft0= ft(t0, s2, xi);
  zeta= sqrt(2*(ftp-ft0));
  if tp < t0 
    zeta= -zeta;
  end  
  y0= (x+2*p)/(x+2*r);
  dzetady= (x*y*(1-y)*(tp-t0)+(2*r+x)*(y-y0))/(2*y*r*(1-y)*zeta);
end
zetaxy=zeta;
end
function y=ft(t, s2, xi)
y= log(t)-s2*log(t-1)+xi*t; 
end