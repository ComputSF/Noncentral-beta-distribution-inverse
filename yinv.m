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
iNewt=1;
r= p+q; s2= q/r; c2= p/r; w= x+2*r; 
y0= (x+2*p)/w; 
rho=x/(2*r);
rho2=rho*rho;
yk(1)=-sqrt(s2)*sqrt(rho2+2*rho+c2)/(rho+1)^2;
if abs(zeta)<0.01
  iNewt=0;
  rho3=rho2*rho;
  rho4=rho3*rho;
  yk(2)=(1/3)*(-rho4+c2^2*rho-3*c2*rho2-4*rho3-2*c2^2-...
      8*c2*rho-3*rho2+c2+3*rho)/((rho2+c2+2*rho)*(rho+1)^3);
  y= y0; zetak= 1; k= 0; d= 1;
  while abs(d) > eps && k < 2
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
if iNewt==1
  y1= y;
  d= 1; k= 0;
  while abs(d) > eps && k < 10
    k= k+1; [zy,zetak]= zetaxypq(x, y1, p, q);
    y= y1-(zetak-zeta)/zy;
    d= y/y1-1; 
    y1= y;
  end
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
  dzetady= -s2/((1-y)*(t0-1)*zeta1);
else
  ftp= ft(tp, s2, xi); 
  ft0= ft(t0, s2, xi);
  zeta= sqrt(2*(ftp-ft0));
  if tp < t0 
    zeta= -zeta;
  end  
  y0= (x+2*p)/(x+2*r);
  dzetady= s2*(t0-tp)/((1-y)*(t0-1)*zeta);
end
zetaxy=zeta;
end
function y=ft(t, s2, xi)
y= log(t)-s2*log(t-1)+xi*t; 
end