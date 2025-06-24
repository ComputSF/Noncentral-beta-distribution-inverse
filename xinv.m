function [x,ier]=xinv(z, y, p, q)
%-----------------------------------------------------------
%Inversion of B_(x,y)(p,q)=z.
%Computes x when z, y, p, q are given.
%Halley's method is used for the main iterations
%--------------------------------------------------------------------
% Outputs:
%     ier: error flag.
%          ier=0,computation succesful
%          ier=1,the maximum number of iterations 
%                has been reached. Possible loss of
%                accuracy in the x-value.
%          ier=2,the z-value is not acceptable for the given
%                y, p, q values: z is larger than I_y(p,q).
%                A negative value for x is returned.
%          ier=3,the initial approximation to x is greater than
%                xlim=1000 (upper limit for the computation of Bpqxy).
%                The initial approximation to x is returned.        
%     x:  approximation found. 
%----------------------------------------------------------------------
%First check: z<I_y(p,q)
g= betaincreg(y, p, q);
if z>g
  ier=2;
  x=-100;
  %It is not possible to find a solution
else  
  eps=1.0e-13;
  ier=0;
  kmax=15;
  r= p+q;
  if p>9 && q>9 
    zeta= sqrt(2/r)*erfcinv(2*z);
    x= xzeta(zeta, y, p, q);
    if x>1000
      ier=3;
    end
    if r<15
      xi=x;
      for j=1:5
        x=gFP(xi,y,p,q);
      end  
    end  
  else
    x0=2*((p+q)*y-p)/(1-y);
    if x0<0
      x=0;
    else
    %fixed point iterations
      xi=x0;
      for j=1:5
        x=gFP(xi,y,p,q);
      end  
    end
  end
  if ier==0
    d= 1; k= 1;
    while abs(d)>eps && k<kmax
      [ib,ier1]=Bpqxy(x,y,p,q);
      [dbdx,ier2]=Bpqxyderx(x,y,p,q);
      if dbdx==0 || ier1>0 || ier2>0
        k=kmax;  
      else    
        f=ib-z;
        Bg=0.5*Bghalf(x, y, p, q);
        coci=f/(dbdx+0.5*Bg*f);
        x1= x-coci;
        if x1<0
          k=kmax;
        else
          d= x/x1-1;
          k= k+1; 
          x=x1;
        end  
      end  
    end  
    if k==kmax
      ier=1; 
      if x<0
        x=0;
      end   
    end
  end
end  
end
function xzet=xzeta(zeta, y, p, q)
eps=1.0e-14;
r= p+q; s2= q/r;
y1= 1.0-y;  
x0= 2*(r*y-p)/y1;
[zx,zeta0]=zetaxypqX(0, y, p, q);
if zeta < zeta0
  x= -100;
  ik=1;
else  
  if abs(zeta) < 0.1
    s4= s2*s2; s6= s4*s2; s8= s4*s4;
	yk(1)= y1;
	for k= 2:12 
	  yk(k)= y1*yk(k-1);
  end  
	w= sqrt(s2-yk(2)); 
	w2= w*w; w4= w2*w2; w5= w*w4; w8= w4*w4; w11= w*w2*w8;
  xk(1)= x0;
	xk(2)= 2*r*w/y1;
	xk(3)= (2.0/3.0)*(s2-yk(3))*r/(w2*y1);
  xk(4)= (1/18)*r*(s4+16*s2*yk(3)+yk(6)-9*s2*yk(2)-9*s2*yk(4))/(y1*w5);
  xk(5)= (1.0/135.0)*(-s6+yk(9)+27*s4*yk(2)-105*s4*yk(3)+135*s4*yk(4)...
		   -54*s4*yk(5)-135*s2*yk(5)+105*s2*yk(6)-27*s2*yk(7)...
           +54*s2*yk(4))*r/(y1*w8);
  xk(6)= -(1/2160)*r*(s8-162*s6*yk(2)-1359*s4*yk(4)...
			-720*s2*yk(6)+yk(12)-720*s6*yk(6)-1359*s4*yk(8)-162*s2*yk(10)...
			-2610*s6*yk(4)-8220*s4*yk(6)-2610*s2*yk(8)+2304*s6*yk(5)...
            +1184*s2*yk(9)+5472*s4*yk(7)...
			+1184*s6*yk(3)+5472*s4*yk(5)+2304*s2*yk(7))/(y1*w11);
	x= xk(1); zetak= 1.0; k= 0; d= 1;
  while abs(d) > eps && k < 5
    k= k+1;  zetak= zeta*zetak; term= xk(k+1)*zetak; 
	  x= x+term; d= term/x;
  end 	
else 
	a= x0;
  c= r/y1; 
	if abs(zeta0)<1.0e-5 
	  b= (4*p*r+x0*(3*p+q))/(2*sqrt(p*q)); 
    else   
	  b=-a/zeta0-c*zeta0;
    end  
	  x= a+zeta*(b+c*zeta);
  end  
  x1= x; 
  d= 1; k=0;
  ik=0; 
  while abs(d) > eps && k < 10 && ik==0
    [zx,zeta0]=zetaxypqX(x1, y, p, q); 
    x= x1-(zeta0-zeta)/zx; 
    if x < 0 
      x= 0;
      ik=1;
    else  
	    d= x/x1-1; 
	    x1= x; k= k + 1;
    end
  end
end 
zeta1=h0(zeta,x,y,p,q);
zet=zeta1/r;
zetan=zeta+zet;
xzetN=xzeta2(zetan, zeta0, y, p, q);
if ik==1
  %Negative x occurred (x=0) or zeta too small (x=-100)
  xzet=0;
else  
  xzet=real(xzetN);
end 
end
function xzet=xzeta2(zeta, zeta0, y, p, q)
eps=1.0e-13;
r= p+q; s2= q/r;
y1= 1.0-y;  
x0= 2*(r*y-p)/y1;
if abs(zeta) < 0.1
  s4= s2*s2; s6= s4*s2; s8= s4*s4;
  yk(1)= y1;
  for k= 2:12 
    yk(k)= y1*yk(k-1);
  end  
  w= sqrt(s2-yk(2)); 
  w2= w*w; w4= w2*w2; w5= w*w4; w8= w4*w4; w11= w*w2*w8;
  xk(1)= x0;
  xk(2)= 2*r*w/y1;
  xk(3)= (2.0/3.0)*(s2-yk(3))*r/(w2*y1);
  xk(4)= (1/18)*r*(s4+16*s2*yk(3)+yk(6)-9*s2*yk(2)-9*s2*yk(4))/(y1*w5);
  xk(5)= (1.0/135.0)*(-s6+yk(9)+27*s4*yk(2)-105*s4*yk(3)+135*s4*yk(4)...
        -54*s4*yk(5)-135*s2*yk(5)+105*s2*yk(6)-27*s2*yk(7)...
        +54*s2*yk(4))*r/(y1*w8);
  xk(6)= -(1/2160)*r*(s8-162*s6*yk(2)-1359*s4*yk(4)...
        -720*s2*yk(6)+yk(12)-720*s6*yk(6)-1359*s4*yk(8)-162*s2*yk(10)...
        -2610*s6*yk(4)-8220*s4*yk(6)-2610*s2*yk(8)+2304*s6*yk(5)...
        +1184*s2*yk(9)+5472*s4*yk(7)...
        +1184*s6*yk(3)+5472*s4*yk(5)+2304*s2*yk(7))/(y1*w11);
  x= xk(1); zetak= 1.0; k= 0; d= 1;
  while abs(d) > eps && k < 5
    k= k+1;  zetak= zeta*zetak; term= xk(k+1)*zetak; 
    x= x+term; d= term/x;
  end   
else 
  a= x0; c= r/y1; 
  if abs(zeta0)<1.0e-5 
    b= (4*p*r+x0*(3*p+q))/(2*sqrt(p*q)); 
  else   
    b= -a/zeta0-c*zeta0;
  end  
  x= a+zeta*(b+c*zeta);
end  
x1= x; 
d= 1; k=0;
ik=0; 
while abs(d) > eps && k < 10 && ik==0
  [zx,zeta0]=zetaxypqX(x1, y, p, q); 
  x= x1-(zeta0-zeta)/zx;
  if x < 0 
    x= 0;
    ik=1;
  else  
    d= x/x1-1; 
    x1= x; k= k + 1;
  end
end 
if ik==1
  %Negative x occurred (x=0) or zeta too small (x=-100)
  xzet=0;
else  
  xzet=x;
end 
end
function [dzetadx, zetaxy]=zetaxypqX(x,y,p,q)
%Calculation of zeta and its derivative wrt x
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
  dzetadx= y/(2*r*zeta1);
else
  ftp= ft(tp, s2, xi); 
  ft0= ft(t0, s2, xi);
  zeta= sqrt(2*(ftp-ft0));
  if tp < t0 
    zeta= -zeta;
  end  
  dzetadx= y*(tp-t0)/(2*r*zeta);
end
zetaxy=zeta;
end
function y=ft(t, s2, xi)
y= log(t)-s2*log(t-1)+xi*t; 
end
function x=gFP(xi, y, p, q)
%---------------------------------------
%Function for the fixed point iteration
%---------------------------------------
a=p+q;
b=p+1;
amb=a-b;
z=0.5*xi*y;
if abs(amb-round(amb))<eps
  c2=Mabx(a+1,b+1,z)/Mabx(a,b,z);
else  
  c2=1/continued_fractionM(a,b,z);
end
Bgh=a*c2*y/b;
x=Bgh*xi;
end

function Bgh=Bghalf(x, y, p, q)
%---------------------------------------
%Function for the fixed point iteration
%---------------------------------------
a=p+q;
b=p+1;
amb=a-b;
z=0.5*x*y;
if abs(amb-round(amb))<eps
  c2=Mabx(a+1,b+1,z)/Mabx(a,b,z);
else  
  c2=1/continued_fractionM(a,b,z);
end
Bgh=1-a*c2*y/b;
end
function [h] = continued_fractionM(a,b,z)
%-----------------------------------------
%Continued fraction for the evaluation of 
% M(a,b,z)/M(a+1,b+1,z)
%-----------------------------------------
fpmin=1e-50;
c=1.0;
h=1;
d=0;
m=0;
del=0.0;
while (abs(del-1.0) > eps) && (m<100)
  m=m+1;
  if mod(m,2)==0
    n=m/2;
    aa=z*(a+n)/((b+2*n-1)*(b+2*n));
  else
    n=(m-1)/2;    
    aa=z*(a-b-n)/((b+2*n)*(b+2*n+1));
  end  
  d=1.0+aa*d;
  if (abs(d)<fpmin)
    d=fpmin;
  end
  c=1.0 +aa/c;
  if (abs(c)<fpmin)
    c=fpmin;
  end
  d=1.0/d;
  del=d*c;
  h=h*del;
end
end
function [Mkummer,ier]=Mabx(a,b,z)
%----------------------------------------------------------
% Computation of the confluent hypergeometric function
% M(a,b,x) for positive values of the parameters a, b and
% argument x
%----------------------------------------------------------
% Outputs:
% Mkummer: M(a,b,x) function value
% ier: error flag
%     ier=0,  computation successful.
%     ier=1,  some loss of accuracy is expected.
%     ier=2,  overflow error
%----------------------------------------------------------
overf=1e+308;
ier=0;
if a>100 && b>100 && z>20
  Mkummer= Mabxlps(a,b,z);
else  
  if a<1 && b<1 && z>30
    [Mkummer,ier]=MabxlargexN(a,b,z);
  else 
    [Mkummer,ier]=Mabxseries(a,b,z);
  end  
end
if ~isfinite(Mkummer)
  Mkummer=overf;
  ier=2;
end   
end
function [m1,ier]=MabxlargexN(a,b,x)
%---------------------------------------------------
% Use of the large x asymptotics to compute M(a,b,x)
% We will use this function for small a, b
%---------------------------------------------------
ier=0;
ba=b-a;
if ba>0.1
  [m1,ier]=Mabxlargex(a,b,x); 
else
  bb1=b+1;
  [m1a,iera]=Mabxlargex(a,b+1,x);
  [m1b,ierb]=Mabxlargex(a,b+2,x);
  m1=(1/(b*bb1))*(bb1*(b+x)*m1a-x*(bb1-a)*m1b);
  if iera>0 || ierb>0
    ier=1;
  end  
end
end
function [m1,ier]=Mabxlargex(a,b,x)
%---------------------------------------------------
% Use of the large x asymptotics to compute M(a,b,x)
%---------------------------------------------------
ier=0;
ba=b-a;
m1=1;
aa=1-a;
v=1;
r=1;
while v>=m1*eps && r<10000
  v=v/(x*r);
  fact=(aa+r-1)*(ba+r-1);
  v=v*fact;
  m1=m1+v;
  r=r+1;
end
m1=gamma(b)*exp(x)*x^(-ba)*m1/gamma(a);
if r==10000
  ier=1;
end
end
function [m1,ier]=Mabxseries(a,b,x)
%----------------------------------------
% Use of series to compute M(a,b,x)
%----------------------------------------
ier=0;
m1=1;
v=1;
r=1;
while v>=m1*eps && r<10000
  v=v*x/r;
  v=v*(a+r-1)/(b+r-1);
  m1=m1+v;
  r=r+1;
end
if r==10000
  ier=1;
end
end
function Mabz=Mabxlps(a,b,z)
%---------------------------------------------
%Use of the expansion for large a, b, z. 
%Computes M(a, b, z) with n=nmax coefficients.
%---------------------------------------------
nmax=5;
alpha=a/z; beta=b/z;
mu=beta-alpha;  w=beta+1;
tau= 2/(w+sqrt(w^2-4*mu)); t0= mu*tau;
fg0= 1/sqrt(beta*mu*tau^2-2*mu*tau+1);
u= beta*tau; v= (1-t0)/(alpha*tau);
Phi= exp(z*(1-t0))*gamstar(b)/gamstar(a)*sqrt(a/b)*fg0*u^b*v^a;
s=1; 
for n=1:nmax 
  fgntilde= pqmutau(mu, tau, n); 
  s= s + fgntilde/z^n;
end  
Mabz=Phi*s; 
end
function [fg]=pqmutau(mu, tau, n) 
% Coefficients for the expansion
fg=0;
mutaum1=(mu*tau-1);
mu2=mu*mu;
tau2=tau*tau;
tau3=tau2*tau;
tau4=tau2*tau2;
tau5=tau4*tau;
mutau2m1=(mu*tau2-1);
mutau2m12=mutau2m1*mutau2m1;
mutau2m13=mutau2m12*mutau2m1;
if n == 1    
  fg=-1/12*mu*tau2*(mu2*tau5-13*mu2*tau4-mu*tau4+21*mu*tau3+4*mu*tau2-...
      9*tau2-2*tau-1)/mutau2m13/mutaum1;
elseif n == 2
  mu3=mu2*mu;
  mu4=mu3*mu;
  mu5=mu4*mu;
  tau6=tau5*tau;
  tau7=tau6*tau;
  tau8=tau7*tau;
  tau9=tau8*tau;
  tau10=tau9*tau;
  mutaum12=mutaum1*mutaum1;
  mutau2m16=mutau2m13*mutau2m13;
  fg=1/288*mu*tau4*(mu5*tau10-26*mu5*tau9+313*mu5*tau8-2*mu4*tau9+...
   68*mu4*tau8-1690*mu4*tau7+mu3*tau8+1048*mu4*tau6-60*mu3*tau7+...
   3255*mu3*tau6-3958*mu3*tau5+18*mu2*tau6+186*mu3*tau4-2678*mu2*tau5+...
   5462*mu2*tau4-490*mu2*tau3+801*mu*tau4-8*mu2*tau2-3276*mu*tau3+...
   454*mu*tau2+4*mu*tau+720*tau2+mu-144*tau)/mutau2m16/mutaum12;
elseif n == 3 
  mu3=mu2*mu;
  mu4=mu3*mu;
  mu5=mu4*mu;
  mu6=mu5*mu;
  mu7=mu6*mu;
  mu8=mu7*mu;
  tau6=tau5*tau;
  tau7=tau6*tau;
  tau8=tau7*tau;
  tau9=tau8*tau;
  tau10=tau9*tau;
  tau11=tau10*tau;
  tau12=tau11*tau;
  tau13=tau12*tau;
  tau14=tau13*tau;
  tau15=tau14*tau;
  tau16=tau15*tau;
  tau17=tau16*tau;
  mutaum12=mutaum1*mutaum1;
  mutaum13=mutaum12*mutaum1; 
  mutau2m16=mutau2m13*mutau2m13;
  mutau2m19=mutau2m16*mutau2m13;
  fg= 1/51840*mu*tau4*(139*mu8*tau17+195*mu8*tau16-4695*mu8*tau15...
      -417*mu7*tau16+56201*mu8*tau14-2001*mu7*tau15+30105*mu7*tau14+...
      417*mu6*tau15-753675*mu7*tau13+4848*mu6*tau14+829668*mu7*tau12...
      -69141*mu6*tau13-139*mu5*tau14+3197085*mu6*tau12-4473*mu5*tau13...
      -6330048*mu6*tau11+73563*mu5*tau12+1797159*mu6*tau10...
      -6394821*mu5*tau11+1431*mu4*tau12+19013532*mu5*tau10...
      -36663*mu4*tau11-10070166*mu5*tau9+6694593*mu4*tau10+...
      700264*mu5*tau8-29304894*mu4*tau9+6831*mu3*tau10+...
      22957125*mu4*tau8-3562623*mu3*tau9-3181371*mu4*tau7...
      +24696258*mu3*tau8+18579*mu4*tau6-27110118*mu3*tau7+...
      763101*mu2*tau8+5732235*mu3*tau6-10855458*mu2*tau7-75291*mu3*tau5...
      +17316315*mu2*tau6+1668*mu3*tau4-5103788*mu2*tau5+1951776*mu*tau6...
      +114675*mu2*tau4-5572800*mu*tau5-5586*mu2*tau3+2244240*mu*tau4...
      -139*mu2*tau2-81648*mu*tau3+680400*tau4+6480*mu*tau2...
      -388800*tau3+432*mu*tau+21600*tau2-1728*tau-432)/mutau2m19/mutaum13;
elseif n==4
  mu3=mu2*mu;
  mu4=mu3*mu;
  mu5=mu4*mu;
  mu6=mu5*mu;
  mu7=mu6*mu;
  mu8=mu7*mu;
  mu9=mu8*mu;
  mu10=mu9*mu;
  mu11=mu10*mu;
  tau5=tau4*tau;
  tau6=tau5*tau;
  tau7=tau6*tau;
  tau8=tau7*tau;
  tau9=tau8*tau;
  tau10=tau9*tau;
  tau11=tau10*tau;
  tau12=tau11*tau;
  tau13=tau12*tau;
  tau14=tau13*tau;
  tau15=tau14*tau;
  tau16=tau15*tau;
  tau17=tau16*tau;
  tau18=tau17*tau;
  tau19=tau18*tau;
  tau20=tau19*tau;
  tau21=tau20*tau;
  tau22=tau21*tau;
  mutaum12=mutaum1*mutaum1;
  mutaum13=mutaum12*mutaum1; 
  mutaum14=mutaum13*mutaum1; 
  mutau2m16=mutau2m13*mutau2m13;
  mutau2m19=mutau2m16*mutau2m13;
  mutau2m112=mutau2m19*mutau2m13;
 fg= -1/2488320*mu*tau6*(571*mu11*tau22-7228*mu11*tau21...
     -9390*mu11*tau20-2284*mu10*tau21+224804*mu11*tau19+...
     28176*mu10*tau20-2697077*mu11*tau18+139336*mu10*tau19...
     +3426*mu9*tau20-3208064*mu10*tau18-40980*mu9*tau19+...
     80235396*mu10*tau17-507368*mu9*tau18-2284*mu8*tau19...
     -112029040*mu10*tau16+15330560*mu9*tau17+26164*mu8*tau18...
     -628961670*mu9*tau16+829908*mu8*tau17+571*mu7*tau18...
     +1589570756*mu9*tau15-36616440*mu8*tau16-5952*mu7*tau17...
     -751445924*mu9*tau14+2295515340*mu8*tau15-698056*mu7*tau16...
     -8478911244*mu8*tau14+49416240*mu7*tau15-180*mu6*tau16...
     +7259504812*mu8*tau13-4662971775*mu7*tau14+295520*mu6*tau15...
     -1338944656*mu8*tau12+23731481796*mu7*tau13-38458640*mu6*tau14...
     -29188066222*mu7*tau12+5619000492*mu6*tau13-49950*mu5*tau14...
     +10084935348*mu7*tau11-39265436532*mu6*tau12+16155720*mu5*tau13...
     -685938350*mu7*tau10+64567434600*mu6*tau11-4009274670*mu5*tau12...
     -32405563128*mu6*tau10+39962916924*mu5*tau11-2844180*mu4*tau12+...
     4322203276*mu6*tau9-86495791886*mu5*tau10+1569566700*mu4*tau11...
     -82229968*mu6*tau8+57846707984*mu5*tau9-24671708676*mu4*tau10...
     -11535226146*mu5*tau8+72070323500*mu4*tau9-260412165*mu3*tau10...
     +453500204*mu5*tau7-62403592332*mu4*tau8+8501653944*mu3*tau9...
     -1278020*mu5*tau6+16856031940*mu4*tau7-36498020364*mu3*tau8...
     -1035350060*mu4*tau6+41264324712*mu3*tau7-1257542496*mu2*tau8...
     +6369444*mu4*tau5-14506908830*mu3*tau6+10254654432*mu2*tau7-...
     9136*mu4*tau4+1250732360*mu3*tau5-16049280384*mu2*tau6...
     -12674316*mu3*tau4+7309198080*mu2*tau5-1218576960*mu*tau6+...
     30488*mu3*tau3-841937760*mu2*tau4+3236526720*mu*tau5+...
     571*mu3*tau2+12632544*mu2*tau3-1977048000*mu*tau4-36288*mu2*tau2...
     +298798848*mu*tau3-235146240*tau4-1728*mu2*tau-6277824*mu*tau2...
     +217728000*tau3+10368*mu*tau-43545600*tau2+1728*mu+1244160*tau)...
      /mutau2m112/mutaum14;
elseif n == 5 
  mu3=mu2*mu;
  mu4=mu3*mu;
  mu5=mu4*mu;
  mu6=mu5*mu;
  mu7=mu6*mu;
  mu8=mu7*mu;
  mu9=mu8*mu;
  mu10=mu9*mu;
  mu11=mu10*mu;
  mu12=mu11*mu;
  mu13=mu12*mu;
  mu14=mu13*mu;
  tau5=tau4*tau;
  tau6=tau5*tau;
  tau7=tau6*tau;
  tau8=tau7*tau;
  tau9=tau8*tau;
  tau10=tau9*tau;
  tau11=tau10*tau;
  tau12=tau11*tau;
  tau13=tau12*tau;
  tau14=tau13*tau;
  tau15=tau14*tau;
  tau16=tau15*tau;
  tau17=tau16*tau;
  tau18=tau17*tau;
  tau19=tau18*tau;
  tau20=tau19*tau;
  tau21=tau20*tau;
  tau22=tau21*tau;
  tau23=tau22*tau;
  tau24=tau23*tau;
  tau25=tau24*tau;
  tau26=tau25*tau;
  tau27=tau26*tau;
  tau28=tau27*tau;
  tau29=tau28*tau;
  mutaum12=mutaum1*mutaum1;
  mutaum13=mutaum12*mutaum1; 
  mutaum14=mutaum13*mutaum1; 
  mutaum15=mutaum14*mutaum1; 
  mutau2m16=mutau2m13*mutau2m13;
  mutau2m19=mutau2m16*mutau2m13;
  mutau2m112=mutau2m19*mutau2m13; 
  mutau2m115=mutau2m112*mutau2m13; 
  fg=-1/209018880*mu*tau6*(163879*mu14*tau29+51961*mu14*tau28...
     -609098*mu14*tau27-819395*mu13*tau28-786814*mu14*tau26...
     -2761957*mu13*tau27+18879539*mu14*tau25+4643870*mu13*tau26...
     +1638790*mu12*tau27-226718347*mu14*tau24+15544662*mu13*tau25...
     +13034367*mu12*tau26-568911959*mu13*tau24+1968834*mu12*tau25...
     -1638790*mu11*tau26+14215984447*mu13*tau23-97150697*mu12*tau24...
     -25557118*mu11*tau25-22442472628*mu13*tau22+4848464754*mu12*tau23...
     -58247532*mu11*tau24+819395*mu10*tau25-190209820899*mu12*tau22+...
     227453744*mu11*tau23+25305307*mu10*tau24+557830157734*mu12*tau21...
     -19951693006*mu11*tau22+141740872*mu10*tau23-163879*mu9*tau24...
     -341662517523*mu12*tau20+1129160678898*mu11*tau21...
     -168506758*mu10*tau22-12578709*mu9*tau23-4810869930640*mu11*tau20...
     +47596236037*mu10*tau21-152292756*mu9*tau22+5349046421180*mu11*tau19...
     -3748733329103*mu10*tau20-157461080*mu9*tau21+...
     2506149*mu8*tau22-1518771856656*mu11*tau18+21445932135974*mu10*tau19...
     -70959965853*mu9*tau20+78909320*mu8*tau21-33983549214008*mu10*tau18...
     +7697540177397*mu9*tau19+361741639*mu8*tau20+...
     17732722429594*mu10*tau17-57554437847498*mu9*tau18+...
     67487840518*mu8*tau19-16113510*mu7*tau20-2415740878110*mu10*tau16...
     +119685779939022*mu9*tau17-10206271195099*mu8*tau18...
     -233955666*mu7*tau19-89169896914134*mu9*tau16+...
     99739751279332*mu8*tau17-39958635430*mu7*tau18+...
     22916849979674*mu9*tau15-263680721374601*mu8*tau16+...
     8788521704858*mu7*tau17+53120970*mu6*tau18-1432757733624*mu9*tau14+...
     256484437391078*mu8*tau15-114802338968730*mu7*tau16+...
     13464713115*mu6*tau17-95430671843570*mu8*tau14+...
     383878606871006*mu7*tau15-4760419029975*mu6*tau16+...
     11622357950056*mu8*tau13-469437600792850*mu7*tau14+...
     87537208811668*mu6*tau15-1976927715*mu5*tau16...
     -291120246102*mu8*tau12+229436532063102*mu7*tau13...
     -376567603133040*mu6*tau14+1476932356935*mu5*tau15...
     -41356829752504*mu7*tau12+572380085309538*mu6*tau13...
     -42610665786192*mu5*tau14+2089045526876*mu7*tau11...
     -351721795944554*mu6*tau12+246854268235548*mu5*tau13...
     -200510972991*mu4*tau14-15515805072*mu7*tau10+...
     84538189051092*mu6*tau11-470778436877402*mu5*tau12...
     +12016334849274*mu4*tau13-6502842188850*mu6*tau10+...
     357518772856658*mu5*tau11-103785546173499*mu4*tau12+...
     100884700315*mu6*tau9-109058409840772*mu5*tau10+...
     257328962259480*mu4*tau11-1496300589504*mu3*tau12...
     -145233639*mu6*tau8+11450092993836*mu5*tau9...
     -242117649171342*mu4*tau10+25313069877120*mu3*tau11...
     -279648512507*mu5*tau8+91671072496444*mu4*tau9-...
     89001855959904*mu3*tau10+847910431*mu5*tau7...
     -12445982186078*mu4*tau8+106429659813792*mu3*tau9...
     -2721696305760*mu2*tau10+3277580*mu5*tau6+428011918936*mu4*tau7...
     -49874490975456*mu3*tau8+17417246973696*mu2*tau9...
     -2063292787*mu4*tau6+8525788578720*mu3*tau7...
     -28351406885760*mu2*tau8-17283238*mu4*tau5-390203087904*mu3*tau6...
     +16755661907712*mu2*tau7-1436872296960*mu*tau8-163879*mu4*tau4...
     +2681952480*mu3*tau5-3578660592192*mu2*tau6+...
     3916478200320*mu*tau7+37161504*mu3*tau4+211650243840*mu2*tau5...
     -3095195120640*mu*tau6+823392*mu3*tau3-1968582528*mu2*tau4+...
     835944883200*mu*tau5-181062604800*tau6-41423616*mu2*tau3...
     -63125733888*mu*tau4+230443315200*tau5-1652832*mu2*tau2+...
     761799168*mu*tau3-82301184000*tau4+24883200*mu*tau2...
     +7965941760*tau3+1658880*mu*tau-121927680*tau2-4976640*tau-829440)...
     /mutau2m115/mutaum15;
 end
end
%% Gamma* function
function [g] = gamstar(x)
giant=realmax/1000;
if (x>=3.0)
    g = exp(stirling(x));
elseif (x>0.0)
    g = gamma(x)/(exp(-x+(x-0.5)*log(x))*sqrt(2*pi));
else
    g = giant;
end
end
%% Stirling
function [s] =  stirling(x)
%Stirling series, function corresponding with asymptotic series for log(gamma(x))}
% that is:  1/(12x)-1/(360x**3)...; x>= 3}
dwarf=realmin*1000.0;
giant=realmax/1000;
lnsqrttwopi=0.9189385332046727418;
if (x<dwarf)
    s =giant;
elseif (x<1.0)
    ln1=log(gamma(1+x));
    s = ln1-(x+0.5)*log(x)+x-lnsqrttwopi;
elseif (x<2.0)
    ln1=log(gamma(x));
    s =ln1-(x-0.5)*log(x)+x-lnsqrttwopi;
elseif (x<3.0)
    ln1=log(gamma(-1+x));
    s =ln1-(x-0.5)*log(x)+x-lnsqrttwopi+log(x-1);
elseif (x<12.0)
    a=[1.996379051590076518221;
        -0.17971032528832887213e-2;
        0.131292857963846713e-4;
        -0.2340875228178749e-6;
        0.72291210671127e-8;
        -0.3280997607821e-9;
        0.198750709010e-10;
        -0.15092141830e-11;
        0.1375340084e-12;
        -0.145728923e-13;
        0.17532367e-14;
        -0.2351465e-15;
        0.346551e-16;
        -0.55471e-17;
        0.9548e-18;
        -0.1748e-18;
        0.332e-19;
        -0.58e-20];
    z=18.0/(x*x)-1.0;
    s =chepolsum(17,z,a)/(12.0*x);
else
    z=1.0/(x*x);
    if (x<1000.0)
        c=[0.25721014990011306473e-1;
            0.82475966166999631057e-1;
            -0.25328157302663562668e-2;
            0.60992926669463371e-3;
            -0.33543297638406e-3;
            0.250505279903e-3;
            0.30865217988013567769];
        s =((((((c(6)*z+c(5))*z+c(4))*z+c(3))*z+c(2))*z+c(1))/(c(7)+z)/x);
    else
        s =(((-z*0.000595238095238095238095238095238+...
            0.000793650793650793650793650793651)*z...
            -0.00277777777777777777777777777778)*z+...
            0.0833333333333333333333333333333)/x;
    end
end
end
%% Function chepolsum
function [chep]=chepolsum(n,t,ak)
u0=0; u1=0; k=n; tt=t+t;
while k>=0
  u2=u1; 
  u1=u0; 
  u0=tt*u1-u2+ ak(k+1); 
  k= k-1; 
end
s=(u0-u2)/2.0;
chep=s;
end
function zeta1=h0(zeta,x,y,p,q)
r= p+q; xi= x*y/(2*r); c2= p/r; s2= q/r;
w= sqrt((xi-c2)*(xi-c2)+4*xi);
tp= 1/y; 
if xi > c2 
  t0= (xi-c2+w)/(2*xi);
else
  t0= 2/(c2-xi+w);
end  
if abs(zeta) < 0.001
  c22=c2*c2;
  s22=s2*s2;
  tp2=tp*tp;
  tp3=tp2*tp;
  tp4=tp3*tp;
  tp5=tp4*tp;
  tp6=tp5*tp;
  tp7=tp6*tp;
  tp8=tp7*tp;
  tp9=tp8*tp;
  tp10=tp9*tp;
  h00=((1/3)*(2*c2*s2*tp5-3*c2*s2*tp4-2*c2*tp5...
     +4*c2*tp4-6*s2*tp4-2*c2*tp3+12*s2*tp3+...
     6*tp4-5*s2*tp2-18*tp3+20*tp2-10*tp+2)...
     /((c2*tp2-2*tp+1)*(s2*tp2-tp2+2*tp-1)^(3/2)));      
  h01=((1/12)*(c22*s22*tp10-12*c22*s22*tp9-2*c22*s2*tp10+...
    12*c22*s22*tp8+16*c22*s2*tp9+c22*tp10-6*c2*s22*tp9-...
    20*c22*s2*tp8-4*c22*tp9+63*c2*s22*tp8+12*c2*s2*tp9+...
    6*c22*tp8-92*c2*s22*tp7-102*c2*s2*tp8-6*c2*tp9+...
    6*c22*s2*tp6-4*c22*tp7+33*c2*s22*tp6+208*c2*s2*tp7+...
    39*c2*tp8-48*s22*tp7+c22*tp6-164*c2*s2*tp6-116*c2*tp7+...
    102*s22*tp6+48*s2*tp7+52*c2*s2*tp5+203*c2*tp6-...
    66*s22*tp5-132*s2*tp6-6*c2*s2*tp4-222*c2*tp5+13*s22*tp4+...
    120*s2*tp5+6*tp6+149*c2*tp4-32*s2*tp4-30*tp5-...
    56*c2*tp3-8*s2*tp3+61*tp4+9*c2*tp2+4*s2*tp2-...
    64*tp3+36*tp2-10*tp+1)/((c2*tp2-2*tp+1)^2*...
    (-c2*tp2+2*tp-1)^(1/2)*(s2*tp2-tp2+2*tp-1)^(5/2)));    
  h0=h00+h01*zeta;
else
  h0=(t0-1)/(1-y*t0)/sqrt(s2*t0^2-(t0-1)*(t0-1))-1/zeta;
end
pp=h0*zeta;
if y==0 
  zeta1= h0; 
elseif abs(y) < 0.1 
  zeta1= 2*atanh(pp/(2+pp))/zeta; 
else 
  zeta1= log(1 + pp)/zeta;
end
end
