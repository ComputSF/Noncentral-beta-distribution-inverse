function [Bdy,ierr]=Bpqxydery(x,y,p,q)
%--------------------------------------------------------
%Noncentral beta: derivative with respect to y
%We use the expression in terms of the Kummer-M function
%--------------------------------------------------------
ierr=0;
if y<1e-300
  Bdy=0;
  ierr=1;
else  
  argu=0.5*x*y;
  a=p+q;
  b=p;
  [m1,ier]=Mabx(a,b,argu);
  if ier>0
    ierr=1;
    Bdy=0;
  else  
    D=Ffactor(y,p,q);
    factor=exp(-0.5*x)*D/(y*(1-y));
    Bdy=factor*m1;
  end
end   
end

%% Factor x^p*(1-x)^q/Beta(p,q)
function [b] = Ffactor (x,p,q)
s=p+q;
xt=p/s;
sigma=(x-xt)/xt;
tau=(xt-x)/(1-xt);
if abs(sigma) <0.01
  sigma2=sigma*sigma;sigma3=sigma2*sigma;
  sigma4=sigma3*sigma;sigma5=sigma4*sigma;
  sigma6=sigma5*sigma;sigma7=sigma6*sigma;
  sigma8=sigma7*sigma;sigma9=sigma8*sigma;
  Ps=-(1/2)*sigma2+(1/3)*sigma3-(1/4)*sigma4...
     +(1/5)*sigma5-(1/6)*sigma6+(1/7)*sigma7...
     -(1/8)*sigma8+(1/9)*sigma9;
else
  Ps=log(1+sigma)-sigma;  
end    
if abs(tau) <0.01
   tau2=tau*tau;tau3=tau2*tau;
   tau4=tau3*tau;tau5=tau4*tau;
   tau6=tau5*tau;tau7=tau6*tau;
   tau8=tau7*tau;tau9=tau8*tau;
   Pt=-(1/2)*tau2+(1/3)*tau3-(1/4)*tau4...
     +(1/5)*tau5-(1/6)*tau6+(1/7)*tau7...
     -(1/8)*tau8+(1/9)*tau9;
else
  Pt=log(1+tau)-tau;  
end  
E=exp(p*Ps+q*Pt);
F=sqrt(p*q/s/(2*pi))*E;
b = F*gamstar(p+q)/(gamstar(p)*gamstar(q));
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
if isfinite(Mkummer)==0
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

function quotgam=quotgamm(x,y) 
%----------------------------------------------------
% Computation of the quotient of two gamma functions
%  gamma(x)/gamma(y)
%----------------------------------------------------
if (x <= 0.0)||(y <= 0.0)
  quotgam=qratio(x,y);
elseif x>y
  quotgam=1.0/quotgamm(y,x);
else
  n=floor(y-x);
  if n == (y - x)
    quotgam=1.0/shiftfact(x, n);
  elseif (n>15) 
    quotgam=qratio(x,y);
  elseif (n>0)
    quotgam=quotgamm(x+n,y)/shiftfact(x,n);
  elseif (x<26.0) 
    quotgam=qratio(x,y);
  else
    quotgam=qasym(x,y);
  end
end
end

function qrat=qratio(x,y)
if (x<=0.0)||(y<=0.0)||(y>2.5*x) 
  qrat=gamma(x)/gamma(y);
else
  b=y-x;
  c=b/x;
  q=lnec(c);
  q=exp(-b*log(x)-(b-0.5)*log(1.0+c)-x*q);
  qrat=q*gamstar(x)/gamstar(y);
end
end

function qas=qasym(x,y)  
w=(x+y-1)/2.0;
w2=1.0/(w*w);
r=(x-y+1)/2.0;
r2=r*r;
r3=r2*r;
r4=r2*r2;
r5=r4*r;
r6=r5*r;
r7=r6*r;
cc(1)=1.0;
cc(2)=r/12.0;
cc(3)=r/1440.0+r2/288.0;
cc(4)=r/90720.0+r2/17280.0+r3/10368.0;
cc(5)=r/4838400.0+101*r2/87091200.0+r3/414720.0+r4/497664.0;
cc(6)=r/239500800.0+13*r2/522547200.0+61*r3/1045094400.0+...
      r4/14929920.0+r5/29859840.0;
cc(7)=691.0*r/7846046208000.0+7999.0*r2/14485008384000.0+...
      59.0*r3/41803776000.0+143.0*r4/75246796800.0...
      +r5/716636160.0+r6/2149908480.0;
cc(8)=r/523069747200.0+ 2357*r2/188305108992000.0...
      +5941.0*r3/173820100608000.0...
      +11*r4/214990848000.0+41*r5/902961561600.0...
      +r6/42998169600.0+r7/180592312320.0;
s=1.0;
k=1;
t=1.0;
u=1.0;
while (abs(u)>eps)&&(k< 8)			
  t=-4.0*w2*(k-r)*(k-r-0.5)*t;
  u=t*cc(k+1);
  s=s+u;
  k=k+1;
  qas=s*exp((x-y)*log(w));
end
end

function shift=shiftfact(x,n)
if n==0
  shift=1.0;
elseif (n<0) 
  shift=1.0/shiftfact(x-n,n);
else
  s=1.0;
  k=0;
  while (k<n)
    s=s*(x+k);
    k=k+1;
  end					
  shift=s;
end
end
    
function y=exmin1minx(x)
%Computes (exp(x)-1-x)/(0.5*x*x) 
if x==0
  y=1.0;
elseif abs(x)>0.9 
  y=(exp(x)-1.0-x)/(x*x*0.5);
else
  t=sinh(x*0.5);
  t2=t*t;
  y=(2*t2+(2.0*t*sqrt(1.0+t2)-x))/(x*x*0.5);
end
end

function ln1=lnec(x)
%x>-1; lnec=ln1=ln(1+x)-x
z=logoneplusx(x);
y0=z-x;
e2=exmin1minx(z);
s=e2*z*z/2;
r=(s+y0)/(s+1+z);
ln1=y0-r*(6.0-r)/(6.0-4.0*r);
end


function s=oddchepolsum(n, x, ak)
if (n==0) 
  s=ak(1)*x; 
elseif (n == 1)
  s=x*(ak(1)+ak(2)*(4*x*x - 3)); 
else
  y=2*(2*x*x - 1); 
  r= ak(n+1); h= ak(n)+ r*y;  
  k=n-2;
  while (k >= 0)
    s=r; r= h; h= ak(k+1)+r*y-s; 
    k= k-1; 
  end
  s=x*(h-r);
end
end

function y=logoneplusx(t)
if (-0.2928<t)&&(t<0.4142)
  p= 1.18920711500272106671749997;
  p=(p-1.0)/(p+1.0);
  pj=p; ck(1)= pj; p2= p*p; j=1; c=1;
  while (abs(c)> 1.0e-20)&&(j<1000)   
    pj=pj*p2; c=pj/(2.0*j+1.0); 
    ck(j+1)=c; j= j+1; 
  end
  x=t/(2.0+t)*(1.0+p2)/(2.0*p);
  y=4*oddchepolsum(j-1, x, ck);
else
  y=log(1.0+t); 
end
end
