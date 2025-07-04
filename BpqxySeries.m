function [betanc,ierr]=BpqxySeries(x,y,p,q)
% Computation of the series of B_(p,q)(x,y)  
% using backward recursion
ierr=0;
z=x/2; 
j0= j0proc(z);
iij=betaincreg(y, p+j0, q);
while iij<realmin && j0>2
  j0=j0-1;
  iij=betaincreg(y, p+j0, q);
end
Ij(j0+1)=iij;
pj=exp(-z+j0*log(z)-gammaln(j0+1));
s= pj*Ij(j0); 
Ij(j0)=betaincreg(y,p+j0-1,q);
iij=Ij(j0);
pj= pj*j0/z; 
s= s + pj*Ij(j0-1);  
j= j0-1;  
k=0;
while j > 0 && k==0
  Ij(j)=((p+j+(p+j+q-1)*y)*Ij(j+1)-(p+j)*Ij(j+2))/(y*(p+j+q-1));
  iij=Ij(j);
  if iij<1
    pj=j*pj/z;
    term= pj*Ij(j);
    s= s+term;
    j= j-1;
  else
    k=1;
  end
end
if k==1
  ierr=1;
end    
betanc=s;
end


function j0=j0proc(z) 
% to compute the number of terms in the series 
e=exp(1); eps100= 10*eps;
r=-log(eps)/(z*e); 
a= e/(e-1)^2; t0=1/e+sqrt(r/a); d= 1;
while d>eps100 
   t1=(e*(r+t0)-1)/(e*(log(t0)+1));
   d= abs(t0/t1-1); 
   t0=t1;
end
j0=1+fix(z*e*t0);
end

function [I] = betaincreg(x,p,q)
%--------------------------------------------------------
%Function implementing the algorithm given in the paper
%"Computation of the regularized incomplete beta function"
% by V. Egorova, A.Gil, J.Segura and N.M. Temme
%---------------------------------------------------------
if x<= 0
  I = 0;
elseif x>=1
  I = 1;
else
  xt=p/(p+q);
  if p>50 && q >50 && p+q>700 && abs(x-xt)<0.2
    % error function approximation
    I = errorfunapprox(x,p,q);
  elseif p>100 && q<10 
    if x<0.65    
    % series expansion 2
      I = series2(x,p,q);
    else
    % Continued fraction  
      I = betainc_contfrac(x,p,q);  
    end    
  else   
    % for q above the line the most precise algorithm is continued fraction
    if q > (1-x)/x*p
      I = betainc_contfrac(x,p,q);
    else
    % for q below the line - series expansion 3
      I = series3(x,p,q);
    end
  end
end
end
%% Error function approximation
function [I] = errorfunapprox(x,p,q)
r= p+q; s=sqrt(p/r); c=sqrt(q/r);
s2=s*s; s3=s2*s; s4=s3*s; s5=s4*s; s6=s5*s; s7=s6*s; s8=s7*s; s9=s8*s;
s10=s9*s; s11=s10*s; s12=s11*s; s13=s12*s; s14=s13*s; s15=s14*s;
s16=s15*s; s17=s16*s; s18=s17*s; s19=s18*s; s20=s19*s; s21=s20*s;
s22=s21*s; s23=s22*s; s24=s23*s; s25=s24*s; s26=s25*s; s27=s26*s;
s28=s27*s; s29=s28*s; s30=s29*s; s31=s30*s; s32=s31*s; s33=s32*s; s34=s33*s;
s35=s34*s; s36=s35*s; s37=s36*s; s38=s37*s; s39=s38*s; s40=s39*s;
s41=s40*s; s42=s41*s; s43=s42*s; s44=s43*s; s45=s44*s; s46=s45*s; s47=s46*s;
s48=s47*s; 
c2=c*c; c3=c2*c; c4=c3*c; c5=c4*c; c6=c5*c; c7=c6*c; c8=c7*c; c9=c8*c; c10=c9*c;
c11=c10*c; c12=c11*c; c13=c12*c; c14=c13*c; c15=c14*c; c16=c15*c; c17=c16*c;
c18=c17*c; c19=c18*c; c20=c19*c; c21=c20*c; c22=c21*c; c23=c22*c; c24=c23*c;
c25=c24*c;
sfactor=1.0+s4-s2;
sfactor2=sfactor*sfactor;
ak(1)= 1.0;
ak(2)= (1.0/3.0)*(2.0*s2-1.0)/(s*c);
ak(3)= (1.0/12.0)*sfactor/(s2*c2);
ak(4)= -(1.0/135.0)*(s2-2.0)*(2.0*s2-1.0)*(s2+1.0)/(s3*c3);
ak(5)= (1.0/864.0)*sfactor2/(s4*c4);          
ak(6)= (1.0/5670.0)*(s2-2.0)*(2.0*s2-1.0)*...
           (s2+1.0)*sfactor/(s5*c5);
ak(7)=-(1.0/777600.0)*(139.0-417.0*s2+402.0*s4...
           +139.0*s12-417.0*s10+402.0*s8-109.0*s6)/(s6*c6);
ak(8)=(1.0/51030.0)*(s2-2.0)*(2.0*s2-1.0)...
           *(s2+1.0)*sfactor2/(s7*c7);
ak(9)=-(1.0/261273600.0)*sfactor*(571.0*s12-1713.0*s10+1698.0*s8...
          -541.0*s6+1698.0*s4-1713.0*s2+571.0)/(s8*c8);
ak(10)=-(1.0/303118200.0)*(s2-2.0)*(2.0*s2-1.0)*(s2+1.0)...
           *(281.0*s12-843.0*s10+1308.0*s8-1211.0*s6...
           +1308.0*s4-843.0*s2+281.0)/(s9*c9);
ak(11)= (1.0/197522841600.0)*(163879.0*s12-491637.0*s10+159882.0*s8...
            +499631.0*s6+159882.0*s4-491637.0*s2+163879.0)*sfactor2...
           /(s10*c10);
ak(12)=-(1.0/59108049000.0)*(s2-2.0)*...
           (2.0*s2-1.0)*(s2+1.0)*sfactor...
           *(5221.0*s12-15663.0*s10+20688.0*s8-15271.0*s6...
            +20688.0*s4-15663.0*s2+5221.0)/(s11*c11);
ak(13)= (1.0/782190452736000.0)*...
          (5246819.0-31480914.0*s2+83618655.0*s4+5246819.0*s24...
          -31480914.0*s22+83618655.0*s20-129518230.0*s18+129487374.0*s16...
          -87130170.0*s14+64799751.0*s12-87130170.0*s10+129487374.0*s8...
          -129518230.0*s6)/(s12*c12);
ak(14)=(1.0/1063944882000.0)*(s2-2.0)*(2.0*s2-1.0)*...
          (s2+1.0)*(5459.0*s12-16377.0*s10+15042.0*s8...
          -2789.0*s6+15042.0*s4-16377.0*s2+5459.0)*...
           sfactor2/(s13*c13);
ak(15)=-(1.0/122021710626816000.0)*sfactor*(534703531.0*s24...
          -3208221186.0*s22+7480602615.0*s20-7994318870.0*s18...
           +3156814926.0*s16+48220470.0*s14+499100559.0*s12...
           +48220470.0*s10+3156814926.0*s8-7994318870.0*s6...
           +7480602615.0*s4-3208221186.0*s2+534703531.0)/(s14*c14);
ak(16)=(1.0/199409869508850000.0)*(s2-2.0)*(2.0*s2-1.0)*...
           (s2+1.0)*(91207079.0*s24-547242474.0*s22+1545209655.0*s20...
           -2709658930.0*s18+3495277134.0*s16-3742822170.0*s14...
           +3827266491.0*s12-3742822170.0*s10+3495277134.0*s8...
           -2709658930.0*s6+1545209655.0*s4-547242474.0*s2+91207079.0)/(s15*c15);
ak(17)=-(1.0/175711263302615040000.0)*(4483131259.0-26898787554.0*s2...
           +62558334615.0*s4-66219453830.0*s6+25800142014.0*s8-1770508170.0*s10...
           +8577414591.0*s12-1770508170.0*s14+25800142014.0*s16-66219453830.0*s18...
           +62558334615.0*s20-26898787554.0*s22+4483131259.0*s24)*sfactor2/(s16*c16);
ak(18)=-(1.0/90930900496035600000.0)*(s2-2.0)*(2.0*s2-1.0)*(s2+1.0)*...
           sfactor*(2650986803.0*s24-15905920818.0*s22+41867798685.0*s20...
           -63534719260.0*s18+68256815688.0*s16-66784076190.0*s14...
           +69549216987.0*s12-66784076190.0*s10+68256815688.0*s8...
           -63534719260.0*s6+41867798685.0*s4-15905920818.0*s2...
           +2650986803.0)/(s17*c17);
ak(19)=(1.0/17743323368298066739200000.0)*(-3890357294511339.0*s2...
           +15559210028072151.0*s4+36125299363900566.0*s24-15384557050481664.0*s22...
           +3759499898757441.0*s20-326522671455545.0*s18+3759499898757441.0*s16...
           -15384557050481664.0*s14+36125299363900566.0*s12-54331224371473014.0*s10...
           +54401508016310970.0*s8-36292248215653524.0*s6+432261921612371.0*s36...
           -3890357294511339.0*s34+15559210028072151.0*s32-36292248215653524.0*s30...
           +54401508016310970.0*s28-54331224371473014.0*s26+432261921612371.0)/(s18*c18);
ak(20)=-(1.0/2455134313392961200000.0)*(s2-2.0)*(2.0*s2-1.0)*(s2+1.0)*...
           (6171801683.0*s24-37030810098.0*s22+92136923385.0*s20-121235524360.0*s18...
           +106807138668.0*s16-107154319590.0*s14+126781382307.0*s12-107154319590.0*s10...
           +106807138668.0*s8-121235524360.0*s6+92136923385.0*s4...
           -37030810098.0*s2+6171801683.0)*sfactor2/(s19*c19);
ak(21)=(1.0/56636688191607429031526400000.0)*sfactor*(6232523202521089.0*s36...
           -56092708822689801.0*s34+224189978990811309.0*s32-522085098612188316.0*s30...
           +780826142447376030.0*s28-779315338018394226.0*s26+523118213261102994.0*s24...
           -233830339096133376.0*s22+70321992581529819.0*s20-20498208665349955.0*s18...
           +70321992581529819.0*s16-233830339096133376.0*s14+523118213261102994.0*s12...
           -779315338018394226.0*s10+780826142447376030.0*s8-522085098612188316.0*s6...
           +224189978990811309.0*s4-56092708822689801.0*s2+6232523202521089.0)/(s20*c20);
ak(22)=(1.0/25410640143617148420000000.0)*(s2-2.0)*(2.0*s2-1.0)*(s2+1.0)...
           *(4283933145517.0*s36-38555398309653.0*s34+161987808938877.0*s32...
            -421980109825548.0*s30+772871761776690.0*s28-1084178676388878.0*s26...
            +1263214939279182.0*s24-1324071137771028.0*s22+1334609313822207.0*s20...
            -1332080936189215.0*s18+1334609313822207.0*s16-1324071137771028.0*s14...
            +1263214939279182.0*s12-1084178676388878.0*s10+772871761776690.0*s8...
            -421980109825548.0*s6+161987808938877.0*s4-38555398309653.0*s2...
            +4283933145517.0)/(s21*c21);
ak(23)=-(1.0/185541790515705937507280486400000.0)*(25834629665134204969.0...
            -232511666986207844721.0*s2+878296385812930194789.0*s4-1756106634816063744636.0*s6...
            +1883776220929228333230.0*s8-900493018129305149346.0*s10...
            +74123143631665666674.0*s12-127743638304015861696.0*s14+...
             77056378023492043299.0*s16+181371030011418519845.0*s18+...
             77056378023492043299.0*s20-127743638304015861696.0*s22+...
             74123143631665666674.0*s24-900493018129305149346.0*s26+...
             1883776220929228333230.0*s28-1756106634816063744636.0*s30...
             +878296385812930194789.0*s32-232511666986207844721.0*s34...
             +25834629665134204969.0*s36)*sfactor2/(s22*c22);
ak(24)=(1.0/838551124739365897860000000.0)*(s2-2.0)*(2.0*s2-1.0)*...
            (s2+1.0)*sfactor*(11963983648109.0*s36-107675852832981.0*s34...
             +441707388667209.0*s32-1093006445123436.0*s30+1850727807014670.0*s28...
             -2369766184192386.0*s26+2570628833094594.0*s24-2620121222397816.0*s22...
             +2630635026351039.0*s20-2618222684809895.0*s18+2630635026351039.0*s16...
             -2620121222397816.0*s14+2570628833094594.0*s12-2369766184192386.0*s10...
             +1850727807014670.0*s8-1093006445123436.0*s6+441707388667209.0*s4...
             -107675852832981.0*s2+11963983648109.0)/(s23*c23);
ak(25)=-(1.0/3072572050940090325120564854784000000.0)*(-18948349666259029037148.0*s2...
            +105766000944332712002454.0*s4+52377964859222906383913.0*s24...
            -122531573397975640756968.0*s22+381108379887241904073102.0*s20...
             -879586600317678654144548.0*s18+1513169866999911172741362.0*s16...
             -1974991189331755995042492.0*s14+1973542318787418583200482.0*s12...
             -1507009328043931734547848.0*s10+866939243266728543070353.0*s8...
             -364437266127070774293920.0*s6-1507009328043931734547848.0*s38...
             +1973542318787418583200482.0*s36-1974991189331755995042492.0*s34...
             +1513169866999911172741362.0*s32-879586600317678654144548.0*s30...
             +381108379887241904073102.0*s28-122531573397975640756968.0*s26...
             +105766000944332712002454.0*s44-364437266127070774293920.0*s42...
             +866939243266728543070353.0*s40+1579029138854919086429.0*s48...
             +1579029138854919086429.0-18948349666259029037148.0*s46)/(s24*c24);
ak(26)= -(1.0/211314883434320206260720000000.0)*(s2-2.0)*(2.0*s2-1.0)...
             *(s2+1.0)*(208697624924077.0*s36-1878278624316693.0*s34...
             +7476558808459167.0*s32-17238154983161628.0*s30+...
              25903980761952330.0*s28-28670257324128798.0*s26+28542389728823772.0*s24...
             -29410943594470128.0*s22+29377831383649017.0*s20-28414949938538155.0*s18...
             +29377831383649017.0*s16-29410943594470128.0*s14+28542389728823772.0*s12...
             -28670257324128798.0*s10+25903980761952330.0*s8-17238154983161628.0*s6...
             +7476558808459167.0*s4-1878278624316693.0*s2+208697624924077.0)*sfactor2...
              /(s25*c25);          
              
dk(26)=0.0; dk(27)=0.0;
for k=24:-1:0 
  dk(k+1)=ak(k+2)+(k+2.0)*dk(k+3)/r;
end
rh=p/r*log(r*x/p)+q/r*log(r*(1.0-x)/q);
if (rh>-5.0e-20) 
  eta=0.0; 
else
  eta=sqrt(-2.0*rh); 
end
if (x-p/r<0) 
  eta=-eta; 
end
seta=dk(1);
etak=1.0;
for k=2:23
  etak=etak*eta;
  seta=seta+dk(k)*etak;
end
Reta=r/(r+dk(2))*exp(r*rh)/sqrt(2.0*pi*r)*seta;
I=0.5*erfc(-eta*sqrt(r/2.0))-Reta;
end 
%% Series 2
function [I] =series2(x,p,q)
s= 1.0; t= 1.0; k= 0; del= 1.0;
xx=x/(x-1.0); 
while del > eps && k<1000 
  k= k+1; t= t*(k-q)*xx/(k+p); 
  term=t;
  s=s+term;
  if s==0
    del=1.0; 
  else
    del=abs(term/s);
  end
end   
F=Ffactor (x,p,q);
I=s*F/p/(1-x);
end
%% Continued fraction
function [I] = betainc_contfrac (x,p,q)
if x<= 0
    I = 0;
elseif x>=1
    I = 1;
else
    pq = p+q;
    r  = p / pq;
    recur = false;
    if x > r
        [p,q] = swap(p,q);
        x = 1-x;
        recur = true;
    end
    f=continued_fraction (x,p,q);
    D=Ffactor (x,p,q);
    I = f*D/p;
    if recur
        I = 1-I;
    end
end
end
%% Lentz-Thompson algorithm for the continued fraction
function [h] = continued_fraction (x,p,q)
fpmin=1e-30;
qab=p+q;
qap=p+1.0;
qam=p-1.0;
c=1.0;
d=1.0-qab*x/qap;
if (abs(d)<fpmin)
    d=fpmin;
end
d=1.0/d;
h=d;
m=0;
del=0.0;
while (abs(del-1.0) > eps) && (m<100)
    m=m+1;
    m2=2.0*m;
    aa=m*x*(q-m)/((qam+m2)*(p+m2));
    d=1.0+aa*d;
    if (abs(d)<fpmin)
        d=fpmin;
    end
    c=1.0 +aa/c;
    if (abs(c)<fpmin)
        c=fpmin;
    end
    d=1.0/d;
    h=h*d*c;
    aa=-(p+m)*(qab+m)*x/((qap+m2)*(p+m2));
    d=1.0 +aa*d;
    if (abs(d)<fpmin)
        d=fpmin;
    end
    c=1.0 +aa/c;
    if (abs(c)<fpmin)
        c=fpmin;
    end
    d=1.0 /d;
    del=d*c;
    h=h*del;
end
end

%% Series expansion
function [I] = series3 (x,p,q)
if x<= 0
    I = 0;
elseif x>=1
    I = 1;
else
    pq = p+q;
    r  = p / pq;
    recur = false;
    if x > r
        [p,q] = swap(p,q);
        x = 1-x;
        recur = true;
    end
    s= 1.0; t= 1.0; k= 0; del= 1.0;
    while del > eps && k<1000
        k= k+1; 
        t= t*(pq+k-1.0)*x/(k+p); 
        term=t;
        s=s+term;
        if (s == 0)
            del=1.0 ;
        else
            del=abs(term/s);
        end
    end
    D=Ffactor (x,p,q);
    I=s*D/p;
    if recur
        I = 1-I;
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
b = F*gammastar(p+q)/(gammastar(p)*gammastar(q));
end

%% Gamma* function
function [g] = gammastar (x)
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
%% Swap function
function [b,a] = swap(a,b)
end
