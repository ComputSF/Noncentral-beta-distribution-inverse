function [Bpq,ierr]=Bpqxy(x,y,p,q)
% Computes the Noncentral Beta Distribution
%
% This function calculates the noncentral beta distribution Bpq(x, y)
% using various numerical methods to ensure accuracy across a wide range of
% input parameters.
%
% INPUT:
%   x - noncentrality parameter (lambda) (x>0)
%   y - quantile (should be in the range [0, 1])
%   p - Shape parameter p (p>0)
%   q - Shape parameter q (q>0)
%
% OUTPUT:
%   Bpq  - Computed value of the noncentral beta distribution
%   ierr - Error flag indicating the success or issues in computation:
%          0: Computation successful
%          1: Computation successful, but some loss of accuracy is expected
%          2: Computation failed due to one or more input values being out of bounds
%
% USAGE:
%   [Bpq, ierr] = Bpqxy(x, y, p, q)
%
% NOTES:
%   - Ensure that the input parameters are within their valid ranges to 
%     avoid computation errors.
%   - This function is designed to handle edge cases, but certain 
%     combinations of inputs may still lead to reduced accuracy.
%
%--------------------------------------------------
global im
im=0;
ierr=0;
maxvalx = 1e3;
maxvalp = 500;
maxvalq = 500;
dwarf=1e-100;
z=0.5*x*y;
if (1-y)<dwarf
  Bpq=1;
elseif y<dwarf
  Bpq=0;
elseif x<dwarf
  Bpq=betaincreg(y,p,q);
elseif x>maxvalx
  Bpq=0;
  ierr=2;
elseif p>maxvalp
  Bpq=0;
  ierr=2;
elseif q>maxvalq
  Bpq=1;
  ierr=2;
else
  y0=(x+2*p)/(x+2*p+2*q);  r=p+q; s=q/r; z=0.5*x*y;
  Ff=z+gammaln(p)-gammaln(p+q)+q*log(z);
  if r>160 && y>0.15
    if s>0.01 && p>50 && q>50
      if y<1.1*y0
        Bpq=BpqxyErFun(x,y,p,q);
      else
        Bpq=1-BpqxyErFunC(x,y,p,q);
      end
    elseif s<0.01 && z<15 && p>100
      if y<0.95
        Bpq=BpqxyLargep(x,y,p,q);
      else
        Bpq=Bpqxyyclose1(x,y,p,q);
      end
    elseif s<0.2 && x>100 && y>0.8
      if y<1.1*y0
        Bpq=BpqxyKummer(x,y,p,q);
      else
        Bpq=1-BpqxyKummerC(x,y,p,q);
      end
    else
      if y<1.1*y0
        [Bpq,ierro]=BpqxySeries(x,y,p,q);
        if ierro==1
          ierr=1;
        end
      else
        [betanc,ier]=BpqxySeriesC(x,y,p,q);
        if ier==1
          betanc=0;
          ierr=1;
        end
        Bpq=1-betanc;
      end
    end
  else
    if (z>250 && q<10 && p<10) 
    %Expansion for z large
      if y<0.95
        Bpq=BpqxyExpz(x,y,p,q);
      else
        Bpq=Bpqxyzlargeyclose1(x,y,p,q);
      end
    else
      if x<200 || Ff>log(1e+240)
        if y<y0
          [Bpq,ierro]=BpqxySeries(x,y,p,q);
          if ierro==1
            ierr=1;
          end
        else
          [betanc,ier]=BpqxySeriesC(x,y,p,q);
          if ier==1
            betanc=0;
            ierr=1;
          end
          Bpq=1-betanc;
        end
      else
        if y<0.95*y0
          [Bpq,ier]=BpqxyKummer(x,y,p,q);
          if ier==1
            Bpq=0;
            ierr=1;
          end
        else
          [BpqC,ier]=BpqxyKummerC(x,y,p,q);
          if ier==1
            BpqC=0;
            ierr=1;
          end
          Bpq=1-BpqC;
        end
      end
    end
    end
    if Bpq>1
      Bpq=1;
      ierr=1;
    end
  end
end



