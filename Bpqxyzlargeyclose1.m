function Bpq=Bpqxyzlargeyclose1(x,y,p,q)
ymn1=1-y;
eta=ymn1/y;
z=0.5*x*y;
X=z*eta;
Q=gammainc(X,q,'upper');
Phi=exp(X)*eta^(-q)*Q;
suma=Phi;
k=0;
term=suma;
while abs(term/suma)>1.e-15
  Phi=eta*Phi-z^(q-k-1)/gamma(q-k);
  k=k+1;
  ck=(-1)^k*pochhammer(1-p-q,k)/factorial(k);
  term=ck*Phi;
  suma=suma+term;
end
Bpq=exp(-X)*y^(p-1)*ymn1^q*suma;
end




