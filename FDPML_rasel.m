clear all
m=28*1.67e-27;
g=44.9;
a=0.272e-9;
n_total=10000;
n_pml=4500;
n_domain=n_total-2*n_pml;
sigma_max=1.68e12;


r=.5;                                        %mass ratio
M=[m*ones(n_total/2,1); r*m*ones(n_total/2,1)];

kna=linspace(0.01,pi,100);

sigma=zeros(n_total,1);

for n=1:1:n_total                          
    if n<n_pml+1
        sigma(n)=(sigma_max*(n_pml+1-n)^2)/n_pml^2;
    elseif n>=n_total-n_pml
        sigma(n)=(sigma_max*(n-(n_total-n_pml))^2)/n_pml^2;
    end 
end
x=[-n_total/2:1:-1];


for j=1:numel(kna)

    U_inc= [exp(1i*kna(j).*x(:)); zeros(n_total/2,1)];
    w = sqrt(2*g/m*(1-cos(kna(j))));              %frequency
    A = [0;-g./(M(2:n_total-1)*w^2)]  ;                                                            %super diagonal element
    b = [1;(-1+2*g./(M(2:n_total-1)*w^2)) + sigma(2:n_total-1).^2/w^2- 2*sigma(2:n_total-1)*1i/w;1];   %diagonal element                                                              
    c = [-g./(M(2:n_total-1)*w^2);0];                                                                %sub diagonal elements                                              
    d = [0;(1 -2*g./(M(2:n_total-1)*w^2)).*U_inc(2:n_total-1)+(g./(M(2:n_total-1)*w^2)).*(U_inc(3:n_total)+U_inc(1:n_total-2));0];                       
%     U_scat=TDMAsolver(A,b,c,d);
%     u_scat=U_scat';
    U=gallery('tridiag',c,b,A);   
    u_scat=U\d;

     
   
    R= M(1:n_pml).*sigma(1:n_pml).*(abs(u_scat(1:n_pml))).^2.*w^2;     %Reflected energy
    T= M(n_total-n_pml+1:n_total).*sigma(n_total-n_pml+1:n_total).*(abs(u_scat(n_total-n_pml+1:n_total))).^2.*w^2; %Transmitted energy

    Re= sum(R);
    Tr= sum(T);
       
Te(j)= Tr/(Re+Tr);   %Reflection coefficient
end

plot(kna,Te)

% end