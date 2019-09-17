clear all
m=28*1.67e-27;
g=44.9;
a=0.272e-9;
n_total=10000;
n_pml=4500;
n_domain=1000;
sigma_max=1.68e12;


r=4;                                        %mass ratio
M=[m*ones(5000,1); r*m*ones(5000,1)];

kna=linspace(0.0001,pi,100);

sigma=zeros(n_total,1);

for n=1:1:n_total                          
    if n<n_pml+1
        sigma(n)=(sigma_max*(n_pml+1-n)^2)/n_pml^2;
    elseif n>=n_total-n_pml
        sigma(n)=(sigma_max*(n-(n_total-n_pml))^2)/n_pml^2;
    end 
end
x=[-5000:1:0];


for j=1:numel(kna)

    U_inc= [exp(1i*kna(j).*x(:)); zeros(5000,1)];
    w = sqrt(2*g/m*(1-cos(kna(j))));              %frequency
    A = [0;0;-g./M(2:n_total-1)]  ;                                                            %super diagonal element
    b = [1;(-w^2+2*g./M(2:n_total-1)) + sigma(2:n_total-1).^2- 2*sigma(2:n_total-1)*1i*w;1];   %diagonal element                                                              
    c = [-g./M(2:n_total-1);0];                                                                %sub diagonal elements                                              
    d = [0;(w^2 -2*g./M(2:n_total-1)).*U_inc(2:n_total-1)+(g./M(2:n_total-1)).*(U_inc(3:n_total)+U_inc(1:n_total-2));0];                       
    U_scat=TDMAsolver(A,b,c,d);
    u_scat=U_scat';
%     U_scat=gallery('tridiag',c,b,A);   %Could not make it work
%     u_scat=U_scat\d;

     
   
    R= M(2:n_pml).*sigma(2:n_pml).*(abs(u_scat(2:n_pml))).^2.*w^2;     %Reflected energy
    T= M(n_total-n_pml+1:n_total).*sigma(n_total-n_pml+1:n_total).*(abs(u_scat(n_total-n_pml+1:n_total))).^2.*w^2; %Transmitted energy

    Re= sum(R);
    Tr= sum(T);
       
Te(j)= Tr/(Re+Tr);   %Reflection coefficient
end

plot(kna,Te)
hold on
% end