close all;
clear all;
clc;

%Assigning parameters from simulation table
snrdb = (-5:35);
snrlin = 10.^(snrdb./10);
r = 45*10^(-9);
d = 500*10^(-9);
D = 4.265*10^(-10);
L =5;
lambda0 = 100;
T = 30*9*10^(-6);
C0 = 2*lambda0*T*10.^(snrdb./10);
K=1;
m=0;
n=0;
P0 = p(r,d,D,T,0);
NTx = 100.* ones(1,length(snrdb))
%NTx = 2*lambda0*T*10.^(snrdb./10)./P0;
P = zeros(1,L-K);
%Vector of P for 1 to L=5 time slots
for jj=1:(L-K)
    tmp = p(r,d,D,T,jj+K);
    P(jj) = tmp;

end 
%For each value of snr first finding op_tou and then using op_tou finding
%BER
for kk=1:length(snrdb)
    ber1(kk)=0;
tou = (0:120);
%Finding Tou minimum for particular SNR
for i= 1:length(tou)
    ber(i) = 0;
     c = NTx(kk)*P;   %For given snr finding C for i=1 to L=5 time slots
     s_t = (0:2^(L-K)-1);    %permutations of 5 bits slot
     s = de2bi(s_t);  %Permutation vector
     ISI = s*c';      %ISI for given SNR 32x1 vector
 
     %Equation number 14 for each tou
     for k=1:length(ISI)
      pe(k) =(0.5*(gammainc(lambda0*T + ISI(k),ceil(tou(i)))) + 1 - gammainc(lambda0*T + C0(kk) + ISI(k) ,ceil(tou(i))));
      ber(i) = ber(i) + pe(k);   %Adding BER of each permutation 
    end
  ber(i) = (1/2^(L-K))*ber(i);
end
[tm,min_in] = min(ber);  %Tou with minimum BER
op_tou = tou(min_in);     %Tou
  s_1 = (0:31);
  s_1_or = de2bi(s_1);
  P1 = [p(r,d,D,T,1) p(r,d,D,T,2) p(r,d,D,T,3) p(r,d,D,T,4) p(r,d,D,T,5)];
  c1 = NTx(kk)*P1;
  ISI1 = s_1_or*c1';
  m1=0.5;
  n1=0.5;
for k=1:length(ISI)
    sel=1;
    sel2=1;
    for j=1:K
     if(s_1_or(k,j)==1)
       sel = sel*n1;
       sel2 = sel2*(1-n1);
    else
        sel = sel*(1-m1);
        sel2 =sel2*(m1);
     end
    end
    ans2 = gammainc(lambda0*T + ISI(k),op_tou);
    ans1= 1 - gammainc(lambda0*T + C0(kk) + ISI(k) ,op_tou)
     % pe(k) =(0.5*( + 1 - gammainc(lambda0*T + C0(kk) + ISI(k) ,op_tou));
      
    m= (ans2*sel + ans2*sel2);
    n = (ans1*sel + ans1*sel2);
   m1= m;
   n1 =n;
    ber1(kk) = (m+n)/2;
      %ber1(kk) = ber1(kk) + pe(k); 
    end
 %ber1(kk) = (1/2^(L-K))*ber1(kk);

end


sum=0;
for ii=1:1:5
    sum = sum +c1(ii);
    
end
avg = sum + lambda0*T;
avg1 = avg + 20;
final_term = 100000;
Co = zeros(1,length(snrdb));

for ii=1:1:length(snrdb)
    %Co(ii) = snr(ii)*2*lambda0*T;
    a = sec_term(C0(ii), avg);
    Pe(ii)= 0.5 - 0.5*(a);   % PLOT FOR SUB OPTIMAL APPROACH
   
end

semilogy(snrdb,Pe,'-.b*');
hold on
semilogy(snrdb,ber1 ,'-.ro');   %Plot of ber



axis([-5 35 10^(-6) 1]);
grid on
legend('Two bit equi probable','Two bit Optimal');
xlabel('SNRdB')
ylabel('BER')

function ans1 = sec_term(Co1, avg1)
    t3 = Co1/avg1;
    t2 = log(1+t3); 
    t1 = Co1/t2;  % tao
    
    ans1 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=0:1:t1
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans1 = ans1 + (d1-d2)/factorial(jj);
    end
end
function p = p(r,d,D,T,i)
p = (r/d).*(erfc((d-r)./sqrt(4*D*(i+1)*T)) - erfc((d-r)./sqrt(4*D*i*T)));

end