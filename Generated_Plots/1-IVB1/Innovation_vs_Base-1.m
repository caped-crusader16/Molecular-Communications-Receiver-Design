clc;
clear;
close all;

snrdb =1:1:60;
snr = 10.^(snrdb/10);

lambda = 100;
r = 45e-9;
d2 = 5e-7;
D2 = 4.265e-10;
rs=4e-8;
rr=1e-5;
n=1250;
D1=7.94e-11;
d1=1e-6;
r0=d1+rr;

delta_T = 9e-6;
T = 30*delta_T;
L = 5;

global sum_Cj1;
global sum_Cj2;
sum1=0;
sum11=0;
sum2=0;
sum12=0;

x = 0:1:60;
sum_Cj1=0;
sum_Cj2=0;
Co = zeros(1,length(snrdb));
data = rand(1,L)>0.5;

for ii=1:1:5
    sum1 = sum1 + Cj_fun1(ii);
    sum11 = sum11 + data(ii)*Cj_fun1(ii);
    
    sum2 = sum2 + Cj_fun2(ii);
    sum12 = sum12 + data(ii)*Cj_fun2(ii);
end

avg1 = sum1 + lambda*T;
avg11 = avg1 + 20;

avg2 = sum2 + lambda*T;
avg12 = avg2 + 20;
final_term1 = 100000;
final_term2 = 100000;

for ii=1:1:length(snrdb)
    Co(ii) = snr(ii)*2*lambda*T;
    a1 = sec_term(Co(ii), avg1);
    a2 = sec_term(Co(ii), avg2);
    Pe1(ii)= 0.5 - 0.5*(a1);   % PLOT FOR SUB OPTIMAL APPROACH
    Pe2(ii)= 0.5 - 0.5*(a2);   % PLOT FOR SUB OPTIMAL APPROACH
    
    for jj=1:1:length(x)      % finding min tao for the optimal approach for a given Co 
        lam11 = lambda*T + sum11;
        lam21 = lam11 + Co(ii);
        first_term1 = Q_fun(lam11, x(jj));
        last_term1 = 1 - Q_fun(lam21, x(jj));
        final_term_next1 = 0.5*(first_term1 + last_term1);

        lam12 = lambda*T + sum12;
        lam22 = lam12 + Co(ii);
        first_term2 = Q_fun(lam12, x(jj));
        last_term2 = 1 - Q_fun(lam22, x(jj));
        final_term_next2 = 0.5*(first_term2 + last_term2);
        
        if(final_term_next1 < final_term1)
            final_term1 = final_term_next1;
            final_tao1 = x(jj);
        end
        
        if(final_term_next2 < final_term2)
            final_term2 = final_term_next2;
            final_tao2 = x(jj);
        end
    end
    
    b1 = sec_term1(Co(ii), avg1, final_tao1);
    b2 = sec_term1(Co(ii), avg2, final_tao2);
    Pe11(ii) = 0.5 - 0.5*(b1);
    Pe12(ii) = 0.5 - 0.5*(b2);
end

x=figure;
plot(snrdb, Pe1, '*-')
hold on
plot(snrdb, Pe11, 'o-')
hold on
plot(snrdb, Pe2, '*--')
hold on
plot(snrdb, Pe12, 'o--')
hold off
grid on
xlabel('SNR (db)')
ylabel('BER')
legend('Optimal zero bit mem. Rx. BER - Inno', 'Equiprobability zer bit mem. Rx. BER - Inno', 'Optimal zero bit mem. Rx. BER - Base', 'Equiprobability zer bit mem. Rx. BER - Base')
savefig(x,'Innovation_vs_Base-1.fig')
%---- functions ------
function ans1 = sec_term(Co1, avg1)
    t3 = Co1/avg1;
    t2 = log(1+t3); 
    t1 = Co1/t2;  % tao
    
    ans1 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=0:1:t1
        tt1 = lg0^(jj) * exp(-lg0);
        tt2 = lg1^(jj) * exp(-lg1);
        ans1 = ans1 + (tt1-tt2)/factorial(jj);
    end
end

function ans2 = sec_term1(Co1, avg1, final_tao)
    ans2 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=0:1:final_tao
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans2 = ans2 + (d1-d2)/factorial(jj);
    end
end

function sum_Cj1 = Cj_fun1(j)
    Ntx = 10^2;
    global sum_Cj1;
    sum_Cj1 = sum_Cj1 + Ntx*prob_j1(j);
    %Ntx*prob_j(j)
end

function sum_Cj2 = Cj_fun2(j)
    Ntx = 10^2;
    global sum_Cj2;
    sum_Cj2 = sum_Cj2 + Ntx*prob_j2(j);
    %Ntx*prob_j(j)
end

function y = prob_j1(j)
    rs=4e-8;
    rr=1e-5;
    n=1250;
    D1=7.94e-11;
    d1=1e-6;
    r0=d1+rr;
    lambda = 100;
    delta_T = 9e-6;
    T = 30*delta_T;
    L = 5;
    
    x = (rr/r0)*((rs*n) / (rs*n + pi*rr));
    b=r0-rr;
    beta = (n*rs + pi*rr)/(pi*rr*rr);
    t12 = D1*beta*beta*T;
    t11 = b*beta;
    t13 = 2*D1*beta*beta*T;
    t141 = sqrt(4*D1*T*j);
    t142 = sqrt(4*D1*T*(j+1));
    t1 = exp(t11 + t12*j)*erfc((b + t13*j)/t141);
    t2 = exp(t11 + t12*(j+1))*erfc((b + t13*(j+1))/t142);
    y = x*(t1 - t2);  
end

function y = prob_j2(j)
    lambda = 100;
    r = 45e-9;
    d2 = 500e-9;
    D2 = 4.265e-10;
    delta_T = 9e-6;
    T = 30*delta_T;
    L = 5;
    
   x = r/d2;
   t1 = ((d2-r)/(sqrt(4*D2*(j+1)*T)));
   t2 = ((d2-r)/(sqrt(4*D2*(j)*T)));
   y = x*(erfc(t1) - erfc(t2));  
end

function ans = Q_fun(lambda, n)
    u1 = factorial(n);
    u3 = lambda.^n;
    u2 = exp(-lambda)*u3;
    ans = u2/u1;
end