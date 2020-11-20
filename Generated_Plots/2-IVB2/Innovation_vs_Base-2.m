clc;
clear;
close all;

% Basic values
rs=4e-8;
rr=1e-5;
n=1250;
D1=7.94e-11;
d1=1e-6;
r0=d1+rr;

lambda_0 = 100;                     % 1/s
r = 45*(10^-9);                     % radius of reciever (m)
d2 = 5*(10^-7);                    % distance (m)
D2 = 4.265*(10^-10);                 % Co-efficient of Diffusion (m*m/s)
del_t = 9*(10^-6);                          % Discrete Time length (s)
T = 30*del_t;               % Slot length (s)
L=5;                                % Channel Length (m)
snr_db = 25;
snr_lin=10.^(snr_db./10);

global sum_Cj1;
global sum_Cj2;
sum1=0;
sum11=0;
sum2=0;
sum12=0;
rec_par = [0:60];
sum_Cj1=0;
sum_Cj2=0;

data = rand(1,L)>0.5;

for ii=1:L
    sum1 = sum1 + Cj_fun1(ii);
    sum2 = sum2 + Cj_fun2(ii);
    sum11 = sum11 + data(ii)*Cj_fun1(ii);
    sum12 = sum12 + data(ii)*Cj_fun2(ii);
end

avg1 = sum1 + lambda_0*T;
avg11 = avg1 + 20;

avg2 = sum2 + lambda_0*T;
avg12 = avg2 + 20;
final_term = 100000;

for jj=1:1:length(rec_par)
    factor1 = avg1.^(rec_par(jj));
    fact = factorial(rec_par(jj));
    final1(jj) = (exp(-avg1)*(factor1))/fact;
    
    
    factor2 = avg2.^(rec_par(jj));
    fact = factorial(rec_par(jj));
    final2(jj) = (exp(-avg2)*(factor2))/fact;
 
    factor11 = avg11.^(rec_par(jj));
    final11(jj) = (exp(-avg11)*(factor11))/fact;
    
    factor12 = avg12.^(rec_par(jj));
    final12(jj) = (exp(-avg12)*(factor12))/fact;
    
    lam11 = lambda_0*T + sum11;
    lam21 = lam11 + 20;
    first_term1 = Q_fun(lam11, rec_par(jj));
    last_term1 = 1 - Q_fun(lam21, rec_par(jj));
    final_term_next1 = 0.5*(first_term1 + last_term1);
    
    lam12 = lambda_0*T + sum12;
    lam22 = lam12 + 20;
    first_term2 = Q_fun(lam12, rec_par(jj));
    last_term2 = 1 - Q_fun(lam22, rec_par(jj));
    final_term_next2 = 0.5*(first_term2 + last_term2);
    
    
    if(final_term_next1 < final_term)
        final_term = final_term_next1;
        final_tao1 = rec_par(jj);
    end
    
    if(final_term_next2 < final_term)
        final_term = final_term_next2;
        final_tao2 = rec_par(jj);
    end
end

t31 = 20/avg1;
t21 = log(1+t31); 
t11 = 20/t21;

t32 = 20/avg2;
t22 = log(1+t32); 
t12 = 20/t22;
final_tao2=36;
xyz=figure;
plot(rec_par,final1,'r-')
hold on
plot(rec_par,final2,'r--')
hold on
plot(rec_par,final11,'k-')
hold on
plot(rec_par,final12,'k--')
stem(t11, 0.1, 'Marker', 'o');
stem(final_tao1, 0.1, 'Marker', 'o')
stem(t12, 0.1, 'Marker', 'o');
stem(final_tao2, 0.1, 'Marker', 'o')
grid on
legend('App. distr. when si=0 - Inno', 'App. distr. when si=0 - Base', 'App. distr. when si=1 - Inno', 'App. distr. when si=1 - Base', 'threshold from equiprobability - Inno', 'optimal threshold - Inno', 'threshold from equiprobability - Base', 'optimal threshold - Base'); 
xlabel('Received Number of Particles')
ylabel('Probability')
savefig(xyz, 'Innovation_vs_Base-2.fig');

%---- functions ------

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
 
    lambda_0 = 100;
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
    t31 = erf(-b/t141);
    t32 = erf(-b/t142);
    t1 = exp(t11 + t12*j)*erfc((b + t13*j)/t141);
    t2 = exp(t11 + t12*(j+1))*erfc((b + t13*(j+1))/t142);
    y = x*(t1 - t2 - t31 + t32); 
end

function y = prob_j2(j)
    lambda_0 = 100;
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