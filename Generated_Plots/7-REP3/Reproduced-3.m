clc;
close all;
clear;

% Basic values
lambda_0 = 100;                     % 1/s
r = 45*(10^-9);                     % radius of reciever (m)
d = 5*(10^-7);                    % distance (m)
D = 4.265*(10^-10);                 % Co-efficient of Diffusion (m*m/s)
del_t = 9*(10^-6);                          % Discrete Time length (s)
T = 30*del_t;               % Slot length (s)
L=5;                                % Channel Length (m)
snr_db = 25;
snr_lin=10.^(snr_db./10);

global sum_Cj;
sum=0;
sum1=0;
rec_par = [0:60];
sum_Cj=0;

data = rand(1,L)>0.5;

for ii=1:L
    sum = sum + Cj_fun(ii);
    sum1 = sum1 + data(ii)*Cj_fun(ii);
end

avg = sum + lambda_0*T;
avg1 = avg + 20;
final_term = 100000;

for jj=1:1:length(rec_par)
    factor = avg.^(rec_par(jj));
    fact = factorial(rec_par(jj));
    final(jj) = (exp(-avg)*(factor))/fact;
 
    factor1 = avg1.^(rec_par(jj));
    final1(jj) = (exp(-avg1)*(factor1))/fact;
    
    lam1 = lambda_0*T + sum1;
    lam2 = lam1 + 20;
    first_term = Q_fun(lam1, rec_par(jj));
    last_term = 1 - Q_fun(lam2, rec_par(jj));
    final_term_next = 0.5*(first_term + last_term);
    
    if(final_term_next < final_term)
        final_term = final_term_next;
        final_tao = rec_par(jj);
    end
end

t3 = 20/avg;
t2 = log(1+t3); 
t1 = 20/t2;

xyz=figure;
plot(rec_par,final,'r-')
hold on
plot(rec_par,final1,'k-')
stem(t1, 0.08, 'Marker', 'o');
stem(final_tao, 0.08, 'Marker', 'o')
grid on
legend('App. distr. when si=0', 'App. distr. when si=1', 'threshold from equiprobability', 'optimal threshold'); 
savefig(xyz,'Reproduce-3')
%---- functions ------

function sum_Cj = Cj_fun(j)
    Ntx = 10^2;
    global sum_Cj;
    sum_Cj = sum_Cj + Ntx*prob_j(j);
    %Ntx*prob_j(j)
end

function y = prob_j(j)
    lambda_0 = 100;
    r = 45e-9;
    d = 500e-9;
    D = 4.265e-10;
    delta_T = 9e-6;
    T = 30*delta_T;
    L = 5;
    
   x = r/d;
   t1 = ((d-r)/(sqrt(4*D*(j+1)*T)));
   t2 = ((d-r)/(sqrt(4*D*(j)*T)));
   y = x*(erfc(t1) - erfc(t2));  
end

function ans = Q_fun(lambda, n)
    u1 = factorial(n);
    u3 = lambda.^n
    u2 = exp(-lambda)*u3;
    ans = u2/u1;
end