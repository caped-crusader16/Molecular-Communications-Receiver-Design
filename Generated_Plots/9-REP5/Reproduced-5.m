%close all;
%clear all;
%clc; 

snrdb =-5:1:55;
snr = 10.^(snrdb/10);

lambda = 100;
r = 45e-9;
d = 500e-9;
D = 4.265e-10;
delta_T = 9e-6;
T = 30*delta_T;
L = 5;
global sum_Cj;
sum=0;
sum1=0;
x =0:1:65;
sum_Cj=0;
Co = zeros(1,length(snrdb));
s=0;
display(snr);
temp=0;
data = rand(1,L)>0.5;
display(data);

data_est = zeros(1,L);
display(data_est);
cj = rand(1,L)>0.5;

for ii=1:1:5
    sum = sum + Cj_fun(ii);
    display(sum);
    sum1 = sum1 + data(ii)*Cj_fun(ii);
    display(sum1);
end
avg = sum + lambda*T;
display(avg);

avgg = sum1 + lambda*T;
display(avgg);
%avg1 = avg + 20;
i=1;
P_0 = (r/d)*(erfc((d-r)/sqrt(4*D*i*T))-erfc((d-r)/sqrt(4*D*(i-1)*T)));
Ntx = 2*lambda*T.*(10.^(snr./10))./P_0;

final_term = 100000;
for ii=1:1:length(snrdb)
    for i = 1:L
        P_i1(i) = (r/d)*(erfc((d-r)/sqrt(4*D*i*T))-erfc((d-r)/sqrt(4*D*(i-1)*T)));
    end
    for j = 1:L
        cj(j) = Ntx(ii).*P_i1(j);
    end
    display(cj);
    Co(ii) = snr(ii)*2*lambda*T;
    average = (lambda*T + sum1);
    average1 = average + data.*Co(ii);
    ri = poissrnd(average1);
    a = sec_term(Co(ii), avgg);
  
  
    Pe(ii)=0.5*a;
    lam1 = lambda*T + sum1;
    lam2 = lam1 + Co(ii);
     for jj=1:1:length(x)     
        
        first_term = (1-Q_fun(lam2, x(jj)));
        last_term = Q_fun(lam1, x(jj));
        final_term_next = 0.5*(first_term + last_term);
        if(final_term_next < final_term)
            final_term = final_term_next;
            final_tao = x(jj);
        end
    end
   if(ri(1)<=final_tao)
            data_est(1)=0;
        else
            data_est(1)=1;
        end
  
    for kk=2:1:L
        if(data(kk-1)==0) && (data_est(kk-1)==1)
            temp=last_term;
        elseif(data(kk-1)==0) && (data_est(kk-1)==0)
            temp=1-last_term;
        elseif(data(kk-1)==1) && (data_est(kk-1)==0)
            temp=first_term;
        else
            temp=1-first_term;
        end
        
        for jj=1:1:length(x)     
       
        first_term = first_term + (1-Q_fun(lam2, x(jj)))*temp;
        last_term = last_term + Q_fun(lam1, x(jj))*temp;
        final_term_next = 0.5*(first_term + last_term);
      
        if(final_term_next < final_term)
            final_term = final_term_next;
            final_tao = x(jj);
        end
        end
        
        if(ri(kk)<=final_tao)
            data_est(kk)=0;
        else
            data_est(kk)=1;
        end
    end
    
    
    b = sec_term1(Co(ii), avgg, final_tao);
    Pe1(ii)=0.5*b;   
end

xyz=figure;
plot(snrdb, Pe, 's-');
hold on
plot(snrdb, Pe1, 'o-');
hold on
plot(snrdb,s,'b-')
ylim([0 0.5])
savefig(xyz,'Reproduced-5.fig')


% functions 
function ans1 = sec_term(Co1, avg1)
    t3 = Co1/avg1;
    t2 = log(1+t3); 
    t1 = Co1/t2; 
    
    ans1 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=1:1:t1
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans1 = ans1 + (-d1+d2)/factorial(jj);
        %ans1=ans1+d1+d2;
    end
    ans1 = 1 + ans1;
end

function ans2 = sec_term1(Co1, avg1, final_tao)
    ans2 = 0;
    lg0 = avg1;
    lg1 = avg1 + Co1;
    
    for jj=1:1:final_tao
        d1 = lg0^(jj) * exp(-lg0);
        d2 = lg1^(jj) * exp(-lg1);
        ans2 = ans2 + (-d1+d2)/factorial(jj);
        %ans2=ans2+d1+d2;
    end
    ans2 = 1 + ans2;
end


function sum_Cj = Cj_fun(j)
    Ntx = 10^2;
    global sum_Cj;
    sum_Cj = sum_Cj + Ntx*prob_j(j);
    %Ntx*prob_j(j)
end

function y = prob_j(j)
    lambda = 100;
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
    ans=0;
    for j=1:1:n
         u1 = factorial(n);
         u3 = lambda.^n;
         u2 = exp(-lambda)*u3;
         ans = ans + u2/u1;
    end
end