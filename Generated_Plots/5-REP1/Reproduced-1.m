clear;
close all;
clc;

% Basic values
lambda_0 = 100;                     % 1/s
r = 45*(10^-9);                     % radius of reciever (m)
d = 5*(10^-7);                    % distance (m)
D = 4.265*(10^-10);                 % Co-efficient of Diffusion (m*m/s)
del_t = 9*(10^-6);                          % Discrete Time length (s)
T1 = 30*del_t;               % Slot length (s)
T2 = 50*del_t;
L=5;                                % Channel Length (m)
snr_db=30;

tau=[20:120]; 
snr_lin=1000;               % 10^(30/10)
C0_1=snr_db*2*lambda_0*T1;
C0_2=snr_db*2*lambda_0*T2;
P_i_0_1 = (r/d)*(erfc((d-r)/sqrt(4*D*T1)));
P_i_0_2 = (r/d)*(erfc((d-r)/sqrt(4*D*T2)));

ntx1 = C0_1/P_i_0_1;
ntx2=C0_2/P_i_0_2;

for j=1:L
    p_1(j)=(r/d)*(erfc((d-r)/sqrt(4*D*(j+1)*T1))-erfc((d-r)/sqrt(4*D*j*T1)));
    p_2(j)=(r/d)*(erfc((d-r)/sqrt(4*D*(j+1)*T2))-erfc((d-r)/sqrt(4*D*j*T2)));
    c_1(j)=p_1(j)*ntx1;
    c_2(j)=p_2(j)*ntx2;
end

sum_sc_1 = zeros(1, 2^L);
sum_sc_2 = zeros(1, 2^L);

for ii=0:((2^L) - 1)
    s=dec2bin(ii, L);
    for j = 1:L
        sum_sc_1(ii+1)=sum_sc_1(ii+1)+ (s(j)*c_1(j));
        sum_sc_2(ii+1)=sum_sc_2(ii+1)+ (s(j)*c_2(j));
    end
end

for i=1:length(tau)
    ber_sum_1=0;
    ber_sum_2=0;
    for ii=1:(2^L)
        ber_sum_1 = ber_sum_1 + (0.5*(1 + gammainc((lambda_0*T1 + sum_sc_1(ii)), ceil(tau(i))) - gammainc((lambda_0*T1 + sum_sc_1(ii) + C0_1), ceil(tau(i)))));
        ber_sum_2 = ber_sum_2 + (0.5*(1 + gammainc((lambda_0*T2 + sum_sc_2(ii)), ceil(tau(i))) - gammainc((lambda_0*T2 + sum_sc_2(ii) + C0_2), ceil(tau(i)))));
    end
    ber_1(i)=(ber_sum_1/(2^L));
    ber_2(i)=(ber_sum_2/(2^L));
end

xyz=figure;
plot(tau, ber_1, 'b*-', tau, ber_2, 'ro-')
legend('slot length 30T', 'slot length 50T')
xlabel('threshold')
ylabel('BER')
grid on;
savefig(xyz, 'Reproduced-1.fig')