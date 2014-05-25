function data = update_plot( data, channel )
%UPDATE_PLOT( data, channel )
%
% Compute water filling and update the plot. Use variables inside data, if
% channel input variable is present, generate new channel.
%
% Reference:
% "MIMO-OFDM Wireless Communications with MATLAB"
% Yong S. Cho, Jaekwon Kim, Won Y. Yang, Chung G. Kang (Korea University)
% Wiley-IEEE Press, 2011

% Giulio Marin
%
% giulio.marin@me.com
% 2013/03/20

N    = 64;        % FFT size
Ptot = data.Ptot; % Maximum power
SNR  = data.SNR;  % sigma_a^2/sigma_w^2
ax1  = data.ax1; 
ax2  = data.ax2;

%% Channel
% h is a random channel of length [2-5] taps and unitary power

if nargin>1
    n = randi([2,5],1);
    h = (rand(1,n)+1i*rand(1,n))/sqrt(2);
    h = h/sqrt(h*h');
    data.h = h;
else
    h = data.h;
    n = length(h);
end

axes(ax2)
stem(0:n-1,abs(h),'^','LineWidth',2)
title('Channel response');
ylabel('h(n)');
xlabel('n')
set(gca,'xTick',0:n-1)
axis([-0.5 n-0.5 0 1])

H = fft(h, N);
Gamma = abs(H).^2*SNR;

axes(ax1);
bar(1./Gamma,'r')
title('Inverse of the frequency channel response');
ylabel('1/Gamma(k)');
xlabel('k')
axis([0 N 0 10])

%% Water filling

Psum  = 0;
c     = 0;      % 1/lambda
delta = 0.01;   % increment of 1/lambda
P = zeros(size(Gamma));

% Start from 1/lambda = 0 and increase the value until the maximum
% available power is reached. Force to 0 the subcarriers where the
% algorithm would have given a negative value.
while Psum < Ptot
    c = c + delta;
    P = c - (1./Gamma);
    P(P<0) = 0;
    Psum = sum(P);
    
end

hold on
plot([0 N],[c c],'--r')
S = 1./Gamma;
S(P==0) = 0;
bar(S,'FaceColor','g')
hold off

end

