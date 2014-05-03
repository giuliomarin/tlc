%% Project Channel Codes and Capacity
%  LDPC - IEEE 802.11n
%

% Giulio Marin
%
% giulio.marin@me.com
% 2013/06/12

clear;
close all;
tic

%% ****************************** Parameter *******************************

SNR = 0; %linspace(-2,0,5);   % SNR values in dB
var_noise = 10.^(-SNR/10);
MAX_iter = 50;  % max number of iteration per message
blocks = 100;   % number of words to decode (depends on BER required)
MAX_R = 1000;

% ---------------- IEEE 802.11n LDPC matrix definition --------------------

Z = 27;
matrix_prototype = ...
    [0 -1 -1 -1  0  0 -1 -1  0 -1 -1  0  1  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    22  0 -1 -1 17 -1  0  0 12 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    6 -1  0 -1 10 -1 -1 -1 24 -1  0 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1 -1;
    2 -1 -1  0 20 -1 -1 -1 25  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1 -1;
    23 -1 -1 -1  3 -1 -1 -1  0 -1  9 11 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1 -1;
    24 -1 23  1 17 -1  3 -1 10 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1 -1;
    25 -1 -1 -1  8 -1 -1 -1  7 18 -1 -1  0 -1 -1 -1 -1 -1  0  0 -1 -1 -1 -1;
    13 24 -1 -1  0 -1  8 -1  6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1 -1;
    7 20 -1 16 22 10 -1 -1 23 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1 -1;
    11 -1 -1 -1 19 -1 -1 -1 13 -1  3 17 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -1;
    25 -1  8 -1 23 18 -1 14  9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0;
    3 -1 -1 -1 16 -1 -1  2 25  5 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0];

H = zeros(size(matrix_prototype)*Z);
P0 = eye(Z);

for r=1:size(matrix_prototype,1)
    for c=1:size(matrix_prototype,2)
        shift = matrix_prototype(r,c);
        if (shift > -1)
            Pi = circshift(P0,[0 shift]);
        else
            Pi = zeros(Z);
        end
        R = (r-1)*Z+1:r*Z;
        C = (c-1)*Z+1:c*Z;
        H(R,C) = Pi;
    end
end

[m, n] = size(H);
% -------------------------------------------------------------------------

%% ****************************** Decoder *********************************

% Based on paper: High-Throughput LDPC Decoders, M. M. Mansour and N. R.
% Shanbhag, 6 December 2003.
% (Algorithm 1 - Memory aware decoder)

MAX_j = max(sum(H,2));  % max degree for check nodes

Q = zeros(1,n);       % total bit-check messages
P = zeros(1,n);       % temporary buffer for total bit-check messages
R = zeros(m,MAX_j);   % check-bit messages
R_1 = zeros(size(R)); % previous check -bit messages

BER = zeros(size(SNR));
PER = zeros(size(SNR));

for snr = 1:length(SNR)
    
    v = 0;
    
    Q = zeros(1,n);       % total bit-check messages
    P = zeros(1,n);       % temporary buffer for total bit-check messages
    R = zeros(m,MAX_j);   % check-bit messages
    R_1 = zeros(size(R)); % previous check -bit messages
    
    bit_errors = 0;
    pkt_errors = 0;
    
    for b = 1:blocks
        
        if mod(b,blocks/10) == 0
            fprintf('%d%%\n', b/blocks * 100)
        end
        
        % Create a random word
        %(use random or all 0 signal to increase speed)
        
        % u =randi([0 1],1,m);
        % c = LDPC_generator(H,u);
        u = zeros(1,m);
        c = zeros(1,n);
        
        s = 2*c-1; % PAM signal
        r = s + randn(size(s))*sqrt(var_noise(snr)); % received signal
        
        % Prior log-likelihood
        Lci = (-2*r./var_noise(snr));
        P = Lci;
        Q = Lci;
        
        k = 0;
        while k < MAX_iter
            
            for i = 1:m % for every check node
                Vi = find(H(i,:)); % Set of incident variables
                z = ones(length(Vi))-eye(length(Vi));
                Rij = Q(Vi)-R_1(i,1:length(Vi)); % Incoming messages
                % avoid NaN values
                Rij(abs(Rij)<1e-8) = 1e-8;
                R(i,1:length(Vi)) = -log(tanh(z*(-log(tanh(abs(Rij)/2)))'/2)).*prod(sign(Rij)).*sign(Rij(1:length(Rij)))'; % outgoing message
                % avoid overflow
                R(i,abs(R(i,:)) > MAX_R) = sign(R(i,abs(R(i,:)) > MAX_R))*MAX_R;
                P(Vi) = P(Vi) + R(i,1:length(Vi));
            end
            R_1 = R;
            Q = P;
            P = Lci;
            v = Q<0; % current word estimation
            
            if ~sum(mod(H*v',2)) % valid code
                break
            end
            
            k = k+1;
        end
        
        errors = sum(u~=v(1:m));
        
        if errors
            bit_errors = bit_errors + errors;
            pkt_errors = pkt_errors + 1;
        end
    end
    
    BER(snr) = bit_errors/(m*blocks);
    PER(snr) = pkt_errors/blocks;
    
    if ~exist('./data', 'dir')
        mkdir('./data');
    end
    
    save(['./data/data' num2str(SNR(snr)) '.mat'])
    fprintf('SNR: %g, BER: %g\n',SNR(snr),BER(snr))
    
end

save('./data/BER.mat','BER','PER','SNR')

toc