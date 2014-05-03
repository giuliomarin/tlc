function c = LDPC_generator(h,u)

% Based on paper: Richardson, T.J.; Urbanke, R.L., "Efficient encoding of
% low-density parity-check codes," Information Theory, IEEE Transactions
% on , vol.47, no.2, pp.638,656, Feb 2001

% Giulio Marin
%
% giulio.marin@me.com
% 2013/06/12

% define sizes of message and code
mlen = size(h,1);
clen = size(h,2);

% rename variables to match those used by the author
n = clen;
m = clen - mlen;

% *******************************************************
% Ttable iii on the paper for details

hrow1 = h(:,end);

% find the 'gap' length
for i=1:clen
    if hrow1(i) == 1
        g = i;
        break;
    end
end
g = mlen - g;

% width and last index of submatrices a, b
wa = n-m;
wb = g;
ea = wa;
eb = wa + wb;

% extract the submatrices a, b, c, d, e and t
a = h(1:m-g,1:ea);
b = h(1:m-g,ea+1:eb);
t = h(1:m-g,eb+1:end);
c = h(m-g+1:end,1:ea);
d = h(m-g+1:end,ea+1:eb);
e = h(m-g+1:end,eb+1:end);

% calculate p1 and p2
invt = (inv(t));
et1 = -(e*invt);

% todo: check singularity
% iup = diag(ones(1,size(et1,2)),0);
% idn = diag(ones(1,size(et1,1)),0);
% x = [iup zeros(size(iup,1),size(idn,2)); et1 idn];
% y = x*h;
% spy(y);

phi = et1*b + d;
xtra = et1*a + c;
p1 = mod(phi*xtra*(u'),2)';
p2 = mod(invt*(a*(u') + b*(p1')),2)';
c = [u p1 p2];

% checking c*h'=0;
zero = mod(c*h',2);
if sum(zero)== 0
    % disp('ok!')
else
    disp('error')
    return
end
end