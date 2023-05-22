close all
clear all

Nx = 2000;
M = 301;
L = 768;
x = rand(1, 2000);
h = rand(1, M);
y = zeros(1, Nx + M - 1);
N = M + L - 1;
y_actual = conv(x, h);

position = 0;
while position + L <= Nx
    y (position + (1:N)) = y(position + (1:N)) + conv(x(position + (1:L)), h);
    position = position + L;
end
figure; plot(x)
figure; plot(h)
figure; plot(y_actual)
hold on; plot(y)

%% Real-time
position = 0;
z = zeros(1, N);
out = [];
while position + L <= Nx 
    z = [z(L+1:end), zeros(1, L)] + conv(x(position + (1:L)), h);
    out = [out, z(1:L)];
    position = position + L;  
end
figure; 
plot(y_actual);hold on; plot(out)

%% Real time with circular convolution
position = 0;
N = 2^(ceil(log2(L + M)));
z = zeros(1, N);
H = fft(h, N);
out = [];
while position + L <= Nx 
    X = fft(x(position+(1:L)), N);
    z = [z(L+1:end), zeros(1, L)] + ifft(X.*H, N);
    out = [out, z(1:L)];
    position = position + L;  
end
figure; 
plot(y_actual);hold on; plot(out)
figure;plot(out)
