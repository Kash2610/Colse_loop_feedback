close all
clear all
clc

Nx = 2000;
M = 301;
L = 16;
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
% figure; plot(x)
% figure; plot(h)
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
% figure;plot(out)

%% With Overlap
% position = 0;
% N = 2^(ceil(log2(L + M)));
% z = zeros(1, N);
% H = fft(h, N);
% out = [];
% while position + L <= Nx 
%     X = fft(x(position+(1:L)), N);
%     z = [z(L+1:end), zeros(1, L)] + ifft(X.*H, N);
%     out = [out, z(1:L)];
%     position = position + L/2;  
% end
% figure; 
% plot(y_actual);hold on; plot(out)
% figure;plot(out)

%% Check overlap
% 
position = 0;
win = hamming(L, 'periodic')';
hop = L/2;
z = zeros(1, N);
prev_recording = zeros(1, L);
prev_out = zeros(1, L);
H = fft(h, N);
COLA_const = sum(win(1:L/2) + win(L/2+1:end))/(L/2);
out = zeros(1, L/2);
out_1 = [];
while position + L <= Nx
   rec = x(position + (1:L)).*win;
   recording = prev_recording(L/2+1:end) + rec(1:L/2);
   prev_recording = rec;
   out_1 = [out_1, recording/COLA_const];
   X = fft(rec, N);
   z = [z(L+1:end), zeros(1, L)] + ifft(X.*H, N);
   overlapped_out = prev_out(L/2+1:end) + z(1:L/2);
   out = [out, overlapped_out/COLA_const];
   prev_out = out(end-L+1:end);
   position = position + L/2; 
end

figure;plot(y_actual*5)  
hold on; plot(out)
% 
% 
% 
% %% Windowing
% h1 = hamming(256);
