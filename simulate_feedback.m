% Simulate the feedback
close all
clear all
clc;

%%
[x1, fs] = audioread('.\examples\ex1.wav');
[x2, ~] = audioread('.\examples\ex2.wav');
[x3, ~] = audioread('.\examples\ex3.wav');
[x4, ~] = audioread('.\examples\ex4.wav');
[x5, fs] = audioread('.\examples\ex5.wav');
x = vertcat(x1);
x = x';
% x = ones(1, length(x))*1e-5;
%% Feedback path
fileID = fopen("bcoff_pixel.txt");
A = fscanf(fileID, '%f');
h = A(1533:end)';
figure;plot(h)
frame_processing_time = 0.008;
dt = fs*frame_processing_time;
% h = [zeros(1, dt), h];
M = length(h);

%% Frame settings
Nx = length(x);
frame_size = 0.016;
L = fs*frame_size;
N = 2^(ceil(log2(L + M)));
G = 1;
position = 0;
z = zeros(1, N);
H = fft(h, N);
out = [];
y = [];
while position + L <= Nx 
    recording = x(position + (1:L)) + z(1:L);
    y = [y, recording];
    %Process
    recording = G*recording;
    X = fft(recording, N);
    z = [z(L+1:end), zeros(1, L)] + ifft(X.*H, N);
    out = [out, z(1:L)];  
    position = position + L;  
end
y1 = y;
out1 = out;
energy = sum((x - y).^2);
%%
%% With overlap
Nx = length(x);
frame_size = 0.016;
L = fs*frame_size;
win = hamming(L, 'periodic')';
hop = L/2;
N = 2^(ceil(log2(L + M)));
G = 1;
position = 0;
z = zeros(1, N);
H = fft(h, N);
out = zeros(1, L/2);
y = [];
z = zeros(1, N);
prev_recording = zeros(1, L);
prev_out = zeros(1, L);  
COLA_const = sum(win(1:L/2) + win(L/2+1:end))/(L/2);

while position + L <= Nx
    recording = x(position + (1:L)) + prev_out;
    recording = recording.*win;
    overlapped_recording = prev_recording(L/2+1:end) + recording(1:L/2);
    y = [y, overlapped_recording/COLA_const];
    prev_recording = recording;
    
    %Process
    recording = G*recording;
    X = fft(recording, N);
    z = [z(L+1:end), zeros(1, L)] + ifft(X.*H, N);
    overlapped_out = prev_out(L/2+1:end) + z(1:L/2);
    out = [out, overlapped_out/COLA_const];
    prev_out = out(end-L+1:end);  
    position = position + L/2;  
end
energy_in = sum((x(L/2+1:Nx-L/2) - y(L/2+1:end)).^2);
energy_out = sum((out1(N:end-N) - out(N:end-N)).^2);

%% 
energy_feedback = sum(y1(L/2:end-L/2)-y(L/2:end)).^2;


%% Spectogram
figure;
spectrogram(out/norm(out), hamming(256), 240, 1024*2, fs,'yaxis','reassigned');title('y');

