clear;clc;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
addpath('sim/');
N = 32;
global M;
M = 16;
R = M/N;
n = log2(N);

info_bits_logical = logical([false;false;false;false;false;false;false;false;false;false;false;true;false;true;true;true;false;false;false;true;false;true;true;true;true;true;true;true;true;true;true;true]).';

%% Add BPSK modulation and calculate the Finite-Word-Length results.
LLR_WIDTH = 8;                          % In fact, the width in FPGA is LLR_WIDTH+n.
Ebn0 = 3.5;                             % in dB.
sigma = 1/sqrt(2*R)*(10^(-Ebn0/20));


%% Control Serial Port.
delete(instrfindall);
app.Com_Obj = serial('COM4','BaudRate', 115200);
set(app.Com_Obj,'BytesAvailableFcnMode','byte');
set(app.Com_Obj,'BytesAvailableFcnCount', 1);

app.ComReceiveCallback = @ComReceiveCallback;
app.Com_Obj.BytesAvailableFcn=@app.ComReceiveCallback;
if(app.Com_Obj.Status == "closed")
	fopen(app.Com_Obj);
	if(app.Com_Obj.Status == "open")
        fprintf(1,'Serial port successfully opened\n'); 
    else
        fprintf(1,'Serial prot open failed\n');
	end  
else
	fprintf(1,'Serial port has been used by other processes.\n'); 
end 

%% Simulation Loop.
global decode_result;
global FPGA_done;

N_err = 0;
min_errors = 50;
cnt_errors = 0;
N_runs = 0;
A = 2.5;        % Amplification factor.

while cnt_errors < min_errors
    random_bits = (rand([1,M])>0.5);
    u = zeros(1,N);
    u(info_bits_logical) = random_bits;
    x = my_polar_encoder(u,n);              % x is the encoded bits.
    
    x_bpsk = 1-2*x;
    
    y = x_bpsk + sigma * randn([1, N]);
    llr = 2*y/(sigma^2);

    % Map llr into finite-word-length numbers of LLR_WIDTH bits.
    % Ensure there is no overflow
    bit_stream = zeros(1, N);
    
    for i = 1:N
        amp = llr(i)*A;
        if amp > 127
            amp=127;
        elseif amp < -128
            amp=-128;
        end

        t = int8(amp);  % Range: -128 ~ 127. Convert it into 0~255 representation.
        if t<0
            bit_stream(i) = double(t) + 256;
        else
            bit_stream(i) = double(t);
        end
    end
    
    % Initialize FPGA.
    fwrite(app.Com_Obj, 0);  % Initialize
    fwrite(app.Com_Obj, N);  % N = 8, trigger FPGA decoder.
    
    FPGA_done = false;
    decode_result = [];         % Initialize FPGA globals.
    for i = 1:N
        fwrite(app.Com_Obj, bit_stream(i));
    end
    while ~FPGA_done
        % do nothing but waiting for the FPGA.
    end
    % FPGA_done is true now.
    
    if any(decode_result ~= random_bits)
        cnt_errors = cnt_errors + 1;
    end
    if mod(N_runs, min_errors) == 1
        fprintf('Complete: %.2f%%\n', cnt_errors/min_errors*100);
    end
    N_runs = N_runs + 1;
end

fclose(app.Com_Obj);
fprintf('FPGA connection closed.\n');
fprintf('BLER = %f\n', cnt_errors/N_runs);


%% Callback functions for FPGA-UART communication.
function ComReceiveCallback(app, src, event)
global decode_result;
global FPGA_done;
global M;
    persistent cnt_bytes;
    N_bytes = M/8;
	Com_DataReceive = fread(app,1); % Read 1 byte from FPGA.
    for k = 1:8
        if mod(Com_DataReceive, 2) == 1
            decode_result = [1, decode_result];
        else
            decode_result = [0, decode_result];
        end
        Com_DataReceive = floor(Com_DataReceive/2);
    end
    
    if isempty(cnt_bytes)
        cnt_bytes = 0;
    end
    
    if cnt_bytes < (N_bytes-1)
        cnt_bytes = cnt_bytes + 1;
    else
        cnt_bytes = 0;
        decode_result = flip(decode_result);
        FPGA_done = true;
    end
    
    return;
end
