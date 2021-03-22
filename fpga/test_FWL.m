clear;clc;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
addpath('sim/');
N = 16;
global M;
M = 8;
R = M/N;
n = log2(N);

info_bits_logical = logical([0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1]);

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

N_sim = 1000;
N_err = 0;

for sim_iter = 1:N_sim
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
    A = 3.5;                      % Amplification factor.

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
    fwrite(app.Com_Obj, 16);  % N = 8, trigger FPGA decoder.
    
    FPGA_done = false;
    for i = 1:N
        fwrite(app.Com_Obj, bit_stream(i));
    end
    while ~FPGA_done
        % do nothing but waiting for the FPGA.
    end
    % FPGA_done is true now.
    
    if any(decode_result ~= random_bits)
        N_err = N_err + 1;
    end
    if mod(sim_iter, N_sim/40) == 1
        fprintf('Complete: %.2f%%\n', sim_iter/N_sim*100);
    end
end

fclose(app.Com_Obj);
fprintf('FPGA connection closed.\n');
fprintf('BLER = %f\n', N_err/N_sim);


%% Callback functions for FPGA-UART communication.
function ComReceiveCallback(app, src, event)
global decode_result;
global FPGA_done;
global M;
	Com_DataReceive = fread(app,1); % Read 1 byte from FPGA.
    decode_result = [];
    for k = 1:M
        if mod(Com_DataReceive, 2) == 1
            decode_result = [1, decode_result];
        else
            decode_result = [0, decode_result];
        end
        Com_DataReceive = floor(Com_DataReceive/2);
    end
    decode_result = flip(decode_result);
    FPGA_done = true;
end
