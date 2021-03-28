function BLER = sim_CA_PCC_SCL_fpga(Ebn0, min_errors, amp_factor)
%% Code parameters are hard-wired inside the FPGA, thus cannot be altered.
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

N = 32;
global M;
M = 16;         % Overall end-to-end Coding rate R = 1/2.
R = M/N;
n = log2(N);

K_CRC = 4;
K_PCC = 2;

[~, ~, g] = get_crc_objective(K_CRC);
[G_crc, ~] = crc_generator_matrix(g, M);
CRC_Qmatrix = G_crc(:, M+1 : end);              % Generator matrix for CRC bits.; G=[I,Q].

% Fixed info_bits.
info_bits_logical = logical([0   0   0   0   0   0   0   1   0   0   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]).';
assert(sum(info_bits_logical) == M + K_CRC + K_PCC, 'info_bits_logical mismatch!');

% Setup PCC config structure, used in PCC-polar encoder.
PCC_conf.info_bits_cnt = M + K_CRC;
PCC_conf.parity_bits_cnt = K_PCC;
PCC_conf.N = N;
PCC_conf.nonfrozen_bits_logical = info_bits_logical;
PCC_conf.parity_bits_index = {[0,2,5]+1, [1,6,10]+1};

%% Add BPSK modulation and calculate the Finite-Word-Length results.
LLR_WIDTH = 8;                                      % In fact, the width in FPGA is LLR_WIDTH+ADDITIONAL.
sigma = 1/sqrt(2*R)*(10^(-Ebn0/20));


%% Control Serial Port.
delete(instrfindall);
app.Com_Obj = serial('COM4','BaudRate', 115200);
set(app.Com_Obj,'BytesAvailableFcnMode','byte');
set(app.Com_Obj,'BytesAvailableFcnCount', 1);

% Setup UART-Receive Callback functions.
app.ComReceiveCallback = @ComReceiveCallback;
app.Com_Obj.BytesAvailableFcn=@app.ComReceiveCallback;  

if(app.Com_Obj.Status == "closed")
	fopen(app.Com_Obj);
	if(app.Com_Obj.Status == "open")
        fprintf(1,'Serial port successfully opened\n'); 
    else
        fprintf(1,'Serial port open failed\n');
	end  
else
	fprintf(1,'Serial port has been used by other processes.\n'); 
end 

%% Simulation Loop.
global decode_result;
global FPGA_done;

if ~exist('min_errors', 'var') || isempty(min_errors)
    min_errors = 50;
end
cnt_errors = 0;
N_runs = 0;

if ~exist('amp_factor', 'var') || isempty(amp_factor)
    amp_factor = 3.5;
end
A = amp_factor;         % LLR Amplification factor in order to increase the accuracy of LLRs in FPGA.
fprintf('Estimating BLER @ Eb/n0=%.2f dB for CRC-PC-Polar SCLF decoder on FPGA.\n', Ebn0);

while cnt_errors < min_errors
    random_bits = (rand([1,M])>0.5);
    CRC_aided_bits = [random_bits, logical(mod(random_bits*CRC_Qmatrix, 2))];  % row vector.
    x = PCC_polar_encoder(CRC_aided_bits, PCC_conf);
    
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
    fwrite(app.Com_Obj, 0);     % Initialize
    fwrite(app.Com_Obj, N);     % Trigger FPGA decoder.
    
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
% fprintf('BLER = %f\n', cnt_errors/N_runs);
BLER = cnt_errors / N_runs;

end


%% Callback functions for FPGA-UART communication.
function ComReceiveCallback(app, ~, ~)
global decode_result;
global FPGA_done;
global M;
    persistent cnt_bytes;
    N_bytes = M/8;
	Com_DataReceive = fread(app,1); % Read 1 byte from FPGA.
    for k = 1:8
        if mod(Com_DataReceive, 2) == 1
            decode_result = [decode_result, 1];
        else
            decode_result = [decode_result, 0];
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
        FPGA_done = true;
    end
    
    return;
end