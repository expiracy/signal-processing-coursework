function Answers = u2227951_lab()
%% ES3C5 lab submission template v2.0
%
% Please DO NOT change this header
% Please DO NOT change any code in this top-level function
% Please DO NOT change function or subfunction arguments (Input OR Output)
% DO Change function name above to u<ID>_lab(), where <ID> is your student #

% ES3C5 2025-2026 Lab Assignment
% Module Leader: Viji Ahanathapillai
% Modify the SUBFUNCTIONS below with the code needed to determine or
% demonstrate the answers requested.

% See the Briefing Sheet for full instructions.

% Initialise answer structure
Answers = [];

%% Template call (dummy)
Q0();

%% Remaining Calls
Answers.Q1 = Q1Fun();
Answers.Q2 = Q2Fun();
Answers.Q3 = Q3Fun();


end

%% Template Question on hypotenuse length
% This subfunction is a sample to demonstrate what is expected for comments
function c0 = Q0()
% Please DO NOT change function arguments (input OR output)
% Assign answer to c0 (double value)

% Define triangle lengths
a0 = 2; % 1st side
b0 = 1; % 2nd side

% Find length of hypothenuse
c0 = sqrt(a0^2 + b0^2); % Pythagorean theorem to find 3rd side

end

%%
%%%%%%%%%%%%%%%%%%% Start Modifying Below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Q1 Audio Signal Processing
% Figures created must be labelled and copied to u<ID>_lab.docx
function Q1 = Q1Fun()
% Please DO NOT change function arguments above (input OR output)
% You must call a plotting function but you can label it manually.

    % Initialising outputs to guarantee they exist (unless you delete them)
    % Please DO NOT modify or move the code below
  

    Q1.audioInput =[];% Q1 part (a)
    Q1.noiseSignal=[];% Q1 part (b)
    Q1.audioNoisy = []; % Q1 part (b)
    Q1.FFTNoisy = []; % Q1 part (b)

    Q1.h = []; % Q1 part (c)
    Q1.filteredAudio = []; % Q1 part (c)
    Q1.audioCustom = []; % Q1 part (d)

    
    % Please DO NOT modify or move the code above
    
    %% Start your Q1 code here
    % DO NOT REMOVE OR MOVE THIS IF STATEMENT
    % WRITE YOUR CODE INSIDE THIS IF STATEMENT
    if exist('u2227951_lab_Audio.mat', 'file') == 2
        load('u2227951_lab_Audio.mat', 'audioRaw');
        
        % 1a) Load audio
        Q1.audioInput = audioRaw;
        
        % 1b) Generate noise and contaminate signal
        N = length(Q1.audioInput);
        Fs = 22050;
        t = (0:N-1)' / Fs;
        n1 = 350;
        
        Q1.noiseSignal = sin(2*pi*n1*t);
        Q1.audioNoisy = Q1.audioInput + Q1.noiseSignal;
        Q1.FFTNoisy = fft(Q1.audioNoisy);
        
        % Compute frequency axis and magnitude
        f = (0:N-1) * (Fs / N);
        mag = abs(Q1.FFTNoisy);
        
        % Plot magnitude spectrum with range [0, 22050) Hz
        figure;
        plot(f, mag);
        xlim([0 Fs]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('Magnitude spectrum of Q1.audioNoisy');
        
        % Find and label the largest peak (excluding DC)
        [peakMag, idx] = max(mag(2:end));
        idx = idx + 1;
        peakFreq = f(idx);
        
        hold on;
        plot(peakFreq, peakMag, 'ro', 'MarkerSize', 10);
        text(peakFreq, peakMag, sprintf('  %.1f Hz', peakFreq), 'VerticalAlignment', 'bottom');
        hold off;

        % 1c) FIR Filter Design
        % Band-stop filter to suppress 350 Hz noise
        % Stopband: 300-400 Hz to remove the 350 Hz noise component
        % Passbands: 0-300 Hz and 400-11025 Hz with unity gain
        
        Fnyq = Fs / 2;
        Wn = [300 400] / Fnyq;
        filterOrder = 1000;
        
        Q1.h = fir1(filterOrder, Wn, 'stop');
        
        % Plot frequency response
        figure;
        freqz(Q1.h, 1, 4096, Fs);
        title('FIR Band-Stop Filter Frequency Response');
        saveas(gcf, 'Q1c_filter_response.png');
        
        % 1c v) Apply filter
        Q1.filteredAudio = filter(Q1.h, 1, Q1.audioNoisy);
        
        % 1c vi) Plot noisy and filtered signals
        figure;
        t_plot = (0:N-1) / Fs;
        
        subplot(2,1,1);
        plot(t_plot, Q1.audioNoisy);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Noisy Audio Signal');
        
        subplot(2,1,2);
        plot(t_plot, Q1.filteredAudio);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Filtered Audio Signal');

        % 1d) Flanger Effect
        % Parameters
        maxDelay = 0.003;       % Maximum delay in seconds (3 ms)
        rate = 0.5;             % LFO frequency in Hz (speed of sweep)
        depth = 0.002;          % Depth of modulation in seconds (2 ms)
        mix = 0.7;              % Mix level (0 to 1)
        
        % Convert to samples
        maxDelaySamp = round(maxDelay * Fs);
        depthSamp = round(depth * Fs);
        
        % Create output array (same length as input)
        Q1.audioCustom = zeros(size(Q1.audioInput));
        
        % Generate LFO (low frequency oscillator)
        t = (0:N-1)' / Fs;
        lfo = depthSamp * sin(2 * pi * rate * t);
        
        % Apply flanger effect
        for n = 1:N
            % Calculate current delay (must be positive integer)
            currentDelay = maxDelaySamp + round(lfo(n));
            
            % Ensure we don't read before the start of the signal
            if (n - currentDelay) > 0
                Q1.audioCustom(n) = Q1.audioInput(n) + mix * Q1.audioInput(n - currentDelay);
            else
                Q1.audioCustom(n) = Q1.audioInput(n);
            end
        end
        
        % Normalise to prevent clipping
        Q1.audioCustom = Q1.audioCustom / max(abs(Q1.audioCustom));
        
        % Plot comparison
        figure;
        t_plot = (0:N-1) / Fs;
        
        subplot(2,1,1);
        plot(t_plot, Q1.audioInput);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Original Audio Signal');
        
        subplot(2,1,2);
        plot(t_plot, Q1.audioCustom);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Audio with Flanger Effect');
        
        saveas(gcf, 'Q1d_flanger_effect.png');
        
        % Play the effect
        disp('Playing original...');
        soundsc(Q1.audioInput, Fs);
        pause(6);
        disp('Playing with flanger...');
        soundsc(Q1.audioCustom, Fs);
    end

end

%% Q2 Filter Design
% Figures created must be labelled and copied to u<ID>_lab.docx
function Q2 = Q2Fun()
% Please DO NOT change function arguments above (input OR output)
% You must call a plotting function but you can label it manually.

    % Initialising outputs to guarantee they exist (unless you delete them)
    % Please DO NOT modify code below
    
    Q2.filterOrder = []; % Q2 part (a)
    Q2.a = []; % Q2 part (a)
    Q2.b = []; % Q2 part (a)
    
    % Please DO NOT modify code above
    
    %% Start your Q2 code here
    % DO NOT REMOVE OR MOVE THIS IF STATEMENT
    % WRITE YOUR CODE INSIDE THIS IF STATEMENT

    % Given specifications
    Fs = 650;           % Sampling frequency (Hz)
    fc = 65;            % Passband edge frequency (Hz)
    fsb = 91;           % Stopband edge frequency (Hz)
    Rp = 0.2;           % Passband ripple (dB)
    As = 60;            % Stopband attenuation (dB)
    
    % Normalise frequencies (divide by Nyquist frequency)
    Fnyq = Fs / 2;
    Wp = fc / Fnyq;     % Normalised passband frequency
    Ws = fsb / Fnyq;    % Normalised stopband frequency
    
    % Q2(a) Calculate minimum filter order and design filter
    % Using Chebyshev Type I
    [n, Wn] = cheb1ord(Wp, Ws, Rp, As);
    [Q2.b, Q2.a] = cheby1(n, Rp, Wn, 'low');
    Q2.filterOrder = n;
    
    % Display results
    fprintf('Filter Order: %d\n', Q2.filterOrder);
    fprintf('Normalised cutoff frequency: %.4f\n', Wn);
    
    % Q2(b) Plot frequency response
    figure;
    [H, f] = freqz(Q2.b, Q2.a, 1024, Fs);
    
    % Magnitude response
    subplot(2,1,1);
    plot(f, 20*log10(abs(H)), 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Chebyshev Type I Lowpass Filter - Magnitude Response');
    grid on;
    hold on;
    
    % Add labels for passband and stopband
    xline(fc, 'g--', 'LineWidth', 1.5);
    xline(fsb, 'r--', 'LineWidth', 1.5);
    yline(-Rp, 'm:', 'LineWidth', 1.5);
    yline(-As, 'b:', 'LineWidth', 1.5);
    
    % Add text labels instead of inline labels
    text(fc+5, -10, 'f_c = 65 Hz', 'Color', 'g');
    text(fsb+5, -30, 'f_{sb} = 91 Hz', 'Color', 'r');
    text(200, -Rp-5, 'Passband Ripple = -0.2 dB', 'Color', 'm');
    text(200, -As+5, 'Stopband Atten. = -60 dB', 'Color', 'b');
    
    hold off;
    ylim([-80 5]);
    xlim([0 Fs/2]);
    
    % Phase response
    subplot(2,1,2);
    plot(f, angle(H) * 180/pi, 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    title('Chebyshev Type I Lowpass Filter - Phase Response');
    grid on;
    xlim([0 Fs/2]);
    
    saveas(gcf, 'Q2_filter_response.png');
       
end

%% Q3 Estimating Unknown Signal
function Q3 = Q3Fun()
% Please DO NOT change function arguments (input OR output)

    % Initialising outputs to guarantee they exist (unless you delete them)
    % Please DO NOT modify code below
    Q3.Obs = []; % Q3 part (a)

    Q3.param = []; % Q3 part (b)

    Q3.yHat = []; % Q3 part (c)

    Q3.mse = []; % Q3 part (d)

    Q3.yFFT = []; % Q3 part (e)
    Q3.fRange = []; % Q3 part (e)
    Q3.yHatFFT = []; % Q3 part (e)
    % Please DO NOT modify code above

    % DO NOT REMOVE THIS IF STATEMENT
    % WRITE YOUR CODE INSIDE THIS IF STATEMENT
    if exist('u2227951_lab_signals.mat', 'file') == 2
        load('u2227951_lab_signals.mat', 'P');
        
        % Sampling parameters
        Ts = 0.045;
        N = length(P);
        t = (0:N-1)' * Ts;
        
        % (a) Construct the observation matrix
        % Model: P(t) = A*(1-exp(-0.14*t)) + B*exp(-0.23*t)*cos(14.74*t) + C*sin(51.97*t)
        % Theta = [basis_A, basis_B, basis_C]
        
        basis_A = 1 - exp(-0.14 * t);
        basis_B = exp(-0.23 * t) .* cos(14.74 * t);
        basis_C = sin(51.97 * t);
        
        Q3.Obs = [basis_A, basis_B, basis_C];
        
        % (b) Parameter estimation using least squares
        % P = Theta * param + noise
        % param = (Theta' * Theta)^(-1) * Theta' * P
        Q3.param = Q3.Obs \ P;
        
        % Display estimated parameters
        fprintf('Estimated parameters:\n');
        fprintf('A = %.4f\n', Q3.param(1));
        fprintf('B = %.4f\n', Q3.param(2));
        fprintf('C = %.4f\n', Q3.param(3));
        
        % (c) Pressure prediction
        Q3.yHat = Q3.Obs * Q3.param;
        
        % (d) Mean squared error
        Q3.mse = mean((P - Q3.yHat).^2);
        fprintf('Mean Squared Error: %.6f\n', Q3.mse);
        
        % (e) Frequency domain analysis
        % i) FFT of noisy pressure data
        Q3.yFFT = fft(P);
        
        % ii) FFT of predicted pressure data
        Q3.yHatFFT = fft(Q3.yHat);
        
        % iii) Frequency vector
        Fs = 1 / Ts;
        Q3.fRange = (0:N-1)' * (Fs / N);
        
        % iv) Plot magnitude spectra
        figure;
        plot(Q3.fRange, abs(Q3.yFFT), 'b', 'DisplayName', 'Noisy Data (P)');
        hold on;
        plot(Q3.fRange, abs(Q3.yHatFFT), 'r--', 'DisplayName', 'Predicted (yHat)');
        hold off;
        
        xlim([0 Fs]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('FFT Magnitude Spectra: Noisy vs Predicted Pressure');
        legend('Location', 'best');
        grid on;
        
        % Label key spectral peaks
        % Find peaks in the predicted signal (cleaner)
        halfN = floor(N/2);
        [pks, locs] = findpeaks(abs(Q3.yHatFFT(1:halfN)), 'MinPeakHeight', max(abs(Q3.yHatFFT))*0.1);
        
        hold on;
        for i = 1:length(locs)
            peakFreq = Q3.fRange(locs(i));
            peakMag = pks(i);
            plot(peakFreq, peakMag, 'ko', 'MarkerSize', 8);
            text(peakFreq, peakMag*1.1, sprintf('%.2f Hz', peakFreq), 'HorizontalAlignment', 'center');
        end
        hold off;
        
        saveas(gcf, 'Q3_fft_comparison.png');
        
        % Optional: Time domain comparison plot
        figure;
        plot(t, P, 'b', 'DisplayName', 'Noisy Data (P)');
        hold on;
        plot(t, Q3.yHat, 'r', 'LineWidth', 1.5, 'DisplayName', 'Predicted (yHat)');
        hold off;
        xlabel('Time (s)');
        ylabel('Pressure (bar)');
        title('Time Domain: Noisy vs Predicted Pressure');
        legend('Location', 'best');
        grid on;
        
        saveas(gcf, 'Q3_time_comparison.png');
    end

end

