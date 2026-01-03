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
        
        % 1a) Loading Recording
        Q1.audioInput = audioRaw;
        
        % 1b) Noise Contamination and FFT Analysis
        % Define parameters for noise signal
        N = length(Q1.audioInput);
        n1 = 350; % From u2227951_lab.txt
        f_s = 22050; % Sampling frequency (used to record audio)
        t = (0:N-1)' / f_s; % N time points in seconds
        
        % Create the noise signal
        Q1.noiseSignal = sin(2*pi*n1*t); 
    
        % Contaminate the recorded signal with the noise signal
        Q1.audioNoisy = Q1.audioInput + Q1.noiseSignal; 
    
        % FFT analysis of the contaminated signal
        Q1.FFTNoisy = fft(Q1.audioNoisy); 
        
        f_k = (0:N-1) * (f_s / N); % Compute frequency axis. Each bin k corresponds to frequency k * (f_s / N)
        mag = abs(Q1.FFTNoisy); 
        
        % Plot magnitude spectrum with range [0, 22050) Hz
        figure;
        plot(f_k, mag);
        xlim([0 f_s]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('Magnitude Spectrum of Q1.audioNoisy');
        grid on;
        
        % Find and label the largest peak (excluding DC)
        [mag_peak, peak_index] = max(mag(2:end));
        peak_index = peak_index + 1; % Compensate for skipping DC component
        f_peak = f_k(peak_index);
        saveas(gcf, 'images/1b_magnitude_spectrum.png');
        
        % Label the peak on the plot
        hold on;
        plot(f_peak, mag_peak, 'ro');
        text(f_peak, mag_peak, sprintf('  %.1f Hz', f_peak));
        hold off;
        saveas(gcf, 'images/1b_magnitude_spectrum_labelled.png');
    
        % 1c) FIR Filter design to isolate the largest noise frequency
        % Define filter parameters
        f_nyq = f_s / 2; % Nyquist frequency needed for filter

        % Need to filter out 350 Hz so I tested lower and upper bounds around that to get these values 
        w_cutoff_lower = 310;
        w_cutoff_upper = 380;
        w_n = [w_cutoff_lower w_cutoff_upper] / f_nyq; % Normalised stopband edges 
        
        % Need to ensure the transition is sharp enough to filter out the 350 Hz component
        filter_order = 2100;

        % Blackman window best suits our needs
        win = blackman(filter_order + 1);
        Q1.h = fir1(filter_order, w_n, 'stop', win);
        Q1.filteredAudio = filter(Q1.h, 1, Q1.audioNoisy);
        %soundsc(Q1.filteredAudio, f_s);
        
        % Plot frequency response
        figure;
        freqz(Q1.h, 1, N, f_s);
        title('FIR Bandstop Filter Frequency Response');
        
        % Access the magnitude subplot (top) to add annotations
        subplot(2,1,1);
        ylim auto;
        hold on;
        xline(w_cutoff_lower/1000, '--r', sprintf('f_{c (lower)} = %d Hz', w_cutoff_lower), 'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'FontSize', 10);
        xline(w_cutoff_upper/1000, '--r', sprintf('f_{c (upper)} = %d Hz', w_cutoff_upper), 'LabelOrientation', 'aligned', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'middle', 'FontSize', 10);
        hold off;
        
        % Add legend or annotation for stopband region
        ylabel('Magnitude (dB)');
        title('FIR Bandstop Filter Frequency Response');
        
        saveas(gcf, 'images/1c_freqz.png');
        
        % Plot noisy and filtered signals for comparison
        figure;
        subplot(2,1,1);
        plot(t, Q1.audioNoisy);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Noisy Audio Signal');
        
        subplot(2,1,2);
        plot(t, Q1.filteredAudio);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Filtered Audio Signal');
        saveas(gcf, 'images/1c_contaminated_vs_filtered.png');
    
        % 1d) Apply an Effect of Your Choice (Echo)
        % Define echo parameters
        delay_time = 0.3; % Delay in seconds (300 ms for audible echo)
        decay = 0.5; % Decay factor for each echo (0 to 1)
        num_echoes = 3; % Number of echo repetitions
        
        % Convert delay time to samples
        delay_samples = round(delay_time * f_s);
        
        % Initialise output with original signal
        Q1.audioCustom = Q1.audioInput;
        
        % Apply multiple echoes
        % y[n] = x[n] + a*x[n-D] + a^2*x[n-2D] + a^3*x[n-3D] + ...
        for echo = 1:num_echoes
            echo_delay = echo * delay_samples;
            echo_amplitude = decay^echo;
            
            for n = (echo_delay + 1):N
                Q1.audioCustom(n) = Q1.audioCustom(n) + echo_amplitude * Q1.audioInput(n - echo_delay);
            end
        end
        
        % Normalise output to prevent clipping
        Q1.audioCustom = Q1.audioCustom / max(abs(Q1.audioCustom));
        
        % Plot original and echo signals for comparison
        figure;
        subplot(2,1,1);
        plot(t, Q1.audioInput);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Original Audio Signal');
        
        subplot(2,1,2);
        plot(t, Q1.audioCustom);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Audio with Echo Effect');

        saveas(gcf, 'images/1d_effect.png')
        %soundsc(Q1.audioCustom, f_s);
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
    f_s = 650; % Sampling frequency (Hz)
    f_c = 65; % Passband edge frequency (Hz)
    f_sb = 91; % Stopband edge frequency (Hz)
    r_p = 0.2; % Passband ripple (dB)
    a_s = 60; % Stopband attenuation (dB)
    
    % Normalise frequencies (divide by Nyquist frequency)
    f_nyq = f_s / 2;
    w_p = f_c / f_nyq; % Normalised passband frequency
    w_s = f_sb / f_nyq; % Normalised stopband frequency
    
    % 2a) Calculate minimum filter order and design filter
    % Using Chebyshev Type I
    [n, w_n] = cheb1ord(w_p, w_s, r_p, a_s);
    [Q2.b, Q2.a] = cheby1(n, r_p, w_n, 'low');
    Q2.filterOrder = n;
    
    % 2b) Plot frequency response
    figure;
    freqz(Q2.b, Q2.a, 1024, f_s);
    
    % Magnitude response
    subplot(2,1,1);
    title('Chebyshev Type I Lowpass Filter - Magnitude Response');
    grid on;

    % Lines to mark cutoffs and stopband attenuation
    hold on;
    xline(f_c, 'g--');
    xline(f_sb, 'r--');
    yline(-a_s, 'r:');
    
    % Add text labels instead of inline labels
    text(f_c+5, -10, 'f_c = 65 Hz', 'Color', 'g');
    text(f_sb+5, -30, 'f_{sb} = 91 Hz', 'Color', 'r');
    text(200, -a_s+5, 'Stopband Attenuation = 60 dB', 'Color', 'r');
    hold off;

    % Phase response
    subplot(2,1,2);
    title('Chebyshev Type I Lowpass Filter - Phase Response');
    grid on;

    saveas(gcf, 'images/2b_freqz.png')
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
        T_s = 0.045;
        n = length(P);
        t = (0:n-1)' * T_s;
    
        % 3a) Construct the Linear Model
        % Model: P(t) = A*(1-exp(-0.14*t)) + B*exp(-0.23*t)*cos(14.74*t) + C*sin(51.97*t)
        basis_a = 1 - exp(-0.14 * t);
        basis_b = exp(-0.23 * t) .* cos(14.74 * t);
        basis_c = sin(51.97 * t);
        Q3.Obs = [basis_a, basis_b, basis_c];
    
        % 3b) Parameter Estimation
        Q3.param = Q3.Obs \ P;
    
        % 3c) Pressure Prediction
        Q3.yHat = Q3.Obs * Q3.param;
    
        % 3d) Model Error Calculation
        Q3.mse = mean((P - Q3.yHat).^2);
    
        % 3e) Frequency Domain Analysis
        % FFT of noisy pressure data
        Q3.yFFT = fft(P);
    
        % FFT of predicted pressure data
        Q3.yHatFFT = fft(Q3.yHat);
    
        % Frequency vector
        f_s = 1 / T_s;
        Q3.fRange = (0:n-1)' * (f_s / n);
    
        % Plot magnitude spectra
        figure;
        plot(Q3.fRange, abs(Q3.yFFT), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Noisy Data (P)');
        hold on;
        plot(Q3.fRange, abs(Q3.yHatFFT), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Predicted (yHat)');
        hold off;
    
        xlim([0 f_s]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('FFT Magnitude Spectra of Noisy vs Predicted Pressure');
        legend('Location', 'best');
        grid on;
    
        % Label key spectral peaks
        half_n = floor(n/2);
        mag = abs(Q3.yFFT(1:half_n));
        [mag_peak, peak_index] = findpeaks(mag, 'MinPeakHeight', max(mag)*0.1);
        
        hold on;
        for i = 1:length(peak_index)
            f_peak = Q3.fRange(peak_index(i));
            plot(f_peak, mag_peak(i), 'ko', 'MarkerSize', 8);
            text(f_peak, mag_peak(i), sprintf('  %.2f Hz', f_peak));
        end
        hold off;
        
        % Save FFT plot
        saveas(gcf, 'images/3e_fft.png');
        
        % Time domain comparison plot for analysis
        figure;
        plot(t, P, 'b-', 'LineWidth', 1, 'DisplayName', 'Noisy Data (P)');
        hold on;
        plot(t, Q3.yHat, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Predicted (yHat)');
        hold off;
    
        xlabel('Time (s)');
        ylabel('Pressure (bar)');
        title('Time Domain: Noisy vs Predicted Pressure');
        legend('Location', 'best');
        grid on;
    end
end

