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
    if exist('u2227951_lab_Audio.mat', 'file') == 2 ... % Update with your student ID      
        load('u2227951_lab_Audio.mat', 'audioRaw'); % Update with your student ID
        
        % 1a)
        Q1.audioInput = audioRaw;

        % 1b)
        N = length(Q1.audioInput);
        Fs = 22050;
        t = [0:N-1]' / Fs; 
        n1 = 350; % Hz
        Q1.noiseSignal = sin(2*pi*n1*t);
        
        Q1.audioNoisy = Q1.audioInput + Q1.noiseSignal;
        
        Q1.FFTNoisy = fft(Q1.audioNoisy);

        f = (0:N-1) * (Fs / N);          % Hz
        mag = abs(Q1.FFTNoisy);          % magnitude (no scaling)
        
        figure;
        plot(f, mag);
        xlim([0 Fs]);                    % visible range [0, 22050)
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('Magnitude spectrum of Q1.audioNoisy');
        
        % (v) Label the largest frequency component (peak)
        % If you want to ignore DC (0 Hz), use 2:N instead of 1:N.
        [~, idx] = max(mag(2:end));
        idx = idx + 1;                   % because we searched from 2:end
        peakFreq = f(idx);
        peakMag  = mag(idx);
        
        hold on;
        plot(peakFreq, peakMag, 'o');    % marker on the peak
        text(peakFreq, peakMag, sprintf('  Peak: %.1f Hz', peakFreq), ...
            'VerticalAlignment', 'bottom');
        hold off;           
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
    if exist('u<ID>_lab_signals.mat', 'file') == 2 % Update with your student ID
        load('u<ID>_lab_signals.mat', 'P') % Update with your student ID
        
        





    end

end

