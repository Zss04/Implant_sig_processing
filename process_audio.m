function process_audio(fileName)
    % Specify the directory containing sound files
    folder = 'Sound Files';
    fullFileName = fullfile(folder, fileName);
    
    % Step 3.1: Read files into Matlab and find their sampling rate.
    [audioData, sampleFrequency] = audioread(fullFileName);
    
    % Step 3.2: Check whether the sound is stereo or mono.
    [numSamples, numChannels] = size(audioData);
    if numChannels == 2
        audioData = sum(audioData, 2);
    end
    
    % Step 3.3: Play the sound in Matlab
    sound(audioData, sampleFrequency);
    
    % Step 3.4: Write the sound to a new file
    outputFileName = ['processed_' fileName];
    audiowrite(outputFileName, audioData, sampleFrequency);
    
    % Step 3.5: Plot the sound waveform as a function of sample number
    figure;
    plot(audioData);
    title(['Waveform of ' fileName]);
    xlabel('Sample Number');
    ylabel('Amplitude');
    
    % Step 3.6: Downsample the audio if it is more than 16kHz
    targetFs = 16000;
    if sampleFrequency > targetFs
        audioData = resample(audioData, targetFs, sampleFrequency);
        sampleFrequency = targetFs;
    elseif sampleFrequency < targetFs 
        error('The sampling rate is less than 16 kHz, please redo Step 3.1 with a higher sample rate file.');
    end

    pause(numSamples / sampleFrequency);
    
    % Step 3.7: Generate a cosine signal at 1 kHz
    duration = numSamples / sampleFrequency;
    t = (0:1/sampleFrequency:(numSamples-1)/sampleFrequency)';
    cosineSignal = cos(2 * pi * 1000 * t);

    % Play the cosine signal 
    % sound(cosineSignal, sampleFrequency);
    
    % % Commented Out: Plot two cycles of the cosine signal 
    % figure;
    % twoCycles = 2 * sampleFrequency / 1000;
    % plot(t(1:twoCycles), cosineSignal(1:twoCycles));
    % title('Cosine Wave');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    
    % --- Phase 2: Task 4: Bandpass Filter Bank ---
% Define the number of channels (frequency bands)
numChannels = 16;
lowerFreq = 100; % Lower cutoff frequency (Hz)
upperFreq = min(8000, sampleFrequency / 2); % Upper cutoff frequency limited to Nyquist

% Check if the frequency range is valid
if lowerFreq >= upperFreq
    error('Frequency range is invalid. Ensure that the Nyquist frequency is higher than 100 Hz.');
end

% Calculate the frequency bands (logarithmic spacing for better perception)
bandEdges = logspace(log10(lowerFreq), log10(upperFreq), numChannels + 1);

% --- Calculate Central Frequencies ---
% Central frequency calculation for each band
centralFrequencies = zeros(numChannels, 1); % Preallocate array for central frequencies
for k = 1:numChannels
    lowCutoff = bandEdges(k);
    highCutoff = bandEdges(k + 1);
    centralFrequencies(k) = (lowCutoff + highCutoff) / 2; % Midpoint of the frequency band
end


% % Debug: Display frequency bands
% disp(['Frequency bands: ', num2str(bandEdges)]);

% Design bandpass filters
filters = cell(numChannels, 1);
for k = 1:numChannels
    lowCutoff = bandEdges(k);
    highCutoff = bandEdges(k + 1);
    
    % Ensure the high cutoff of the last band is below Nyquist
    if k == numChannels
        highCutoff = highCutoff * 0.99; % Reduce slightly to ensure it is strictly less than 1
    end
    
    % Normalize the cutoff frequencies to the Nyquist frequency
    normalizedLowCutoff = lowCutoff / (sampleFrequency / 2);
    normalizedHighCutoff = highCutoff / (sampleFrequency / 2);
    
    % % Debug: Display normalized frequencies
    % disp(['Channel ', num2str(k), ': Low = ', num2str(normalizedLowCutoff), ', High = ', num2str(normalizedHighCutoff)]);
    % 
    % Ensure the normalized frequencies are within (0, 1)
    if normalizedLowCutoff <= 0 || normalizedHighCutoff >= 1
        error(['Cutoff frequencies for channel ', num2str(k), ...
               ' must be within the range (0, Nyquist frequency).']);
    end
    
    % Design the bandpass filter using a Butterworth filter
    [b, a] = butter(4, [normalizedLowCutoff, normalizedHighCutoff], 'bandpass');
    filters{k} = {b, a};
end
% Step 4: Apply bandpass filters three times to each signal
filteredSignals = cell(numChannels, 1); % Initialize cell array for storing filtered outputs
for k = 1:numChannels
    % Retrieve filter coefficients for the current channel
    [b, a] = filters{k}{:};
    
    % Apply the filter three times
    tempSignal = audioData; % Start with the original audio data
    for i = 1:3 % Apply the bandpass filter 3 times
        tempSignal = filter(b, a, tempSignal); % Apply bandpass filter
    end
    
    % Store the final filtered signal
    filteredSignals{k} = tempSignal;
end

%     % Apply the bandpass filters
% filteredSignals = cell(numChannels, 1); % Initialize cell array for storing filtered outputs
% for k = 1:numChannels
%     [b, a] = filters{k}{:}; % Retrieve filter coefficients
%     filteredSignals{k} = filter(b, a, audioData); % Apply filter to the input audio
% end

% Step 7: Rectify the output signals
rectifiedSignals = cell(numChannels, 1); % Initialize cell array for rectified signals
for k = 1:numChannels
    rectifiedSignals{k} = abs(filteredSignals{k}); % Rectify by taking the absolute value
end

% % Debug: Check rectified signals
% disp('First few samples of the lowest frequency rectified signal:');
% disp(rectifiedSignals{1}(1:10)); % Display first 10 samples of the lowest frequency channel
% disp('First few samples of the highest frequency rectified signal:');
% disp(rectifiedSignals{end}(1:10)); % Display first 10 samples of the highest frequency channel


% % Play the lowest and highest frequency channels
% disp('Playing the lowest frequency channel...');
% sound(filteredSignals{1}, sampleFrequency); % Play the lowest frequency channel
% pause(length(filteredSignals{1}) / sampleFrequency); % Wait for playback to finish
% 
% disp('Playing the highest frequency channel...');
% sound(filteredSignals{end}, sampleFrequency); % Play the highest frequency channel
% pause(length(filteredSignals{end}) / sampleFrequency); % Wait for playback to finish

% Debug: Check rectifiedSignals content
disp(['Rectified signals created for ', num2str(numChannels), ' channels.']);
for k = 1:numChannels
    disp(['Channel ', num2str(k), ' Signal Size: ', num2str(size(rectifiedSignals{k}))]);
end

% % Play the lowest and highest frequency rectified signals
% disp('Playing the lowest frequency rectified channel...');
% sound(rectifiedSignals{1}, sampleFrequency); % Play the lowest frequency rectified channel
% pause(length(rectifiedSignals{1}) / sampleFrequency); % Wait for playback to finish
% 
% disp('Playing the highest frequency rectified channel...');
% sound(rectifiedSignals{end}, sampleFrequency); % Play the highest frequency rectified channel
% pause(length(rectifiedSignals{end}) / sampleFrequency); % Wait for playback to finish

% Debug: Check filteredSignals content
disp(['Number of bandpass channels processed: ', num2str(length(filteredSignals))]);
for k = 1:numChannels
    disp(['Channel ', num2str(k), ' Signal Size: ', num2str(size(filteredSignals{k}))]);
end



% % Step 8: Envelope detection with triple Butterworth filtering
% % Design a low-pass filter with a cutoff frequency of 400 Hz
% lowPassCutoff = 400; % Cutoff frequency in Hz
% nyquist = sampleFrequency / 2; % Nyquist frequency
% normalizedCutoff = lowPassCutoff / nyquist; % Normalize cutoff frequency
% [lp_b, lp_a] = butter(4, normalizedCutoff, 'low'); % 4th-order Butterworth low-pass filter

% % Apply the low-pass filter three times to the rectified signals
% envelopeSignals = cell(numChannels, 1); % Initialize cell array for envelope signals
% for k = 1:numChannels
%     tempSignal = rectifiedSignals{k}; % Start with the rectified signal
%     for i = 1:3 % Apply the filter 3 times
%         tempSignal = filter(lp_b, lp_a, tempSignal); % Apply low-pass filter
%     end
%     envelopeSignals{k} = tempSignal; % Store the final result
% end

% % Apply the low-pass filter to the rectified signals
% envelopeSignals = cell(numChannels, 1); % Initialize cell array for envelope signals
% for k = 1:numChannels
%     envelopeSignals{k} = filter(lp_b, lp_a, rectifiedSignals{k}); % Apply low-pass filter
% end
% 
% % Debug: Check envelopeSignals content
% disp(['Envelope signals created for ', num2str(numChannels), ' channels.']);
% for k = 1:numChannels
%     disp(['Channel ', num2str(k), ' Signal Size: ', num2str(size(envelopeSignals{k}))]);
% end



% Step 8: Envelope extraction using manually designed FIR low-pass filter
% Design an FIR low-pass filter manually
fc = 400; % Cutoff frequency (Hz)
N = 101; % Filter order (larger N results in sharper cutoff but increases computation)
nyquist = sampleFrequency / 2; % Nyquist frequency
normalizedCutoff = fc / nyquist; % Normalize cutoff frequency to (0, 1)

% Create the FIR filter coefficients
h = sinc(2 * normalizedCutoff * ((0:N-1) - (N-1)/2)); % Sinc function for FIR design
h = h .* hamming(N)'; % Apply Hamming window for smooth transition
h = h / sum(h); % Normalize the filter coefficients to ensure unity gain

% Apply the FIR low-pass filter to the rectified signals
envelopeSignals = cell(numChannels, 1); % Initialize cell array for envelope signals
for k = 1:numChannels
    rectifiedSignal = rectifiedSignals{k}; % Get the rectified signal from Task 7
    % Apply the FIR low-pass filter using convolution
    envelopeSignals{k} = conv(rectifiedSignal, h, 'same'); % 'same' keeps the original length
end

% Debug: Check envelope signals
disp(['Envelope signals created for ', num2str(numChannels), ' channels using FIR low-pass filter.']);
for k = 1:numChannels
    disp(['Channel ', num2str(k), ' Signal Size: ', num2str(size(envelopeSignals{k}))]);
end

% % Play the envelope signals of the lowest and highest frequency channels
% disp('Playing the envelope of the lowest frequency channel...');
% sound(envelopeSignals{1}, sampleFrequency); % Play the lowest frequency envelope
% pause(length(envelopeSignals{1}) / sampleFrequency); % Wait for playback to finish
% 
% disp('Playing the envelope of the highest frequency channel...');
% sound(envelopeSignals{end}, sampleFrequency); % Play the highest frequency envelope
% pause(length(envelopeSignals{end}) / sampleFrequency); % Wait for playback to finish

% % Plot the envelope signals of the lowest and highest frequency channels
% figure;
% subplot(2, 2, 3);
% plot(envelopeSignals{1}); % Plot the lowest frequency channel envelope
% title('Envelope Signal: Lowest Frequency Channel (FIR Low-Pass Filter)');
% xlabel('Sample Number');
% ylabel('Amplitude');
% 
% subplot(2, 2, 4);
% plot(envelopeSignals{end}); % Plot the highest frequency channel envelope
% title('Envelope Signal: Highest Frequency Channel (FIR Low-Pass Filter)');
% xlabel('Sample Number');
% ylabel('Amplitude');

    % Step 6: Plot the output signals of the lowest and highest frequency channels
    figure;
    subplot(2, 2, 1);
    plot(filteredSignals{1});
    title('Lowest Frequency Channel Output');
    xlabel('Sample Number');
    ylabel('Amplitude');

    subplot(2, 2, 2);
    plot(filteredSignals{end});
    title('Highest Frequency Channel Output');
    xlabel('Sample Number');
    ylabel('Amplitude');

    % Envelope extraction will be implemented in subsequent tasks.

% % Play the lowest and highest frequency envelope signals
% disp('Playing the envelope of the lowest frequency channel...');
% sound(envelopeSignals{1}, sampleFrequency); % Play the lowest frequency envelope
% pause(length(envelopeSignals{1}) / sampleFrequency); % Wait for playback to finish
% 
% disp('Playing the envelope of the highest frequency channel...');
% sound(envelopeSignals{end}, sampleFrequency); % Play the highest frequency envelope
% pause(length(envelopeSignals{end}) / sampleFrequency); % Wait for playback to finish

% Plot the envelope signals of the lowest and
% highest frequency channels
subplot(2, 2, 3);
plot(envelopeSignals{1}); % First channel (lowest frequency)
title('Envelope Signal: Lowest Frequency Channel');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(2, 2, 4);
plot(envelopeSignals{end}); % Last channel (highest frequency)
title('Envelope Signal: Highest Frequency Channel');
xlabel('Sample Number');
ylabel('Amplitude');


% --- Phase 3: Task 10 ---
    % Step 10: Generate Cosine Signals for Each Channel
    cosineSignals = cell(numChannels, 1);
    for k = 1:numChannels
        numSamples = length(rectifiedSignals{k});
        t = (0:numSamples-1)' / sampleFrequency; % Time vector
        cosineSignals{k} = cos(2 * pi * centralFrequencies(k) * t); % Cosine at central frequency
    end

    % % Debug: Display sizes of generated cosine signals
    % disp(['Generated cosine signals for ', num2str(numChannels), ' channels.']);
    % for k = 1:numChannels
    %     disp(['Channel ', num2str(k), ': Signal Length = ', num2str(length(cosineSignals{k}))]);
    % end

    % % Debug: Play one of the cosine signals (e.g., first channel)
    % disp('Playing a cosine signal for the first channel...');
    % sound(cosineSignals{7}, sampleFrequency);
    % pause(numSamples / sampleFrequency);


% % --- Plot Cosine Wave for Channel 1 ---
% % Retrieve the first cosine signal and its corresponding time vector
% channel1Cosine = cosineSignals{1};
% numSamples = length(channel1Cosine);
% t = (0:numSamples-1) / sampleFrequency; % Time vector for plotting
% 
% % Plot two cycles of the cosine wave
% cyclesToPlot = 2; % Number of cycles to plot
% timePerCycle = 1 / centralFrequencies(1); % Time for one cycle at central frequency
% samplesToPlot = floor(cyclesToPlot * timePerCycle * sampleFrequency); % Number of samples for two cycles
% 
% % Plot the cosine wave
% figure;
% plot(t(1:samplesToPlot), channel1Cosine(1:samplesToPlot));
% title('Cosine Wave for Channel 1 (First Band)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% end


% Task 11: Amplitude Modulation
modulatedSignals = cell(numChannels, 1); % Initialize cell array for modulated signals
for k = 1:numChannels
    modulatedSignals{k} = rectifiedSignals{k} .* cosineSignals{k};
end
% Debug: Display information about the modulated signals
disp(['Amplitude modulated signals generated for ', num2str(numChannels), ' channels.']);
for k = 1:numChannels
    disp(['Channel ', num2str(k), ': Signal Length = ', num2str(length(modulatedSignals{k}))]);
end


% Task 12: Summation and Normalization
finalOutputSignal = zeros(size(modulatedSignals{1}));
for k = 1:numChannels
    finalOutputSignal = finalOutputSignal + modulatedSignals{k};
end

finalOutputSignal = finalOutputSignal / max(abs(finalOutputSignal));

% Play the final output signal
disp('Playing the final output signal...');
sound(finalOutputSignal, sampleFrequency);
pause(length(finalOutputSignal) / sampleFrequency);

% Save the final output signal
outputFileName = ['final_output_', fileName];
audiowrite(outputFileName, finalOutputSignal, sampleFrequency);
disp(['Final output signal saved to file: ', outputFileName]);