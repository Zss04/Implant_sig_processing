fileName = 'MaleSpeaking.opus';
folder = 'audio_samples\';
fullfilepath = fullfile(folder, fileName);
[audioData, sampleFrequency] = audioread(fullfilepath);

% Convert stereo to mono if needed
if size(audioData, 2) == 2
    audioData = sum(audioData, 2);
end

% Resample to 16 kHz target sampling rate
targetFs = 16000;
if sampleFrequency > targetFs
    audioData = resample(audioData, targetFs, sampleFrequency);
    sampleFrequency = targetFs;
elseif sampleFrequency < targetFs 
    error('The sampling rate is less than 16 kHz. Use a higher rate input.');
end

% Apply 22-channel filterbank and extract fine structure
numChannels = 22;
[filteredSignals, centralFreqs] = bandpass_filter(audioData, sampleFrequency, numChannels);

% Extract envelope using FIR convolution
rectified = cellfun(@abs, filteredSignals, 'UniformOutput', false);
envelopes = envelope_extraction(rectified, sampleFrequency);

% Modulate with original filtered signal instead of synthetic carrier
modulated = amplitude_modulation(envelopes, filteredSignals);

% Combine channels and normalize
finalSignal = synthesize_output(modulated, filteredSignals, envelopes);

% Output final signal
sound(finalSignal, sampleFrequency);
outputFolder = 'audio_samples';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
outputFileName = fullfile(outputFolder, ['final_output_' fileName]);
audiowrite(outputFileName, finalSignal, sampleFrequency);

