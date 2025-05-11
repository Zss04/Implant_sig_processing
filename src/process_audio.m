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