function [filteredSignals, centralFrequencies] = bandpass_filter(audioData, fs, numChannels)
    lowerFreq = 100;
    upperFreq = min(8000, fs / 2);
    if lowerFreq >= upperFreq
        error('Invalid frequency range.');
    end

    % Logarithmically spaced band edges
    bandEdges = logspace(log10(lowerFreq), log10(upperFreq), numChannels + 1);
    centralFrequencies = zeros(numChannels, 1);
    filters = cell(numChannels, 1);

    for k = 1:numChannels
        low = bandEdges(k); high = bandEdges(k+1);
        if k == numChannels
            high = high * 0.99;  % stay under Nyquist
        end
        normLow = low / (fs/2);
        normHigh = high / (fs/2);

        [b, a] = butter(4, [normLow, normHigh], 'bandpass');
        filters{k} = {b, a};
        centralFrequencies(k) = (low + high) / 2;
    end

    % Apply each bandpass filter 3 times
    filteredSignals = cell(numChannels, 1);
    for k = 1:numChannels
        [b, a] = filters{k}{:};
        temp = audioData;
        for i = 1:3
            temp = filter(b, a, temp);
        end
        filteredSignals{k} = temp;
    end
end
