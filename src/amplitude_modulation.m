function modulatedSignals = amplitude_modulation(envelopes, filteredSignals)
    numChannels = length(envelopes);
    modulatedSignals = cell(numChannels, 1);

    for k = 1:numChannels
        % Multiply envelope with original band signal (preserves fine structure)
        modulatedSignals{k} = envelopes{k} .* filteredSignals{k};
    end
end
