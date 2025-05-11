function outputSignal = synthesize_output(modulatedSignals)
    outputSignal = zeros(size(modulatedSignals{1}));
    for k = 1:length(modulatedSignals)
        outputSignal = outputSignal + modulatedSignals{k};
    end
    % Normalize amplitude to avoid clipping
    outputSignal = outputSignal / max(abs(outputSignal));
end
