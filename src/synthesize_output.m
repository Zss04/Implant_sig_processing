function outputSignal = synthesize_output(modulatedSignals, filteredSignals, envelopes)
    outputSignal = zeros(size(modulatedSignals{1}));
    for k = 1:length(modulatedSignals)
        outputSignal = outputSignal + modulatedSignals{k};
    end
    % Normalize amplitude to avoid clipping
    outputSignal = outputSignal / max(abs(outputSignal));

    % Plot filtered outputs and envelopes for lowest & highest bands
    figure;
    subplot(2, 2, 1);
    plot(filteredSignals{1});
    title('Filtered Output - Lowest Frequency Channel');
    xlabel('Sample Number'); ylabel('Amplitude');

    subplot(2, 2, 2);
    plot(filteredSignals{end});
    title('Filtered Output - Highest Frequency Channel');
    xlabel('Sample Number'); ylabel('Amplitude');

    subplot(2, 2, 3);
    plot(envelopes{1});
    title('Envelope Signal - Lowest Frequency Channel');
    xlabel('Sample Number'); ylabel('Amplitude');

    subplot(2, 2, 4);
    plot(envelopes{end});
    title('Envelope Signal - Highest Frequency Channel');
    xlabel('Sample Number'); ylabel('Amplitude');
end
