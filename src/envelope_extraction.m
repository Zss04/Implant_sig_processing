function envelopeSignals = envelope_extraction(rectifiedSignals, fs)
    % FIR-based envelope smoothing
    fc = 400; N = 101;
    nyquist = fs / 2;
    normCutoff = fc / nyquist;

    h = sinc(2 * normCutoff * ((0:N-1) - (N-1)/2));
    h = h .* hamming(N)';
    h = h / sum(h);  % Normalize filter gain

    numChannels = length(rectifiedSignals);
    envelopeSignals = cell(numChannels, 1);
    for k = 1:numChannels
        envelopeSignals{k} = conv(rectifiedSignals{k}, h, 'same');
    end
end
