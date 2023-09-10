function estimated_dc = estimateDC(dc1, dc2, vtg1, vtg2, targetvtg)
    % Given two points (dc1, vtg1) and (dc2, vtg2), this function fits a linear
    % relationship between dc and vtg. It then estimates the dc value that
    % would yield the targetvtg using the linear fit.

    % Validate the input
    if dc1 == dc2 || vtg1 == vtg2
        error('Input data points cannot have the same DC or VTG values.');
    end
    
    % Fit the linear relationship: vtg = m*dc + c
    m = (vtg2 - vtg1) / (dc2 - dc1);
    c = vtg1 - m*dc1;

    % Estimate the dc for the target vtg
    estimated_dc = (targetvtg - c) / m;
end
