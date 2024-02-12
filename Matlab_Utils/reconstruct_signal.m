function [output_signal, G1, G2] = reconstruct_signal(input_signal, first_order_kernel, second_order_kernel)
% Adapted from a script by Ayad Al-Rumaithi from the Mathworks File Exchange 
% https://uk.mathworks.com/matlabcentral/fileexchange/71799-volterra-wiener-characterization

    % Calculate the length of input signal
    N_input = length(input_signal);
    
    % Calculate the length of first-order kernel
    N_first_order_kernel = length(first_order_kernel);
    
    % Calculate the power of the input signal
    input_power = sum(input_signal.^2) / N_input;
    
    % Calculate G1 (First-order component)
    % ---------------------------------------------------------------------
    % Convolve the input signal with the first-order kernel
    G1 = conv(first_order_kernel, input_signal);
    
    % Take only the first N_input elements of the convolution result
    G1 = G1(1:N_input);
    
    % Calculate G2 (Second-order component)
    % ---------------------------------------------------------------------
    
    % Initialize the second-order component
    G2 = zeros(N_input, 1);
    
    % Maximum lag for second-order convolution
    max_lag_second_order = N_first_order_kernel - 1;
    
    % Iterate through different lags for second-order convolution
    for lag = 0:max_lag_second_order
        % Calculate the product of delayed and non-delayed input signals
        delayed_input_product = input_signal((lag + 1):N_input) .* input_signal(1:(N_input - lag));
        
        % Extract the relevant part of the second-order kernel
        second_order_kernel_part = zeros(N_first_order_kernel - lag, 1);
        for kk = 1:(N_first_order_kernel - lag)
            second_order_kernel_part(kk, 1) = second_order_kernel(kk, kk + lag);
        end
        
        % Perform convolution between the extracted kernel part and the product
        convolution_result = conv(second_order_kernel_part, delayed_input_product);
        
        % Take only the first [N_input - lag] elements of the convolution result
        convolution_result = convolution_result(1:(N_input - lag));
        
        % Adjust the convolution result based on the lag
        if lag == 0
            scaling_factor = 1;
        else
            scaling_factor = 2;
        end
        
        % Accumulate the convolution result into the second-order component
        G2 = G2 + scaling_factor * [zeros(lag, 1); convolution_result];
    end
    
    % Subtract the contribution from the input power multiplied by the diagonal elements of the second-order kernel
    G2 = G2 - input_power * sum(diag(second_order_kernel));
    
    % Combine G1 and G2 to reconstruct the output signal
    output_signal = G1 + G2;
end