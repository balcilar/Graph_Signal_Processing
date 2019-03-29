function y = activation_tanh(alpha) 
    % activation of tanh function
    
    y = (exp(alpha)-exp(-alpha))./(exp(alpha)+exp(-alpha));
end