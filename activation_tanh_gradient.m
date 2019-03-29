function gradient = activation_tanh_gradient(y)
    % derivation of tanh function    
    % https://theclevermachine.wordpress.com/tag/tanh-function/
    gradient = 1-y.^2;
end