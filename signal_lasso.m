function w_pred = signal_lasso(X, y, alpha, alpha_2, num_iters, w0, intercept)
%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file is the MATLAB implementation of the signal lasso algorithm
%%% based on the paper L. Shi, C. Shen, L. Jin, Q. Shi, Z. Wang, 
%%% and S. Boccaletti, Physical Review Research 3, 043210 (2021)

    [m,n] = size(X);   
    w = w0(:);

    res = 1;
    iters = 0;
    while iters < num_iters && res > 1e-11
        w1 = w;


        for j = 1:n
            y_pred = X*w; 
            y_diff = (y - y_pred);
            
            X_j = X(:, j);
            rho = y_diff + (w(j) * X_j);

            if dot(X_j', X_j) == 0
                z = 0;
            else
                z = dot(rho, X_j) / dot(X_j', X_j);
            end
            %disp(z);%(y_pred));%dot(X_j', X_j))

            if intercept == 0
                w(j) = threshold_signal(z, X_j, alpha, alpha_2);
            end
        end

        res = sum(abs(w1 - w));
        iters = iters + 1;
    end

    w_pred = w';
end

function out = threshold_signal(z, X_j, alpha, alpha_2)
    if dot(X_j', X_j) == 0
        delta_1 = 0;
        delta_2 = 0;
    else
        delta_1 = (alpha + alpha_2) / dot(X_j', X_j);
        delta_2 = (alpha - alpha_2) / dot(X_j', X_j);
    end

    if z <= 0
        out = rcpp_min(0, z + delta_1);
    elseif z > 0 && z <= 1 + delta_2
        out = rcpp_max(0, z - delta_2);
    elseif z > 1 + delta_2
        out = rcpp_max(1, z - delta_1);
    end
end

function result = rcpp_min(a, b)
    if a < b
        result = a;
    else
        result = b;
    end
end

function result = rcpp_max(a, b)
    if a < b
        result = b;
    else
        result = a;
    end
end

