function r = So(tau, S)
    % shrinkage operator
    r = sign(S) .* max(abs(S) - tau, 0);
end