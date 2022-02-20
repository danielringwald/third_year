% Task 2.3 
function hnew = newstep(tol, err, errold, h0ld, k)
hnew = ((tol/err)^(2/(3*k)))*((tol/errold)^(-1/(3*k)))*h0ld;
end