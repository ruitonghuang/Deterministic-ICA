function error = calError2(A1,A2,y)

sig1 = A1\y;
sig2 = A2\y;
[error,~] = calError(sig1', sig2');