function x = getdata()
seed = 6;
rngs = rng(seed);   onCleanupHandle = onCleanup(@()rng(rngs));

N = 18;
tht = 2*(-N/2:1:N/2-1)/N;
x = (1+cos(pi*tht))/2;
x = x + 0.1*randn(size(x));
x = abs(x);
