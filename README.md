# CommOT framework details
A new framework named Communication Optimal Transport (CommOT) for computing the rate distortion (RD) function.

The Alternating Sinkhorn (AS) algorithm is presented in “A Communication Optimal Transport Approach to the Computation of Rate Distortion Functions”
https://ieeexplore.ieee.org/abstract/document/10161675

# The code for the AS algorithm
The latest code is in the file "RD_CBA.m". 

This implement is more efficient than the version in the files "BlahutA.m" and "sinkhorn_admm.m".

Note: the major modification in the file "RD_CBA.m" is the way of matrix multiplication in matlab.

The code in "RD_CBA.m" is of $O(N^2)$ complexity for each multiplication.
