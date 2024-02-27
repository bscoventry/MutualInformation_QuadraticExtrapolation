# Calculation of mutual information with quadratic estimation bias correction

This repository implements mutual information calculation with quadratic bias correction. The goal of this repository is to demonstrate how QE might be performed for bias correction. Implemented in Matlab with Python coming very soon.
# Calculation of MI
This was originally done using neuroscientific data, with Shannon MI calculations made in the form of Borst and Thuessen Nature Neuroscience.
$$$
I(R,S_x) = \sum_i(p(r_i|s_x)*log_2(\frac{p(r_i|s_x)}{p(r_i)})$$
$$$
'''math
$$ I(R,S) = \sum_i\sum_j(p(s_j)p(r_i|s_j)*log_2(\frac{p(r_i|s_j)}{p(r_i)})$$
'''

#Please Cite
If this is used, we would appreciated a reference to Coventry et al, https://doi.org/10.1093/pnasnexus/pgae082 where this was implemented and use.
