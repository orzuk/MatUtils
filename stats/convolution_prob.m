% Return the convolution of two probabilities: p+q-2pq
function conv_p = convolution_prob(p, q)

conv_p = p+q-2.*p.*q;
