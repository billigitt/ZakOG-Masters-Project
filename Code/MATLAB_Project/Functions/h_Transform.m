function Output = h_Transform(c,delta)
%Function that computes the  dot product of v(c) and delta

N = length(delta);

v = exp(-c*(1:N));

Output = dot(delta, v);

end

