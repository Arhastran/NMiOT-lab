function [out] = Analitical(lambda,z,a)

k = 2*pi/lambda;
Func = @(z) z.*(exp(i*k*z)./z-((exp(i.*k.*sqrt(z.^2+a.^2)))/(sqrt(z.^2+a^2))));
out = Func(z);
end