%{
function dydt = vdp(t, y)
	dydt = zeros([numel(y) 1]);
	dydt(1) = y(2);
	dydt(2) = (1.0 - y(1)^2) * y(2) - y(1);
end
%}

function dydt = odesrhs(t, y, coeff)
	dydt = zeros([numel(y) 1]);
	dydt(1) = -coeff * y(1);
end

