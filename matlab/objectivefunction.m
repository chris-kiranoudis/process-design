%{
function [y] = objectivefunction (parameters, texp, yexp)
	for i = 1:numel(yexp)	
		y(i) = parameters(1) * exp(-parameters(2) * texp(i)) - yexp(i);
	end
end
%}


function [y] = objectivefunction (parameters, tmpexp, texp, yexp)
	for i = 1:numel(yexp)
		coeff = parameters(2) * 1.e-2 * exp(-parameters(3) / tmpexp(i));
		ymodel = parameters(1) * exp(-coeff * texp(i));
		%y(i) = (ymodel - yexp(i)) / yexp(i);
		y(i) = ymodel - yexp(i);
	end
end

