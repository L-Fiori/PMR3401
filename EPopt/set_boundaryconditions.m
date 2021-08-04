function [KGM, FM] = set_boundaryconditions(KG, V_cont)

Ngdl = size(KG, 1);
FM = zeros(Ngdl, 1);
KGM = KG;

for i=1:size(V_cont, 1)
	id = V_cont(i, 1);
	val = V_cont(i, 2);
	FM = FM - KG(:, id)*val;
	KGM(:, id) = 0; 
	KGM(id, 0) = 0;
	KGM(id, id) = 1;
end

FM(V_cont(:, 1)) = V_cont(:, 2);

end