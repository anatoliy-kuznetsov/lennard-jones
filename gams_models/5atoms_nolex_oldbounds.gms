Set i atoms / 1*5 /;
Alias(i,j);
Parameters
    max_sq_distance / 25 /
    min_sq_distance / 0.8008854 /
;
Parameter box_width;
box_width = sqrt(max_sq_distance);
Positive Variables sq_distance(i,j);
sq_distance.lo(i,j) = min_sq_distance;
sq_distance.up(i,j) = max_sq_distance;
Variables x(i), y(i), z(i), v, individual_potential(i,j);
x.up(i) = box_width; y.up(i) = box_width; z.up(i) = box_width;
x.lo(i) = 0; y.lo(i) = -box_width; z.lo(i) = -box_width;
x.fx('1') = 0; y.fx('1') = 0; z.fx('1') = 0;
x.fx('2') = 0; y.fx('2') = 0;
x.fx('3') = 0; 
Equations potential, sq_distance_def(i,j), individual_potential_def(i,j);
potential.. v =e= sum(i,
	sum(j$(ord(j)>ord(i)),
		individual_potential(i,j)
	)
);
sq_distance_def(i,j)$(ord(j)>ord(i)).. sq_distance(i,j) =e= power(x(i)-x(j),2) + power(y(i)-y(j),2) + power(z(i)-z(j),2);
individual_potential_def(i,j)$(ord(j)>ord(i)).. individual_potential(i,j) =e= power(sq_distance(i,j),-3) * (power(sq_distance(i,j),-3) - 2);
Option optcr = 0;
Model m / all /;
Option reslim = 3600;
Solve m using nlp minimizing v;
