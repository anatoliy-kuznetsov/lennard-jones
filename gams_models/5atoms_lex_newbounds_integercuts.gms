Set i atoms / 1*5 /;
Alias(i,j,k);
Parameters
    max_sq_distance / 2.6563849 /
    min_sq_distance / 0.8008854 /
;
Parameter box_width;
box_width = sqrt(max_sq_distance);
Positive Variables sq_distance(i,j);
sq_distance.lo(i,j) = min_sq_distance;
sq_distance.up(i,j) = max_sq_distance;
Variables x(i), y(i), z(i), p, individual_potential(i,j);
x.up(i) = box_width; y.up(i) = box_width; z.up(i) = box_width;
x.lo(i) = 0; y.lo(i) = -box_width; z.lo(i) = -box_width;
x.fx('1') = 0; y.fx('1') = 0; z.fx('1') = 0;
x.fx('2') = 0; y.fx('2') = 0;
x.fx('3') = 0; 
Equations potential, sq_distance_def(i,j), individual_potential_def(i,j);
potential.. p =e= sum(i,
	sum(j$(ord(j)>ord(i)),
		individual_potential(i,j)
	)
);
sq_distance_def(i,j)$(ord(j)>ord(i)).. sq_distance(i,j) =e= power(x(i)-x(j),2) + power(y(i)-y(j),2) + power(z(i)-z(j),2);
individual_potential_def(i,j)$(ord(j)>ord(i)).. individual_potential(i,j) =e= power(sq_distance(i,j),-3) * (power(sq_distance(i,j),-3) - 2);
Binary Variables u(i,j), v(i,j);
Parameter bigM;
bigM = box_width;
Equations uij_bigM_1(i,j), uij_bigM_2(i,j), uij_implication(i,j), vij_bigM_1(i,j), vij_bigM_2(i,j), vij_implication(i,j);
uij_bigM_1(i,j)$(ord(i)<ord(j)).. -bigM * u(i,j) =l= x(j) - x(i);
uij_bigM_2(i,j)$(ord(i)<ord(j)).. x(j) - x(i) =l= bigM * (1 - u(i,j));
uij_implication(i,j)$(ord(i)<ord(j)).. y(i) - y(j) =l= bigM * (1 - u(i,j));
vij_bigM_1(i,j)$(ord(i)<ord(j)).. -bigM * v(i,j) =l= y(j) - y(i);
vij_bigM_2(i,j)$(ord(i)<ord(j)).. y(j) - y(i) =l= bigM * (1 - v(i,j));
vij_implication(i,j)$(ord(i)<ord(j)).. z(i) - z(j) =l= bigM * (2 - u(i,j) - v(i,j));
Equations ic1(i,j,k), ic2(i,j,k), ic3(i,j,k), ic4(i,j,k), ic5(i,j,k), ic6(i,j,k), ic7(i,j,k);
ic1(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. u(i,j) + u(j,k) - u(i,k) =l= 1;
ic2(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. u(i,j) =l= u(i,j);
ic3(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. u(i,k) =l= u(j,k);
ic4(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. v(i,j) + v(j,k) - v(i,k) =l= 1;
ic5(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. v(i,j) + v(j,k) =g= v(i,k);
ic6(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. u(j,k) + v(i,k) - v(i,j) =l= 1;
ic7(i,j,k)$((ord(i)<ord(j)) and (ord(j)<ord(k))).. u(i,j) + v(i,k) - v(j,k) =l= 1;
Option optcr = 0;
Model m / all /;
Option reslim = 3600;
Solve m using minlp minimizing p;
