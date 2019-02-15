% main file that runs analysis for sisters, randomised and mother daughters

% This analysis was run on a server, and 10,000 MCMC runs takes approx 10
% days. The number of MCMC runs can be adjusted in each of the individual
% files
% 
% Figures 4,5 and 6 can then be recreated, and the expected output
% corresponds to the figures in the paper

% analysis for sister cells
% more detail can be found in Mapfinder1

a = batch('Mapfinder1')
b = batch('Mapfinder2')
c = batch('Mapfinder3')
d = batch('Mapfinder4')
e = batch('Mapfinder5')
f = batch('Mapfinder6')
g = batch('Mapfinder7')
h = batch('Mapfinder8')
i = batch('Mapfinder9')
j = batch('Mapfinder10')
k = batch('Mapfinder11')
l = batch('Mapfinder12')

% analysis of randomised pairings of cells
a1 = batch('Copy_of_Mapfinder1')
b1 = batch('Copy_of_Mapfinder2')
c1 = batch('Copy_of_Mapfinder3')
d1 = batch('Copy_of_Mapfinder4')
e1 = batch('Copy_of_Mapfinder5')
f1 = batch('Copy_of_Mapfinder6')
g1 = batch('Copy_of_Mapfinder7')
h1 = batch('Copy_of_Mapfinder8')
i1 = batch('Copy_of_Mapfinder9')

% analysis of cells in same micro-environment
a2 = batch('Micro_1')
b2 = batch('Micro_2')
c2 = batch('Micro_3')
d2 = batch('Micro_4')

wait(a)
wait(b)
wait(c)
wait(d)
wait(e)
wait(f)
wait(g)
wait(h)
wait(i)
wait(j)
wait(k)
wait(l)

wait(a1)
wait(b1)
wait(c1)
wait(d1)
wait(e1)
wait(f1)
wait(g1)
wait(h1)
wait(i1)

wait(a2)
wait(b2)
wait(c2)
wait(d2)

wait(a3)


quit