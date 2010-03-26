L3a_coll_m1 = load('L3_aff01a_coll-1int.txt');
L3a_coll_m2 = load('L3_aff01a_coll-2int.txt');
L3a_coll_p1 = load('L3_aff01a_coll1int.txt');
L3a_coll_p2 = load('L3_aff01a_coll2int.txt');
L3a_coll_p3 = load('L3_aff01a_coll3int.txt');
L3a_coll_p4 = load('L3_aff01a_coll4int.txt');
L3a_coll_p5 = load('L3_aff01a_coll5int.txt');
L3a_coll_p6 = load('L3_aff01a_coll6int.txt');
L3a = load('../L3_aff01adcint.txt');

plot3(L3a_coll_m1(:,1),L3a_coll_m1(:,2),L3a_coll_m1(:,3)); hold on;
plot3(L3a_coll_m2(:,1),L3a_coll_m2(:,2),L3a_coll_m2(:,3));
plot3(L3a_coll_p1(:,1),L3a_coll_p1(:,2),L3a_coll_p1(:,3));
plot3(L3a_coll_p2(:,1),L3a_coll_p2(:,2),L3a_coll_p2(:,3));
plot3(L3a_coll_p3(:,1),L3a_coll_p3(:,2),L3a_coll_p3(:,3));
plot3(L3a_coll_p4(:,1),L3a_coll_p4(:,2),L3a_coll_p4(:,3));
plot3(L3a_coll_p5(:,1),L3a_coll_p5(:,2),L3a_coll_p5(:,3));
plot3(L3a_coll_p6(:,1),L3a_coll_p6(:,2),L3a_coll_p6(:,3));

plot3(L3a(:,1),L3a(:,2),L3a(:,3));