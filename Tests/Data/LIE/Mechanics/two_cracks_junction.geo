Mesh.Algorithm = 1;
Mesh.Optimize = 1;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.CharacteristicLengthMin = 0;
Mesh.CharacteristicLengthMax = 1e+22;
Mesh.RecombinationAlgorithm = 0; 
Mesh.SubdivisionAlgorithm = 1;
Mesh.RecombineAll = 1; 
lc=.01;
lc2=.002;
Point(0)={0,0,0,lc};
Point(1)={.05,0,0,lc};
Point(2)={.05,.1,0,lc};
Point(3)={0,.1,0,lc};
Point(4)={.015,.05,0,lc2};
Point(5)={.025,.05,0,lc2};
Point(6)={.035,.05,0,lc2};
Point(7)={.025,.06,0,lc2};
Point(8)={.025,.04,0,lc2};
Line (1) ={0,1};
Line (2) ={1,2};
Line (3) ={2,3};
Line (4) ={3,0};
Line (5) ={4,5};
Line (6) ={5,6};
Line (7) ={7,5};
Line (8) ={5,8};
Line Loop(1) ={1,2,3,4};
Plane Surface(1) = {1};
Line {5} In Surface {1};
Line {6} In Surface {1};
Line {7} In Surface {1};
Line {8} In Surface {1};
Point {5} In Surface {1};
Physical Surface(1) = {1};
Physical Line(2) = {5,6};
Physical Line(3) = {7,8};
