// Gmsh project created on Fri Apr 1 9:13:32 2019
lc = 0.01;

Point(1) = {-0.05,-0.03, 0, lc};
Point(2) = { 0.05,-0.03, 0, lc};
Point(3) = { 0.05, 0.03, 0, lc};
Point(4) = {-0.05, 0.03, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

// Set mesh Physical parameters
Physical Line("free_rec") = {1,2,3,4};
Physical Surface("fluid") = {6};

Mesh.MshFileVersion = 2.2;
