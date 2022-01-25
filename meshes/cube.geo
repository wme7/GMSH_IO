// Gmsh project created on Fri Apr 1 9:13:32 2019
lc = 0.02;
Point(1) = {-0.05,-0.03,-0.03, lc};
Point(2) = { 0.05,-0.03,-0.03, lc};
Point(3) = { 0.05, 0.03,-0.03, lc};
Point(4) = {-0.05, 0.03,-0.03, lc};
Point(5) = {-0.05,-0.03, 0.03, lc};
Point(6) = { 0.05,-0.03, 0.03, lc};
Point(7) = { 0.05, 0.03, 0.03, lc};
Point(8) = {-0.05, 0.03, 0.03, lc};
Line(1) = {8, 7};
Line(2) = {7, 6};
Line(3) = {6, 5};
Line(4) = {5, 8};
Line(5) = {3, 2};
Line(6) = {2, 1};
Line(7) = {1, 4};
Line(8) = {4, 3};
Line(9) = {3, 7};
Line(10) = {2, 6};
Line(11) = {8, 4};
Line(12) = {5, 1};
Line Loop(13) = {9, 2, -10, -5};  Plane Surface(14) = {13};
Line Loop(15) = {1, -9, -8, -11}; Plane Surface(16) = {15};
Line Loop(17) = {8, 5, 6, 7};     Plane Surface(18) = {17};
Line Loop(19) = {3, 12, -6, 10};  Plane Surface(20) = {19};
Line Loop(21) = {12, 7, -11, -4}; Plane Surface(22) = {21};
Line Loop(23) = {2, 3, 4, 1};     Plane Surface(24) = {-23};
Surface Loop(25) = {24, 14, 16, 18, 20, 22};
Volume(26) = {25};

// Set mesh Physical parameters
Physical Surface("free") = {14,16,18,20,22,24};
Physical Volume("fluid") = {26};