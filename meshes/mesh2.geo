//SetFactory("OpenCASCADE");

myh=h;
myhf=myh/2.;
myhfine=myh/4.;
myhsuperfine=myhfine/4.;

H = 0.5;
R = 1.0;
ratio = 10.;
Rf = R/ratio;
Hext = 0.2*H;

Point(1) = {0., 0., 0., myhf};
Point(2) = {Rf + 0.2*R, 0., 0., myh};
Point(3) = {R, 0., 0., myh};
Point(4) = {R, H, 0., myh};
Point(5) = {Rf + 0.2*R, H, 0., myh};
Point(6) = {Rf, H, 0., myhsuperfine};
Point(7) = {Rf, H + Hext, 0., myhfine};
Point(8) = {0, H + Hext, 0., myhfine};
Point(9) = {0, 0.5*H, 0., myhfine};

Point(10) = {Rf, 0.5*H, 0., myhsuperfine};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};

Line(10) = {9, 10};
Line(11) = {10, 6};

Line Loop(12) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Line Loop(13) = {10, 11, 6, 7, 8};

// surface fiber
Plane Surface(1) = {12, 13};

// surface matrix
Plane Surface(2) = {13};

Physical Surface(1) = {1};     
Physical Surface(2) = {2};     
