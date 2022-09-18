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

Point(1) = {0., -Hext, 0., myhf};
Point(2) = {Rf, -Hext, 0., myhf};
Point(3) = {Rf, H + Hext, 0., myhfine};
Point(4) = {0., H + Hext, 0., myhfine};

Point(5) = {Rf, 0., 0., myhfine};
Point(6) = {R, 0., 0., myh};
Point(7) = {R, H, 0., myh};
Point(8) = {Rf, H , 0., myhsuperfine};
Point(9) = {Rf+0.2*R, 0 , 0., myh};
Point(10) = {Rf+0.2*R, H , 0., myh};
Point(11) = {Rf, 0.5*H , 0., myhsuperfine};


Line(1) = {1, 2};
Line(2) = {2, 5};
Line(5) = {8, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(6) = {5, 9};
Line(7) = {9, 6};
Line(8) = {6, 7};
Line(9) = {7, 10};
Line(10) = {10, 8};
Line(11) = {8, 11};
Line(12) = {11, 5};

Line Loop(12) = {1, 2, 6, 7, 8, 9, 10, 5, 3, 4};
Line Loop(13) = {6, 7, 8, 9, 10, 11, 12};

// surface fiber
Plane Surface(1) = {12, 13};

// surface matrix
Plane Surface(2) = {13};

Physical Surface(1) = {1};     
Physical Surface(2) = {2};     
