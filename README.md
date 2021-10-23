# Mixed-lubrication-Helical-gear
Han's Model(Template: A-A)

SetDirectory[NotebookDirectory[]];

(*Get["1023 A-A.mx"]*)

TypeG1 = 2;(*A*)
TypeG2 = 2;(*A*)

Basic Parameters

Setting Parameters

(*數字1、2分別表示gear1和gear2*)
mn = 2.5;(*Normal module,mm*)
\[Alpha]n = 20 Degree;(*Normal pressure angle*)
T1 = 35;(*tooth number of G1(Pinion)*)
T2 = 35;(*tooth number of G2(Gear)*)
B = 30.0;(*gear width,mm*)
\[Beta]op = 25 Degree;(*helix angle at pitch circle*)
xn1 = 0.2;(*shifted coefficient of G1(Pinion)*)
xn2 = 0.2;(*shifted coefficient of G2(Gear)*)
rotav = 800;(*輸入轉速，單位rpm*)
\[Omega]1 = rotav*2*Pi/60;(*角速度，單位rad/s*)
\[Omega]2 = \[Omega]1*T1/T2;

Normal Parameters

involuten = 
  2*Tan[\[Alpha]n]*(xn1 + xn2)/(T1 + T2) + Tan[\[Alpha]n] - \[Alpha]n;
\[Alpha]wn = \[Alpha]wn /. 
   FindRoot[Tan[\[Alpha]wn] - \[Alpha]wn == involuten, {\[Alpha]wn, 
     0.04}];(*normal mesh pressure angle*)
ycorn = (T1 + T2)/
    2*(Cos[\[Alpha]n]/Cos[\[Alpha]wn] - 
     1);(*correction coefficient of centre distance*)
(*radius of Normal pitch circle,mm*)
rpn1 = N[T1*mn/2, 16];
rpn2 = N[T2*mn/2, 16];
(*radius of Normal base circle,mm*)
rbn1 = N[rpn1*Cos[\[Alpha]n], 16];
rbn2 = N[rpn2*Cos[\[Alpha]n], 16];
(*radius of Normal mesh pitch circle,mm*)
rwn1 = N[ rbn1/Cos[\[Alpha]wn], 16 ];
rwn2 = N[ rbn2/Cos[\[Alpha]wn], 16 ];
(*centre distance*)
O1O2n = rwn1 + rwn2;
(*height of addendum*)
han1 = N[(1 + ycorn - xn2)*mn, 16 ];
han2 = N[(1 + ycorn - xn1)*mn, 16 ];
(*total height of tooth groove*)
hn = N[(2.25 + ycorn - (xn1 + xn2))*mn, 16 ];
(*radius of addendum circle,mm*)
ran1 = N[rpn1 + han1, 16];
ran2 = N[rpn2 + han2, 16];
(*radius of dedendum circle,mm*)
rfn1 = N[ran1 - hn, 16];
rfn2 = N[ran2 - hn, 16];
(*Normal base pitch,mm*)
Pbn = N[2*Pi*rbn1/T1, 
   16];(*Notice:Normal base pitches of two gears is equal to each other*)

Transverse Parameters

mt = mn/Cos[\[Beta]op];(*Transverse module,mm*)
\[Alpha]t = ArcTan[ Tan[\[Alpha]n]/Cos[\[Beta]op]];(*Normal pressure angle*)
involutet = 
  2*Tan[\[Alpha]n]*(xn1 + xn2)/(T1 + T2) + Tan[\[Alpha]t] - \[Alpha]t;
(*咬合壓力角*)
\[Alpha]wt = \[Alpha]wt /. 
   FindRoot[Tan[\[Alpha]wt] - \[Alpha]wt == involutet, {\[Alpha]wt, 
     0.04}];(*transverse mesh pressure angle*)
ycort = (T1 + T2)/(2*Cos[\[Beta]op])*(Cos[\[Alpha]t]/Cos[\[Alpha]wt] - 
     1);(*correction coefficient of centre distance*)
(*radius of pitch cylinder,mm*)
rpt1 = N[T1*mt/2, 16];
rpt2 = N[T2*mt/2, 16];
(*radius of base cylinder,mm*)
rbt1 = N[rpt1*Cos[\[Alpha]t], 16];
rbt2 = N[rpt2*Cos[\[Alpha]t], 16];
(*radius of mesh pitch cylinder,mm*)
rwt1 = N[ rbt1/Cos[\[Alpha]wt], 16 ];
rwt2 = N[ rbt2/Cos[\[Alpha]wt], 16 ];
O1O2t = rwt1 + rwt2;(*centre distance*)
(*height of addendum*)
hat1 = N[(1 + ycort - xn2)*mn, 16 ];
hat2 = N[(1 + ycort - xn1)*mn, 16 ];
(*total height of tooth groove*)
ht = N[(2.25 + ycort - (xn1 + xn2))*mn, 16 ];
(*radius of addendum cylinder,mm*)
rat1 = N[rpt1 + hat1, 16];
rat2 = N[rpt2 + hat2, 16];
(*radius of dedendum cylinder,mm*)
rft1 = N[rat1 - ht, 16];
rft2 = N[rat2 - ht, 16];
(*transverse base pitch,mm*)
Pbt = N[2*Pi*rbt1/T1, 
   16];(*Notice:Transverse base pitches of two gears is equal to each other*)
\

(*length of the "line of action" at transverse plane, i.e. A1A2*)
Lat = N[Sqrt[rat1^2 - rbt1^2] + Sqrt[rat2^2 - rbt2^2] - O1O2t*Sin[\[Alpha]t], 
   16];
\[Beta]b = 
  ArcTan[ Tan[\[Beta]op]/rpt1*rbt1 ];(*helix angle at base cylinder*)
\[Epsilon]\[Alpha] = N[ Lat/Pbt, 16 ];(*transverse contact ratio*)
\[Epsilon]\[Beta] = N[ B*Tan[\[Beta]b]/Pbt, 16 ];(*axial contact ratio*)
\[Epsilon]\[Gamma] = \[Epsilon]\[Alpha] + \[Epsilon]\[Beta];(*total contact \
ratio*)
(*端面作用線上各線段長度*)
N1N2 = N[O1O2t*Sin[\[Alpha]t], 16];
N2A1 = N[Sqrt[rat2^2 - rbt2^2], 16];
N1A2 = N[Sqrt[rat1^2 - rbt1^2], 16];
N1A1 = N[N1N2 - N2A1, 16];(*式[1f]Subscript[中之lN, 1]A*)
N2A2 = N[N1N2 - N1A2, 16];
O1A1 = N[  Sqrt[rbt1^2 + N1A1^2], 16  ];
O1A2 = rat1;
O2A1 = rat2;
O2A2 = N[  Sqrt[rbt2^2 + N2A2^2], 16  ];
(*Fn=P/(\[Omega]1*rbt1*10^-3);(*各嚙合齒齒面正向力之和，N*)*)
(*FT=T/rpt1;*)

Texture

Reishauer coordinate system (+2\[Degree])

tdata[2] = Transpose[ Flatten[Import["A 2kong.mat"], 1] ];
{Min[tdata[2]], Max[tdata[2]]};
colorrange = {-2.0, 2.0};
colorbar = BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row"];
topview[2] = 
  ReliefPlot[tdata[2], 
   ColorFunction -> ColorData[{"Rainbow", colorrange}],(*AspectRatio\[Rule]2,*)
   Frame -> False, ColorFunctionScaling -> False, LightingAngle -> None, 
   BoxRatios -> {1, 1, 1}];
(*Export["A +2.wmf",topview[2](*,ImageResolution\[Rule] 1000*)];*)

Reishauer coordinate system (+1\[Degree])

tdata[1] = Transpose[ Flatten[Import["B 1kong.mat"], 1] ];
{Min[tdata[1]], Max[tdata[1]]};
colorrange = {-2.0, 2.0};
colorbar = BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row"];
topview[1] = 
  ReliefPlot[tdata[1], 
   ColorFunction -> ColorData[{"Rainbow", colorrange}],(*AspectRatio\[Rule]2,*)
   Frame -> False, ColorFunctionScaling -> False, LightingAngle -> None, 
   BoxRatios -> {1, 1, 1}];
(*Export["B +1.wmf",topview[1]];*)

Reishauer coordinate system (0\[Degree])

tdata[0] = Transpose[ Flatten[Import["C 0kong.mat"], 1] ];
{Min[tdata[0]], Max[tdata[0]]};
colorrange = {-2.0, 2.0};
colorbar = BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row"];
topview[0] = 
  ReliefPlot[tdata[0], 
   ColorFunction -> ColorData[{"Rainbow", colorrange}],(*AspectRatio\[Rule]2,*)
   Frame -> False, ColorFunctionScaling -> False, LightingAngle -> None, 
   BoxRatios -> {1, 1, 1}];
(*Export["C 0.wmf",topview[0]];*)

Reishauer coordinate system (-1\[Degree])

tdata[-1] = Transpose[ Flatten[Import["D -1kong.mat"], 1] ];
{Min[tdata[-1]], Max[tdata[-1]]};
colorrange = {-2.0, 2.0};
colorbar = BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row"];
topview[-1] = 
  ReliefPlot[tdata[-1], 
   ColorFunction -> ColorData[{"Rainbow", colorrange}],(*AspectRatio\[Rule]2,*)
   Frame -> False, ColorFunctionScaling -> False, LightingAngle -> None, 
   BoxRatios -> {1, 1, 1}];
(*Export["D -1.wmf",topview[-1]];*)

Reishauer coordinate system (-2\[Degree])

tdata[-2] = Transpose[ Flatten[Import["E -2kong.mat"], 1] ];
{Min[tdata[-2]], Max[tdata[-2]]};
colorrange = {-2.0, 2.0};
colorbar = BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row"];
topview[-2] = 
  ReliefPlot[tdata[-2], 
   ColorFunction -> ColorData[{"Rainbow", colorrange}],(*AspectRatio\[Rule]2,*)
   Frame -> False, ColorFunctionScaling -> False, LightingAngle -> None, 
   BoxRatios -> {1, 1, 1}];
(*Export["E -2.wmf",topview[-2]];*)

Datum surface of Gear1 and Gear2 (右齒面)

Datum surface

(*齒輪1基準齒面及基圓柱*)
(*parameter equation of involute segment, \[Phi] is action angle*)
x1[\[Theta]_] := rbt1*(Cos[\[Theta] + Pi/2] + \[Theta]*Sin[\[Theta] + Pi/2]);
y1[\[Theta]_] := rbt1*(Sin[\[Theta] + Pi/2] - \[Theta]*Cos[\[Theta] + Pi/2]);
\[Theta]1[z_] := z*Tan[\[Beta]op]/rpt1;(*deflection angle*)
\[Theta]1min = N1A1/rbt1;
\[Theta]1max = N1A2/rbt1;
SurG1[\[Theta]_, z_] := ({
     {Cos[\[Theta]1[z]], -Sin[\[Theta]1[z]], 0},
     {Sin[\[Theta]1[z]], Cos[\[Theta]1[z]], 0},
     {0, 0, 1}
    }).({
     {x1[\[Theta]]},
     {y1[\[Theta]]},
     {z}
    });


(*齒輪2基準齒面及基圓柱*)
x2[\[Theta]_] := rbt2*(Cos[\[Theta] + Pi/2] + \[Theta]*Sin[\[Theta] + Pi/2]) ;
y2[\[Theta]_] := rbt2*(Sin[\[Theta] + Pi/2] - \[Theta]*Cos[\[Theta] + Pi/2]);
\[Theta]2[z_] := z*Tan[-\[Beta]op]/rpt2;(*deflection angle*)
\[Theta]2min = N2A2/rbt2;
\[Theta]2max = N2A1/rbt2;
SurG2[\[Theta]_, z_] := ({
     {Cos[\[Theta]2[z]], -Sin[\[Theta]2[z]], 0},
     {Sin[\[Theta]2[z]], Cos[\[Theta]2[z]], 0},
     {0, 0, 1}
    }).({
     {x2[\[Theta]]},
     {y2[\[Theta]]},
     {z}
    });


Check Coordinate Transform

Tz = 200;(*Reishauer坐標系原始數據中，z軸之無因次長度（取樣點數）*)
dr = 1/Tz;(*Reishauer坐標系中r軸單位間隔，單位mm*)
dz = 1/Tz;(*Reishauer坐標系中z軸單位間隔，單位mm*)
\[CapitalDelta]r = Tz*dr; \[CapitalDelta]z = Tz*dr;(*下圖中取樣間隔*)
(*漸開線作用角\[Theta] \[Rule] 連心線長r\[Theta]1*)
r\[Theta]1[\[Theta]_] := Sqrt[(x1[\[Theta]]^2 + y1[\[Theta]]^2)];
(*連心線長r\[Theta]1 \[Rule] 漸開線作用角\[Theta]r1*)
\[Theta]r1[r_] := \[Theta] /. 
   FindRoot[r\[Theta]1[\[Theta]] == 
     r, {\[Theta], (\[Theta]1min + \[Theta]1max)/2}];
rclt1 = r\[Theta]1[\[Theta]1min];(*形圓(底隙圓)半徑*)
(*在直角坐標係中，離散參數點(r,z)對應的坐標*)
DisSurG1[r_, 
   z_] := { { Cos[\[Theta]1[z]], -Sin[\[Theta]1[z]], 0}, {Sin[\[Theta]1[z]], 
     Cos[\[Theta]1[z]], 0}, {0, 0, 1} }.{x1[\[Theta]r1[r]], y1[\[Theta]r1[r]],
     z};


(*漸開線作用角\[Theta] \[Rule] 連心線長r\[Theta]2*)
r\[Theta]2[\[Theta]_] := Sqrt[(x2[\[Theta]]^2 + y2[\[Theta]]^2)];
(*連心線長r\[Theta]2 \[Rule] 漸開線作用角\[Theta]r2*)
\[Theta]r2[r_] := \[Theta] /. 
   FindRoot[r\[Theta]2[\[Theta]] == 
     r, {\[Theta], (\[Theta]2min + \[Theta]2max)/2}];
rclt2 = r\[Theta]2[\[Theta]2min];
DisSurG2[r_, 
   z_] := { { Cos[\[Theta]2[z]], -Sin[\[Theta]2[z]], 0}, {Sin[\[Theta]2[z]], 
     Cos[\[Theta]2[z]], 0}, {0, 0, 1} }.{x2[\[Theta]r2[r]], y2[\[Theta]r2[r]],
     z};


Real Surface in Space Rectangular Coordinate System

Subroutine

(*內插法求距離r0最近的Reishauer無因次坐標（矩陣序數）*)
ReishauerMatNum[rinput_, rtop_] := Module[{},
   n = (rtop - rinput)/dr;
   If[n - Floor[n] < 0.5, output = Floor[n], output = Ceiling[n]]
   ];
(*距離rinput最近的Reishauer半徑，rtop為齒頂圓半徑*)
ReishauerRaius[rinput_, rtop_] := Module[{},
   n = (rtop - rinput)/dr;
   rr1 = N[rtop - Floor[n]*dr];
   rr2 = N[rtop - Ceiling[n]*dr];
   If[ rinput > (rr1 + rr2)/2, output = rr1, output = rr2 ]
   ];
(*距離zinput最近的Reishauer齒長，B為總齒長*)
ReishauerZ[zinput_] := Module[{},
   n = zinput/dz;
   zz1 = N[Floor[n]*dz];
   zz2 = N[Ceiling[n]*dz];
   If[ zinput < (zz1 + zz2)/2, output = zz1, output = zz2 ]
   ];
(*空間坐標平面化。說明：
1.OriPoint:新平面直角坐標系之原點，在原3D空間中之坐標；
2.PointOnX:原空間中一點，此點與原點之連線，為新平面直角坐標之x軸；
3.PointP：原空間中一點，欲轉換入新平面直角坐標中的任意點*)
Planarity[OriPoint_, PointOnX_, PointP_] := Module[{},
   \[Rho] = Norm[PointP - OriPoint];
   \[Theta] = VectorAngle[PointP - OriPoint, PointOnX - OriPoint];
   {xp, yp} = 
    N[CoordinateTransform[ "Polar" -> "Cartesian", {\[Rho], \[Theta]}]];
   TargetPoint = {xp, yp};
   ];

Gear1

取樣區域圖示(Gear1)

(*利用Reishauer坐標系中的坐標參數r,z，表示直角坐標系中的齒面方程，mm*)
xG1[r_, z_] := 
  x1[\[Theta]r1[r]]*Cos[\[Theta]1[z]] - y1[\[Theta]r1[r]]*Sin[\[Theta]1[z]];
yG1[r_, z_] := 
  x1[\[Theta]r1[r]]*Sin[\[Theta]1[z]] + y1[\[Theta]r1[r]]*Cos[\[Theta]1[z]];
zG1[z_] := z;
Hws = 50;
rmin = rwt1 - Hws*dr;
rmax = rwt1 + Hws*dr;
zmin = B/2 - Hws*dz;
zmax = B/2 + Hws*dz;
SB = ListPointPlot3D[{{xG1[rmin, zmax], yG1[rmin, zmax], zmax}}, 
   PlotStyle -> {PointSize[0.02], Orange}];(*原點空間坐標*)
SA = ListPointPlot3D[{{xG1[rmin, zmin], yG1[rmin, zmin], zmin}}, 
   PlotStyle -> {PointSize[0.02], Black}];(*參考點空間坐標*)
figSampleG1 = 
  ParametricPlot3D[
   SurG1[\[Theta], z], {\[Theta], \[Theta]r1[rmin], \[Theta]r1[rmax]}, {z, 
    zmin, zmax}, PlotStyle -> {Blue, Opacity[0.8]}, Mesh -> None];
figrwt1cyl2 = 
  ParametricPlot3D[{rmin*Cos[\[Theta]], rmin*Sin[\[Theta]], z}, {\[Theta], 
    Pi/2 + 0.048*Pi, Pi/2 + 0.058*Pi}, {z, zmin, zmax}, 
   PlotStyle -> Opacity[0.5], Mesh -> None];(*半徑為取樣下限rmin的圓*)
figrwt1cyl3 = 
  ParametricPlot3D[{rmax*Cos[\[Theta]], rmax*Sin[\[Theta]], z}, {\[Theta], 
    Pi/2 + 0.048*Pi, Pi/2 + 0.058*Pi}, {z, zmin, zmax}, 
   PlotStyle -> {Opacity[0.5], Purple}, Mesh -> None];(*半徑為取樣上限rmax的圓*)
(*Show[figSampleG1,SA,SB,PlotRange->All,Boxed->False]*)

平面化(Gear1)

(*取樣函數FunRealSur：
0.磨紋種類Type。
1.在Reishauer坐標系中取樣後，將取樣點轉換至直角坐標系。
2.取樣間隔放大率為Amp（即無因次取樣週期），實際取樣週期為\[CapitalDelta]r,\[CapitalDelta]z。
3.取樣半寬：(half width of sample)為Hws，即某方向取樣點數的一半。
4.取樣中心：參數坐標(r,z)；取樣區中任意點符號(r0,z0)。
5.範圍：rmin<r0<rmax; zmin<B<zmax*)
FunRealSur1[Type_, Hws_, r_, z_] := Module[{},
   RealSurG1 = {};
   rmin = r - Hws*dr;
   rmax = r + Hws*dr;
   zmin = z - Hws*dz;
   zmax = z + Hws*dz;
   rmin = ReishauerRaius[rmin, rat1];
   rmax = ReishauerRaius[rmax, rat1];
   zmin = ReishauerZ[zmin];
   zmax = ReishauerZ[zmax];
   SB = {xG1[rmin, zmax], yG1[rmin, zmax], zmax};(*原點*)
   SA = {xG1[rmin, zmin], yG1[rmin, zmin], zmin};
   For[z1 = zmin, z1 <= zmax, z1 += dz,
    numz = Floor[Mod[ z1/dz, Tz ]] + 1;
    For[r1 = rmax, r1 >= rmin, r1 -= dr,
     numr = Floor[(rat1 - r1)/dr] + 1;
     Nd = tdata[Type][[numr, numz]]*10^-3;(*取樣點P(r1,z1)的法向偏距Nd，單位：mm*)
     PointP = {xG1[r1, z1], yG1[r1, z1], z1};(*取樣點PointP的空間坐標*)
     Planarity[SB, SA, PointP];(*以SB為原點，SBSA為x軸平面化*)
     AppendTo[RealSurG1, {xp, yp, Nd}];
     ];
    ];
   ];

FunRealSur1[TypeG1, 30, rwt1, B/2]

G1 = ListPlot3D[RealSurG1, Mesh -> None, 
   ColorFunction -> Function[{x, y, z}, ColorData["Rainbow"][z]], 
   PlotRange -> All, Boxed -> False, 
   AxesLabel -> {"\!\(\*SubscriptBox[\(x\), \(1\)]\)", 
     "\!\(\*SubscriptBox[\(y\), \(1\)]\)"}, AxesStyle -> Thick, 
   LabelStyle -> {FontFamily -> "Times New Roman", 20, GrayLevel[0]}, 
   ViewPoint -> Above];
G1ForCom = 
  ListPlot3D[RealSurG1, Mesh -> None, 
   ColorFunction -> Function[{x, y, z}, ColorData["Rainbow"][z]], 
   PlotRange -> All, Boxed -> False, 
   AxesLabel -> {"\!\(\*SubscriptBox[\(x\), \(1\)]\)", 
     "\!\(\*SubscriptBox[\(y\), \(1\)]\)"}, AxesStyle -> Thick, 
   LabelStyle -> {FontFamily -> "Times New Roman", 20, GrayLevel[0]}, 
   ViewPoint -> Above];


Gear2

取樣區域圖示(Gear2)

(*利用Reishauer坐標系中的坐標參數r,z，表示直角坐標系中的齒面方程，mm*)
xG2[r_, z_] := 
  x2[\[Theta]r2[r]]*Cos[\[Theta]2[z]] - y2[\[Theta]r2[r]]*Sin[\[Theta]2[z]];
yG2[r_, z_] := 
  x2[\[Theta]r2[r]]*Sin[\[Theta]2[z]] + y2[\[Theta]r2[r]]*Cos[\[Theta]2[z]];
zG2[z_] := z;
Hws = 50;
rmin = rwt2 - Hws*dr;
rmax = rwt2 + Hws*dr;
zmin = B/2 - Hws*dz;
zmax = B/2 + Hws*dz;
SB = ListPointPlot3D[{{xG2[rmax, zmax], yG2[rmax, zmax], zmax}}, 
   PlotStyle -> {PointSize[0.02], Orange}];(*原點空間坐標*)
SA = ListPointPlot3D[{{xG2[rmax, zmin], yG2[rmax, zmin], zmin}}, 
   PlotStyle -> {PointSize[0.02], Black}];(*參考點空間坐標*)
figSampleG2 = 
  ParametricPlot3D[
   SurG2[\[Theta], z], {\[Theta], \[Theta]r2[rmin], \[Theta]r2[rmax]}, {z, 
    zmin, zmax}, PlotStyle -> {Red, Opacity[0.8]}, Mesh -> None];
figrwt2cyl2 = 
  ParametricPlot3D[{rmin*Cos[\[Theta]], rmin*Sin[\[Theta]], z}, {\[Theta], 
    Pi/2 - 0.045*Pi, Pi/2 - 0.03*Pi}, {z, zmin, zmax}, 
   PlotStyle -> Opacity[0.5], Mesh -> None];(*半徑為取樣下限rmin的圓*)
figrwt2cyl3 = 
  ParametricPlot3D[{rmax*Cos[\[Theta]], rmax*Sin[\[Theta]], z}, {\[Theta], 
    Pi/2 - 0.045*Pi, Pi/2 - 0.03*Pi}, {z, zmin, zmax}, 
   PlotStyle -> {Opacity[0.5], Purple}, Mesh -> None];(*半徑為取樣上限rmax的圓*)


平面化(Gear2)

FunRealSur2[Type_, Hws_, r_, z_] := Module[{},
   RealSurG2 = {};
   RealSurForComG2 = {};
   rmin = r - Hws*dr;
   rmax = r + Hws*dr;
   zmin = z - Hws*dz;
   zmax = z + Hws*dz;
   rmin = ReishauerRaius[rmin, rat2];
   rmax = ReishauerRaius[rmax, rat2];
   zmin = ReishauerZ[zmin];
   zmax = ReishauerZ[zmax];
   SB = {xG2[rmin, zmax], yG2[rmin, zmax], zmax};(*原點*)
   SA = {xG2[rmin, zmin], yG2[rmin, zmin], zmin};
   For[z2 = zmax, z2 >= zmin, z2 -= dz,
    numz = Floor[Mod[ z2/dz, Tz ]];
    For[r2 = rmax, r2 >= rmin, r2 -= dr,
     numr = Floor[(rat2 - r2)/dr] + 1;
     Nd = tdata[Type][[numr, numz]]*10^-3;(*取樣點P(r2,z2)的法向偏距Nd，單位：mm*)
     PointP = {xG2[r2, z2], yG2[r2, z2], z2};(*取樣點PointP的空間坐標*)
     Planarity[SB, SA, PointP];(*以SB為原點，SBSA為x軸平面化*)
     AppendTo[RealSurG2, {xp, yp, Nd}];
     AppendTo[RealSurForComG2, {xp, 0.5 - yp, 0.006 - Nd}];
     ];
    ];
   ];

FunRealSur2[TypeG2, 30, rwt1, B/2]

G2 = ListPlot3D[RealSurG2, Mesh -> None, 
   ColorFunction -> Function[{x, y, z}, ColorData["Rainbow"][z]], 
   PlotRange -> All, Boxed -> False, 
   AxesLabel -> {"\!\(\*SubscriptBox[\(x\), \(2\)]\)", 
     "\!\(\*SubscriptBox[\(y\), \(2\)]\)"}, AxesStyle -> Thick, 
   LabelStyle -> {FontFamily -> "Times New Roman", 20, GrayLevel[0]}, 
   ViewPoint -> Above];
G2ForCom = 
  ListPlot3D[RealSurForComG2, Mesh -> None, 
   ColorFunction -> Function[{x, y, z}, ColorData["Rainbow"][z]], 
   PlotRange -> All, Boxed -> False, 
   AxesLabel -> {"\!\(\*SubscriptBox[\(x\), \(1\)]\)", 
     "\!\(\*SubscriptBox[\(y\), \(1\)]\)"}, AxesStyle -> Thick, 
   LabelStyle -> {FontFamily -> "Times New Roman", 20, GrayLevel[0]}, 
   ViewPoint -> Above];

取樣區域合併展示

Show[G1ForCom, G2ForCom]

赫茲接觸半寬a(Semi-contact Width)

Constant

FT = 500;(*transmission force,N*)
Eeq = 2.2831*10^11/10^6;(*equivalent Young's modulus,N/mm^2*)
t1 = B*Tan[\[Beta]b]/(rbt1*\[Omega]1);(*time,s*)
t2 = Lat/(rbt1*\[Omega]1);(*time,s*)
t3 = (B*Tan[\[Beta]b] + Lat)/(rbt1*\[Omega]1);(*time,s*)
E\[Gamma] = Floor[\[Epsilon]\[Gamma]];(*總接觸率的整數部分*)
b0 = (0.5*(1 + \[Epsilon]\[Alpha]/2)^2 - 1)^-0.5;(*負載函數中的常數*)
stepnum = 30;(*步階數*)
dt = t3/stepnum;(*單位步階時間間隔,s*)


Subroutine

AuxiLoad[] := Module[{},
   zsvalue[step] = {};
   (*auxilary function to calculate load*)
   zetasup = \[Zeta]0[t] + i + \[Epsilon]\[Alpha] - Min[\[Zeta]0[t] + i, 0] - 
     Max[\[Zeta]0[t] + i, \[Epsilon]\[Alpha]];
   zetainf = \[Zeta]0[t] + i - \[Epsilon]\[Beta] + \[Epsilon]\[Alpha] - 
     Min[\[Zeta]0[t] + i - \[Epsilon]\[Beta], 0] - 
     Max[\[Zeta]0[t] + i - \[Epsilon]\[Beta], \[Epsilon]\[Alpha]];
   \[CapitalIota]\[Nu][t] = 1/b0*\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 
        0\), \(E\[Gamma]\)]\((Sin[\ 
         b0*\((zetasup - \[Epsilon]\[Alpha]/2)\)\ ] - 
        Sin[\ b0*\((zetainf - \[Epsilon]\[Alpha]/2)\)\ ])\)\);
   ];
DynaCharAlongLoC[] := Module[{},
   AppendTo[zsvalue[step], zs];
   R2[t, zs] = N1N2 - R1[t, zs];(*mm*)
   R[t, zs] = R1[t, zs]*R2[t, zs]/(R1[t, zs] + R2[t, zs]);(*mm*)
   Rstep[step, zs] = R[t, zs];
   AppendTo[eqRoC, {t, zs, R[t, zs]}];
   (*Surface velocity*) 
   U1[t, zs] = \[Omega]1*R1[t, zs];(*unit:mm/s*)
   U2[t, zs] = \[Omega]2*R2[t, zs];(*unit:mm/s*)
   (*rolling velocity*)
   U[t, zs] = U1[t, zs] + U2[t, zs];(*unit:mm/s*)
   Ustep[step, zs] = U[t, zs];
   (*sliding velocity*)
   Us[t, zs] = Abs[ U1[t, zs] - U2[t, zs] ];(*unit:mm/s*)
   Usstep[step, zs] = Us[t, zs];
   \[Nu][t, zs] = Cos[b0*(\[Zeta][t, zs] - \[Epsilon]\[Alpha]/2)];
   f[t, zs] = (\[Epsilon]\[Beta]*Cos[\[Beta]b])/
     B*\[Nu][t, zs]/\[CapitalIota]\[Nu][t]*FT;(*unit:N/mm*)
   fstep[step, zs] = f[t, zs];
   AppendTo[Load, {t, zs, f[t, zs]}];
   (*contact stress, N/mm^2*)
   ph[t, zs] = Sqrt[ (Eeq*f[t, zs])/(Pi*R[t, zs])];
   phstep[step, zs] = ph[t, zs];
   (*Hertzian contact half width, mm*)
   a[t, zs] = 2*Sqrt[ f[t, zs]*R[t, zs]/(Pi*Eeq) ];
   A[step, zs] = a[t, zs];
   AppendTo[HalfWidth, {t, zs, a[t, zs]}];
   AppendTo[SemiWidth, {t, zs, A[step, zs]}];
   ];

Main Code

eqRoC = {};
Load = {};
HalfWidth = {};
SemiWidth = {};
step = 1;
\[CapitalDelta]z = Tz*dz;(*接觸線上取樣中心點的間距在齒長方向(z軸)的投影,mm*)
For[ t = dt, t <= t1, t += dt,
  zmin = 0;
  zmax = rbt1*\[Omega]1*t/Tan[\[Beta]b];(*接觸線端點的齒長(z軸)坐標,mm*)
  \[Zeta]0[t] = rbt1*\[Omega]1*t/Pbt;
  AuxiLoad[];
  For[ zs = zmin, zs <= zmax, zs += \[CapitalDelta]z,
   \[Zeta][t, zs] = (rbt1*\[Omega]1*t - zs*Tan[\[Beta]b])/Pbt;
   R1[t, zs] = N1A1 + rbt1*\[Omega]1*t - zs*Tan[\[Beta]b];(*齒輪1嚙合曲率半徑,mm*)
   DynaCharAlongLoC[];
   ];
  step = step + 1;
  ];
tlast = t;
For[ t = tlast, t <= t2, t += dt,
  zmin = 0;
  zmax = B;
  \[Zeta]0[t] = rbt1*\[Omega]1*t/Pbt;
  AuxiLoad[];
  For[ zs = zmin, zs <= zmax, zs += \[CapitalDelta]z,
   \[Zeta][t, zs] = (rbt1*\[Omega]1*t - zs*Tan[\[Beta]b])/Pbt;
   R1[t, zs] = N1A1 + rbt1*\[Omega]1*t - zs*Tan[\[Beta]b];(*mm*)
   DynaCharAlongLoC[];
   ];
  step = step + 1;
  ];
tlast = t;
For[ t = tlast, t <= t3, t += dt,
  zmin = B - (Lat - rbt1*\[Omega]1*(t - t1))/Tan[\[Beta]b];
  zmax = B;
  \[Zeta]0[t] = Lat/Pbt;
  AuxiLoad[];
  For[ zs = zmin, zs <= zmax, zs += \[CapitalDelta]z,
   \[Zeta][t, zs] = (Lat - (zs - zmin)*Tan[\[Beta]b])/Pbt;
   R1[t, zs] = N1A1 + Lat - (zs - zmin)*Tan[\[Beta]b];(*mm*)
   DynaCharAlongLoC[];
   ];
  step = step + 1;
  ];
tlast = t;

取樣截線(RealSampleLine)

Constant

Tan[\[Beta]b] < Lat/B

True

結論：嚙合類型為次臨界嚙合(Subcritical Meshing)

t1 = B*Tan[\[Beta]b]/(rbt1*\[Omega]1);
t2 = Lat/(rbt1*\[Omega]1);
t3 = (B*Tan[\[Beta]b] + Lat)/(rbt1*\[Omega]1);

Gear1(RealSampleLine)

Subroutine(Gear1)

LineFunction[x1_, y1_, z1_, x2_, y2_, 
   z2_] := {x1 + (x2 - x1)*(z - z1)/(z2 - z1), 
   y1 + (y2 - y1)*(z - z1)/(z2 - z1), z};
(*Reishauer參數形式的空間節點坐標。(r1,z1),(r2,z2)為兩個已知參數點，z為此直線上其他任意點之參數坐標*)
ReToSpa[r1_, z1_, r2_, z2_, z_] :=
  {xG1[r1, z1] + (xG1[r2, z2] - xG1[r1, z1])*(z - z1)/(z2 - z1), 
   yG1[r1, z1] + (yG1[r2, z2] - yG1[r1, z1])*(z - z1)/(z2 - z1), z};
(*將空間坐標轉換為Reishauer參數坐標*)
SpaToRe[x_, y_, z_] := {Sqrt[x^2 + y^2], z};
(*在(r,z)=(r0,z0)處，對變數r的離散偏微分的函數值*)
xG1r[r0_, z0_] := (xG1[r0 + 10^-6, z0] - xG1[r0 - 10^-6, z0])/(2*10^-6);
yG1r[r0_, z0_] := (yG1[r0 + 10^-6, z0] - yG1[r0 - 10^-6, z0])/(2*10^-6);
(*在(r,z)=(r0,z0)處，對變數z的離散偏微分的函數值*)
xG1z[r0_, z0_] := (xG1[r0, z0 + 10^-6] - xG1[r0, z0 - 10^-6])/(2*10^-6);
yG1z[r0_, z0_] := (yG1[r0, z0 + 10^-6] - yG1[r0, z0 - 10^-6])/(2*10^-6);
(*在(r0,z0)處的曲面法向量*)
nG1[r0_, z0_] := {yG1r[r0, z0], -xG1r[r0, z0], 
   xG1r[r0, z0]*yG1z[r0, z0] - xG1z[r0, z0]*yG1r[r0, z0]};
(*在參考端面的作用角*)
\[Theta]1t[t_] := \[Theta]1min + \[Omega]1*t;
(*空間中兩點直線方程*)

TakeSampleG1[Type_, rat_, Amp_, \[CapitalDelta]z_, rcR_, zcR_, rcL_, zcL_] := 
  Module[{},
   zsvalue1[step] = {};
   K1 = {xG1[rcR, zcR], yG1[rcR, zcR], zcR};(*接觸線右端點*)
   K2 = {xG1[rcL, zcL], yG1[rcL, zcL], zcL};(*接觸線左端點*)
   
   For[zs = zcR, zs <= zcL, zs += \[CapitalDelta]z,(*接觸線上的位置參數zs(zSample)*)
    (*定義取樣截線：以幾何取樣中心(Geometrical Sample Center)為中點，垂直於接觸線之空間直線。
    利用接觸線兩端點，得接觸線直線方程，進而得幾何取樣中心之空間坐標Gsc*)
    Gsc = ReToSpa[rcR, zcR, rcL, zcL, zs];
    rs = Sqrt[Gsc[[1]]^2 + Gsc[[2]]^2];(*幾何取樣中心之參數坐標中的半徑坐標rs*)
    VecAB = Cross[nG1[rs, zs], K2 - K1];(*利用該點之平面法向量nG1，與接觸線方向向量(K2-
    K1)，求得取樣線方向向量VecAB*)
    u = Amp*a[t, zs]/Norm[VecAB];(*求解取樣截線AB所在直線之方程中的參數u，a[t,zs]為接觸半寬函數*)
    SpaA = Gsc + u*VecAB;(*上端點空間坐標*)   
    SpaB = Gsc - u*VecAB;(*下端點空間坐標*)
    (*轉換為參數坐標*)
    {rA, zA} = SpaToRe[  SpaA[[1]], SpaA[[2]], SpaA[[3]]  ];  
    {rB, zB} = SpaToRe[  SpaB[[1]], SpaB[[2]], SpaB[[3]]  ];  
    (*節點化*)
    rmin = ReishauerRaius[Min[rA, rB], rat];
    rmax = ReishauerRaius[Max[rA, rB], rat];
    zmin = ReishauerZ[Min[zA, zB]];
    zmax = ReishauerZ[Max[zA, zB]];
    (*判定該取樣截線是否有效：若取樣截線端點超出齒面範圍，則此點作廢；若未超出，則記錄此step下，參數zs的取值*)
    If[rmin < Min[rcR, rcL] \[Or] zmin < Min[zcR, zcL] \[Or] 
      rmax > Max[rcR, rcL] \[Or] zmax > Max[zcR, zcL],  Continue[] ];
    AppendTo[zsvalue1[step], zs];
    appz = ReishauerZ[zs];(*取樣中心的近似z坐標*)
    SpaA = {xG1[rmax, appz], yG1[rmax, appz], appz};
    SpaB = {xG1[rmin, appz], yG1[rmin, appz], appz};(*重新定義取樣截線端點*)
    (*以上端點SpaA作為坐標原點，並以取樣截線AB作為x軸，以曲面法向量方向為z軸，建立新坐標系。在取樣截線上取點，z坐標近似等於取樣中心之值zs。\
i: r On Sample Line*)
    numz = Tz - Round[Mod[appz/dz, Tz]];(*無因次齒長矩陣坐標*)
    RealSampleLineG1[step, zs] = {};
    prex = -1;(*在取樣節點賦值之前，先對輔助變數賦值，以確保第一個點不遺漏*)
    For[i = rmax, i >= rmin, i -= dr,(*以點SpaA為起點，SpaB為終點設置取樣節點NodeSi*)
     numr = Floor[(rat - i)/dr] + 1;(*無因次半徑坐標*)
     If[numr > Length[tdata[Type]], numr = Length[tdata[Type]]];
     If[numr < 1, numr = 1];
     Nd = tdata[Type][[numr, numz]]*10^-3;(*法向偏距，單位：mm*)
     NodeSi = {xG1[i, appz], yG1[i, appz], appz};(*NodeSi的空間坐標*)
     xsi = Norm[NodeSi - SpaA];(*取樣節點在新坐標系中橫坐標x的值*)
     If[xsi == prex, Continue[]];(*若本次賦值的x坐標和上次相同，則跳過*)
     AppendTo[RealSampleLineG1[step, zs], {xsi, Nd}];
     prex = xsi;(*記錄本次x坐標*)
     ];
    ];
   ];

Main code : Take Sample on Sample Line(Gear1)

Amp = 3;(*接觸寬度放大係數*)
step = 1;
\[CapitalDelta]z = dz*Tz;
For[t = dt, t <= t1, t += dt,
  {rcR, zcR} = {rbt1*Sqrt[1 + \[Theta]1t[t]^2], 0};
  {rcL, zcL} = {rclt1, rbt1*\[Omega]1*t/Tan[\[Beta]b]};
  TakeSampleG1[TypeG1, rat1, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];
tlast = t;
For[t = tlast, t <= t2, t += dt,
  {rcR, zcR} = {rbt1*Sqrt[1 + \[Theta]1t[t]^2], 0};
  {rcL, zcL} = {rbt1*Sqrt[1 + \[Theta]1t[t - t1]^2], B};
  TakeSampleG1[TypeG1, rat1, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];
tlast = t;
For[t = tlast, t <= t3, t += dt,
  {rcR, zcR} = {rat1, B - (Lat - rbt1*\[Omega]1*(t - t1))/Tan[\[Beta]b]};
  {rcL, zcL} = {rbt1*Sqrt[1 + \[Theta]1t[t - t1]^2], B};
  TakeSampleG1[TypeG1, rat1, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];

Gear2(RealSampleLine)

Subroutine(Gear2)

xG2r[r0_, z0_] := (xG2[r0 + 10^-6, z0] - xG2[r0 - 10^-6, z0])/(2*10^-6);
yG2r[r0_, z0_] := (yG2[r0 + 10^-6, z0] - yG2[r0 - 10^-6, z0])/(2*10^-6);
(*在(r,z)=(r0,z0)處，對變數z的離散偏微分的函數值*)
xG2z[r0_, z0_] := (xG2[r0, z0 + 10^-6] - xG2[r0, z0 - 10^-6])/(2*10^-6);
yG2z[r0_, z0_] := (yG2[r0, z0 + 10^-6] - yG2[r0, z0 - 10^-6])/(2*10^-6);
(*在(r0,z0)處的曲面法向量*)
nG2[r0_, z0_] := {yG2r[r0, z0], -xG2r[r0, z0], 
   xG2r[r0, z0]*yG2z[r0, z0] - xG2z[r0, z0]*yG2r[r0, z0]};
\[Theta]2t[t_] := \[Theta]2max - \[Omega]2*t;
(*空間中兩點直線方程*)

TakeSampleG2[Type_, rat_, Amp_, \[CapitalDelta]z_, rcR_, zcR_, rcL_, zcL_] := 
  Module[{},
   zsvalue2[step] = {};
   K1 = {xG2[rcR, zcR], yG2[rcR, zcR], zcR};(*接觸線右端點*)
   K2 = {xG2[rcL, zcL], yG2[rcL, zcL], zcL};(*接觸線左端點*)
   
   For[zs = zcR, zs <= zcL, zs += \[CapitalDelta]z,
    Gsc = ReToSpa[rcR, zcR, rcL, zcL, zs];
    rs = Sqrt[Gsc[[1]]^2 + Gsc[[2]]^2];
    VecAB = Cross[nG2[rs, zs], K2 - K1];
    u = Amp*a[t, zs]/Norm[VecAB];
    SpaA = Gsc + u*VecAB;(*上端點空間坐標*)   
    SpaB = Gsc - u*VecAB;(*下端點空間坐標*)
    {rA, zA} = SpaToRe[  SpaA[[1]], SpaA[[2]], SpaA[[3]]  ];  
    {rB, zB} = SpaToRe[  SpaB[[1]], SpaB[[2]], SpaB[[3]]  ];  
    rmin = ReishauerRaius[Min[rA, rB], rat];
    rmax = ReishauerRaius[Max[rA, rB], rat];
    zmin = ReishauerZ[Min[zA, zB]];
    zmax = ReishauerZ[Max[zA, zB]];
    If[rmin < Min[rcR, rcL] \[Or] zmin < Min[zcR, zcL] \[Or] 
      rmax > Max[rcR, rcL] \[Or] zmax > Max[zcR, zcL],  Continue[] ];
    AppendTo[zsvalue2[step], zs];
    appz = ReishauerZ[zs];
    SpaA = {xG2[rmax, appz], yG2[rmax, appz], appz};
    SpaB = {xG2[rmin, appz], yG2[rmin, appz], appz};
    numz = Tz - Round[Mod[appz/dz, Tz]];
    RealSampleLineG2[step, zs] = {};
    prex = -1;
    (*以點SpaB為起點，SpaA為終點設置取樣節點NodeSi*)
    For[i = rmin, i <= rmax, i += dr,
     numr = Floor[(rat - i)/dr] + 1;
     If[numr > Length[tdata[Type]], numr = Length[tdata[Type]]];
     If[numr < 1, numr = 1];
     Nd = tdata[Type][[numr, numz]]*10^-3;
     NodeSi = {xG2[i, appz], yG2[i, appz], appz};
     xsi = Norm[NodeSi - SpaB];
     If[xsi == prex, Continue[]];
     AppendTo[RealSampleLineG2[step, zs], {xsi, Nd}];
     prex = xsi;
     ];
    ];
   ];

Main code : Take Sample on Sample Line(Gear2)

Amp = 3;(*接觸寬度放大係數*)
step = 1;
\[CapitalDelta]z = dz*Tz;
For[t = dt, t <= t1, t += dt,
  {rcR, zcR} = {rbt2*Sqrt[1 + \[Theta]2t[t]^2], 0};
  {rcL, zcL} = {rat2, rbt2*\[Omega]2*t/Tan[\[Beta]b]};
  TakeSampleG2[TypeG2, rat2, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];
tlast = t;
For[t = tlast, t <= t2, t += dt,
  {rcR, zcR} = {rbt2*Sqrt[1 + \[Theta]2t[t]^2], 0};
  {rcL, zcL} = {rbt2*Sqrt[1 + \[Theta]2t[t - t1]^2], B};
  TakeSampleG2[TypeG2, rat2, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];
tlast = t;
For[t = tlast, t <= t3, t += dt,
  {rcR, zcR} = {rclt2, B - (Lat - rbt2*\[Omega]2*(t - t1))/Tan[\[Beta]b]};
  {rcL, zcL} = {rbt2*Sqrt[1 + \[Theta]2t[t - t1]^2], B};
  TakeSampleG2[TypeG2, rat2, Amp, \[CapitalDelta]z, rcR, zcR, rcL, zcL];
  step = step + 1;
  ];

Topography  Parameters

均方根粗糙度(RMS roughness)

Constant

\[CapitalDelta]z = dz*Tz;

Clear[m1, m2, Ra1, Ra2, \[Sigma]1, \[Sigma]2, \[Sigma]]

Equivalent Rough Surface(RMS roughness, mm)

test1 = {};
test2 = {};
test3 = {};
For[step = 1, step <= stepnum, step++,
  If[ zsvalue1[step] == {} \[Or] zsvalue2[step] == {}, 
   Continue[] ];(*若此步中不包含有效取樣中心，則跳過此步*)
  minzs = Max[  zsvalue1[step][[1]], zsvalue2[step][[1]]  ];(*zs的最小值，取較大者*)
  zs1 = Length[ zsvalue1[step] ];
  zs2 = Length[ zsvalue2[step] ];
  If[  zs1 > zs2, 
   maxzs = zsvalue2[step][[zs2]], 
   maxzs = zsvalue1[step][[zs1]]  ];(*zs的最大值，取較小者*)
  
  For[zs = minzs, zs <= maxzs , zs += \[CapitalDelta]z,
   (*第step步，取樣中心z坐標等於zs之取樣截面上的節點數Nodenum[step,zs]*)
   NodeNum[step, zs] = 
    Min[    Length[ RealSampleLineG1[step, zs] ], 
     Length[ RealSampleLineG2[step, zs] ]    ];
   num = NodeNum[step, zs];
   (*平均線高度，mm*)
   m1 = Mean[  Table[RealSampleLineG1[step, zs][[i, 2]], {  i, 1, num  }]  ];
   m2 = Mean[  Table[RealSampleLineG2[step, zs][[i, 2]], {  i, 1, num  }]  ];
   (*平均線粗糙度，mm*)
   Ra1[step, zs] = 
    Mean[  Table[RealSampleLineG1[step, zs][[i, 2]] - m1, {  i, 1, num  }]  ];
   Ra2[step, zs] = 
    Mean[  Table[RealSampleLineG2[step, zs][[i, 2]] - m2, {  i, 1, num  }]  ];
   (*均方根粗糙度，mm*)
   \[Sigma]1[step, zs] = 
    Sqrt[Sum[(RealSampleLineG1[step, zs][[i, 2]] - m1)^2, {  i, 1, num  }]/
      num];
   \[Sigma]2[step, zs] = 
    Sqrt[Sum[(RealSampleLineG1[step, zs][[i, 2]] - m2)^2, {  i, 1, num  }]/
      num];
   (*等效曲面的均方根粗糙度，mm*)
   \[Sigma][step, zs] = Sqrt[\[Sigma]1[step, zs]^2 + \[Sigma]2[step, zs]^2];
   AppendTo[test1, {step, zs, \[Sigma]1[step, zs]}];
   AppendTo[test2, {step, zs, \[Sigma]2[step, zs]}];
   AppendTo[test3, {step, zs, \[Sigma][step, zs]}];
   ];
  
  ];

峰點分布密度(Density of Asperity)

Equivalent Rough Surface(Density of Asperity)

(*兩基準齒面之間距(Distance Between Datum Surface)，單位：mm*)
DBDS = 0;
\[CapitalDelta]z = Tz*dz;

matDA = {};(*頂點密度的值*)
matDAdis = {};(*頂點密度的時空分佈*)
For[step = 1, step <= stepnum, step++,
  If[ zsvalue1[step] == {} \[Or] zsvalue2[step] == {}, Continue[] ];
  minzs = Max[  zsvalue1[step][[1]], zsvalue2[step][[1]]  ];
  zs1 = Length[ zsvalue1[step] ];
  zs2 = Length[ zsvalue2[step] ];
  If[  zs1 > zs2, 
   maxzs = zsvalue2[step][[zs2]], 
   maxzs = zsvalue1[step][[zs1]]  ];
  
  For[zs = minzs, zs <= maxzs , zs += \[CapitalDelta]z,
   (*Approximate Sample Line*)
   AppSL1 = Interpolation[RealSampleLineG1[step, zs]];
   AppSL2 = Interpolation[RealSampleLineG2[step, zs]];
   X1max = RealSampleLineG1[step, zs][[Length[RealSampleLineG1[step, zs]], 1]];
   X2max = RealSampleLineG2[step, zs][[Length[RealSampleLineG2[step, zs]], 1]];
   Xmax = Min[X1max, X2max];(*最右端取樣節點x值*)
   dX1 = X1max/(Length[RealSampleLineG1[step, zs]] - 1);
   dX2 = X2max/(Length[RealSampleLineG2[step, zs]] - 1);
   dX = Min[dX1, dX2];(*節點間隔*)
   num = Floor[Xmax/dX];(*節點數-1*)
   AsperityNum = 0;
   
   For[i = 1, i <= num - 1, i++,
    z0 = AppSL1[(i - 1)*dX] + AppSL2[(i - 1)*dX] - DBDS;(*等效曲面的左鄰點高度*)
    z1 = AppSL1[i*dX] + AppSL2[i*dX] - DBDS;(*等效曲面的取樣點高度*)
    z2 = AppSL1[(i + 1)*dX] + AppSL2[(i + 1)*dX] - DBDS;(*等效曲面的左鄰點高度*)
    If[z1 > z0 \[And] z1 > z2, AsperityNum = AsperityNum + 1];
    (*若取樣點高度等於左鄰點，則擴大比較範圍*)
    If[  z1 == z0,
     If[z0 > zpre \[And] z1 > z2, AsperityNum = AsperityNum + 1]   ];
    (*記錄本次循環左鄰點高度*)
    zpre = z0;
     ];
   
   DA[step, zs] = (AsperityNum/(2*A[step, zs]))^2;(*等效曲面的峰點密度，1/mm^2*)
   AppendTo[matDA, DA[step, zs]];
   AppendTo[matDAdis, {step, zs, DA[step, zs]}];
   ];
  
  ];

峰點半徑

Subroutine

(*三點法求頂點圓半徑\[Rho]i*)
Rhoi[x0_, z0_, x1_, z1_, x2_, z2_] := Module[{},
   sol = {dd, ee, ff} /. Solve[x0^2 + z0^2 + dd*x0 + ee*z0 + ff == 0
       && x1^2 + z1^2 + dd*x1 + ee*z1 + ff == 0
       && x2^2 + z2^2 + dd*x2 + ee*z2 + ff == 0, {dd, ee, ff}, Reals];
   d1 = sol[[1, 1]];
   e1 = sol[[1, 2]];
   f1 = sol[[1, 3]];
   0.5*(d1^2 + e1^2 - 4*f1)^0.5
   ];

Gear1

testRho1 = {};
For[step = 1, step <= stepnum, step++,(*循環1：確定時間步*)
  If[ zsvalue1[step] == {}, Continue[] ];
  minzs = zsvalue1[step][[1]];
  zs1 = Length[ zsvalue1[step] ];
  maxzs = zsvalue1[step][[zs1]];
  
  For[zs = minzs, zs <= maxzs , zs += \[CapitalDelta]z,(*循環2：確定取樣中心*)
   num = Length[ RealSampleLineG1[step, zs] ];
   AsperityNum = 0;
   \[Rho] = 0;
   For[i = 2, i <= num - 1, i++,(*循環3：在取樣線上取樣*)
    x0 = RealSampleLineG1[step, zs][[i - 1, 1]];(*左鄰點橫坐標*)
    x1 = RealSampleLineG1[step, zs][[i, 1]];(*取樣節點橫坐標*)
    x2 = RealSampleLineG1[step, zs][[i + 1, 1]];(*右鄰點橫坐標*)
    z0 = RealSampleLineG1[step, zs][[i - 1, 2]];(*左鄰點縱坐標*)
    z1 = RealSampleLineG1[step, zs][[i, 2]];(*取樣節點縱坐標*)
    z2 = RealSampleLineG1[step, zs][[i + 1, 2]];(*右鄰點縱坐標*)
    If[  z1 > z0 \[And] z1 > z2,
     AsperityNum = AsperityNum + 1;(*峰點數累加*)
     \[Rho] = \[Rho] + Rhoi[x0, z0, x1, z1, x2, z2](*半徑累加*)  ];
    If[  z1 == z0,
     If[z0 > zpre \[And] z1 > z2,
      AsperityNum = AsperityNum + 1;(*峰點數累加*)
      \[Rho] = \[Rho] + Rhoi[x0, z0, x1, z1, x2, z2],(*半徑累加*)
      ]  ];
    zpre = z0;
     ];
   If[  AsperityNum == 0,
    Rho1[step, zs] = 0.064,
    Rho1[step, zs] = \[Rho]/AsperityNum  ];
   (*Average radius of asperity tips of each Sample Line of G1, 1/mm*)
   AppendTo[testRho1, Rho1[step, zs]];
   ];
  
  ];

Gear2

testRho2 = {};
For[step = 1, step <= stepnum, step++,(*循環1：確定時間步*)
  If[ zsvalue2[step] == {}, Continue[] ];
  minzs = zsvalue2[step][[1]];(*取樣中心最小zs*)
  zs2 = Length[ zsvalue2[step] ];
  maxzs = zsvalue2[step][[zs2]];(*取樣中心最大zs*)
  
  For[zs = minzs, zs <= maxzs , zs += \[CapitalDelta]z,(*循環2：確定取樣中心*)
   num = Length[ RealSampleLineG2[step, zs] ];
   AsperityNum = 0;
   \[Rho] = 0;
   For[i = 2, i <= num - 1, i++,(*循環3：在取樣線上取樣*)
    x0 = RealSampleLineG2[step, zs][[i - 1, 1]];(*左鄰點橫坐標*)
    x1 = RealSampleLineG2[step, zs][[i, 1]];(*取樣節點橫坐標*)
    x2 = RealSampleLineG2[step, zs][[i + 1, 1]];(*右鄰點橫坐標*)
    z0 = RealSampleLineG2[step, zs][[i - 1, 2]];(*左鄰點縱坐標*)
    z1 = RealSampleLineG2[step, zs][[i, 2]];(*取樣節點縱坐標*)
    z2 = RealSampleLineG2[step, zs][[i + 1, 2]];(*右鄰點縱坐標*)
    If[  z1 > z0 \[And] z1 > z2,
     AsperityNum = AsperityNum + 1;(*峰點數累加*)
     \[Rho] = \[Rho] + Rhoi[x0, z0, x1, z1, x2, z2](*半徑累加*)  ];
    If[  z1 == z0,
     If[z0 > zpre \[And] z1 > z2,
      AsperityNum = AsperityNum + 1;(*峰點數累加*)
      \[Rho] = \[Rho] + Rhoi[x0, z0, x1, z1, x2, z2],(*半徑累加*)
      ]  ];
    zpre = z0;(*記錄本次循環左鄰點高度*)
     ];
   (*Average radius of asperity tips of each Sample Line, 1/mm*)
   If[  AsperityNum == 0,
    Rho2[step, zs] = 0.064,
    Rho2[step, zs] = \[Rho]/AsperityNum ];
   AppendTo[testRho2, Rho2[step, zs]];
   ];
  
  ];

Equivalent Rough Surface

matRho = {};(*等效曲面平均峰點半徑數值*)
matRhoDis = {};(*等效曲面平均峰點半徑之時空分佈*)
For[step = 1, step <= stepnum, step++,
  If[ zsvalue1[step] == {} \[Or] zsvalue2[step] == {}, Continue[] ];
  minzs = Max[  zsvalue1[step][[1]], zsvalue2[step][[1]]  ];
  zs1 = Length[ zsvalue1[step] ];
  zs2 = Length[ zsvalue2[step] ];
  If[  zs1 > zs2, 
   maxzs = zsvalue2[step][[zs2]], 
   maxzs = zsvalue1[step][[zs1]]  ];
  
  For[zs = minzs, zs <= maxzs , zs += \[CapitalDelta]z,
   Rho[step, zs] = 
    Rho1[step, zs]*Rho2[step, zs]/(Rho1[step, zs] + Rho2[step, zs]);(*mm*)
   AppendTo[matRho, Rho[step, zs]];
   AppendTo[matRhoDis, {step, zs, Rho[step, zs]}];
   ];
  
  ];

聯立求解壓力方程

Subroutine and Constant

F52[H_] := If[ H < 4,  4.4086*10^-5*(4 - H)^6.804  ,  0];
\[Eta]0 = 0.08;(*  在一定工作溫度下的粘度，單位：Pa*s  *)
mEeq = 2.2831*10^11;(*equivalent elastic module,N/m^2*)
\[Alpha]EHL = 2.19*10^-8;(*,單位：(Pa^-1)*)
z = 0.6;
\[CapitalDelta]z = Tz*dz;
fC = 0.1;

Solving Pressure Equation

hcvalue = {};
hcmat = {};
\[Mu]value = {};
\[Mu]mat = {};
For[step = 1, step <= stepnum, step++,
  If[ zsvalue1[step] == {} \[Or] zsvalue2[step] == {}, Continue[] ];
  minzs = Max[  zsvalue1[step][[1]], zsvalue2[step][[1]]  ];
  zs1 = Length[ zsvalue1[step] ];
  zs2 = Length[ zsvalue2[step] ];
  If[  zs1 > zs2, 
   maxzs = zsvalue2[step][[zs2]], 
   maxzs = zsvalue1[step][[zs1]]  ];
  
  For[zs = minzs, zs <= maxzs, zs += \[CapitalDelta]z,
   (*單位轉換，mm\[Rule]m*)
   mU = Ustep[step, zs]*10^-3;(*滾動速度, m*)
   mUs = Usstep[step, zs]*10^-3;(*滑動速度, m/s*)
   mR = Rstep[step, zs]*10^-3;(*等效嚙合曲率半徑, m*)
   mph = phstep[step, zs]*10^6;(*赫茲接觸壓力, N/m^2即Pa*)
   mf = fstep[step, zs]*10^3;(*負載, N/m*)
   mA = A[step, zs]*10^-3;(*接觸半寬, m*)
   Usig = \[Eta]0*mU/(mEeq*mR);
   W = mf/(mEeq*mR);
   G = Round[\[Alpha]EHL*mEeq];
   M = W*Usig^-0.5;
   L = G*Usig^0.25;
   RI = N[ 3/M, 16 ];
   EI = N[ 2.621*M^-0.2, 16 ];
   RP = N[ 1.287*L^(2/3), 16  ];
   EP = N[ 1.311*M^(-1/8)*L^(3/4), 16  ];
   s = (7 + 8*Exp[-2*\[Gamma]1^-0.4*EI/RI])/5;
   hc = (   \[Gamma]1^(s/2)*(   RI^(7/3) + \[Gamma]1^(-14/15)*EI^(7/3) )^(
        3*s/7)  +   \[Gamma]1^(-s/2)*(  RP^(-7/2)  +  EP^(-7/2)  )^(-2*
         s/7)  )^(1/s)*  \[Gamma]1^(1/2)*mR*Usig^0.5;
   (*單位轉換，mm\[Rule]m*)
   m\[Sigma] = \[Sigma][step, zs]*10^-3;(*等效RMS, m*)
   mDA = DA[step, zs]*10^6;(*等效曲面頂點分佈密度, 1/m^2*)
   mRho = Rho[step, zs]*10^-3;(*等效曲面平均頂點半徑, m*)
   a1 = 1.558;  a2 = 0.0337;  a3 = -0.442;  a4 = -1.70;(*無因次係數*)
   pc1 = (1 - 1/\[Gamma]1)*
     mph*(  1 + ( 
        a1*(mDA*\[Gamma]1/(\[Gamma]1 - 1)*mR^1.5*mRho^0.5)^a2*(m\[Sigma]/mR)^
         a3*W^(a2 - a3) )^a4  )^(1/a4);(*壓力恆等式左端*)
   pc2 = 8*2^0.5/15*Pi*mDA^2*mRho^1.5*m\[Sigma]^2.5*mEeq*
     F52[(hc - 1.15*m\[Sigma])/m\[Sigma]];(*壓力恆等式右端*)
   (*二分法解壓力方程*)
   solA = 1.000001;
   solB = 10^6;
   error = 10^-3;
   For[i = 1, i <= 1000, i++,
    EqpcA = (pc1 - pc2) /. \[Gamma]1 -> solA;
    EqpcB = (pc1 - pc2) /. \[Gamma]1 -> solB;
    EqpcMid = (pc1 - pc2) /. \[Gamma]1 -> (solA + solB)/2;
    If[Abs[EqpcMid] < error, Break[]];
    If[EqpcA*EqpcMid < 0, solB = (solA + solB)/2, solA = (solA + solB)/2];
    ];
   If[ i > 1000, solA = 1.49632; solB = 1.50001 ];(*防呆*)
   gamma1[step, zs] = (solA + solB)/2;(*濕摩擦比例因子*)
   gamma2[step, zs] = gamma1[step, zs]/(gamma1[step, zs] - 1);(*乾摩擦比例因子*)
   SOLhc[step, zs] = (hc /. \[Gamma]1 -> gamma1[step, zs])*10^6;(*中心膜厚，\[Mu]m*)

   
   pc[step, zs] = pc1 /. \[Gamma]1 -> (solA + solB)/2;(*乾摩擦壓力，Pa*)
   phyd[step, zs] = mph - pc[step, zs];(*濕摩擦壓力，Pa*)
   \[Eta] = \[Eta]0*
     Exp[(Log[\[Eta]0] + 9.67)*((1 + 5.1*10^-9*phyd[step, zs])^z - 
         1)];(*黏度，Pa*s*)
   FfH[step, zs] = 
    2*mA*\[Eta]*mUs/SOLhc[step, zs]*\[CapitalDelta]z;(*濕摩擦力，N*)
   
   FfC[step, zs] = fC*FT/gamma2[step, zs];(*乾摩擦力，N*)
   \[Mu][step, zs] = (FfH[step, zs] + FfC[step, zs])/FT;(*摩擦係數，無因次*)
   AppendTo[hcvalue, SOLhc[step, zs]];
   AppendTo[hcmat, {step, zs, SOLhc[step, zs]}];
   AppendTo[\[Mu]value, \[Mu][step, zs]];
   AppendTo[\[Mu]mat, {step, zs, \[Mu][step, zs]}];
   ];
  ];


輸出

fighcmat = 
 ListPlot3D[hcmat, 
  ColorFunction -> Function[{x, y, z}, ColorData["Rainbow"][z]], 
  PlotRange -> All, BoxRatios -> {1, 1, 0.2}, Boxed -> False, 
  AxesStyle -> Thick, 
  LabelStyle -> {FontFamily -> "Times New Roman", 12, GrayLevel[0]}]
Export["A-A hc.wmf", fighcmat];
colorrange = {Floor[Min[hcvalue]], Ceiling[Max[hcvalue]]};
colorbar = 
 BarLegend[{"Rainbow", colorrange}, 50, LegendLayout -> "Row", 
  LabelStyle -> {FontFamily -> "Times New Roman", 12, GrayLevel[0]}]
Export["A-A hc Legend.wmf", colorbar];
(*平均膜厚，\[Mu]m*)
Mean[hcvalue]

