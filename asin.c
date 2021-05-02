#define MLA mla
#define C2V(x) (x)
#define mla(a,b,c) (a*b+c)

#define POLY2(x, c1, c0) MLA(x, C2V(c1), C2V(c0))
#define POLY3(x, x2, c2, c1, c0) MLA(x2, C2V(c2), MLA(x, C2V(c1), C2V(c0)))
#define POLY4(x, x2, c3, c2, c1, c0) MLA(x2, MLA(x, C2V(c3), C2V(c2)), MLA(x, C2V(c1), C2V(c0)))
#define POLY5(x, x2, x4, c4, c3, c2, c1, c0) MLA(x4, C2V(c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY6(x, x2, x4, c5, c4, c3, c2, c1, c0) MLA(x4, POLY2(x, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY7(x, x2, x4, c6, c5, c4, c3, c2, c1, c0) MLA(x4, POLY3(x, x2, c6, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0) MLA(x4, POLY4(x, x2, c7, c6, c5, c4), POLY4(x, x2, c3, c2, c1, c0))
#define POLY9(x, x2, x4, x8, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, C2V(c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY10(x, x2, x4, x8, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY2(x, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY11(x, x2, x4, x8, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY3(x, x2, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY12(x, x2, x4, x8, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY4(x, x2, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
  

double xasin(double d) 
{
  int o = fabsk(d) < 0.5;
  double x2 = o ? (d*d) : ((1-fabsk(d))*0.5), x = o ? fabsk(d) : SQRT(x2), u;

  double x4 = x2 * x2, x8 = x4 * x4, x16 = x8 * x8;
  u = POLY12(x2, x4, x8, x16,
	     +0.3161587650653934628e-1,
	     -0.1581918243329996643e-1,
	     +0.1929045477267910674e-1,
	     +0.6606077476277170610e-2,
	     +0.1215360525577377331e-1,
	     +0.1388715184501609218e-1,
	     +0.1735956991223614604e-1,
	     +0.2237176181932048341e-1,
	     +0.3038195928038132237e-1,
	     +0.4464285681377102438e-1,
	     +0.7500000000378581611e-1,
	     +0.1666666666666497543e+0);

  u = mla(u, x * x2, x);
  
  double r = o ? u : (M_PI/2 - 2*u);
  r = mulsign(r, d);

  return r;
}
