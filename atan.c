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
#define POLY13(x, x2, x4, x8, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY5(x, x2, x4, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY14(x, x2, x4, x8, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY6(x, x2, x4, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY15(x, x2, x4, x8, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY7(x, x2, x4, ce, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, POLY8(x, x2, x4, cf, ce, cd, cc, cb, ca, c9, c8), POLY8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY17(x, x2, x4, x8, x16, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, C2V(d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY18(x, x2, x4, x8, x16, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, POLY2(x, d1, d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define POLY19(x, x2, x4, x8, x16, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, POLY3(x, x2, d2, d1, d0), POLY16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))

double atan(double s) {
  double t, u;
  int q = 0;

  if (sign(s) == -1) { s = -s; q = 2; }
  if (s > 1) { s = 1.0 / s; q |= 1; }

  t = s * s;

  double t2 = t * t, t4 = t2 * t2, t8 = t4 * t4, t16 = t8 * t8;
  u = POLY19(t, t2, t4, t8, t16,
	     -1.88796008463073496563746e-05,
	     0.000209850076645816976906797,
	     -0.00110611831486672482563471,
	     0.00370026744188713119232403,
	     -0.00889896195887655491740809,
	     0.016599329773529201970117,
	     -0.0254517624932312641616861,
	     0.0337852580001353069993897,
	     -0.0407629191276836500001934,
	     0.0466667150077840625632675,
	     -0.0523674852303482457616113,
	     0.0587666392926673580854313,
	     -0.0666573579361080525984562,
	     0.0769219538311769618355029,
	     -0.090908995008245008229153,
	     0.111111105648261418443745,
	     -0.14285714266771329383765,
	     0.199999999996591265594148,
	     -0.333333333333311110369124);

  t = s + s * (t * u);

  if ((q & 1) != 0) t = 1.570796326794896557998982 - t;
  if ((q & 2) != 0) t = -t;

  return t;
}
