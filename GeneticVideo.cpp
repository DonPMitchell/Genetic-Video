// GeneticVideo.cpp : Defines the entry point for the console application.
//
#pragma intrinsic(sqrt, sin, cos, exp, log, fabs, atan2)

#include "stdafx.h"
#include "texture.h"
#include "Image.h"
#include "Video.h"
#include "MonteCarlo.h"

#define D_PI    3.1415926535897932384626433832795028842
#define WIDE    1920
#define HIGH    1080
#define NSMALL  6
#define HORZ_IMAGES	NSMALL
#define VERT_IMAGES	NSMALL
#define SMALL_WIDE (WIDE/HORZ_IMAGES)
#define SMALL_HIGH (HIGH/VERT_IMAGES)

#define NUM 0
#define INT 1
#define ABS 2
#define SQT 3
#define EXP 4
#define LOG 5
#define SIN 6
#define COS 7
#define INV 8
#define NEG 9
#define BLR 10
#define GRD 11
#define DIR 12
#define WRP 13
#define TRB 14
#define IFS 15

#define DIV 16
#define MOD 17
#define MIN 18
#define MAX 19
#define AND 20
#define ORR 21
#define XOR 22
#define WP2 23
#define XXX 24
#define YYY 25
#define FRC 26
#define BMP 27
#define ADD 28
#define SUB 29
#define RAD 30
#define ANG 31
#define MUL 32
#define IF2 33
#define TTT 34

#define NTYPES  35
#define NSINGLE 16

struct Expr {
    int             nType;
    int             nConst;
    DisplayRGB      rgb;
    Expr            *pn1;
    Expr            *pn2;

            Expr() : pn1(0), pn2(0), nConst(0) {}
};

struct VectorXYT {
    double x, y, t;

            VectorXYT() {}
            VectorXYT(double a, double b, double c) : x(a), y(b), t(c) {}
    VectorXYT operator+ (VectorXYT v) { return VectorXYT(x+v.x, y+v.y, t+v.t); }
    VectorXYT operator* (double f)    { return VectorXYT(x*f,   y*f,   t*f  ); }
};

inline VectorXYT
Cross(VectorXYT &vec1, VectorXYT &vec2)
{
    return VectorXYT((vec1.y)*(vec2.t)-(vec1.t)*(vec2.y),
                     (vec1.t)*(vec2.x)-(vec1.x)*(vec2.t),
                     (vec1.x)*(vec2.y)-(vec1.y)*(vec2.x));
}

inline VectorXYT
Normalize(VectorXYT &v)
{
    double f;

    f = 1.0/sqrt(v.x*v.x + v.y*v.y + v.t*v.t);
    return VectorXYT(v.x*f, v.y*f, v.t*f);
}

static double   H = 1.0/WIDE;
static unsigned  s_nMax;
static Expr     *s_eMax;
static double   s_fHi, s_fLo;

static Expr *
RandomExpression()
{
    Expr *e;
    
    e = new Expr;
    e->rgb = DisplayRGB(ML_RandomFloat(), ML_RandomFloat(), ML_RandomFloat());
    if (ML_RandomFloat() < 0.4) {
        e->nType = NUM;
    } else if (ML_RandomFloat() < 0.1) {
        e->nType = XXX;
    } else if (ML_RandomFloat() < 0.1) {
        e->nType = YYY;
    } else if (ML_RandomFloat() < 0.1) {
        e->nType = TTT;
    } else {
        e->nType = ML_RandomLessThanN(NTYPES);
        e->pn1 = RandomExpression();
        if (e->nType >= NSINGLE)
            e->pn2 = RandomExpression();
    }
    return e;
}

static void
GetRandomSubExpression(Expr *e)
{
    int n;

    if (e == 0)
        return;
    n = ML_RandomUnsigned();
    if (n > s_nMax) {
        s_nMax = n;
        s_eMax = e;
    }
    GetRandomSubExpression(e->pn1);
    GetRandomSubExpression(e->pn2);
}

static Expr *
RandomSubExpression(Expr *e)
{
    s_nMax = 0;
    s_eMax = e;
    GetRandomSubExpression(e);
    return s_eMax;
}

static Expr *
CopyExpression(Expr *e)
{
    Expr *pn;

    if (e == 0)
        return 0;
    pn = new Expr;
    *pn = *e;
    pn->pn1 = CopyExpression(e->pn1);
    pn->pn2 = CopyExpression(e->pn2);
    return pn;
}

static Expr *
MutateExpression(Expr *e)
{
    int nMutationType, n;
    Expr *eSub, *eNew;

    if (e == 0)
        return 0;
    e = CopyExpression(e);
    eSub = RandomSubExpression(e);
    if (eSub->nType == NUM) {
        eSub->rgb = DisplayRGB(ML_RandomFloat(), ML_RandomFloat(), ML_RandomFloat());
        return e;
    }
    if (eSub->nType == XXX || eSub->nType == YYY || eSub->nType == TTT)
        return e;
    nMutationType = ML_RandomLessThanN(4);
    switch (nMutationType) {
            //
            //  Replace a random subexpression
            //
    case 0: eNew = RandomExpression();
            *eSub = *eNew;
            return e;
            //
            //  Randomly change an op code
            //
    case 1: if (eSub->nType < NSINGLE)
                eSub->nType = ML_RandomLessThanN(NSINGLE);
            else
                eSub->nType = NSINGLE + ML_RandomLessThanN(NTYPES - NSINGLE);
            return e;
            //
            //  Make e into a subexpression of a new operation
            //
    case 2: eNew = new Expr;
            do {
                n = ML_RandomLessThanN(NTYPES);
            } while (n != NUM && n != XXX && n != YYY && n != TTT);
            eNew->nType = n;
            eNew->pn1 = CopyExpression(eSub);
            if (n >= NSINGLE)
                eNew->pn2 = RandomExpression();
            *eSub = *eNew;
            return e;
            //
            //  Replace expression with a subexpression
            //
    case 3: if (eSub->nType < NSINGLE || ML_RandomFloat() < 0.5)
                *eSub = *eSub->pn1;
            else
                *eSub = *eSub->pn2;
            return e;
    }
    return e;
}

static Expr*
BreedExpressions(Expr *e1, Expr *e2)
{
    Expr *eSub1, *eSub2;

    if (ML_RandomFloat() < 0.33333)
        e1 = MutateExpression(e1);
    if (ML_RandomFloat() < 0.33333)
        e2 = MutateExpression(e2);
    e1 = CopyExpression(e1);
    eSub1 = RandomSubExpression(e1);
    eSub2 = RandomSubExpression(e2);
    *eSub1 = *eSub2;
    return e1;
}

static int s_nConst= 0;
static int s_nConstCount = 1;

static DisplayRGB
Eval(Expr *e, double x, double y, double t)
{
	DisplayRGB a, b, c, d, g1, g2, g3, g4, g5, g6;
	union {
		long l;
		float f;
	} u, v;
    float rgf[3];
    double r, theta;
    int nConst;

    s_nConst = 0;
    nConst = 0;
    if (e->nConst == s_nConstCount) {               // constant expression value
        s_nConst = s_nConstCount;
        return e->rgb;
    }
	if (e->nType == NUM) {
        s_nConst = s_nConstCount;
		return e->rgb;
    }
	if (e->nType == XXX)
		return DisplayRGB(x);
	if (e->nType == YYY)
		return DisplayRGB(y);
    if (e->nType == TTT)
        return DisplayRGB(t);
    if (e->nType != WRP && e->nType != TRB) {
	    a = Eval(e->pn1, x, y, t);
        if (s_nConst == s_nConstCount)
            nConst = s_nConst;
	    if (e->nType >= NSINGLE) {
		    b = Eval(e->pn2, x, y, t);
            if (nConst && (s_nConst == s_nConstCount))
                nConst = s_nConst;
            else
                nConst = 0;
        }
    }
    c = ML_Black;
	switch (e->nType) {

	case SQT:
		c.red = sqrt(fabs(a.red));
		c.grn = sqrt(fabs(a.grn));
		c.blu = sqrt(fabs(a.blu));
        e->nConst = s_nConst = nConst; e->rgb = c;
		return c;
	case SIN:
		c.red = sin(10.0*D_PI*a.red);
		c.grn = sin(10.0*D_PI*a.grn);
		c.blu = sin(10.0*D_PI*a.blu);
		return c;
		c.red = sin(a.red);
		c.grn = sin(a.grn);
		c.blu = sin(a.blu);
		return c;
	case COS:
		c = c * cos(10.0*D_PI*x);
		return c;
		c.red = cos(a.red);
		c.grn = cos(a.grn);
		c.blu = cos(a.blu);
		return c;
	case EXP:
        r = hypot(x, y);
        theta = atan2(x, y);
        c = Eval(e->pn1, r, theta, t);
        return c;
		c.red = exp(a.red);
		c.grn = exp(a.grn);
		c.blu = exp(a.blu);
		return c;
	case LOG:
        r = x;
        theta = y;
        c = Eval(e->pn1, r*sin(theta), r*cos(theta), t);
        return c;
		c.red = (a.red == 0.0) ? -10.0 : log(fabs(a.red));
		c.grn = (a.grn == 0.0) ? -10.0 : log(fabs(a.grn));
		c.blu = (a.blu == 0.0) ? -10.0 : log(fabs(a.blu));
		return c;
	case ABS:
		c.red = fabs(a.red);
		c.grn = fabs(a.grn);
		c.blu = fabs(a.blu);
		return c;
	case INT:
		c.red = floor(a.red);
		c.grn = floor(a.grn);
		c.blu = floor(a.blu);
		return c;
	case NEG:
		c.red = -(a.red);
		c.grn = -(a.grn);
		c.blu = -(a.blu);
		return c;
	case INV:
		c.red = (a.red == 0.0) ? 0.0 : 1.0/(fabs(a.red));
		c.grn = (a.grn == 0.0) ? 0.0 : 1.0/(fabs(a.grn));
		c.blu = (a.blu == 0.0) ? 0.0 : 1.0/(fabs(a.blu));
		return c;
    case GRD:
        g1 = Eval(e->pn1, x + H, y, t);
        g2 = Eval(e->pn1, x - H, y, t);
        g3 = Eval(e->pn1, x, y + H, t);
        g4 = Eval(e->pn1, x, y - H, t);
        g5 = Eval(e->pn1, x, y, t + H);
        g6 = Eval(e->pn1, x, y, t - H);
        a = (g1 - g2)/(2.0*H);
        b = (g3 - g4)/(2.0*H);
        c = (g5 - g6)/(2.0*H);
        d.red = hypot(a.red, b.red);
        d.grn = hypot(b.grn, c.grn);
        d.blu = hypot(c.blu, a.blu);
        return c;
    case DIR:
        g1 = Eval(e->pn1, x + H, y, t);
        g2 = Eval(e->pn1, x - H, y, t);
        g3 = Eval(e->pn1, x, y + H, t);
        g4 = Eval(e->pn1, x, y - H, t);
        g5 = Eval(e->pn1, x, y, t + H);
        g6 = Eval(e->pn1, x, y, t - H);
        a = (g1 - g2)/(2.0*H);
        b = (g3 - g4)/(2.0*H);
        c = (g5 - g6)/(2.0*H);
        c.red = atan2(a.red, b.red)/D_PI;
        c.grn = atan2(b.grn, c.grn)/D_PI;
        c.blu = atan2(c.blu, a.blu)/D_PI;
        return c;
    case BLR:
        return Eval(e->pn1, a.blu, a.grn, a.red);       // Blur is not interesting, do more interated function
        g1 = Eval(e->pn1, x + H, y, t);
        g2 = Eval(e->pn1, x - H, y, t);
        g3 = Eval(e->pn1, x, y + H, t);
        g4 = Eval(e->pn1, x, y - H, t);
        g5 = Eval(e->pn1, x, y, t + H);
        g6 = Eval(e->pn1, x, y, t - H);
        c = (a + g1 + g2 + g3 + g4 + g5 + g6)/7.0;
        return c;
    case WRP:
        rgf[0] = x;
        rgf[1] = y;
        rgf[2] = t;
        x += 0.1*BandNoise(rgf, 3);
        c = Eval(e->pn1, x, y, t);
        return c;
    case TRB:
        rgf[0] = x;
        rgf[1] = y;
        rgf[2] = t;
        x += 0.1*FractalNoise(rgf, 3, 3, 1.0, 2.1, 1);
        c = Eval(e->pn1, x, y, t); 
        return c;
	case AND:
		u.f = a.red; v.f = b.red; u.l &= v.l; c.red = u.f;
		u.f = a.grn; v.f = b.grn; u.l &= v.l; c.grn = u.f;
		u.f = a.blu; v.f = b.blu; u.l &= v.l; c.blu = u.f;
		return c;
	case XOR:
		u.f = a.red; v.f = b.red; u.l ^= v.l; c.red = u.f;
		u.f = a.grn; v.f = b.grn; u.l ^= v.l; c.grn = u.f;
		u.f = a.blu; v.f = b.blu; u.l ^= v.l; c.blu = u.f;
		return c;
	case ORR:
		u.f = a.red; v.f = b.red; u.l |= v.l; c.red = u.f;
		u.f = a.grn; v.f = b.grn; u.l |= v.l; c.grn = u.f;
		u.f = a.blu; v.f = b.blu; u.l |= v.l; c.blu = u.f;
		return c;
    case IFS:
        c = Eval(e->pn1, a.red, a.grn, a.blu);
        return c;
	case ADD:
		c.red = a.red + b.red; c.grn = a.grn + b.grn; c.blu = a.blu + b.blu;
		return c;
	case SUB:
		c.red = a.red - b.red; c.grn = a.grn - b.grn; c.blu = a.blu - b.blu;
		return c;
	case MUL:
		c.red = a.red * b.red; c.grn = a.grn * b.grn; c.blu = a.blu * b.blu;
		return c;
	case DIV:
		c.red = (b.red == 0.0) ? 0.0 : a.red / b.red;
		c.grn = (b.grn == 0.0) ? 0.0 : a.grn / b.grn;
		c.blu = (b.blu == 0.0) ? 0.0 : a.blu / b.blu;
		return c;
	case MOD:
		c.red = (b.red == 0.0) ? 0.0 : fmod(a.red, b.red);
		c.grn = (b.grn == 0.0) ? 0.0 : fmod(a.grn, b.grn);
		c.blu = (b.blu == 0.0) ? 0.0 : fmod(a.blu, b.blu);
		return c;
	case RAD:
		c.red = hypot(a.red, b.red);
		c.grn = hypot(a.grn, b.grn);
		c.blu = hypot(a.blu, b.blu);
		return c;
	case ANG:
		c.red = (a.red == 0.0 && b.red == 0.0) ? 0.0 : atan2(a.red, b.red)/D_PI;
		c.grn = (a.grn == 0.0 && b.grn == 0.0) ? 0.0 : atan2(a.grn, b.grn)/D_PI;
		c.blu = (a.blu == 0.0 && b.blu == 0.0) ? 0.0 : atan2(a.blu, b.blu)/D_PI;
		return c;
	case MAX:
		c.red = (a.red > b.red) ? a.red : b.red;
		c.grn = (a.grn > b.grn) ? a.grn : b.grn;
		c.blu = (a.blu > b.blu) ? a.blu : b.blu;
		return c;
	case MIN:
		c.red = (a.red < b.red) ? a.red : b.red;
		c.grn = (a.grn < b.grn) ? a.grn : b.grn;
		c.blu = (a.blu < b.blu) ? a.blu : b.blu;
		return c;
    case BMP:
        rgf[0] = x; rgf[1] = y; rgf[2] = t;
        c = a + (b - a)*BandNoise(rgf, 3);
        return c;
        rgf[0] = a.red; rgf[1] = b.red; c.red = BandNoise(rgf, 2);
        rgf[0] = a.grn; rgf[1] = b.grn; c.grn = BandNoise(rgf, 2);
        rgf[0] = a.blu; rgf[1] = b.blu; c.blu = BandNoise(rgf, 2);
        return c;
    case FRC:
        rgf[0] = x; rgf[1] = y; rgf[2] = t;
        c = a + (b - a)*FractalNoise(rgf, 3);
        return c;
        rgf[0] = a.red; rgf[1] = b.red; c.red = FractalNoise(rgf, 2);
        rgf[0] = a.grn; rgf[1] = b.grn; c.grn = FractalNoise(rgf, 2);
        rgf[0] = a.blu; rgf[1] = b.blu; c.blu = FractalNoise(rgf, 2);
        return c;
    case WP2:
        rgf[0] = a.red;
        rgf[1] = a.grn;
        y += 0.1*BandNoise(rgf, 2);
        c = Eval(e->pn2, x, y, t);
        return c;
    case IF2:
        c = Eval(e->pn2, a.red+x, a.grn+y, a.blu+t);
        return c;
    }
	return c;
}

void
PrintExpr(Expr *e)
{
	if (e->nType == NUM) {
        printf("(%d,%d,%d)", int(255.0*e->rgb.red), int(255.0*e->rgb.grn), int(255.0*e->rgb.blu));
		return;
    }
	if (e->nType == XXX) {
		printf("X");
        return;
    }
	if (e->nType == YYY) {
		printf("Y");
        return;
    }
    if (e->nType == TTT) {
        printf("T");
        return;
    }
	switch (e->nType) {

	case SQT:
        printf("sqrt("); PrintExpr(e->pn1); printf(")");
		return;
	case SIN:
        printf("sin("); PrintExpr(e->pn1); printf(")");
		return;
	case COS:
        printf("cos("); PrintExpr(e->pn1); printf(")");
		return;
	case EXP:
        printf("polar("); PrintExpr(e->pn1); printf(")");
		return;
	case LOG:
        printf("invpolar("); PrintExpr(e->pn1); printf(")");
		return;
	case ABS:
        printf("abs("); PrintExpr(e->pn1); printf(")");
		return;
	case INT:
        printf("floor("); PrintExpr(e->pn1); printf(")");
		return;
	case NEG:
        printf("-("); PrintExpr(e->pn1); printf(")");
		return;
	case INV:
        printf("1.0/("); PrintExpr(e->pn1); printf(")");
		return;
    case GRD:
        printf("grad("); PrintExpr(e->pn1); printf(")");
		return;
    case DIR:
        printf("dir("); PrintExpr(e->pn1); printf(")");
		return;
    case BLR:
        printf("IFS0("); PrintExpr(e->pn1); printf(")");
		return;
    case WRP:
        printf("warp("); PrintExpr(e->pn1); printf(")");
		return;
    case TRB:
        printf("fractal("); PrintExpr(e->pn1); printf(")");
		return;
	case AND:
        printf("("); PrintExpr(e->pn1); printf(")&("); PrintExpr(e->pn2); printf(")");
        return;
	case XOR:
        printf("("); PrintExpr(e->pn1); printf(")^("); PrintExpr(e->pn2); printf(")");
        return;
	case ORR:
        printf("("); PrintExpr(e->pn1); printf(")|("); PrintExpr(e->pn2); printf(")");
        return;
    case IFS:
        printf("IFS1("); PrintExpr(e->pn1); printf(")");
		return;
	case ADD:
        printf("("); PrintExpr(e->pn1); printf(")+("); PrintExpr(e->pn2); printf(")");
        return;
	case SUB:
        printf("("); PrintExpr(e->pn1); printf(")-("); PrintExpr(e->pn2); printf(")");
        return;
	case MUL:
        printf("("); PrintExpr(e->pn1); printf(")*("); PrintExpr(e->pn2); printf(")");
        return;
	case DIV:
        printf("("); PrintExpr(e->pn1); printf(")/("); PrintExpr(e->pn2); printf(")");
        return;
	case MOD:
        printf("("); PrintExpr(e->pn1); printf(")%("); PrintExpr(e->pn2); printf(")");
        return;
	case RAD:
        printf("hypot("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
	case ANG:
        printf("atan2("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
	case MAX:
        printf("max("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
	case MIN:
        printf("min("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
    case BMP:
        printf("bump("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
    case FRC:
        printf("fractal("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
    case WP2:
        printf("warp2("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
    case IF2:
        printf("IF2("); PrintExpr(e->pn1); printf(","); PrintExpr(e->pn2); printf(")");
        return;
    }
 
}

static void
Label(int i, int j, Image &im)
{
    char sz[16];

    sprintf_s(sz, sizeof(sz), "%d %d", i, j);
    im.DrawTextRGB(ML_Black, double(i*SMALL_WIDE + 20), double(HIGH - j*SMALL_HIGH - 20), sz, 1.0);
}

static double tFrame = 0.0;

static int
Plot(int i, int j, Expr *e, Image &im)
{
	int ix, iy;
	double red, grn, blu;
	double x, y, hi, lo;
    static DisplayRGB image[SMALL_WIDE][SMALL_HIGH];

    // printf(".");
    H = 1.0/SMALL_WIDE;
	if (e->nType == XXX || e->nType == YYY || e->nType == NUM || e->nType == TTT)
		return 0;
	for (ix = -SMALL_WIDE/2; ix < SMALL_WIDE/2; ix++) {
		x = ((double)ix + ML_RandomFloat() - 0.5)/(SMALL_WIDE/2);
		for (iy = -SMALL_HIGH/2; iy < SMALL_HIGH/2; iy++) {
			y = ((double)iy + ML_RandomFloat() - 0.5)/(SMALL_HIGH/2);
			image[ix+SMALL_WIDE/2][iy+SMALL_HIGH/2] = Eval(e, x, y, tFrame);
		}
	}
	hi = lo = image[0][0].red;
	for (ix = 0; ix < SMALL_WIDE; ix++)
		for (iy = 0; iy < SMALL_HIGH; iy++) {
			if(image[ix][iy].red > hi) hi = image[ix][iy].red;
			if(image[ix][iy].grn > hi) hi = image[ix][iy].grn;
			if(image[ix][iy].blu > hi) hi = image[ix][iy].blu;
			if(image[ix][iy].red < lo) lo = image[ix][iy].red;
			if(image[ix][iy].grn < lo) lo = image[ix][iy].grn;
			if(image[ix][iy].blu < lo) lo = image[ix][iy].blu;
		}
	if (lo < -100000000.0 || hi > 100000000 || hi == lo)
			return 0;
	for (iy = 0; iy < SMALL_HIGH; iy++) {
		for (ix = 0; ix < SMALL_WIDE; ix++) {
			red = (image[ix][iy].red-lo)/(hi-lo);
			grn = (image[ix][iy].grn-lo)/(hi-lo);
			blu = (image[ix][iy].blu-lo)/(hi-lo);
			im.SetRGB(DisplayRGB(red, grn, blu), ix+SMALL_WIDE*i, iy+SMALL_HIGH*j);
		}
	}
    Label(i, j, im);
    s_fLo = lo;     // for use by BigPlot
    s_fHi = hi;
	return 1;
}

static void
VideoSamples(Expr *rgeTrials[HORZ_IMAGES][VERT_IMAGES])
{
    Image im;
    Video vid;
    int i, j;

    im.NewImage(WIDE, HIGH, 3);
    if (vid.NewVideo("VideoSample.mpeg", im) == 0) {
        printf("cannot open file\n");
        return;
    }
    for (tFrame = -1.0; tFrame <= 1.0; tFrame += 0.1) 
        printf("-");
    printf("\n");
    for (tFrame = -1.0; tFrame <= 1.0; tFrame += 0.1) {
        for (i = 0; i < HORZ_IMAGES; i++) {
            for (j = 0; j < VERT_IMAGES; j++) {
                Plot(i, j, rgeTrials[i][j], im);
            }
        }
        vid.WriteFrame(im);
        vid.WriteFrame(im);
        vid.WriteFrame(im);
        vid.WriteFrame(im);
        vid.WriteFrame(im);
        vid.WriteFrame(im);
        printf(".");
    }
    vid.Close();
}

static void
BigPlot(Expr *e, Image &im)
{
    double x, y;
    int ix, iy;
    DisplayRGB rgb;
    char szFile[64];
    static int nBig;

    tFrame = 0.0;
    Plot(0, 0, e, im);
    //H = 1.0/(IMAGE_VERT*VERT_IMAGES);
	for (ix = -WIDE/2; ix < WIDE/2; ix++) {
		x = (double)ix/(HIGH/2);
		for (iy = -HIGH/2; iy < HIGH/2; iy++) {
			y = (double)iy/(HIGH/2);
            rgb = Eval(e, x, y, tFrame);
            rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
            rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
            rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
            im.SetRGB(rgb, ix + WIDE/2, iy + HIGH/2);
        }
        if (ix % WIDE == 0)
            printf(".");
    }
    sprintf_s(szFile, sizeof(szFile), "Big%02d.bmp", nBig++);
    im.WriteBMP(szFile);
}

#define NSUPER 2

static void
AntiAliasedPlot(Expr *e, Image &im)
{
    double x, y;
    int ix, iy;
    DisplayRGB rgb;
    char szFile[64];
    static int nBig;
    unsigned nC, nX;

    tFrame = 0.0;
    Plot(0, 0, e, im);
    //H = 1.0/(IMAGE_VERT*VERT_IMAGES);
    im.FillRGB(ML_Black);
    ML_SaveStateRandom(nX, nC);
	for (ix = -NSUPER*WIDE/2; ix < NSUPER*WIDE/2; ix++) {
		for (iy = -NSUPER*HIGH/2; iy < NSUPER*HIGH/2; iy++) {
		    x = ((double)ix + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
			y = ((double)iy + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
            rgb = Eval(e, x, y, tFrame);
            rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
            rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
            rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
            im.SplatRGB(rgb/float(NSUPER*NSUPER), (ix + NSUPER*WIDE/2)/NSUPER,
                                                  (iy + NSUPER*HIGH/2)/NSUPER);
        }
        if (ix % WIDE == 0)
            printf(".");
    }
    ML_InitializeRandom(nX, nC);
    sprintf_s(szFile, sizeof(szFile), "AABig%02d.bmp", nBig++);
    im.WriteBMP(szFile);
}

static void
AntiAliasedPlot4K(Expr *e, Image &im, int nSuper = NSUPER)
{
    double x, y;
    int ix, iy;
    DisplayRGB rgb;
    char szFile[64];
    static int nBig;
    unsigned nC, nX;

    tFrame = 0.0;
    Plot(0, 0, e, im);
    //H = 1.0/(IMAGE_VERT*VERT_IMAGES);
    im.FillRGB(ML_Black);
    ML_SaveStateRandom(nX, nC);
	for (ix = -nSuper*WIDE; ix < nSuper*WIDE; ix++) {
		for (iy = -nSuper*HIGH; iy < nSuper*HIGH; iy++) {
		    x = ((double)ix + ML_RandomFloat() - 0.5)/(nSuper*HIGH);
			y = ((double)iy + ML_RandomFloat() - 0.5)/(nSuper*HIGH);
            rgb = Eval(e, x, y, tFrame);
            rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
            rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
            rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
            im.SplatRGB(rgb/float(nSuper*nSuper), (ix + nSuper*WIDE)/nSuper,
                                                  (iy + nSuper*HIGH)/nSuper);
        }
        if (ix % WIDE == 0)
            printf(".");
    }
    ML_InitializeRandom(nX, nC);
    sprintf_s(szFile, sizeof(szFile), "AA_4K_%02d.bmp", nBig++);
    im.WriteBMP(szFile);
}


static void
AntiAliasedVector(Expr *e, Image &im, VectorXYT vOrg, VectorXYT vRight)
{
    VectorXYT vUp, v;
    double x, y, t;
    int ix, iy;
    DisplayRGB rgb;
    unsigned nC, nX;

    vUp = VectorXYT(0.0, 1.0, 0.0);
    vRight = Normalize(vRight);
    tFrame = 0.0;
    Plot(0, 0, e, im);
    ML_SaveStateRandom(nX, nC);
    im.FillRGB(ML_Black);
	for (ix = -NSUPER*WIDE/2; ix < NSUPER*WIDE/2; ix++) {
		for (iy = -NSUPER*HIGH/2; iy < NSUPER*HIGH/2; iy++) {
		    x = ((double)ix + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
			y = ((double)iy + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
            v = vOrg + vRight*x + vUp*y;
            rgb = Eval(e, v.x, v.y, v.t);
            rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
            rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
            rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
            im.SplatRGB(rgb/float(NSUPER*NSUPER), (ix + NSUPER*WIDE/2)/NSUPER,
                                                  (iy + NSUPER*HIGH/2)/NSUPER);
        }
    }
    printf(".");
    ML_InitializeRandom(nX, nC);
}

static void
VideoCycle(Expr *e)
{
    VectorXYT vOrg, vRight;
    Image im;
    Video vid;
    double theta, rad;

    im.NewImage(WIDE, HIGH, 3);
    vid.NewVideo("Cycle.mpeg", im);
    for (theta = 0.0; theta < 360.0; theta += 1.0)
        printf("-");
    printf("\n");
    for (theta = 0.0; theta < 360.0; theta += 1.0) {
        rad = D_PI*theta/180.0;
        vOrg = VectorXYT(cos(rad), 0.0, sin(rad))*3.0 + VectorXYT(2.5, -0.5, 0.0);
        vRight = VectorXYT(cos(rad), 0.0, sin(rad));
        AntiAliasedVector(e, im, vOrg, vRight);
        vid.WriteFrame(im);
    }
    vid.Close();
}

void
VideoAntiAliased(Expr *e, double tLimit, double framerate = 1.0/30.0)
{
    double x, y;
    int ix, iy;
    DisplayRGB rgb;
    char szFile[64];
    static int nBig;
    unsigned nC, nX;
    Image im;
    Video vid;

    im.NewImage(WIDE, HIGH, 3);
    if (vid.NewVideo("AAvideo.mpeg", im) == 0) {
        printf("cannot open file\n");
        return;
    }
    // tFrame = 0.0;
    // Plot(0, 0, e, im);
    //H = 1.0/(IMAGE_VERT*VERT_IMAGES);
    ML_SaveStateRandom(nX, nC);
    for (tFrame = -tLimit; tFrame <= tLimit; tFrame += framerate)
        printf("-");
    printf("\n");
    for (tFrame = -tLimit; tFrame <= tLimit; tFrame += framerate) {
        Plot(0, 0, e, im);      // sets hi and lo brightness
        im.FillRGB(ML_Black);
	    for (ix = -NSUPER*WIDE/2; ix < NSUPER*WIDE/2; ix++) {
		    for (iy = -NSUPER*HIGH/2; iy < NSUPER*HIGH/2; iy++) {
		        x = ((double)ix + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
			    y = ((double)iy + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
                rgb = Eval(e, x, y, tFrame);
                rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
                rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
                rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
                im.SplatRGB(rgb/float(NSUPER*NSUPER), (ix + NSUPER*WIDE/2)/NSUPER,
                                                      (iy + NSUPER*HIGH/2)/NSUPER);
            }
        }
        vid.WriteFrame(im);
        printf(".");
    }
    vid.Close();
    ML_InitializeRandom(nX, nC);
}
/*
void
VideoSlowMotion(Expr *e, double tLimit)
{
    double x, y;
    int ix, iy;
    DisplayRGB rgb;
    char szFile[64];
    static int nBig;
    unsigned nC, nX;
    Image im;
    Video vid;

    im.NewImage(WIDE, HIGH, 3);
    vid.NewVideo("SAvideo.mpeg", im);
    //tFrame = 0.0;
    //Plot(0, 0, e, im);
    //H = 1.0/(IMAGE_VERT*VERT_IMAGES);
    ML_SaveStateRandom(nX, nC);
    for (tFrame = -tLimit; tFrame <= tLimit; tFrame += 1.0/120.0)
        printf("-");
    printf("\n");
    for (tFrame = -tLimit; tFrame <= tLimit; tFrame += 1.0/120.0) {
        Plot(0, 0, e, im);      // sets hi and lo brightness
        im.FillRGB(ML_Black);
	    for (ix = -NSUPER*WIDE/2; ix < NSUPER*WIDE/2; ix++) {
		    for (iy = -NSUPER*HIGH/2; iy < NSUPER*HIGH/2; iy++) {
		        x = ((double)ix + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
			    y = ((double)iy + ML_RandomFloat() - 0.5)/(NSUPER*HIGH/2);
                rgb = Eval(e, x, y, tFrame);
                rgb.red = (rgb.red - s_fLo)/(s_fHi - s_fLo);
                rgb.grn = (rgb.grn - s_fLo)/(s_fHi - s_fLo);
                rgb.blu = (rgb.blu - s_fLo)/(s_fHi - s_fLo);
                im.SplatRGB(rgb/float(NSUPER*NSUPER), (ix + NSUPER*WIDE/2)/NSUPER,
                                                      (iy + NSUPER*HIGH/2)/NSUPER);
            }
        }
        vid.WriteFrame(im);
        printf(".");
    }
    vid.Close();
    ML_InitializeRandom(nX, nC);
}
*/
inline double
Lum(DisplayRGB rgb)
{
    return (rgb.red + rgb.grn + rgb.blu)/3.0;
}

double
Correlation(Image &im1, Image &im2, int iPic, int jPic)
{
    double fCorr, fN, fMean1, fMean2, fVar1, fVar2, f1, f2;
    int i, j;

    fN = fMean1 = fMean2 = fVar1 = fVar2 = fCorr = 0.0;
    //
    //  Accurate two-pass cross correlation
    //
    for (j = 1; j < im1.m_nHeight-1; j++) {
        for (i = 1; i < im1.m_nWidth-1; i++) {
                fN += 1.0;
                f2 = Lum(im2.GetRGB(i + iPic, j + jPic));
                f1 = Lum(im1.GetRGB(i, j));
                fMean1 += f1;
                fMean2 += f2;
        }
    }
    if (fN < 2.0)
        return 0.0;     // No information
    fMean1 /= fN;
    fMean2 /= fN;
    for (j = 1; j < im1.m_nHeight-1; j++) {
        for (i = 1; i < im1.m_nWidth-1; i++) {
                fN += 1.0;
                f2 = Lum(im2.GetRGB(i + iPic, j + jPic)) - fMean2;
                f1 = Lum(im1.GetRGB(i, j)) - fMean1;
                fVar1 += f1*f1;
                fVar2 += f2*f2;
                fCorr += f1*f2;
        }
    }
    fVar1 /= fN - 1.0;
    fVar2 /= fN - 1.0;
    fCorr /= fN - 1.0;
    return fCorr/sqrt(fVar1*fVar2);
}

struct ExprCorr {
    Expr    *pe;
    double  fCorrelation;
};

void
MatchImage(Image &imGoal)
{
    ExprCorr rgEC[HORZ_IMAGES*VERT_IMAGES];
}
//
//  q - quit (exit)
//  e - evolve
//  c - continue 
//  m - video sample of trials
//  r u1 u2 - reset random generator
//  v i j - plot a full size frame
//  P - print expression
//  a - 1080 antialised plot
//  k - 4k plot
//  s - survival of i j
//  M - 4 second video
//  L - 10 second video (longer)
//  S - 8 seconds of 4x slower sample rate
//  X - 8 seconds of 16x slower rate
//  K - 24 seconds of 4x slower rate
//  C - Toroidal sampling of a cycle
//
void
Evolve()
{
    Image im;
    Expr *pe, *rgeTrials[HORZ_IMAGES][VERT_IMAGES], *rgeSurvive[HORZ_IMAGES*VERT_IMAGES];
    int i, j, m, n, nSurvivors, nGeneration;
    unsigned nSeed1, nSeed2;
    FILE *psLogFile;
    char szCommand[64];

    nSeed1 = ML_TrueRandomUnsigned();
    nSeed2 = ML_TrueRandomUnsigned();
    sprintf_s(szCommand, sizeof(szCommand), "%u_%u", nSeed1, nSeed2);
    CreateDirectory(szCommand, NULL);
    SetCurrentDirectory(szCommand);
    ML_InitializeRandom(nSeed1, nSeed2);
	fopen_s(&psLogFile, "logfile.txt", "a");
	sprintf_s(szCommand, sizeof(szCommand), "r %u %u", nSeed1, nSeed2);
	fprintf(psLogFile, "%s\n", szCommand);
	printf("rseed, select, view, antialiased-view, evolve, multi-video, quit, comment\n");
    printf("Motion-picture, Long-motion-picture, Slow-motion, Cylindrial-animation\n");
	printf("%s\n", szCommand);
    im.NewImage(WIDE, HIGH, 3);
    nSurvivors = 0;
    for (nGeneration = 0; nGeneration < 100; nGeneration++) {
        s_nConstCount++;
        tFrame = 0.0;
        for (i = 0; i < HORZ_IMAGES; i++) {
            for (j = 0; j < VERT_IMAGES; j++) {
                printf("-");
            }
        }
        printf("\n");
        if (nSurvivors == 0) {
            for (i = 0; i < HORZ_IMAGES; i++) {
                for (j = 0; j < VERT_IMAGES; j++) {
                    do {
                        pe = RandomExpression();
                    } while (Plot(i, j, pe, im) == 0);
                    rgeTrials[i][j] = pe;
                    printf(".");
                }
            }
        } else if (nSurvivors == 1) {
            for (i = 0; i < HORZ_IMAGES; i++) {
                for (j = 0; j < VERT_IMAGES; j++) {
                    do {
                        pe = MutateExpression(rgeSurvive[0]);
                    } while (Plot(i, j, pe, im) == 0);
                    rgeTrials[i][j] = pe;
                    printf(".");
                }
            }
        } else {
            for (i = 0; i < HORZ_IMAGES; i++) {
                for (j = 0; j < VERT_IMAGES; j++) {
                    m = ML_RandomLessThanN(nSurvivors);
                    n = ML_RandomLessThanN(nSurvivors);
                    do {
                        pe = BreedExpressions(rgeSurvive[m], rgeSurvive[n]);
                    } while (Plot(i, j, pe, im) == 0);
                    rgeTrials[i][j] = pe;
                    printf(".");
                }
            }
        }
        sprintf_s(szCommand, sizeof(szCommand), "Generation%02d.bmp", nGeneration);
        im.WriteBMP(szCommand);
        nSurvivors = 0;
	    for (;;) {
		    printf("\nrsveq: ");
		    gets_s(szCommand, sizeof(szCommand));
		    fprintf(psLogFile, "%s\n", szCommand); fflush(psLogFile);
		    if (szCommand[0] == 'q')
			    exit(0);
		    if (szCommand[0] == 'c')
			    continue;
            if (szCommand[0] == 'e') {
			    break;
            }
            if (szCommand[0] == 'm') {
                VideoSamples(rgeTrials);
                continue;
            }
		    if (szCommand[0] == 'r') {
			    sscanf_s(szCommand+1, "%u %u", &nSeed1, &nSeed2);
			    ML_InitializeRandom(nSeed1, nSeed2);
			    nSurvivors = 0;		/* mass extinction */
			    break;
		    }
		    sscanf_s(szCommand+1, "%d%d", &i, &j);
		    if (szCommand[0] == 'v')
			    BigPlot(rgeTrials[i][j], im);
            if (szCommand[0] == 'P') {
                PrintExpr(rgeTrials[i][j]);
                printf("\n");
            }
		    if (szCommand[0] == 'a')
			    AntiAliasedPlot(rgeTrials[i][j], im);
            if (szCommand[0] == 'k') {
                im.NewImage(2*WIDE, 2*HIGH, 3);
                AntiAliasedPlot4K(rgeTrials[i][j], im, 4);
                im.NewImage(WIDE, HIGH, 3);
            }
            if (szCommand[0] == 's')
			    rgeSurvive[nSurvivors++] = rgeTrials[i][j];
            if (szCommand[0] == 'M')
                VideoAntiAliased(rgeTrials[i][j], 2.0);             //  4 seconds
            if (szCommand[0] == 'L')
                VideoAntiAliased(rgeTrials[i][j], 5.0);             // 10 seconds
            if (szCommand[0] == 'S')
                VideoAntiAliased(rgeTrials[i][j], 1.0, 1.0/120.0);  //  8 seconds slow
            if (szCommand[0] == 'X')
                VideoAntiAliased(rgeTrials[i][j], 0.25, 1.0/480.0); //  8 seconds super slow
            if (szCommand[0] == 'Y')
                VideoAntiAliased(rgeTrials[i][j], 0.25*0.25, 1.0/1920.0); //  8 seconds super super slow
            if (szCommand[0] == 'K')
                VideoAntiAliased(rgeTrials[i][j], 3.0, 1.0/120.0);  // 24 seconds slow
            if (szCommand[0] == 'C')
                VideoCycle(rgeTrials[i][j]);                        // 12 seconds
        }
    }
}

extern int RandomnessTests(unsigned (*RandomSource)(), double fTries, char *szName, char *szLabel);

int _tmain(int argc, _TCHAR* argv[])
{
    // RandomnessTests(ML_RandomUnsigned, 1000000.0, "Performance Counter", " ");    return 0;
    Evolve();
	return 0;
}