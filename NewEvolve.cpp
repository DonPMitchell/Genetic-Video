#include "stdafx.h"
#include "texture.h"
#include "Image.h"
#include "Video.h"
#include "MonteCarlo.h"
#pragma intrinsic(sqrt, sin, cos, exp, log, fabs, atan2)

#define ARG_XXX     1
#define ARG_YYY     2
#define ARG_TTT     3
#define ARG_XYT     4
#define ARG_ARG1    5 
#define ARG_ARG2    6

struct Expression {
    int         nOperation;
    DisplayRGB  rgb;
    Expression  *pArg1;
    Expression  *pArg2;

    DisplayRGB  operator() (double x, double y, double t);
};

