#pragma once
class Diff_Parameters_test0
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        switch (wi)
        {
        case 1:
            return 2;
        case 2:
            return 1;
        case 3:
            return 0;
        }
    }

    double lambda(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return 1;
        case 2:
            return 10;
        case 3:
            return 10;
        }
    }

    double f(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return 2 * x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return 0;
        }
    }

    double u_g(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return 2;
        case 2:
            return 1.8 + 0.1 * x;
        }
    }

    double theta(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};


class Diff_Parameters_test00
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 3;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y) {
        return exp(x + y);
    }

    double u_g(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return exp(x);
        case 2:
            return exp(y + 4);
        case 3:
            return exp(x + 4);
        case 4:
            return exp(y);
        }
    }

    double theta(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};


class Diff_Parameters_test1
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 0;
    }

    double sigma(int wi) {
        return 0;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 0;
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return y + 4;
        case 3:
            return x + 4;
        case 4:
            return y;
        }
    }

    double u0(int wi, double x, double y) {
        return x + y;
    }

    double u1(int wi, double x, double y, double t1) {
        return x + y;
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};


class Diff_Parameters_test2
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 0;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 2;
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x + 2 * t;
        case 2:
            return y + 4 + 2 * t;
        case 3:
            return x + 4 + 2 * t;
        case 4:
            return y + 2 * t;
        }
    }

    double u0(int wi, double x, double y) {
        return x + y;
    }

    double u1(int wi, double x, double y, double t1) {
        return x + y + 1;
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_test3
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 0;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 8;
    }

    double u0(int wi, double x, double y) {
        return x + y;
    }

    double u1(int wi, double x, double y, double t1) {
        return x + y + 2 * 0.5 + 4 * 0.5 * 0.5;
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x + 2*t + 4 * t * t;
        case 2:
            return y + 2 * t + 4 + 4 * t * t;
        case 3:
            return x + 4 + 2 * t + 4 * t * t;
        case 4:
            return y + 2 * t + 4 * t * t;
        }
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_test4
{
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 8*t + 10;
    }

    double u0(int wi, double x, double y) {
        return x + y;
    }

    double u1(int wi, double x, double y, double t1) {
        return x + y + 2 * 0.5 + 4 * 0.5 * 0.5;
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x + 2 * t + 4 * t * t;
        case 2:
            return y + 2 * t + 4 + 4 * t * t;
        case 3:
            return x + 4 + 2 * t + 4 * t * t;
        case 4:
            return y + 2 * t + 4 * t * t;
        }
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_test5
{
protected:
    double chi(int wi) {// греческая буква Хи
        switch (wi)
        {
        case 1:
            return 3;
        case 2:
            return 2;
        case 3:
            return 1;
        }
    }

    double sigma(int wi) {
        return 0;
    }

    double lambda(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return 1;
        case 2:
            return 10;
        case 3:
            return 10;
        }
    }

    double f(int wi, double x, double y, double t) {
        switch (wi)
        {
        case 1:
            return 24;
        case 2:
            return 16;
        case 3:
            return 8;
        }
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 2 + 4 * t * t;
        case 2:
            return 1.8 + 0.1 * x + 4 * t * t;
        }
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x + 4 * t * t;
        case 2:
            return 1.8 + 0.1 * x + 4 * t * t;
        case 3:
            return -1 + 4 * t * t;
        }
    }

    double u0(int wi, double x, double y) {
        switch (wi)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return 1.8 + 0.1 * x;
        }
    }

    double u1(int wi, double x, double y, double t1) {
        switch (wi)
        {
        case 1:
            return x + 1;
        case 2:
            return 1.8 + 0.1 * x + 1;
        case 3:
            return 1.8 + 0.1 * x + 1;
        }
    }
};

class Diff_Parameters_research_uniform_x {
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 4 * t * t * exp(t * t) + 2 * t * exp(t * t) + 2 * exp(t * t);
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x + exp(t * t);
        case 2:
            return y + 4 + exp(t * t);
        case 3:
            return x + 4 + exp(t * t);
        case 4:
            return y + exp(t * t);
        }
    }

    double u0(int wi, double x, double y) {
        return x + y + 1;
    }

    double u1(int wi, double x, double y, double t1) {
        return x + y + exp(t1 * t1);
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_research_uniform_t {
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return -2 * exp(x + y) + 2 * t + 2;
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return t * t + exp(x);
        case 2:
            return t * t + exp(y+4);
        case 3:
            return t * t + exp(x+4);
        case 4:
            return t * t + exp(y);
        }
    }

    double u0(int wi, double x, double y) {
        return exp(x + y);
    }

    double u1(int wi, double x, double y, double t1) {
        return t1*t1 + exp(x+y);
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_research_uniform_xt {
protected:
    double chi(int wi) {// греческая буква Хи
        return 1;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return 4 * t * t * exp(x + y + t * t) + 2 * t * exp(x + y + t * t);
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return exp(x + t * t);
        case 2:
            return exp(y + 4 + t * t);
        case 3:
            return exp(x + 4 + t * t);
        case 4:
            return exp(y + t * t);
        }
    }

    double u0(int wi, double x, double y) {
        return exp(x + y);
    }

    double u1(int wi, double x, double y, double t1) {
        return exp(x + y + t1 * t1);
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};

class Diff_Parameters_research_uniform_tt{
protected:
    double chi(int wi) {// греческая буква Хи
        return 0;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double x, double y) {
        return 1;
    }

    double f(int wi, double x, double y, double t) {
        return exp(t);
    }

    double u_g(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return exp(t);
        case 2:
            return exp(t);
        case 3:
            return exp(t);
        case 4:
            return exp(t);
        }
    }

    double u0(int wi, double x, double y) {
        return 1;
    }

    double u1(int wi, double x, double y, double t1) {
        return exp(t1);
    }

    double theta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 0;
        }
    }

    double beta(int si) {
        switch (si)
        {
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 0.5;
        }
    }

    double u_beta(int si, double x, double y, double t) {
        switch (si)
        {
        case 1:
            return x;
        case 2:
            return 1.8 + 0.1 * x;
        case 3:
            return -1;
        }
    }
};