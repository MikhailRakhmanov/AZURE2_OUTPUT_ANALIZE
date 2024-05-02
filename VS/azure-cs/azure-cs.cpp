// azure-cs.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex>

#include "..\Constants.h"
#include "..\my_services.h"

double Clebsh_Gordon(int JX1, int JX2, int JX3, int MX1, int MX2);
double YXFCT(int M, int N);
void doAngelarDistribution();
void doPhaseCalculation();

using namespace std;

int main()
{
    //doAngelarDistribution();
    doPhaseCalculation();

    return 0;
}
//------------------------------------------------------------------------------
void doPhaseCalculation()
{
    FILE
        *fin1, *fin2, *fin3, *fout;

    fopen_s(&fin1, "L=0_J=0_5.extrap", "r");
    fopen_s(&fin2, "L=1_J=0_5.extrap", "r");
    fopen_s(&fin3, "L=1_J=1_5.extrap", "r");
    fopen_s(&fout, "phase.csv", "w");

    char
        str[500],
        stmp[50];
    int
        nth = 91,
        ne = 499,
        _2j = 3,
        ie, il, ith, itmp, iM, _iM;
    double
        ** cs = new double* [ne],
        ** phase = new double* [ne],
        * e = new double[ne],
        csect,
        ecm,
        theta, cth,
        pene,
        re, im, x;
    complex <double>
        I = complex <double>(0.0, 1.0),
        ** fMM = new complex <double> * [_2j + 1];

    for (iM = 0; iM < _2j + 1; iM++)
    {
        fMM[iM] = new complex <double> [_2j + 1];
    }

    for (ie = 0; ie < ne; ie++)
    {
        phase[ie] = new double[3];

        fgets(str, 500, fin1);
        sscanf_s(str, "%lf %lf %lf %lf %lf", &ecm, &x, &x, &re, &x);
        e[ie] = ecm;
        phase[ie][0] = re * Constants::degTorad;

        fgets(str, 500, fin2);
        sscanf_s(str, "%lf %lf %lf %lf %lf", &ecm, &x, &x, &re, &x);
        phase[ie][1] = re * Constants::degTorad;

        fgets(str, 500, fin3);
        sscanf_s(str, "%lf %lf %lf %lf %lf", &ecm, &x, &x, &re, &x);
        phase[ie][2] = re * Constants::degTorad;
    }
    fclose(fin1); fclose(fin2); fclose(fin3);

    for (ie = 0; ie < ne; ie++)
    {
        fprintf(fout, "%f;%e;%e;%e\n", e[ie], phase[ie][0], phase[ie][1], phase[ie][2]);
    }
    fclose(fout);

    fopen_s(&fout,"cs.csv", "w");
    // amplitude
    int
        Lmax = 1,
        i2MJ, _i2MJ,
        iJL, _2J, _2S = 1,
        _2L, L, _2mL, mL;
    double
        h2_mn = 41.8015876123,
        A1 = 1.0, A2 = 6.0,
        mu = A1 * A2 / (A1 + A2),
        k,
        cg1, cg2, Y, Plm, SF, CS;
    complex <double>
        fJL;

    ie = 44;//for(ie=0;ie<ne;ie++)
    {

        k = sqrt(2.0 * mu * e[ie] / h2_mn);
        for (ith = 0; ith < nth; ith++)
        {
            theta = (2.0 * ith + 0.0001);
            cth = cos(theta * Constants::degTorad);
            for (iM = 0; iM < _2j + 1; iM++) for (_iM = 0; _iM < _2j + 1; _iM++) fMM[iM][_iM] = 0.0;

            for (iM = 0, i2MJ = -3; iM < _2j + 1; iM++, i2MJ += 2)
            {
                for (_iM = 0, _i2MJ = -3; _iM < _2j + 1; _iM++, _i2MJ += 2)
                {
                    _2mL = i2MJ - _i2MJ; if (abs(_2mL) > 2 * Lmax) continue;
                    mL = _2mL / 2;
                    for (iJL = 0; iJL < _2j; iJL++)
                    {
                        if (iJL == 0) // l=0, j=1/2
                        {
                            SF = 1.0; //continue;
                            _2L = 0; L = 0;
                            _2J = 1;
                        }
                        else if (iJL == 1) // l=1, j=1/2
                        {
                            SF = 1.0;
                            _2L = 2; L = 1;
                            _2J = 1;
                        }
                        else // l=1, j=3/2
                        {
                            SF = 1.0; //continue;
                            _2L = 2; L = 1;
                            _2J = 3;
                        }
                        if (abs(_2mL) > _2L) continue;
                        if (abs(i2MJ) > _2S) continue;
                        if (abs(i2MJ) > _2J || abs(_i2MJ) > _2J) continue;
                        if (abs(_i2MJ) > _2S) continue;


                        cg1 = Clebsh_Gordon(_2L, _2S, _2J, 0, i2MJ);
                        cg2 = Clebsh_Gordon(_2L, _2S, _2J, _2mL, _i2MJ);
                        Y = (_2L + 1.0) * sqrt(YXFCT(L + abs(mL), L - abs(mL)));

                        if (L == 0) Plm = 1.0;
                        else if (L == 1 && abs(mL) == 0) Plm = cth;
                        else if (L == 1 && abs(mL) == 1) Plm = -sqrt(1.0 - cth * cth);
                        else Plm = 0.0;

                        fJL = (exp(2.0 * I * phase[ie][iJL]) - 1.0) / (2.0 * I * k);

                        fMM[iM][_iM] += SF * cg1 * cg2 * Y * Plm * fJL;
                    }
                }
            }

            // cross-section
            CS = 0.0;
            for (iM = 0, i2MJ = -3; iM < _2j + 1; iM++, i2MJ += 2)
            {
                for (_iM = 0, _i2MJ = -3; _iM < _2j + 1; _iM++, _i2MJ += 2)
                {
                    CS += norm(fMM[iM][_iM]);
                }
            }
            CS = 10.0 * Constants::M_PI * CS / (_2S + 1.0);

            fprintf(fout, "%f;%f;%e\n", e[ie], theta, CS);
        }
    }
    fclose(fout);
}
//------------------------------------------------------------------------------
void doAngelarDistribution()
{
    FILE
        *fin,
        *fout;
    fopen_s(&fin, "AZUREOut_aa=1_R=1.extrap", "r");
    fopen_s(&fout, "out.dat", "w");

    char
        str[200];
    int
        nth = 91,
        ne = 499,
        ie, ith;
    double
        ** cs = new double* [ne],
        * e = new double[ne],
        * th = new double[nth],
        csect,
        ecm,
        theta,
        x;

    for (ith = 0; ith < nth; ith++)
    {
        cs[ith] = new double[ne];
        for (ie = 0; ie < ne; ie++)
        {
            fgets(str, 200, fin);
            sscanf_s(str, "%lf %lf %lf %lf %lf", &ecm, &x, &theta, &csect, &x);

            if (ith == 0) e[ie] = ecm;
            if (ie == 0) th[ith] = theta;

            cs[ith][ie] = csect;
        }
    }
    //---------------------------------------------------------------------------
    fprintf(fout, "# Ntheta = %d, NE = %d\n", nth, ne);
    fprintf(fout, "# Ecm(MeV)  theta_cm(deg) ds/dW(b/sr)\n");
    for (ie = 0; ie < ne; ie++)
    {
        for (ith = 0; ith < nth; ith++)
        {
            fprintf(fout, "%13.6e %13.6e %13.6e\n", e[ie], th[ith], cs[ith][ie]);
        }
    }

    fclose(fin);
    fclose(fout);
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
