#pragma once
#include "Area.h"
#include <iomanip>

class Grid : public Area
{
private:
    std::vector <double> X, Y;  // Êîîðäèíàòû óçëîâ ñåòêè
    std::vector <int> IXw, IYw;  // Õðàíÿò ïîçèöèè êîîðäèíàò ãðàíèö â âåêòîðàõ X è Y
    std::vector <int> nx, ny; // Ðàçáèåíèÿ ïî X è Y
    std::vector <double> qx, qy;  // Êîýôôèöèåíòû ñãóùåíèÿ
    int Nx, Ny;  // Êîëè÷åñòâî óçëîâ ïî X è ïî Y

    double t0, tn;
    int nt;
    double qt;
    std::vector <double> T;

    void read_grid(std::string test_folder) {
        std::string path = "Input_Data/" + test_folder + "/";
        std::ifstream GridX, GridY;
        GridX.open(path + "GridX.txt");
        GridY.open(path + "GridY.txt");
        nx.resize(nXw - 1);
        ny.resize(nYw - 1);
        qx.resize(nXw - 1);
        qy.resize(nYw - 1);
        for (int i = 0; i < nXw - 1; i++) {
            GridX >> nx[i] >> qx[i];
        }
        for (int i = 0; i < nYw - 1; i++) {
            GridY >> ny[i] >> qy[i];
        }
        GridX.close();
        GridY.close();

        std::ifstream Time;
        Time.open(path + "Time.txt");
        Time >> t0 >> tn >> nt >> qt;
        Time.close();
    }

    void makeGrid(std::vector <double>& XY, double left, double right, int n, double qxy, int& i0, int& j, std::vector <int>& IXYw) {
        double h0;
        if (qxy - 1 < 1E-16)
            h0 = (right - left) / n;
        else
            h0 = (right - left) * (1 - qxy) / (1 - pow(qxy, n));

        XY[i0] = left;
        IXYw[j] = i0; j++;
        for (int i = i0 + 1; i < n + i0; i++) {
            XY[i] = XY[i - 1] + h0;
            h0 *= qxy;
        }
        i0 = n + i0;
    };

    void makeGrid(std::vector <double>& T, double left, double right, int n, double q) {
        double h0;
        if (q - 1 < 1E-16)
            h0 = (right - left) / n;
        else
            h0 = (right - left) * (1 - q) / (1 - pow(q, n));

        T[0] = left;
        for (int i = 1; i < n; i++) {
            T[i] = T[i - 1] + h0;
            h0 *= q;
        }
        T[n] = right;
    };

public:
    Grid(std::string &test_folder) {
        read_area(test_folder);
        read_grid(test_folder);
        Nx = 0;
        Ny = 0;
        for (int i = 0; i < nXw - 1; i++)
            Nx += nx[i];
        for (int i = 0; i < nYw - 1; i++)
            Ny += ny[i];
        Nx++; Ny++;
        X.resize(Nx);
        Y.resize(Ny);
        IXw.resize(nXw);
        IYw.resize(nYw);
        int ix0 = 0, iy0 = 0;
        int jx = 0, jy = 0;
        for (int i = 0; i < nXw - 1; i++)
            makeGrid(X, Xw[i], Xw[i + 1], nx[i], qx[i], ix0, jx, IXw);
        for (int i = 0; i < nYw - 1; i++)
            makeGrid(Y, Yw[i], Yw[i + 1], ny[i], qy[i], iy0, jy, IYw);
        X[Nx - 1] = Xw[nXw - 1];
        Y[Ny - 1] = Yw[nYw - 1];
        IXw[nXw - 1] = ix0;
        IYw[nYw - 1] = iy0;

        T.resize(nt + 1);
        makeGrid(T, t0, tn, nt, qt);
    }


    // Îïðåäåëÿåò íîìåð ïîäîáëàñòè, â êîòîðîé ðàñïîëîæåí êîíå÷íûé ýëåìåíò
    int inSubArea(int px, int py) {
        int ixw1, ixw2, iyw1, iyw2;
        for (int i = 0; i < nw; i++) {
            ixw1 = IXw[Mw[i].nx1];
            ixw2 = IXw[Mw[i].nx2];
            iyw1 = IYw[Mw[i].ny1];
            iyw2 = IYw[Mw[i].ny2];
            bool flag1, flag2;
            flag1 = px >= ixw1 && px <= ixw2 && px + 1 >= ixw1 && px + 1 <= ixw2;
            flag2 = py >= iyw1 && py <= iyw2 && py + 1 >= iyw1 && py + 1 <= iyw2;
            if (flag1 && flag2)
                return Mw[i].ni;
        }
        return -1;  // Ôèêòèâíûé óçåë
    }


    double getX(int i) { return X[i]; }
    double getY(int i) { return Y[i]; }
    double getT(int i) { return T[i]; }
    int getIXw(int i) { return IXw[i]; }
    int getIYw(int i) { return IYw[i]; }
    int getNx() { return Nx; }
    int getNy() { return Ny; }
    int getNt() { return nt + 1; }
    void print_grid(std::ofstream &ofile) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                ofile << std::setprecision(13) << std::scientific << X[i] << " " << Y[j] << '\n';
            }
        }
    }
};
