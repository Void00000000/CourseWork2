#pragma once
#include "Grid.h"
#include <algorithm>
#include "Func.h"
#include "SolveSoLE_LOS.h"
#include <iomanip>

class FEM : Diff_Parameters_research_uniform_tt
{
public:
    FEM(Grid& grid, std::string &test_folder) {
        this->test_folder = test_folder;
        Nx = grid.getNx();
        Ny = grid.getNy();
        N = Nx * Ny;
        di.resize(N);
        q.resize(N);
        q_1.resize(N);
        q_2.resize(N);
        b.resize(N);
        generate_portrait();
        read_boundary_conditions();
        make_prev_q(grid, q_2, true);
        make_prev_q(grid, q_1, false);
    }

    void solve(Grid& grid, double t_2, double t_1, double t) {
        double delta_t = t - t_2;
        double delta_t1 = t_1 - t_2;
        double delta_t0 = t - t_1;
        G_matrix_assembly(grid);
        M_matrix_assembly(grid, 2 / (delta_t * delta_t0), false);
        M_matrix_assembly(grid, (delta_t + delta_t0) / (delta_t * delta_t0), true);
        b_vector_assembly(grid, t, 1, true);
        b_vector_assembly(grid, t, -2 / (delta_t1 * delta_t), false, false, false);
        b_vector_assembly(grid, t, 2 / (delta_t1 * delta_t0), false, true, false);
        b_vector_assembly(grid, t, -delta_t0 / (delta_t1 * delta_t), false, false, true);
        b_vector_assembly(grid, t, delta_t / (delta_t1 * delta_t0), false, true, true);
        use_bc2(grid, t);
        use_bc3(grid, t);
        use_bc1(grid, t);
        process_fict_nodes();
        LOS_solve(q, di, ggl, ggu, b, ig, jg);
        q_2 = q_1;
        q_1 = q;
    }

    void print_vector(std::ofstream &ofile, double t) {
        ofile << std::setprecision(13) << std::scientific << "t = " << t << '\n';
        ofile << "--------------------------------------\n";
        for (int i = 0; i < N; i++) {
            //ofile.setf(std::ios::scientific);
            ofile << std::setprecision(13) <<std::scientific << q[i] << '\n';
        }
        ofile << "--------------------------------------\n\n";
    }


    void update() {
        for (int i = 0; i < N; i++) {
            di[i] = 0;
            b[i] = 0;
            q[i] = 0;
        }
        for (int i = 0; i < n_jgg; i++) {
            ggl[i] = 0;
            ggu[i] = 0;
        }
    }

private:
    std::string test_folder;
    int N;  // Количество узлов
    int Nx;
    int Ny;
    std::vector <double> q;
    std::vector <double> q_1;
    std::vector <double> q_2;
    std::vector <double> b;
    std::vector <int> fict_nodes;
    int n_jgg;            // Размерность векторов ggl, ggu и jg
    // Матрица A--------------
    std::vector <double> di;
    std::vector <double> ggl;
    std::vector <double> ggu;
    std::vector <int> ig;
    std::vector <int> jg;
    // -----------------------

    // Краевые условия-----------------
    struct bc {
        int ni;  // Номер границы
        int nx1;  // Номер элемента в векторе Xw, с которого начинается начало фрагмента границы
        int nx2;
        int ny1;
        int ny2;
    };
    std::vector <bc> bc1;  // Данные о границе s1
    std::vector <bc> bc2;  // Данные о границе s2
    std::vector <bc> bc3;  // Данные о границе s3
    int ns1, ns2, ns3; // Количество частей границы Si
    std::vector <std::vector<double>> As3 = { {0, 0},
                                              {0, 0} };  // Локальная матрица для учёта третьих краевых условий
    //----------------------------------

    int L[4];    // Вектор для хранения номеров глобальных базисных функций
    std::vector<std::vector<double>> G = { {0 , 0, 0, 0},
                                           {0, 0, 0, 0},
                                           {0, 0, 0, 0},
                                           {0, 0, 0, 0} };  // Локальная матрица жёсткости
    std::vector <std::vector <double>> G1 = { {1.0, -1.0, 1.0 / 3, -1.0 / 3},
                                    {-1.0, 1.0, -1.0 / 3, 1.0 / 3},
                                    {1.0 / 3, -1.0 / 3, 1.0 / 3, -1.0 / 3},
                                    {-1.0 / 3, 1.0 / 3, -1.0 / 3, 1.0 / 3} };

    std::vector <std::vector <double>> G2 = { {1.0, 1.0 / 3, -1.0, -1.0 / 3},
                                    {1.0 / 3, 1.0 / 3, -1.0 / 3, -1.0 / 3},
                                    {-1.0, -1.0 / 3, 1.0, 1.0 / 3},
                                    {-1.0 / 3, -1.0 / 3, 1.0 / 3, 1.0 / 3} };

    std::vector <std::vector <double>> G3 = { {1.0 / 3, 1.0 / 3, -1.0 / 3, -1.0 / 3},
                                    {1.0 / 3, 1.0, -1.0 / 3, -1.0},
                                    {-1.0 / 3, -1.0 / 3, 1.0 / 3, 1.0 / 3},
                                    {-1.0 / 3, -1.0, 1.0 / 3, 1.0} };

    std::vector <std::vector <double>> G4 = { {1.0 / 3, -1.0 / 3, 1.0 / 3, -1.0 / 3},
                                    {-1.0 / 3, 1.0 / 3, -1.0 / 3, 1.0 / 3},
                                    {1.0 / 3, -1.0 / 3, 1.0, -1.0},
                                    {-1.0 / 3, 1.0 / 3, -1.0, 1.0} };

    std::vector <std::vector <double>> M = { {0 , 0, 0, 0},
                                             {0, 0, 0, 0},
                                             {0, 0, 0, 0},
                                             {0, 0, 0, 0} };  // локальная матрица массы
    std::vector <std::vector <int>> C = { {4, 2, 2, 1},
                                {2, 4, 1, 2},
                                {2, 1, 4, 2},
                                {1, 2, 2, 4} };


    // Возвращает глоабальный номер базисной функции
    int global_num(int i, int j) {
        return j * Nx + i;
    }

    void generate_portrait() {
        std::vector <std::vector <int>> list(N);
        list[0].push_back(0);
        int g1, g2;  // Глобальные номера базисных функций
        bool not_in;
        // Цикл по конечным элементам
        for (int j = 0; j < Ny - 1; j++)
            for (int i = 0; i < Nx - 1; i++) {
                L[0] = global_num(i, j);
                L[1] = global_num(i + 1, j);
                L[2] = global_num(i, j + 1);
                L[3] = global_num(i + 1, j + 1);
                // Цикл по ненулевым базисным функциям
                for (int in = 0; in < 4; in++) {
                    g1 = L[in];
                    for (int jn = in + 1; jn < 4; jn++) {
                        // g2 > g1
                        g2 = L[jn];
                        // Перед добавлением проверяем наличие элемента в списке
                        not_in = true;
                        for (int l = 0; l < list[g2].size() && not_in; l++)
                            if (g1 == list[g2][l])
                                not_in = false;

                        // Добавляем
                        if (not_in)
                            list[g2].push_back(g1);
                    }
                }
            }

            
        // Сортировка списков по возрастанию
        for (int i = 0; i < N; i++)
            sort(list[i].begin(), list[i].end());

        // Формирование вектора ig
        ig.resize(N + 1);
        ig[0] = 0;
        for (int i = 0; i < list.size(); i++)
            ig[i + 1] = ig[i] + list[i].size();

        for (int i = 1; i < N + 1; i++)
            ig[i] -= 1;

        n_jgg = ig[N];
        jg.resize(n_jgg);
        ggl.resize(n_jgg);
        ggu.resize(n_jgg);
        // Формирование вектора jg
        for (int i = 1, j = 0; i < N; i++)
            for (int k = 0; k < list[i].size(); k++, j++)
                jg[j] = list[i][k];
    }




    // Занесение локальной матрицы размером k на k в глобальную матрицу A
    void add_local_matrix(std::vector<std::vector <double>> &local_matrix, int k) {
        int ibeg, iend, med;
        for (int i = 0; i < k; i++)
            di[L[i]] = di[L[i]] + local_matrix[i][i];

        for (int i = 0; i < k; i++) {
            ibeg = ig[L[i]];
            for (int j = 0; j <= i - 1; j++) {
                iend = ig[L[i] + 1] - 1;
                while (jg[ibeg] != L[j]) {
                    med = (ibeg + iend) / 2;
                    if (jg[med] < L[j])
                        ibeg = med + 1;
                    else
                        iend = med;
                }
                ggl[ibeg] += local_matrix[i][j];
                ggu[ibeg] += local_matrix[j][i];
                ibeg++;
            }
        }

    }


    void G_matrix_assembly(Grid& grid) {
        double lambda1, lambda2, lambda3, lambda4;
        double h_x, h_y;
        double xmin, xmax, ymin, ymax;
        int wi;  // Номер подобласти
        for (int j = 0; j < Ny - 1; j++)
        {
            ymax = grid.getY(j + 1);
            ymin = grid.getY(j);
            h_y = ymax - ymin;
            for (int i = 0; i < Nx - 1; i++)
            {
                xmax = grid.getX(i + 1);
                xmin = grid.getX(i);
                h_x = xmax - xmin;
                L[0] = global_num(i, j);
                L[1] = global_num(i + 1, j);
                L[2] = global_num(i, j + 1);
                L[3] = global_num(i + 1, j + 1);
                wi = grid.inSubArea(i, j);
                if (wi < 0) {
                    fict_nodes.push_back(L[1]);
                    continue;
                };
                lambda1 = lambda(wi, xmin, ymin);
                lambda2 = lambda(wi, xmax, ymin);
                lambda3 = lambda(wi, xmin, ymax);
                lambda4 = lambda(wi, xmax, ymax);

                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < i; j++) {
                        G[i][j] = lambda1 * h_y / (8 * h_x) * G1[i][j] + lambda1 * h_x / (8 * h_y) * G2[i][j] +
                            lambda2 * h_y / (8 * h_x) * G1[i][j] + lambda2 * h_x / (8 * h_y) * G3[i][j] +
                            lambda3 * h_y / (8 * h_x) * G4[i][j] + lambda3 * h_x / (8 * h_y) * G2[i][j] +
                            lambda4 * h_y / (8 * h_x) * G4[i][j] + lambda4 * h_x / (8 * h_y) * G3[i][j];
                        G[j][i] = G[i][j];
                    }
                    G[i][i] = lambda1 * h_y / (8 * h_x) * G1[i][i] + lambda1 * h_x / (8 * h_y) * G2[i][i] +
                        lambda2 * h_y / (8 * h_x) * G1[i][i] + lambda2 * h_x / (8 * h_y) * G3[i][i] +
                        lambda3 * h_y / (8 * h_x) * G4[i][i] + lambda3 * h_x / (8 * h_y) * G2[i][i] +
                        lambda4 * h_y / (8 * h_x) * G4[i][i] + lambda4 * h_x / (8 * h_y) * G3[i][i];

                }
                add_local_matrix(G, 4);
            }
        }
    }


    void M_matrix_assembly(Grid& grid, double m, bool isSigma) {
        double h_x, h_y;
        double gamma;
        double xmin, xmax, ymin, ymax;
        int wi;  // Номер подобласти
        for (int j = 0; j < Ny - 1; j++)
        {
            ymax = grid.getY(j + 1);
            ymin = grid.getY(j);
            h_y = ymax - ymin;
            for (int i = 0; i < Nx - 1; i++)
            {
                xmax = grid.getX(i + 1);
                xmin = grid.getX(i);
                h_x = xmax - xmin;
                L[0] = global_num(i, j);
                L[1] = global_num(i + 1, j);
                L[2] = global_num(i, j + 1);
                L[3] = global_num(i + 1, j + 1);
                wi = grid.inSubArea(i, j);
                if (wi < 0) {
                    continue;
                };
                isSigma ? gamma = sigma(wi) : gamma = chi(wi);
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < i; j++) {
                        M[i][j] = m * gamma * h_x * h_y * C[i][j] / 36.0;
                        M[j][i] = M[i][j];

                    }
                    M[i][i] = m * gamma * h_x * h_y * C[i][i] / 36.0;

                }
                add_local_matrix(M, 4);
            }
        }
    }

    // whichQ = true -> q_1 else q_2
    void b_vector_assembly(Grid& grid, double t, double m, bool isF, bool whichQ = false, bool isSigma = false) {
        double f1, f2, f3, f4;
        double gamma;
        double h_x, h_y;
        double xmin, xmax, ymin, ymax;
        int wi;  // Номер подобласти
        for (int j = 0; j < Ny - 1; j++)
        {
            ymax = grid.getY(j + 1);
            ymin = grid.getY(j);
            h_y = ymax - ymin;
            for (int i = 0; i < Nx - 1; i++)
            {
                xmax = grid.getX(i + 1);
                xmin = grid.getX(i);
                h_x = xmax - xmin;
                L[0] = global_num(i, j);
                L[1] = global_num(i + 1, j);
                L[2] = global_num(i, j + 1);
                L[3] = global_num(i + 1, j + 1);
                wi = grid.inSubArea(i, j);
                if (wi < 0) {
                    continue;
                };
                if (isF) {
                    f1 = f(wi, xmin, ymin, t);
                    f2 = f(wi, xmax, ymin, t);
                    f3 = f(wi, xmin, ymax, t);
                    f4 = f(wi, xmax, ymax, t);
                    gamma = 1;
                }
                else
                {
                    isSigma ? gamma = sigma(wi) : gamma = chi(wi);
                    if (whichQ) {
                        f1 = q_1[L[0]];
                        f2 = q_1[L[1]];
                        f3 = q_1[L[2]];
                        f4 = q_1[L[3]];
                    }
                    else
                    {
                        f1 = q_2[L[0]];
                        f2 = q_2[L[1]];
                        f3 = q_2[L[2]];
                        f4 = q_2[L[3]];
                    }
                }
                for (int i = 0; i < 4; i++)
                    b[L[i]] += gamma * m * h_x * h_y * (C[i][0] * f1 + C[i][1] * f2 + C[i][2] * f3 + C[i][3] * f4) / 36.0;
            }
        }
    }


    void read_boundary_conditions() {
        std::string path = "Input_Data/" + test_folder + "/";
        std::ifstream bc_file;
        bc_file.open(path + "bc.txt");
        bc_file >> ns1;
        bc1.resize(ns1);
        for (int i = 0; i < ns1; i++) {
            bc_file >> bc1[i].ni >> bc1[i].nx1 >> bc1[i].nx2 >> bc1[i].ny1 >> bc1[i].ny2;
            --bc1[i].nx1;
            --bc1[i].nx2;
            --bc1[i].ny1;
            --bc1[i].ny2;
        }
        bc_file >> ns2;
        bc2.resize(ns2);
        for (int i = 0; i < ns2; i++) {
            bc_file >> bc2[i].ni >> bc2[i].nx1 >> bc2[i].nx2 >> bc2[i].ny1 >> bc2[i].ny2;
            --bc2[i].nx1;
            --bc2[i].nx2;
            --bc2[i].ny1;
            --bc2[i].ny2;
        }
        bc_file >> ns3;
        bc3.resize(ns3);
        for (int i = 0; i < ns3; i++) {
            bc_file >> bc3[i].ni >> bc3[i].nx1 >> bc3[i].nx2 >> bc3[i].ny1 >> bc3[i].ny2;
            --bc3[i].nx1;
            --bc3[i].nx2;
            --bc3[i].ny1;
            --bc3[i].ny2;
        }
        bc_file.close();
    }


    void use_bc1(Grid& grid, double t) {
        int p;
        int i_beg, i_end, j_beg, j_end;
        double x1, y1;
        int l;
        for (int s = 0; s < ns1; s++) {
            i_beg = grid.getIXw(bc1[s].nx1);
            i_end = grid.getIXw(bc1[s].nx2);
            j_beg = grid.getIYw(bc1[s].ny1);
            j_end = grid.getIYw(bc1[s].ny2);
            p = bc1[s].ni;
            if (i_beg == i_end) {
                int i = i_beg;
                x1 = grid.getX(i);
                for (int j = j_beg; j <= j_end; j++) {
                    y1 = grid.getY(j);
                    l = global_num(i, j);
                    di[l] = 1;
                    for (int i = ig[l]; i < ig[l + 1]; i++)
                      ggl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                      if (jg[i] == l)
                          ggu[i] = 0;
                    b[l] = u_g(p, x1, y1, t);
                }
            }
            else {
                int j = j_beg;
                y1 = grid.getY(j);
                for (int i = i_beg; i <= i_end; i++) {
                    x1 = grid.getX(i);
                    l = global_num(i, j);
                    di[l] = 1;
                    for (int i = ig[l]; i < ig[l + 1]; i++)
                        ggl[i] = 0;
                    for (int i = 0; i < n_jgg; i++)
                        if (jg[i] == l)
                            ggu[i] = 0;
                    b[l] = u_g(p, x1, y1, t);
                }
            }
        }
    }


    void use_bc2(Grid& grid, double t) {
        int p;
        int i_beg, i_end, j_beg, j_end;
        double x1, x2, y1, y2, h;
        double theta1, theta2;
        int l1, l2;
        for (int s = 0; s < ns2; s++) {
            i_beg = grid.getIXw(bc2[s].nx1);
            i_end = grid.getIXw(bc2[s].nx2);
            j_beg = grid.getIYw(bc2[s].ny1);
            j_end = grid.getIYw(bc2[s].ny2);
            p = bc2[s].ni;
            if (i_beg == i_end) {
                int i = i_beg;
                x1 = grid.getX(i);
                for (int j = j_beg; j < j_end; j++) {
                    y1 = grid.getY(j);
                    y2 = grid.getY(j + 1);
                    l1 = global_num(i, j);
                    l2 = global_num(i, j + 1);
                    theta1 = theta(p, x1, y1, t);
                    theta2 = theta(p, x1, y2, t);
                    h = y2 - y1;
                    b[l1] += h * (2 * theta1 + theta2) / 6.0;
                    b[l2] += h * (theta1 + 2 * theta2) / 6.0;

                }
            }
            else {
                int j = j_beg;
                y1 = grid.getY(j);
                for (int i = i_beg; i < i_end; i++) {
                    x1 = grid.getX(i);
                    x2 = grid.getX(i + 1);
                    l1 = global_num(i, j);
                    l2 = global_num(i + 1, j);
                    theta1 = theta(p, x1, y1, t);
                    theta2 = theta(p, x2, y1, t);
                    h = x2 - x1;
                    b[l1] += h * (2 * theta1 + theta2) / 6.0;
                    b[l2] += h * (theta1 + 2 * theta2) / 6.0;
                }
            }
        }
    }


    void use_bc3(Grid& grid, double t) {
        int p;
        int i_beg, i_end, j_beg, j_end;
        double x1, x2, y1, y2, h;
        double beta_const, u_beta1, u_beta2;
        int l1, l2;
        for (int s = 0; s < ns3; s++) {
            i_beg = grid.getIXw(bc3[s].nx1);
            i_end = grid.getIXw(bc3[s].nx2);
            j_beg = grid.getIYw(bc3[s].ny1);
            j_end = grid.getIYw(bc3[s].ny2);
            p = bc3[s].ni;
            if (i_beg == i_end) {
                int i = i_beg;
                x1 = grid.getX(i);
                for (int j = j_beg; j < j_end; j++) {
                    y1 = grid.getY(j);
                    y2 = grid.getY(j + 1);
                    h = y2 - y1;
                    L[0] = global_num(i, j);
                    L[1] = global_num(i, j + 1);
                    beta_const = beta(p);
                    u_beta1 = u_beta(p, x1, y1, t);
                    u_beta2 = u_beta(p, x1, y2, t);

                    As3[0][1] = beta_const * h / 6.0;
                    As3[0][0] = 2 * As3[0][1];
                    As3[1][0] = As3[0][1];
                    As3[1][1] = As3[0][0];
                    add_local_matrix(As3, 2);

                    b[L[0]] += beta_const * h * (2 * u_beta1 + u_beta2) / 6.0;
                    b[L[1]] += beta_const * h * (u_beta1 + 2 * u_beta2) / 6.0;
                }
            }
            else {
                int j = j_beg;
                y1 = grid.getY(j);
                for (int i = i_beg; i < i_end; i++) {
                    x1 = grid.getX(i);
                    x2 = grid.getX(i + 1);
                    h = x2 - x1;
                    L[0] = global_num(i, j);
                    L[1] = global_num(i + 1, j);
                    beta_const = beta(p);
                    u_beta1 = u_beta(p, x1, y1, t);
                    u_beta2 = u_beta(p, x2, y1, t);

                    As3[0][1] = beta_const * h / 6.0;
                    As3[0][0] = 2 * As3[0][1];
                    As3[1][0] = As3[0][1];
                    As3[1][1] = As3[0][0];
                    add_local_matrix(As3, 2);

                    b[L[0]] += beta_const * h * (2 * u_beta1 + u_beta2) / 6.0;
                    b[L[1]] += beta_const * h * (u_beta1 + 2 * u_beta2) / 6.0;
                }
            }
        }
    }

    void make_prev_q(Grid& grid, std::vector<double>& q_k, bool isU0) {
        double t_u1 = grid.getT(1);
        double x, y;
        int wi, l;
        for (int j = 0; j < Ny; j++)
        {
            y = grid.getY(j);
            for (int i = 0; i < Nx; i++)
            {
                x = grid.getX(i);
                l = global_num(i, j);
                if (i + 1 != Nx && j + 1 != Ny) {
                    wi = grid.inSubArea(i, j);
                    if (wi < 0) {
                        q_k[l] = 0;
                    }
                    else
                        isU0 ? q_k[l] = u0(wi, x, y) : q_k[l] = u1(wi, x, y, t_u1);
                }
                else
                    isU0 ? q_k[l] = u0(wi, x, y) : q_k[l] = u1(wi, x, y, t_u1);
            }
        }
    }

    void process_fict_nodes() {
        for (int i : fict_nodes) {
            di[i] = 1;
        }
    }
};