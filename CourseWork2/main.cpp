#include "FEM.h"
#include "replace_dots.h"  // Для замены точек на запятые

int main() {
	std::string fname_res = "Result.txt";
	std::string fname_grid = "Grid.txt";  // Файл в котором выводится координаты сетки
	std::string test_folder = "research_uniform_1";

	std::ofstream ofile_grid;
	ofile_grid.open(fname_grid);
	Grid grid(test_folder);
	grid.print_grid(ofile_grid);
	ofile_grid.close();
	replace_dots(fname_grid);

	std::ofstream ofile;
	ofile.open(fname_res);
	FEM fem(grid, test_folder);

	int Nt = grid.getNt();
	double t_2, t_1, t;
	for (int j = 2; j < Nt; j++) {
		t = grid.getT(j);
		t_1 = grid.getT(j - 1);
		t_2 = grid.getT(j - 2);
		fem.solve(grid, t_2, t_1, t);
		fem.print_vector(ofile, t);
		fem.update();
	}
	ofile.close();
	replace_dots(fname_res);
}