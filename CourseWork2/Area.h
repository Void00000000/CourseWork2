#pragma once
#include <vector>
#include <fstream>
#include <string>

// ����������
struct sub_area {
    // ������� � �������� Xw � Yw
    int nx1;  
    int ny1;
    int nx2;
    int ny2;
    int ni;  // ����� ����������
};


// ����� ��������� ��������� �������
class Area
{
protected:
    int nw;  // ���������� �����������
    int nXw; // ����� ������� Xw
    int nYw; // ����� ������� Yw;
    std::vector <double> Xw, Yw;  //���������� �����������
    std::vector <sub_area> Mw;  // ������, ���������� ����������
    void read_area(std::string &test_folder) {
        std::string path = "Input_Data/" + test_folder + "/";
        std::ifstream area_file;
        area_file.open(path + "Area.txt");
        area_file >> nXw;
        Xw.resize(nXw);
        for (int i = 0; i < nXw; i++)
            area_file >> Xw[i];
        area_file >> nYw;
        Yw.resize(nYw);
        for (int i = 0; i < nYw; i++)
            area_file >> Yw[i];
        area_file >> nw;
        Mw.resize(nw);
        for (int i = 0; i < nw; i++) {
            area_file >> Mw[i].ni >> Mw[i].nx1 >> Mw[i].nx2 >> Mw[i].ny1 >> Mw[i].ny2;
            Mw[i].nx1 -= 1;
            Mw[i].nx2 -= 1;
            Mw[i].ny1 -= 1;
            Mw[i].ny2 -= 1;
        }
        area_file.close();
    }
};