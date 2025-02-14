#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

void plot_poincare() {
    std::string energy, path;
    std::cout << "Enter the  name and energy value: ";
    std::cin >> energy;
    std::cout << "Enter the path of the .dat file: ";
    std::cin >> path;

    std::ifstream infile(path);
    if (!infile) {
        std::cerr << "Could not open " << path << "\n";
        return;
    }

    int step;
    std::cout << "Enter the step value to control the number of points: ";
    std::cin >> step;

    std::map<int, std::vector<std::pair<double, double>>> data;
    double y_val, py_val;
    int y0_index;

    std::cout << "Reading data from: " << path << std::endl;
    while (infile >> y_val >> py_val >> y0_index) {
        data[y0_index].push_back({y_val, py_val});
    }
    infile.close();

    TCanvas *c1 = new TCanvas("c1", "Poincare Section", 1600, 800);
    
    // 🔹 Usar una paleta de colores variada
    gStyle->SetPalette(kBird);  // Otras opciones: kRainBow, kCool, kDarkBodyRadiator

    double y_min = 9999, y_max = -9999, py_min = 9999, py_max = -9999;
    for (const auto &entry : data) {
        for (const auto &point : entry.second) {
            if (point.first < y_min) y_min = point.first;
            if (point.first > y_max) y_max = point.first;
            if (point.second < py_min) py_min = point.second;
            if (point.second > py_max) py_max = point.second;
        }
    }
    y_min -= 0.05; y_max += 0.05;
    py_min -= 0.05; py_max += 0.05;

    std::cout << "Enter x-axis minimum value (default " << y_min << "): ";
    std::cin >> y_min;
    std::cout << "Enter x-axis maximum value (default " << y_max << "): ";
    std::cin >> y_max;
    std::cout << "Enter y-axis minimum value (default " << py_min << "): ";
    std::cin >> py_min;
    std::cout << "Enter y-axis maximum value (default " << py_max << "): ";
    std::cin >> py_max;

    std::string title = "Poincare Section-" + energy ;
    TH1F *frame = c1->DrawFrame(y_min, py_min, y_max, py_max);
    frame->SetTitle(title.c_str());
    frame->GetXaxis()->SetTitle("y");
    frame->GetYaxis()->SetTitle("p_{y}");
    frame->Draw();

    int color_index = 0;
    bool first_graph = true;
    for (const auto &entry : data) {
        std::vector<double> y_reduced, py_reduced;
        for (size_t i = 0; i < entry.second.size(); i += step) {
            y_reduced.push_back(entry.second[i].first);
            py_reduced.push_back(entry.second[i].second);
        }

        if (y_reduced.size() < 2) {
            std::cout << "Skipping empty dataset for y0_index=" << entry.first << std::endl;
            continue;
        }

        std::cout << "Plotting y0_index=" << entry.first << " with " << y_reduced.size() << " points." << std::endl;

        TGraph *graph = new TGraph(y_reduced.size(), &y_reduced[0], &py_reduced[0]);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.2);  // Aumentar tamaño de los puntos

        // varius colors

        int r = 50 + (color_index * 100) % 205;  
        int g = 50 + (color_index * 150) % 205;
        int b = 50 + (color_index * 200) % 205;
        int color = TColor::GetColor(r, g, b);
        graph->SetMarkerColor(color);

        // 🔹 Aplicar el título solo al primer gráfico
        if (first_graph) {
            graph->SetTitle(title.c_str());
            first_graph = false;
        }

        graph->Draw("P SAME");
        color_index++;
    }

    gStyle->SetOptStat(0);
    gPad->SetGrid();

    std::string filename = "poincare-" + energy + ".png";
    c1->SaveAs(filename.c_str());
    std::cout << "Graph saved as " << filename << std::endl;
}
