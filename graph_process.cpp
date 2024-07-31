#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

void processFiles(const std::string& fileA, const std::string& fileB, const std::string& outputFile) {
    std::ifstream infileA(fileA);
    std::ifstream infileB(fileB);
    std::ofstream outfile(outputFile);

    if (!infileA.is_open() || !infileB.is_open()) {
        std::cerr << "can not open infile" << std::endl;
        return;
    }

    std::string lineA, lineB;

    long long int num = 1;

    while (std::getline(infileA, lineA) && std::getline(infileB, lineB)) {
        if(num++ == 1) continue;
        std::istringstream streamA(lineA);
        std::istringstream streamB(lineB);

        std::string colA1, colA2, colA3, colA4;
        std::string colB1, colB2, colB3, colB4;

        streamA >> colA1 >> colA2 >> colA3 >> colA4;
        streamB >> colB1 >> colB2 >> colB3 >> colB4;

        if (colA2 == colB2 && colA3 == colB3) {
            outfile << colA2 << " " << colA3 << " " << colA4 << " " << colB4 << std::endl;
        } else {
            std::cerr << "data not match" << lineA << " 与 " << lineB << std::endl;
        }
    }

    infileA.close();
    infileB.close();
    outfile.close();
}

int main() {
    std::string fileA = "dataset/graph/USA-road-t.FLA.gr";
    std::string fileB = "dataset/graph/USA-road-d.FLA.gr";
    std::string outputFile = "USA-road-t-d.FLA.gr";

    processFiles(fileA, fileB, outputFile);

    std::cout << "process complete" << outputFile << std::endl;

    return 0;
}
