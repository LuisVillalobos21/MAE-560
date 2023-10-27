#include "ReadData.hpp"

Eigen::VectorXd readCSVtoEigenVector(const std::string& filename) {
    std::ifstream file(filename);

    // Check if file was opened correctly
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    std::vector<double> vectorEntries;

    // Reading only the first line assuming a single row of data
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        // Split values based on space
        while (ss >> value) {
            vectorEntries.push_back(std::stod(value));
        }
    }
    else {
        throw std::runtime_error("File is empty");
    }

    Eigen::VectorXd vec(vectorEntries.size());

    for (int i = 0; i < vectorEntries.size(); ++i) {
        vec(i) = vectorEntries[i];
    }

    return vec;
}

