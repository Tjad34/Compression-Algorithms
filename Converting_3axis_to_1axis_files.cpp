#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

int main() {
    // Open input CSV file
    std::ifstream infile("data0.txt");
    if (!infile.is_open()) {
        std::cerr << "Error opening input file!" << std::endl;
        return 1;
    }

    // Create three output files
    std::ofstream file1("data0_x.txt");
    std::ofstream file2("data0_y.txt");
    std::ofstream file3("data0_z.txt");
    
    if (!file1 || !file2 || !file3) {
        std::cerr << "Error creating output files!" << std::endl;
        return 1;
    }

    int counter = 0;
    std::string line;

    // Read each line from input CSV
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string value;

        // Split line by commas
        while (std::getline(ss, value, ',')) {
            // Write to appropriate file based on counter
            switch (counter % 3) {
                case 0:
                    file1 << value << "\n";
                    break;
                case 1:
                    file2 << value << "\n";
                    break;
                case 2:
                    file3 << value << "\n";
                    break;
            }
            counter++;
        }
    }

    // Close all files
    infile.close();
    file1.close();
    file2.close();
    file3.close();

    std::cout << "Distribution complete!\n"
              << "Total values processed: " << counter << std::endl;
    return 0;
}