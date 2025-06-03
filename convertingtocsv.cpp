#include <iostream>
#include <fstream>

int main() {
    const char* input_filename = "VM3-00260_662_1677764577.raw";    // Your raw file
    const char* output_filename = "VM3-00260_662_1677764577.csv";  // Output CSV file
    const int values_per_line = 10;               // How many values per CSV line

    std::ifstream infile(input_filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error opening input file: " << input_filename << std::endl;
        return 1;
    }

    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Error creating output file: " << output_filename << std::endl;
        return 1;
    }

    unsigned char byte;
    int count = 0;

    while (infile.read(reinterpret_cast<char*>(&byte), sizeof(byte))) {
        outfile << static_cast<int>(byte);
        count++;

        if (count % values_per_line == 0) {
            outfile << '\n';  // New line after fixed number of values
        } else {
            outfile << ',';   // Comma between values
        }
    }

    // If last line didn't end with newline, add it
    if (count % values_per_line != 0) {
        outfile << '\n';
    }

    infile.close();
    outfile.close();

    std::cout << "Conversion complete. Total bytes processed: " << count << std::endl;
    return 0;
}