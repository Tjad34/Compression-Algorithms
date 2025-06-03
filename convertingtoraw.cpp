#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.csv output.raw\n", argv[0]);
        return 1;
    }
    std::ifstream csv(argv[1]);
    std::ofstream raw(argv[2], std::ios::binary);
    std::string line;
    while (std::getline(csv, line)) {
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ',')) {
            int value = std::stoi(item);
            int16_t sample = static_cast<int16_t>(value);
            raw.write(reinterpret_cast<const char*>(&sample), sizeof(sample));
        }
    }
    return 0;
}