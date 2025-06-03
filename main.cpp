#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <windows.h>  // For GetModuleFileName
#include <string>
#include <queue>
#include <unordered_map>
#include <memory>
#include <bitset>
#include <algorithm>
#include <iomanip>

// Huffman Tree Node
struct HuffmanNode {
    int16_t value;
    size_t frequency;
    std::shared_ptr<HuffmanNode> left, right;
    
    HuffmanNode(int16_t val, size_t freq) : value(val), frequency(freq), left(nullptr), right(nullptr) {}
};

// Custom comparator for priority queue
struct CompareNodes {
    bool operator()(const std::shared_ptr<HuffmanNode>& a, const std::shared_ptr<HuffmanNode>& b) {
        return a->frequency > b->frequency;
    }
};

class DeltaHuffmanCompressor {
private:
    // Convert raw data to delta values
    std::vector<int16_t> computeDeltaValues(const std::vector<uint16_t>& data) {
        std::vector<int16_t> deltaValues(data.size());
        deltaValues[0] = static_cast<int16_t>(data[0]);
        
        for (size_t i = 1; i < data.size(); i++) {
            deltaValues[i] = static_cast<int16_t>(data[i]) - static_cast<int16_t>(data[i-1]);
        }
        return deltaValues;
    }

    // Build Huffman tree and get codes
    std::unordered_map<int16_t, std::string> buildHuffmanCodes(const std::vector<int16_t>& deltaValues) {
        // Calculate frequencies
        std::unordered_map<int16_t, size_t> frequencies;
        for (const auto& value : deltaValues) {
            frequencies[value]++;
        }

        // Create priority queue
        std::priority_queue<std::shared_ptr<HuffmanNode>, 
                          std::vector<std::shared_ptr<HuffmanNode>>, 
                          CompareNodes> pq;

        // Add nodes to priority queue
        for (const auto& pair : frequencies) {
            pq.push(std::make_shared<HuffmanNode>(pair.first, pair.second));
        }

        // Build Huffman tree
        while (pq.size() > 1) {
            auto left = pq.top(); pq.pop();
            auto right = pq.top(); pq.pop();
            
            auto parent = std::make_shared<HuffmanNode>(0, left->frequency + right->frequency);
            parent->left = left;
            parent->right = right;
            
            pq.push(parent);
        }

        // Generate Huffman codes
        std::unordered_map<int16_t, std::string> huffmanCodes;
        if (!pq.empty()) {
            generateCodes(pq.top(), "", huffmanCodes);
        }
        
        return huffmanCodes;
    }

    void generateCodes(const std::shared_ptr<HuffmanNode>& node, 
                      std::string code,
                      std::unordered_map<int16_t, std::string>& codes) {
        if (!node) return;
        
        if (!node->left && !node->right) {
            codes[node->value] = code;
        }
        
        generateCodes(node->left, code + "0", codes);
        generateCodes(node->right, code + "1", codes);
    }

public:
    // Compress data using Delta-Huffman coding
    std::vector<uint8_t> compress(const std::vector<uint16_t>& data, size_t& compressedSize) {
        // Step 1: Compute delta values
        std::vector<int16_t> deltaValues = computeDeltaValues(data);
        
        // Step 2: Build Huffman codes
        auto huffmanCodes = buildHuffmanCodes(deltaValues);
        
        // Step 3: Encode data
        std::string encodedBits;
        for (const auto& value : deltaValues) {
            encodedBits += huffmanCodes[value];
        }
        
        // Convert bit string to bytes
        std::vector<uint8_t> compressedData;
        for (size_t i = 0; i < encodedBits.length(); i += 8) {
            std::string byte = encodedBits.substr(i, 8);
            byte.resize(8, '0');  // Pad with zeros if needed
            compressedData.push_back(static_cast<uint8_t>(std::bitset<8>(byte).to_ulong()));
        }
        
        // Store Huffman table size and codes
        size_t tableSize = 0;
        for (const auto& pair : huffmanCodes) {
            tableSize += sizeof(int16_t) + 1 + pair.second.length() / 8 + 1;
        }
        
        compressedSize = compressedData.size() + tableSize;
        return compressedData;
    }
};

// LZ77 Compression Implementation
class LZ77Compressor {
private:
    struct Match {
        size_t offset;
        size_t length;
        uint16_t nextSymbol;

        Match(size_t o, size_t l, uint16_t n) : offset(o), length(l), nextSymbol(n) {}
    };

    static const size_t WINDOW_SIZE = 4096;  // Search buffer size
    static const size_t LOOKAHEAD_SIZE = 16; // Look-ahead buffer size

    Match findLongestMatch(const std::vector<uint16_t>& data, size_t currentPos) {
        size_t maxLength = 0;
        size_t maxOffset = 0;
        
        // Define search and look-ahead buffer boundaries
        size_t searchStart = (currentPos >= WINDOW_SIZE) ? currentPos - WINDOW_SIZE : 0;
        size_t lookaheadEnd = std::min(currentPos + LOOKAHEAD_SIZE, data.size());
        
        // Search for matches
        for (size_t i = searchStart; i < currentPos; ++i) {
            size_t length = 0;
            while (currentPos + length < lookaheadEnd && 
                   data[i + length] == data[currentPos + length] &&
                   length < LOOKAHEAD_SIZE) {
                length++;
            }
            
            if (length > maxLength) {
                maxLength = length;
                maxOffset = currentPos - i;
            }
        }
        
        uint16_t nextSymbol = (currentPos + maxLength < data.size()) ? 
                             data[currentPos + maxLength] : 0;
        
        return Match(maxOffset, maxLength, nextSymbol);
    }

public:
    std::vector<uint8_t> compress(const std::vector<uint16_t>& data, size_t& compressedSize) {
        std::vector<Match> matches;
        size_t pos = 0;
        
        while (pos < data.size()) {
            Match match = findLongestMatch(data, pos);
            matches.push_back(match);
            pos += match.length + 1;
        }
        
        // Calculate compressed size and create output
        // Format: [2 bytes offset][1 byte length][2 bytes next symbol]
        compressedSize = matches.size() * 5;  // 5 bytes per match
        std::vector<uint8_t> compressed;
        compressed.reserve(compressedSize);
        
        for (const auto& match : matches) {
            // Store offset (2 bytes)
            compressed.push_back(static_cast<uint8_t>(match.offset >> 8));
            compressed.push_back(static_cast<uint8_t>(match.offset & 0xFF));
            
            // Store length (1 byte)
            compressed.push_back(static_cast<uint8_t>(match.length));
            
            // Store next symbol (2 bytes)
            compressed.push_back(static_cast<uint8_t>(match.nextSymbol >> 8));
            compressed.push_back(static_cast<uint8_t>(match.nextSymbol & 0xFF));
        }
        
        return compressed;
    }
};

// Function to run FLAC compression
size_t runFlac(const std::vector<uint16_t>& data, int channels) {
    // Create temporary raw input file
    const std::string inputPath = "temp_flac_input.raw";
    {
        std::ofstream outFile(inputPath, std::ios::binary);
        if (!outFile) {
            throw std::runtime_error("Failed to create temporary raw file");
        }
        outFile.write(reinterpret_cast<const char*>(data.data()), 
                     data.size() * sizeof(uint16_t));
    }

    // Get current executable directory
    char exePath[MAX_PATH];
    GetModuleFileNameA(NULL, exePath, MAX_PATH);
    std::string exeDir(exePath);
    exeDir = exeDir.substr(0, exeDir.find_last_of("\\/"));

    // Create output FLAC file path
    const std::string outputPath = exeDir + "\\temp_flac_output.flac";

    // Build FLAC command with optimized parameters
    std::string cmd = "flac --force-raw-format "
                      "--endian=little "
                      "--sign=unsigned "  // Changed to unsigned
                      "--channels=" + std::to_string(channels) + " "
                      "--bps=16 "
                      "--sample-rate=16000 "
                      "--compression-level-8 "  // Maximum compression
                      "--lax "  // Allow more flexible encoding
                      "-s -f "
                      "\"" + inputPath + "\" "
                      "-o \"" + outputPath + "\"";

    // Execute FLAC command
    int result = std::system(cmd.c_str());
    
    // Cleanup input file
    std::remove(inputPath.c_str());

    // Check result
    if (result != 0) {
        std::string errorMsg = "FLAC compression failed. Possible reasons:\n";
        errorMsg += "1. FLAC not installed or not in PATH\n";
        errorMsg += "2. Missing libFLAC.dll\n";
        errorMsg += "3. Invalid parameters\n";
        errorMsg += "4. File permission issues\n";
        errorMsg += "Return code: " + std::to_string(result);
        throw std::runtime_error(errorMsg);
    }

    // Get compressed file size
    size_t compressedSize = std::filesystem::file_size(outputPath);
    
    // Cleanup output file
    std::remove(outputPath.c_str());
    
    return compressedSize;
}

// Helper functions for string operations (C++17 compatible)
bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) return false;
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

bool starts_with(const std::string& str, const std::string& prefix) {
    if (str.length() < prefix.length()) return false;
    return str.compare(0, prefix.length(), prefix) == 0;
}

int main(int argc, char* argv[]) {
    try {
        // Check if directory path is provided
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <input_directory>\n";
            std::cerr << "Note: Files starting with 'VM' will be processed as 3-axis data, others as 1-axis\n";
            return 1;
        }

        // Configuration
        const std::string inputDir = argv[1];
        
        // Create output CSV file
        std::string csvFile = "compression_results.csv";
        std::ofstream csv(csvFile);
        if (!csv) {
            throw std::runtime_error("Could not create output CSV file");
        }
        
        // Write CSV header
        csv << "Filename,Data Type,File Size (bytes),Samples per Channel,Channels,FLAC Ratio,Delta-Huffman Ratio,LZ77 Ratio\n";
        
        // Get total number of files for progress tracking
        size_t totalFiles = 0;
        for (const auto& entry : std::filesystem::directory_iterator(inputDir)) {
            if (entry.is_regular_file() && !ends_with(entry.path().string(), ".results")) {
                totalFiles++;
            }
        }
        
        size_t currentFile = 0;

        // Process each file in the directory
        for (const auto& entry : std::filesystem::directory_iterator(inputDir)) {
            if (!entry.is_regular_file()) continue;
            
            const std::string originalFile = entry.path().string();
            // Skip .results files
            if (ends_with(originalFile, ".results")) continue;
            
            currentFile++;
            
            // Determine if file is 3-axis based on filename
            std::string filename = entry.path().filename().string();
            bool is3Axis = starts_with(filename, "VM");
            const int CHANNELS = is3Axis ? 3 : 1;
            
            std::cout << "\nProcessing file " << currentFile << " of " << totalFiles << ": " 
                      << entry.path().filename().string() << "\n\n";
            std::cout << "Compressing " << originalFile << "...\n";
            std::cout << "Data type: " << (CHANNELS == 1 ? "1-axis" : "3-axis") << " data\n\n";

            const std::string inputFilePath = "input.raw";
            const size_t BYTES_PER_SAMPLE = CHANNELS * sizeof(uint16_t);

            // Copy input file to standardized name
            std::ifstream src(originalFile, std::ios::binary);
            if (!src) {
                std::cerr << "Could not open input file: " << originalFile << ". Skipping...\n";
                continue;
            }

            // Check file size
            src.seekg(0, std::ios::end);
            const size_t actualSize = src.tellg();
            src.seekg(0, std::ios::beg);

            if (actualSize == 0) {
                std::cerr << "File is empty: " << originalFile << ". Skipping...\n";
                continue;
            }

            // Validate file size based on number of channels
            if (actualSize % sizeof(uint16_t) != 0) {
                std::cerr << "Invalid file size for " << originalFile << ": " << actualSize << 
                         " bytes. File size must be a multiple of 2 bytes. Skipping...\n";
                continue;
            }

            if (CHANNELS == 3 && actualSize % BYTES_PER_SAMPLE != 0) {
                std::cerr << "Invalid file size for three-axis data in " << originalFile << ": " << actualSize << 
                         " bytes. File size must be a multiple of 6 bytes. Skipping...\n";
                continue;
            }

            std::ofstream dst(inputFilePath, std::ios::binary);
            if (!dst) {
                std::cerr << "Could not create temporary file. Skipping " << originalFile << "...\n";
                continue;
            }
            dst << src.rdbuf();
            src.close();
            dst.close();
            
            // Calculate total samples (per channel)
            const int totalSamples = static_cast<int>(actualSize / BYTES_PER_SAMPLE);

            // Read raw data
            std::vector<uint16_t> rawData(totalSamples * CHANNELS);
            std::ifstream file(inputFilePath, std::ios::binary);
            file.read(reinterpret_cast<char*>(rawData.data()), actualSize);
            file.close();

            // Run compressions and store results
            double flacRatio = 0.0, dhRatio = 0.0, lz77Ratio = 0.0;
            
            // Run FLAC compression
            size_t flacCompressedSize = runFlac(rawData, CHANNELS);
            flacRatio = static_cast<double>(actualSize) / flacCompressedSize;

            // Run Delta-Huffman compression
            DeltaHuffmanCompressor dhCompressor;
            size_t dhCompressedSize;
            auto dhCompressedData = dhCompressor.compress(rawData, dhCompressedSize);
            dhRatio = static_cast<double>(actualSize) / dhCompressedSize;

            // Run LZ77 compression
            LZ77Compressor lz77Compressor;
            size_t lz77CompressedSize;
            auto lz77CompressedData = lz77Compressor.compress(rawData, lz77CompressedSize);
            lz77Ratio = static_cast<double>(actualSize) / lz77CompressedSize;

            // Write to CSV file
            csv << std::fixed << std::setprecision(3);
            csv << filename << ","
                << (CHANNELS == 1 ? "Single-axis" : "Three-axis") << ","
                << actualSize << ","
                << totalSamples << ","
                << CHANNELS << ","
                << flacRatio << ","
                << dhRatio << ","
                << lz77Ratio << "\n";

            // Display results to console
            std::cout << "Input file: " << originalFile << "\n";
            std::cout << "File size: " << actualSize << " bytes (" << std::fixed << std::setprecision(2) 
                      << (actualSize / 1024.0) << " KB)\n";
            std::cout << "Data type: " << (CHANNELS == 1 ? "Single-axis" : "Three-axis") << "\n";
            std::cout << "Number of samples per channel: " << totalSamples << "\n";
            std::cout << "Number of channels: " << CHANNELS << "\n";
            std::cout << "Bytes per sample point: " << BYTES_PER_SAMPLE << "\n\n";
            
            std::cout << "Compression Results:\n";
            std::cout << "FLAC Ratio: " << std::fixed << std::setprecision(3) << flacRatio << "\n";
            std::cout << "Delta-Huffman Ratio: " << dhRatio << "\n";
            std::cout << "LZ77 Ratio: " << lz77Ratio << "\n";
            std::cout << "----------------------------------------\n";

            // Cleanup temporary file
            std::remove(inputFilePath.c_str());
        }
        
        csv.close();
        std::cout << "\nResults have been saved to " << csvFile << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}