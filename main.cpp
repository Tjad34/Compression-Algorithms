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

    // Convert delta values back to original data
    std::vector<uint16_t> reconstructFromDelta(const std::vector<int16_t>& deltaValues) {
        std::vector<uint16_t> originalData(deltaValues.size());
        int32_t accumulator = deltaValues[0];  // Use 32-bit to handle potential overflow
        originalData[0] = static_cast<uint16_t>(accumulator);
        
        for (size_t i = 1; i < deltaValues.size(); i++) {
            accumulator += deltaValues[i];
            originalData[i] = static_cast<uint16_t>(accumulator);
        }
        return originalData;
    }

    // Build Huffman tree and get codes
    std::unordered_map<int16_t, std::string> buildHuffmanCodes(const std::vector<int16_t>& deltaValues) {
        // Calculate frequencies
        std::unordered_map<int16_t, size_t> frequencies;
        for (const auto& value : deltaValues) {
            frequencies[value]++;
        }

        // Handle case where all values are the same
        if (frequencies.size() == 1) {
            auto value = frequencies.begin()->first;
            return {{value, "0"}};  // Use single bit for repeated value
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
        try {
            if (data.empty()) {
                throw std::runtime_error("Input data is empty");
            }

            std::cout << "Computing delta values...\n";
            std::vector<int16_t> deltaValues = computeDeltaValues(data);
            
            std::cout << "Building Huffman codes...\n";
            auto huffmanCodes = buildHuffmanCodes(deltaValues);
            
            std::cout << "Creating output buffer...\n";
            std::vector<uint8_t> compressedData;
            
            // Store magic number to identify DHC format
            const uint16_t MAGIC = 0x4448; // "DH"
            compressedData.push_back(static_cast<uint8_t>(MAGIC >> 8));
            compressedData.push_back(static_cast<uint8_t>(MAGIC & 0xFF));
            
            // Store original data size
            uint32_t originalSize = static_cast<uint32_t>(data.size());
            compressedData.push_back(static_cast<uint8_t>((originalSize >> 24) & 0xFF));
            compressedData.push_back(static_cast<uint8_t>((originalSize >> 16) & 0xFF));
            compressedData.push_back(static_cast<uint8_t>((originalSize >> 8) & 0xFF));
            compressedData.push_back(static_cast<uint8_t>(originalSize & 0xFF));
            
            // Store number of Huffman table entries
            uint16_t numEntries = static_cast<uint16_t>(huffmanCodes.size());
            compressedData.push_back(static_cast<uint8_t>(numEntries >> 8));
            compressedData.push_back(static_cast<uint8_t>(numEntries & 0xFF));
            
            std::cout << "Storing Huffman table (" << numEntries << " entries)...\n";
            // Store Huffman table
            for (const auto& pair : huffmanCodes) {
                // Store value
                compressedData.push_back(static_cast<uint8_t>(pair.first >> 8));
                compressedData.push_back(static_cast<uint8_t>(pair.first & 0xFF));
                
                // Store code length
                uint8_t codeLength = static_cast<uint8_t>(pair.second.length());
                compressedData.push_back(codeLength);
                
                // Store code bits
                std::string codeBits = pair.second;
                while (codeBits.length() % 8 != 0) {
                    codeBits += '0';  // Pad with zeros
                }
                
                for (size_t i = 0; i < codeBits.length(); i += 8) {
                    std::string byte = codeBits.substr(i, 8);
                    compressedData.push_back(static_cast<uint8_t>(std::bitset<8>(byte).to_ulong()));
                }
            }
            
            std::cout << "Encoding data...\n";
            // Encode data
            std::string encodedBits;
            for (const auto& value : deltaValues) {
                auto it = huffmanCodes.find(value);
                if (it == huffmanCodes.end()) {
                    throw std::runtime_error("Value not found in Huffman table: " + std::to_string(value));
                }
                encodedBits += it->second;
            }
            
            // Convert encoded data bit string to bytes
            std::cout << "Converting encoded bits to bytes...\n";
            std::cout << "Total encoded bits: " << encodedBits.length() << "\n";
            
            // Store number of valid bits in the last byte
            uint8_t lastByteBits = encodedBits.length() % 8;
            if (lastByteBits == 0) lastByteBits = 8;
            compressedData.push_back(lastByteBits);
            
            // Convert bit string to bytes
            for (size_t i = 0; i < encodedBits.length(); i += 8) {
                std::string byte = encodedBits.substr(i, std::min(size_t(8), encodedBits.length() - i));
                while (byte.length() < 8) {
                    byte += '0';  // Pad with zeros
                }
                compressedData.push_back(static_cast<uint8_t>(std::bitset<8>(byte).to_ulong()));
            }
            
            compressedSize = compressedData.size();
            std::cout << "Compression complete. Original size: " << data.size() * sizeof(uint16_t) 
                      << " bytes, Compressed size: " << compressedSize << " bytes\n";
            return compressedData;
            
        } catch (const std::exception& e) {
            std::cerr << "Error during compression: " << e.what() << "\n";
            throw;
        }
    }

    // Decompress data using Delta-Huffman coding
    std::vector<uint16_t> decompress(const std::vector<uint8_t>& compressedData) {
        try {
            if (compressedData.size() < 9) {  // Minimum size for header
                throw std::runtime_error("Compressed data too small");
            }

            size_t pos = 0;
            
            // Read and verify magic number
            uint16_t magic = (static_cast<uint16_t>(compressedData[pos]) << 8) | 
                             static_cast<uint16_t>(compressedData[pos + 1]);
            pos += 2;
            
            if (magic != 0x4448) {
                throw std::runtime_error("Invalid DHC file format");
            }
            
            // Read original data size
            uint32_t originalSize = (static_cast<uint32_t>(compressedData[pos]) << 24) |
                                  (static_cast<uint32_t>(compressedData[pos + 1]) << 16) |
                                  (static_cast<uint32_t>(compressedData[pos + 2]) << 8) |
                                   static_cast<uint32_t>(compressedData[pos + 3]);
            pos += 4;
            
            // Read number of Huffman table entries
            uint16_t numEntries = (static_cast<uint16_t>(compressedData[pos]) << 8) | 
                                 static_cast<uint16_t>(compressedData[pos + 1]);
            pos += 2;
            
            if (numEntries == 0) {
                throw std::runtime_error("Invalid Huffman table: no entries");
            }
            
            std::cout << "Reading Huffman table (" << numEntries << " entries)...\n";
            
            // Reconstruct Huffman codes
            std::unordered_map<std::string, int16_t> huffmanDecodeMap;
            for (uint16_t i = 0; i < numEntries; i++) {
                if (pos + 3 > compressedData.size()) {
                    throw std::runtime_error("Unexpected end of compressed data in Huffman table");
                }

                // Read value
                int16_t value = (static_cast<int16_t>(compressedData[pos]) << 8) | 
                               static_cast<int16_t>(compressedData[pos + 1]);
                pos += 2;
                
                // Read code length
                uint8_t codeLength = compressedData[pos++];
                if (codeLength == 0) {
                    throw std::runtime_error("Invalid Huffman code length: 0");
                }
                
                // Read code bits
                size_t bytesToRead = (codeLength + 7) / 8;  // Round up to nearest byte
                if (pos + bytesToRead > compressedData.size()) {
                    throw std::runtime_error("Unexpected end of compressed data in Huffman codes");
                }

                std::string code;
                for (size_t j = 0; j < bytesToRead; j++) {
                    std::bitset<8> bits(compressedData[pos++]);
                    code += bits.to_string();
                }
                code = code.substr(0, codeLength);  // Trim padding
                
                huffmanDecodeMap[code] = value;
            }
            
            if (pos >= compressedData.size()) {
                throw std::runtime_error("Unexpected end of compressed data before encoded data");
            }

            // Read number of valid bits in last byte
            uint8_t lastByteBits = compressedData[pos++];
            if (lastByteBits == 0 || lastByteBits > 8) {
                throw std::runtime_error("Invalid number of bits in last byte: " + std::to_string(lastByteBits));
            }
            
            std::cout << "Decoding compressed data...\n";
            
            // Decode compressed data
            std::vector<int16_t> deltaValues;
            std::string currentCode;
            
            size_t totalBytes = compressedData.size() - pos;
            for (size_t i = 0; i < totalBytes; i++) {
                std::bitset<8> bits(compressedData[pos + i]);
                std::string bitString = bits.to_string();
                
                // For the last byte, only use valid bits
                if (i == totalBytes - 1) {
                    bitString = bitString.substr(8 - lastByteBits);
                }
                
                for (char bit : bitString) {
                    currentCode += bit;
                    auto it = huffmanDecodeMap.find(currentCode);
                    if (it != huffmanDecodeMap.end()) {
                        deltaValues.push_back(it->second);
                        currentCode.clear();
                        if (deltaValues.size() == originalSize) break;
                    }
                    
                    // Check for invalid codes
                    if (currentCode.length() > 32) {  // Reasonable maximum code length
                        throw std::runtime_error("Invalid Huffman code encountered");
                    }
                }
            }
            
            if (deltaValues.size() != originalSize) {
                throw std::runtime_error("Decompression failed: got " + 
                                       std::to_string(deltaValues.size()) + 
                                       " values, expected " + 
                                       std::to_string(originalSize));
            }
            
            std::cout << "Converting delta values back to original data...\n";
            return reconstructFromDelta(deltaValues);
            
        } catch (const std::exception& e) {
            std::cerr << "Error during decompression: " << e.what() << "\n";
            throw;
        }
    }

    // Save compressed data to a DHC file
    void saveToFile(const std::vector<uint8_t>& compressedData, const std::string& outputPath) {
        std::cout << "Creating DHC file at: " << outputPath << "\n";
        
        // Create directories if they don't exist
        std::filesystem::path filePath(outputPath);
        std::filesystem::create_directories(filePath.parent_path());
        
        std::ofstream outFile(outputPath, std::ios::binary);
        if (!outFile) {
            throw std::runtime_error("Failed to create output file: " + outputPath);
        }
        
        // Write the compressed data
        outFile.write(reinterpret_cast<const char*>(compressedData.data()), 
                     compressedData.size());
        
        if (!outFile) {
            outFile.close();
            throw std::runtime_error("Failed to write data to file: " + outputPath);
        }
        
        outFile.close();
        
        // Verify file was created and has correct size
        if (!std::filesystem::exists(outputPath)) {
            throw std::runtime_error("File was not created: " + outputPath);
        }
        
        auto fileSize = std::filesystem::file_size(outputPath);
        if (fileSize != compressedData.size()) {
            throw std::runtime_error("File size mismatch. Expected: " + 
                                   std::to_string(compressedData.size()) + 
                                   ", Got: " + std::to_string(fileSize));
        }
        
        std::cout << "Successfully wrote " << fileSize << " bytes to " << outputPath << "\n";
    }

    // Read and decompress from a DHC file
    std::vector<uint16_t> decompressFromFile(const std::string& inputPath) {
        std::ifstream inFile(inputPath, std::ios::binary);
        if (!inFile) {
            throw std::runtime_error("Failed to open input file: " + inputPath);
        }

        // Read the entire file into a vector
        std::vector<uint8_t> compressedData;
        inFile.seekg(0, std::ios::end);
        size_t fileSize = inFile.tellg();
        inFile.seekg(0, std::ios::beg);
        
        compressedData.resize(fileSize);
        inFile.read(reinterpret_cast<char*>(compressedData.data()), fileSize);
        inFile.close();

        // Use existing decompress function
        return decompress(compressedData);
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
        std::string csvFile = "compression_results_new.csv";
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
            try {
                if (!entry.is_regular_file()) continue;
                
                const std::string originalFile = entry.path().string();
                // Skip .results files
                if (ends_with(originalFile, ".results")) continue;
                
                currentFile++;
                
                std::cout << "\nProcessing file " << currentFile << " of " << totalFiles << ": " 
                          << entry.path().filename().string() << "\n\n";
                std::cout << "Compressing " << originalFile << "...\n";

                // Determine if file is 3-axis based on filename
                std::string filename = entry.path().filename().string();
                bool is3Axis = starts_with(filename, "VM");
                const int CHANNELS = is3Axis ? 3 : 1;
                
                std::cout << "Data type: " << (CHANNELS == 1 ? "1-axis" : "3-axis") << " data\n\n";

                const std::string inputFilePath = "input.raw";
                const size_t BYTES_PER_SAMPLE = CHANNELS * sizeof(uint16_t);

                // Copy input file to standardized name
                std::cout << "Opening input file...\n";
                std::ifstream src(originalFile, std::ios::binary);
                if (!src) {
                    throw std::runtime_error("Could not open input file: " + originalFile);
                }

                // Check file size
                src.seekg(0, std::ios::end);
                const size_t actualSize = src.tellg();
                src.seekg(0, std::ios::beg);

                if (actualSize == 0) {
                    throw std::runtime_error("File is empty: " + originalFile);
                }

                // Validate file size
                if (actualSize % sizeof(uint16_t) != 0) {
                    throw std::runtime_error("Invalid file size: must be multiple of 2 bytes");
                }

                if (CHANNELS == 3 && actualSize % BYTES_PER_SAMPLE != 0) {
                    throw std::runtime_error("Invalid file size for three-axis data: must be multiple of 6 bytes");
                }

                std::cout << "Creating temporary file...\n";
                std::ofstream dst(inputFilePath, std::ios::binary);
                if (!dst) {
                    throw std::runtime_error("Could not create temporary file");
                }
                dst << src.rdbuf();
                src.close();
                dst.close();

                // Calculate total samples
                const int totalSamples = static_cast<int>(actualSize / BYTES_PER_SAMPLE);

                // Read raw data
                std::cout << "Reading data...\n";
                std::vector<uint16_t> rawData(totalSamples * CHANNELS);
                std::ifstream file(inputFilePath, std::ios::binary);
                if (!file) {
                    throw std::runtime_error("Failed to open temporary file for reading");
                }
                file.read(reinterpret_cast<char*>(rawData.data()), actualSize);
                if (!file) {
                    throw std::runtime_error("Failed to read data from temporary file");
                }
                file.close();

                // Run compressions
                std::cout << "Running compressions...\n";
                double flacRatio = 0.0, dhRatio = 0.0, lz77Ratio = 0.0;

                // Run FLAC compression
                size_t flacCompressedSize = runFlac(rawData, CHANNELS);
                flacRatio = static_cast<double>(actualSize) / flacCompressedSize;

                // Run Delta-Huffman compression
                std::cout << "Running Delta-Huffman compression...\n";
                DeltaHuffmanCompressor dhCompressor;
                size_t dhCompressedSize;
                auto dhCompressedData = dhCompressor.compress(rawData, dhCompressedSize);
                dhRatio = static_cast<double>(actualSize) / dhCompressedSize;

                // Save DHC compressed data to file
                std::cout << "Saving DHC file...\n";
                std::string dhcOutputPath = std::filesystem::path(originalFile).string() + ".dhc";
                bool dhcVerified = false;
                try {
                    dhCompressor.saveToFile(dhCompressedData, dhcOutputPath);
                    std::cout << "DHC file saved successfully: " << dhcOutputPath << "\n";
                    
                    // Verify the file exists and has correct size
                    if (!std::filesystem::exists(dhcOutputPath)) {
                        throw std::runtime_error("DHC file was not created: " + dhcOutputPath);
                    }
                    
                    auto dhcFileSize = std::filesystem::file_size(dhcOutputPath);
                    std::cout << "DHC file size: " << dhcFileSize << " bytes\n";
                    
                    // Verify compression
                    std::cout << "Verifying compression...\n";
                    auto decompressedData = dhCompressor.decompressFromFile(dhcOutputPath);
                    dhcVerified = (decompressedData == rawData);
                    size_t decompressedSize = decompressedData.size() * sizeof(uint16_t);
                    
                    if (!dhcVerified) {
                        std::cout << "\nWarning: Decompression verification failed!\n";
                        std::cout << "Original size: " << rawData.size() << " samples\n";
                        std::cout << "Decompressed size: " << decompressedData.size() << " samples\n";
                        // Print first few values for comparison
                        std::cout << "First few values comparison:\n";
                        for (size_t i = 0; i < std::min(size_t(10), rawData.size()); ++i) {
                            std::cout << "Original[" << i << "]: " << rawData[i] 
                                    << ", Decompressed[" << i << "]: " 
                                    << (i < decompressedData.size() ? decompressedData[i] : 0) << "\n";
                        }
                    } else {
                        std::cout << "Decompression verification successful!\n";
                        std::cout << "All " << decompressedData.size() << " samples match the original data.\n";
                    }

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
                    std::cout << "\nResults for " << filename << ":\n";
                    std::cout << "File size: " << actualSize << " bytes (" << std::fixed << std::setprecision(2) 
                              << (actualSize / 1024.0) << " KB)\n";
                    std::cout << "Data type: " << (CHANNELS == 1 ? "Single-axis" : "Three-axis") << "\n";
                    std::cout << "Number of samples per channel: " << totalSamples << "\n";
                    std::cout << "Number of channels: " << CHANNELS << "\n";
                    std::cout << "Bytes per sample point: " << BYTES_PER_SAMPLE << "\n\n";
                    
                    std::cout << "Compression Results:\n";
                    std::cout << "FLAC Ratio: " << std::fixed << std::setprecision(3) << flacRatio << "\n";
                    std::cout << "Delta-Huffman:\n";
                    std::cout << "  Compression Ratio: " << dhRatio << "\n";
                    std::cout << "  Original Size: " << actualSize << " bytes\n";
                    std::cout << "  Compressed Size: " << dhCompressedSize << " bytes\n";
                    std::cout << "  Decompressed Size: " << decompressedSize << " bytes\n";
                    std::cout << "  Decompression Verified: " << (dhcVerified ? "Yes" : "No") << "\n";
                    std::cout << "  DHC File: " << dhcOutputPath << "\n";
                    std::cout << "LZ77 Ratio: " << lz77Ratio << "\n";
                    std::cout << "----------------------------------------\n";
                    
                } catch (const std::exception& e) {
                    std::cout << "Error during DHC compression/verification: " << e.what() << "\n";
                    
                    // Still write to CSV even if there's a decompression error
                    csv << std::fixed << std::setprecision(3);
                    csv << filename << ","
                        << (CHANNELS == 1 ? "Single-axis" : "Three-axis") << ","
                        << actualSize << ","
                        << totalSamples << ","
                        << CHANNELS << ","
                        << flacRatio << ","
                        << dhRatio << ","  // Use the compression ratio even if decompression failed
                        << lz77Ratio << "\n";
                    
                    std::cout << "Added to CSV with compression ratios.\n";
                    std::cout << "Continuing with next file...\n\n";
                    continue;
                }

                // Cleanup temporary file
                std::remove(inputFilePath.c_str());

            } catch (const std::exception& e) {
                std::cerr << "Error processing file: " << e.what() << "\n";
                std::cerr << "Continuing with next file...\n\n";
                continue;
            }
        }
        
        csv.close();
        std::cout << "\nResults have been saved to " << csvFile << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}