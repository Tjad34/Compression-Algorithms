# Sensor Data Compression Analysis

A tool to compare different compression algorithms (FLAC, Delta-Huffman, and LZ77) for sensor data analysis. Handles both single-axis and three-axis sensor data files.

## Project Structure

```
AllAlgo/
├── main.cpp        # Main implementation file containing all compression algorithms
└── README.md       # This file
```

## Features

- Supports both single-axis and three-axis sensor data
- Automatically detects data type based on filename (VM* files are treated as three-axis)
- Implements three compression algorithms:
  - FLAC (Free Lossless Audio Codec)
  - Delta-Huffman encoding
  - LZ77 compression
- Generates comprehensive CSV report with compression ratios
- Console output for real-time progress monitoring

## Prerequisites

1. C++ Compiler with C++17 support (g++ recommended)
2. FLAC codec installed and available in system PATH
3. Windows operating system

## Setup and Running

1. Install FLAC codec:
   - Download from: https://xiph.org/flac/download.html
   - Install and ensure it's added to system PATH
   - Verify installation by running: `flac --version`

2. Compile the program:
```bash
g++ -o compression_program AllAlgo/main.cpp -std=c++17
```

3. Run the program:
```bash
./compression_program "path/to/your/data/directory"
```

Example:
```bash
./compression_program "Data_Input(raw)"
```

## Input Requirements

- Files must contain 16-bit unsigned integer data
- Single-axis files:
  - File size must be multiple of 2 bytes
  - Named without "VM" prefix
- Three-axis files:
  - File size must be multiple of 6 bytes
  - Must start with "VM" prefix

## Output Format

1. Console Output:
   - Real-time progress updates
   - File-by-file compression results
   - Data type and size information

2. CSV Report (`compression_results.csv`):
   - Filename
   - Data Type (Single-axis/Three-axis)
   - File Size (bytes)
   - Samples per Channel
   - Number of Channels
   - FLAC Compression Ratio
   - Delta-Huffman Compression Ratio
   - LZ77 Compression Ratio

## Implementation Notes

- FLAC: Uses 16-bit depth, 16kHz sample rate
- Delta-Huffman: Optimized for sensor data with small variations
- LZ77: Uses 4096-byte window size with 16-byte look-ahead buffer
- Error handling for invalid files and missing dependencies 