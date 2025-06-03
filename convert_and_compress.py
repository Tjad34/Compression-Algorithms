import os
import sys
import csv
from datetime import datetime

def run_compression(raw_file, is_three_axis):
    """Run all compression algorithms on a raw file."""
    print(f"\nCompressing {raw_file}...")
    num_channels = 3 if is_three_axis else 1
    print(f"Data type: {num_channels}-axis data")
    
    file_size = os.path.getsize(raw_file)
    results = {
        'File': os.path.basename(raw_file),
        'Size (bytes)': file_size,
        'Type': f'{num_channels}-axis',
        'FLAC Ratio': 0,
        'Delta-Huffman Ratio': 0,
        'LZ77 Ratio': 0,
        'Status': 'Success'
    }
    
    try:
        # Run compression
        cmd = f"flac_compressor.exe {raw_file} {num_channels}"
        result = os.system(cmd)
        
        if result != 0:
            results['Status'] = 'Failed'
            print(f"Warning: Compression failed for {raw_file}")
        else:
            # Read compression ratios from results file
            results_file = raw_file + ".results"
            if os.path.exists(results_file):
                with open(results_file, 'r') as f:
                    results['FLAC Ratio'] = float(f.readline().strip())
                    results['Delta-Huffman Ratio'] = float(f.readline().strip())
                    results['LZ77 Ratio'] = float(f.readline().strip())
                os.remove(results_file)  # Clean up results file
            else:
                results['Status'] = 'Failed - No results file'
                
    except Exception as e:
        results['Status'] = f'Error: {str(e)}'
        print(f"Error processing {raw_file}: {str(e)}")
    
    print("----------------------------------------")
    return results

def process_directory(directory):
    """Process all raw files in the directory."""
    # Get all .raw files
    raw_files = [f for f in os.listdir(directory) if f.lower().endswith('.raw')]
    
    if not raw_files:
        print("No .raw files found in directory")
        return
    
    print(f"\nFound {len(raw_files)} raw files to process")
    print("----------------------------------------")
    
    # Prepare results storage
    results = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_file = os.path.join(directory, f"compression_results_{timestamp}.csv")
    
    # Process each raw file
    for i, raw_file in enumerate(raw_files, 1):
        full_path = os.path.join(directory, raw_file)
        print(f"\nProcessing file {i} of {len(raw_files)}: {raw_file}")
        
        # Check if it's three-axis data (VM prefix)
        is_three_axis = raw_file.startswith('VM')
        
        try:
            # Run compression and collect results
            result = run_compression(full_path, is_three_axis)
            results.append(result)
        except Exception as e:
            print(f"Error processing {raw_file}: {str(e)}")
            results.append({
                'File': raw_file,
                'Size (bytes)': os.path.getsize(full_path),
                'Type': '3-axis' if is_three_axis else '1-axis',
                'FLAC Ratio': 0,
                'Delta-Huffman Ratio': 0,
                'LZ77 Ratio': 0,
                'Status': f'Error: {str(e)}'
            })
            continue
    
    # Save results to CSV
    if results:
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults saved to: {csv_file}")
    
    # Print summary
    print("\nCompression Summary:")
    print("----------------------------------------")
    print(f"Total files processed: {len(raw_files)}")
    successful = sum(1 for r in results if r['Status'] == 'Success')
    print(f"Successfully compressed: {successful}")
    print(f"Failed: {len(raw_files) - successful}")
    
    # Print average ratios for successful compressions
    if successful > 0:
        avg_flac = sum(r['FLAC Ratio'] for r in results if r['Status'] == 'Success') / successful
        avg_dh = sum(r['Delta-Huffman Ratio'] for r in results if r['Status'] == 'Success') / successful
        avg_lz77 = sum(r['LZ77 Ratio'] for r in results if r['Status'] == 'Success') / successful
        
        print("\nAverage Compression Ratios:")
        print(f"FLAC: {avg_flac:.3f}")
        print(f"Delta-Huffman: {avg_dh:.3f}")
        print(f"LZ77: {avg_lz77:.3f}")

if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            print("Usage: python convert_and_compress.py <directory>")
            print("Example: python convert_and_compress.py path/to/your/files")
            print("\nThis will:")
            print("1. Find all .raw files in the directory")
            print("2. Automatically detect three-axis data (VM prefix)")
            print("3. Run all compression algorithms")
            print("4. Save results to CSV file")
            sys.exit(1)
            
        directory = sys.argv[1]
        
        if not os.path.exists(directory):
            print(f"Error: Directory '{directory}' not found")
            sys.exit(1)
            
        print(f"\nProcessing raw files in: {directory}")
        process_directory(directory)
        print("\nAll processing complete!")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)