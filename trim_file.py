import os
import sys

def trim_file(input_file, bytes_to_trim=2):
    """
    Trim specified number of bytes from the end of a file.
    
    Args:
        input_file (str): Path to the input file
        bytes_to_trim (int): Number of bytes to remove from the end
    """
    # Get original file size
    original_size = os.path.getsize(input_file)
    print(f"Original file size: {original_size} bytes")
    
    if original_size <= bytes_to_trim:
        raise ValueError(f"File is too small to trim {bytes_to_trim} bytes")
    
    # Create backup
    backup_file = input_file + '.backup'
    with open(input_file, 'rb') as src, open(backup_file, 'wb') as dst:
        dst.write(src.read())
    print(f"Created backup at: {backup_file}")
    
    # Read file content except last 2 bytes
    with open(input_file, 'rb') as f:
        content = f.read()[:-bytes_to_trim]
    
    # Write trimmed content back
    with open(input_file, 'wb') as f:
        f.write(content)
    
    # Verify new size
    new_size = os.path.getsize(input_file)
    print(f"New file size: {new_size} bytes")
    print(f"Removed {original_size - new_size} bytes")

if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            print("Usage: python trim_file.py <input_raw_file>")
            print("Example: python trim_file.py data.raw")
            sys.exit(1)
            
        input_file = sys.argv[1]
        
        # Verify input file exists
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' not found")
            sys.exit(1)
            
        print(f"\nTrimming 2 bytes from: {input_file}")
        trim_file(input_file)
        print("\nTrim complete!")
        print("A backup of the original file has been created with .backup extension")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1) 