import struct
import numpy as np

def csv_to_raw(csv_file):
    """Convert CSV file back to RAW format"""
    print(f"Reading CSV file: {csv_file}")
    raw_file = csv_file.rsplit('.', 1)[0] + '.raw'
    
    # Read CSV values
    values = []
    with open(csv_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            try:
                value = int(float(line.strip()))
                values.append(value)
            except ValueError as e:
                print(f"Warning: Skipping invalid value at line {line_num}: {line.strip()}")
                continue
    
    print(f"Read {len(values)} values from CSV")
    
    if not values:
        raise ValueError("No valid values found in CSV file")
    
    # Convert to numpy array for statistics
    data = np.array(values, dtype=np.uint16)
    print(f"Data statistics:")
    print(f"Min: {np.min(data)}")
    print(f"Max: {np.max(data)}")
    print(f"Mean: {np.mean(data)}")
    
    # Write to RAW file
    with open(raw_file, 'wb') as f:
        for value in values:
            f.write(struct.pack('<H', value))  # little-endian 16-bit unsigned
    
    raw_size = len(values) * 2  # Each value is 2 bytes
    print(f"Created RAW file: {raw_file}")
    print(f"RAW file size: {raw_size} bytes")
    return raw_file

if __name__ == "__main__":
    try:
        csv_file = "TX-01333_8209_1673371806.csv"
        raw_file = csv_to_raw(csv_file)
    except Exception as e:
        print(f"Error: {str(e)}") 