def add_headers_to_epitopes(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    with open(filename, 'w') as file:
        for line in lines:
            stripped_line = line.strip()
            file.write(f">{stripped_line}\n")
            file.write(f"{stripped_line}\n")

# Change 'epitopes.txt' to your file path if necessary
add_headers_to_epitopes(r"/Users/cliffordkim/Library/CloudStorage/OneDrive-Personal/Education/Georgetown University/Georgetown University School of Medicine/Research/Targeted Personalized Cancer Vaccine Project/Datasets/Manual/GSM6996090 - Cancer Stage IA1/MHC2 Epitopes - Cancer Stage IA1.txt")