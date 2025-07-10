def process_antigenicity_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        
        for line in lines:
            if line.startswith('>'):
                parts = line.strip().split(' ')
                sequence = parts[0][1:]
                prediction_part = line.split("Overall Protective Antigen Prediction = ")[1]
                prediction = prediction_part.split(' ')[0]
                antigen_status = prediction_part.split("Probable ")[1].split(").")[0]
                outfile.write(f"{sequence},{prediction},{antigen_status}\n")

# Usage
input_file = r"/Users/cliffordkim/Library/CloudStorage/OneDrive-Personal/Education/Georgetown University/Georgetown University School of Medicine/Research/Targeted Personalized Cancer Vaccine Project/Datasets/Manual/GSM6996090 - Cancer Stage IA1/MHC2 Antigens - Cancer Stage IA1.txt"
output_file = r"/Users/cliffordkim/Library/CloudStorage/OneDrive-Personal/Education/Georgetown University/Georgetown University School of Medicine/Research/Targeted Personalized Cancer Vaccine Project/Datasets/Manual/GSM6996090 - Cancer Stage IA1/MHC2 Processed Antigens - Cancer Stage IA1.csv"
process_antigenicity_file(input_file, output_file)
