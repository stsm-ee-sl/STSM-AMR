import os
import pandas as pd

# Define the directory where your AMRFinderPlus output files are stored
amrfinder_directory_path = '/home/user/Desktop/stsm/amrfinderplus'

# Specify the columns to keep from the original .tsv file
columns_to_keep = ['Contig id', 'Start', 'Stop', 'Strand', 'Gene symbol',
                   'Sequence name', 'Class', '% Coverage of reference sequence', 'Method', 'Accession of closest sequence']

# New names for the columns to be renamed
columns_rename = {
    'Class': 'Drug class',
    '% Coverage of reference sequence': 'Coverage',
    'Method': 'Genetic variation',
    'Accession of closest sequence': 'Accession nr.'
}



for file_name in os.listdir(amrfinder_directory_path):
    if file_name.endswith('results.tsv'):
        file_path = os.path.join(amrfinder_directory_path, file_name)
        sample_id = file_name.split('_')[0]
        df = pd.read_csv(file_path, sep='\t')
        df = df[columns_to_keep]
        df = df.rename(columns=columns_rename)
        df.insert(0, 'Sample_ID', sample_id)
        df.insert(1, 'Software', 'AMRFinderPlus')

        new_file_name = f'{sample_id}_amrfinder_edited.tsv'
        new_file_path = os.path.join(amrfinder_directory_path, new_file_name)

df['Genetic variation'] = df['Genetic variation'].replace({'POINTX': 'mutation', 'BLASTX': 'gene'})
def extract_mutation(row):
    if row['Genetic variation'] == 'mutation':
        parts = row['Gene symbol'].split('_')
        return parts[-1] if parts else ''
    else:
        return ''

df['Amino acid mutation'] = df.apply(extract_mutation, axis=1)

df.to_csv(new_file_path, sep='\t', index=False)

print("AMRFinderPlus files have been processed and columns renamed.")










rgi_directory_path = '/home/user/Desktop/stsm/rgi'


columns_keep = ['Contig', 'Start', 'Stop', 'Orientation', 'Best_Hit_ARO',
                   'AMR Gene Family', 'Drug Class', 'Percentage Length of Reference Sequence',
                   'SNPs_in_Best_Hit_ARO', 'Model_type']


for rgi_file_name in os.listdir(rgi_directory_path):
    if rgi_file_name.endswith('.tsv'):
        rgi_file_path = os.path.join(rgi_directory_path, rgi_file_name)
        sample_id_rgi = rgi_file_name.split('_')[0]
        df_rgi = pd.read_csv(rgi_file_path, sep='\t')
        df_rgi = df_rgi[columns_keep]

        df_rgi.insert(0, 'Sample_id', sample_id_rgi)
        df_rgi.insert(1, 'rgi', 'rgi')


        new_name = f'{sample_id_rgi}_rgi_edited.tsv'
        new_path = os.path.join(rgi_directory_path, new_name)


        df_rgi.to_csv(new_path, sep='\t', index=False)

print("RGI files have been processed.")




def convert_to_tsv(input_file, output_file):

    input_delimiter = '\t'  # Adjust if your file uses a different delimiter
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Assuming the file is already tab-separated. If not, you might need to split and join by '\t'
            outfile.write(line)



base_directory = '/home/regina/Desktop/stsm/resfinder'

for subdir, dirs, files in os.walk(base_directory):
    for file in files:
        # Check if the file is one of the target files for conversion
        if file in ['PointFinder_results.txt', 'ResFinder_results_tab.txt']:
            input_file_path = os.path.join(subdir, file)
            output_file_path = os.path.join(subdir, file.replace('.txt', '.tsv'))

            # Convert the file to TSV
            convert_to_tsv(input_file_path, output_file_path)
            print(f'Converted {input_file_path} to {output_file_path}')

import os
import pandas as pd

# Define the directory where your ResFinder output directories are stored
resfinder_directory_path = '/home/regina/Desktop/stsm/resfinder'

# Specify the columns to keep from the original .tsv file
columns_keep_resfinder = ['Contig', 'Position in contig', 'Resistance gene',
                          'Phenotype', 'Accession no.', 'Coverage']


# Function to safely split 'Position in contig' and return start, stop, and orientation
def split_position(row):
    try:
        start, stop = map(int, row['Position in contig'].split('..'))
        orientation = "-" if start > stop else "+"
    except ValueError:  # Handles cases where split fails or conversion to int fails
        start, stop, orientation = None, None, None
    return pd.Series([start, stop, orientation])



for sample_dir in os.listdir(resfinder_directory_path):
    sample_dir_path = os.path.join(resfinder_directory_path, sample_dir)


    if os.path.isdir(sample_dir_path):
 
        sample_id_resfinder = sample_dir

        for resfinder_file_name in os.listdir(sample_dir_path):
            if resfinder_file_name.endswith('tab.tsv'):

                resfinder_file_path = os.path.join(sample_dir_path, resfinder_file_name)

              
                df_resfinder = pd.read_csv(resfinder_file_path, sep='\t')

                df_resfinder = df_resfinder[columns_keep_resfinder]

                df_resfinder[['start', 'stop', 'orientation']] = df_resfinder.apply(split_position, axis=1)

                df_resfinder.insert(0, 'Sample_id', sample_id_resfinder)
                df_resfinder.insert(1, 'Tool', 'resfinder')

                new_name_resfinder = f'{sample_id_resfinder}_resfinder_edited.tsv'
                new_path_resfinder = os.path.join(sample_dir_path, new_name_resfinder)

                df_resfinder.to_csv(new_path_resfinder, sep='\t', index=False)

print("ResFinder files have been processed and updated with start, stop, and orientation columns.")
