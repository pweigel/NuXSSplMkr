import numpy as np
import click, tqdm
import os, glob, shutil

"""
This script modified a pdf to remove the c/cbar component. This is used for computing q->c processes.
"""
    
@click.command()
@click.option('--pdfset', required=True, help='Name of the pdf set/folder name')
@click.option('--nmembers', default=1, help='Number of members')
@click.option('--overwrite', is_flag=True, help='Overwrite existing files')
def main(pdfset, nmembers, overwrite):
    lhapdf_path = os.environ['LHAPDF_DATA_PATH']
    print('Running pdf set: ' + pdfset)
    print('Number of members: ' + str(nmembers))
    print('Overwrite: ' + str(overwrite))
    for member in tqdm.tqdm(range(nmembers)):
    
        possible_files = glob.glob(lhapdf_path + '/' + pdfset + '/*.dat')
        pdf_file = lhapdf_path + '/' + pdfset + '/' + pdfset + '_' + str(member).zfill(4) + '.dat'
        info_file = lhapdf_path + '/' + pdfset + '/' + pdfset + '.info'
        if pdf_file not in possible_files:
            print('Could not find file!')
            return 1
        
        z = '0.00000000E+000' # Zero entry
        
        with open(pdf_file, 'r') as infile:
            data = infile.readlines()
            
        preamble = data[:6]
        mapping = [q for q in preamble[-1].rstrip('\n').split(' ') if len(q) > 0]
        # print(mapping)
        cbar, c = mapping.index('-4'), mapping.index('4')
        # print(cbar, c)
        
        _data    = data[6:-1]
        suffix   = data[-1]
        
        # TODO: I think the numbering will depend on nf!
        #                      d  u  s  c  b  g
        # [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
        # [0,   1,  2,  3,  4, 5, 6, 7, 8, 9, 10]
        new_data = []
        for entry in _data:
            if len(entry) == 0:
                continue
            # print(entry)
            split_entry = [q for q in entry.lstrip(' ').rstrip('\n').split(' ') if len(q) > 0]
            # print(split_entry)
            
            # Remove the charm
            for i in [cbar, c]:
                split_entry[i] = z
            
            # Remove all but the charm
            # for i in [0, 2, 3, 4, 5, 6, 7, 9, 10]:
            #     split_entry[i] = z
            
            # Rebuild the string
            new_entry = ''
            for val in split_entry:
                new_entry += '  ' + val
            new_entry += '\n'
            
            new_data.append(new_entry)

        new_pdfset = pdfset + '_nocharm'
        new_pdf_file = lhapdf_path + '/' + new_pdfset + '/' + new_pdfset + '_' + str(member).zfill(4) + '.dat'
        new_info_file = lhapdf_path + '/' + new_pdfset + '/' + new_pdfset + '.info'
        # Make the output folder
        os.makedirs(lhapdf_path + '/' + new_pdfset, exist_ok=True)
        
        if os.path.exists(new_pdf_file) and not overwrite:
            print('File exists!')
            print('The --overwrite flag can be used to write over existing files. Exiting.')
            return 1
        
        with open(new_pdf_file, 'w') as outfile:
            
            for line in preamble:
                outfile.write(line)
                
            for line in new_data:
                outfile.write(line)
                
            outfile.write(suffix)

        shutil.copy(info_file, new_info_file)

    print('Done!')

if __name__ == "__main__":
    main()