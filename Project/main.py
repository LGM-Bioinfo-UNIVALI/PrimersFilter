import pandas as pd
from utils import read_config_file
from PrimerFilter import PrimerFilter
from ExcelBuilder import ExcelBuilder

if __name__ == '__main__':
	config = read_config_file('config.yaml')
	primers = pd.read_csv(config['IN_FILE_NAME'], sep=';')
	
	primer_filter = PrimerFilter(primers, config)
	filtered_primers = primer_filter.filter_primers()

	filtered_primers.to_csv(config['OUT_FLAT_FILE_NAME'], sep='\t', index=False)

	excel_builder = ExcelBuilder(filtered_primers, config['OUT_EXCEL_FILE_NAME'], config['SHEET_NAME'])
	excel_builder.build_excel_file()
