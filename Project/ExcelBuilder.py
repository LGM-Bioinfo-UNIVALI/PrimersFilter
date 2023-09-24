import pandas as pd


class ExcelBuilder():
	def __init__(self, df, file_name, sheet_name):
		super(ExcelBuilder, self).__init__()
		self.df = df
		self.file_name = file_name
		self.sheet_name = sheet_name


	def add_table(self, writer):
		self.df.to_excel(writer, sheet_name=self.sheet_name, startrow=1, header=False, index=False)
		worksheet = writer.sheets[self.sheet_name]
		
		(max_row, max_col) = self.df.shape
		column_settings = [{'header': column} for column in self.df.columns]
		worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings, 'style': 'Table Style Medium 9'})

		return worksheet


	def format_table(self, workbook, worksheet):
		centered_cell = workbook.add_format({'align': 'center'})
		for pos, col in enumerate(self.df.columns):
			if col == 'id':
				worksheet.set_column(pos, pos, 8, centered_cell)
			elif col in ['cocktail', 'direction', 'High']:
				worksheet.set_column(pos, pos, 10, centered_cell)
			elif col in ['notes', 'reference']:
				worksheet.set_column(pos, pos, 24)
			elif col == 'nuc':
				worksheet.set_column(pos, pos, 35)
			elif col in ['code', 'name', 'marker']:
				worksheet.set_column(pos, pos, 12, centered_cell)
			elif 'High' in col and 'Qtd.' in col or 'Low' in col and 'Qtd.' in col:
				worksheet.set_column(pos, pos, 17, centered_cell)
			elif 'Medium' in col and 'Qtd.' in col:
				worksheet.set_column(pos, pos, 21, centered_cell)
			elif 'High Pair' in col:
				worksheet.set_column(pos, pos, 14, centered_cell)
			elif 'Medium Pair' in col:
				worksheet.set_column(pos, pos, 16, centered_cell)
			else:
				worksheet.set_column(pos, pos, 12, centered_cell)

		return workbook, worksheet


	def build_excel_file(self):
		writer = pd.ExcelWriter(self.file_name, engine='xlsxwriter')
		
		workbook = writer.book
		
		worksheet = self.add_table(writer)
		workbook, worksheet = self.format_table(workbook, worksheet)

		writer.close()
