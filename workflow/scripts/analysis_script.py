import sys
import pandas as pd

def process_data(file_name):
	data = pd.read_csv(file_name, sep = "\t", header = None)
	data.columns = ["symbol",'value']
	return data.groupby(["symbol"]).mean()

if __name__ == "__main__":
	input_file_name = sys.argv[1]
	output_file_name = sys.argv[2]
	res = process_data(input_file_name)
	res.to_csv(output_file_name, sep = "\t")
