import pandas as pd
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--pickle-path",
        type=str,
        help="Path to pickle files folder",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=str,
        help="Path to combined pickle file",
        required=True,
    )
    args = parser.parse_args()

    list_of_df = []
    input_path = args.pickle_path
    for a_file in os.listdir(input_path):
        if a_file.endswith('.pickle.gz'):
            list_of_df.append(pd.read_pickle(os.path.join(input_path , a_file)))

    b = pd.concat(list_of_df)
    pd.to_pickle(b, args.output_path)
    b.to_csv(args.output_path+'.csv')
    with open(args.output_path+'_column_names.txt', 'w') as f:
        f.write(str(b.columns))
    # combined pickle stored.



