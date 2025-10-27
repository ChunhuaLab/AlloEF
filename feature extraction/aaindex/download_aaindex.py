from ftplib import FTP
import pandas as pd

def download_aaindex_ftp(host, path, filename):
    ftp = FTP(host)
    ftp.login()
    with open(filename, 'wb') as file:
        ftp.retrbinary('RETR ' + path, file.write)
    ftp.quit()


host = 'ftp.genome.jp'
path = 'pub/db/community/aaindex/aaindex2'
aaindex2_file = 'aaindex2.txt'


download_aaindex_ftp(host, path, aaindex2_file)



with open(aaindex2_file, 'r') as file:
    lines = file.readlines()


aaindex_dict = {}
current_index = None
allowed_keys = set("ARNDCQEGHILKMFPSTWYV")

for line in lines:
    if line.startswith("H "):
        current_index = line[2:].strip()
        aaindex_dict[current_index] = {}
    elif line.startswith("I "):
        if current_index is not None:
            values = line[2:].strip().split()
            if set(values).issubset(allowed_keys) or len(values) == 20:
                try:
                    aaindex_dict[current_index]['values'] = {aa: float(value) for aa, value in zip("ARNDCQEGHILKMFPSTWYV", values)}
                except ValueError as e:
                    print(f"Error converting line: {line.strip()}")
                    print(e)
            else:
                print(f"Ignoring invalid line: {line.strip()}")


data = []
for feature, content in aaindex_dict.items():
    if 'values' in content:
        for aa, value in content['values'].items():
            data.append((feature, aa, value))


aaindex_df = pd.DataFrame(data, columns=['Feature', 'AminoAcid', 'Value'])


aaindex_wide_df = aaindex_df.pivot(index='AminoAcid', columns='Feature', values='Value')


print(aaindex_wide_df.head())
