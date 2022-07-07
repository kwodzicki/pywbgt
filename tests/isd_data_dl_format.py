import csv
from seus_hvi_wbgt import get_isd_data


def main():
  usaf = wban = year = None
  with open('/Data/NCICS/Southeash_HVI/ASOS_CSV/2947451.csv', newline='') as csvfile:
    reader = csv.reader(csvfile )
    for row in reader:
      if not row[0].isdigit(): continue
      usaf1 = row[0][:6]
      wban1 = row[0][6:]
      year1 = row[1][:4]
      if (usaf1 != usaf) and (wban1 != wban) and (year1 != year):
        usaf, wban, year = usaf1, wban1, year1
        x = get_isd_data.getYear( year, usaf, wban )

if __name__ == "__main__":
  main()  
