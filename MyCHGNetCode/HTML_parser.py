import requests
from bs4 import BeautifulSoup

def get_melting_point(halide):
    base_url = "https://www.chem.msu.su/cgi-bin/tkv.pl?brutto=&show=search&joules=0"  # Replace "example.com" with the actual website URL.
    search_url = f"{base_url}search?q={halide}"
    response = requests.get(search_url)
    
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, "html.parser")
        table = soup.find("table", class_="halide-table")  # Replace "halide-table" with the actual table class name.
        
        if table:
            rows = table.find_all("tr")
            for row in rows:
                data = row.find_all("td")
                if len(data) >= 2 and data[0].text.strip() == halide:
                    return data[1].text.strip()
    
    return "N/A"

def main():
    input_file = "/home/joachim/Documents/Stage/CH4Net/chgnet/MyCHGNetCode/cif_files/cif_library.txt"
    output_file = "output_experimental_melting_points.txt"

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            halide = line.strip()
            melting_point = get_melting_point(halide)
            outfile.write(f"{halide}: {melting_point}\n")
            print(f"Processed {halide}: {melting_point}")

if __name__ == "__main__":
    main()
