import requests
import re
from sys import argv

def parse_translation_tables(url):
    """ Parses translation tables from the given URL and returns them as a dictionary. """ 
    table_regex = re.compile(r'^Genetic-code-table\s+(\d+).*\n([\s\S]+?)(?=^Genetic-code-table|\Z)', re.MULTILINE) 
    codon_regex = re.compile(r'^\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)', re.MULTILINE) 
    response = requests.get(url) 
    tables = {} 
    for match in table_regex.finditer(response.text):
        genetic_code = int(match.group(1)) 
        codons = {} 
        for codon_match in codon_regex.finditer(match.group(2)):
            codon = codon_match.group(1) 
            amino_acids = [codon_match.group(i) for i in range(2, 14)]
            codons[codon] = amino_acids
        tables[genetic_code] = codons
    return tables


if __name__ == '__main__':
    print(parse_translation_tables(argv[1]))
