import configparser
import sys
from urllib.parse import urlparse

def read_input_file(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    catalog_names = [line.strip() for line in lines if not line.startswith('Warning') and line.strip()]
    return catalog_names

def get_catalog_urls(catalog_names, config_files):
    urls = {}
    for config_file in config_files:
        config = configparser.ConfigParser(interpolation=None)
        config.read(config_file)
        for catalog_name in catalog_names:
            if catalog_name in config:
                urls[catalog_name] = config[catalog_name]['url']
    return urls

def truncate_url(url):
    parsed_url = urlparse(url)
    return f"{parsed_url.scheme}://{parsed_url.netloc}"

def print_catalog_urls(urls):
    for catalog_name, url in urls.items():
        truncated_url = truncate_url(url)
        print(f"***Warning: the following catalog: {catalog_name:<10} from site {truncated_url} was not queried successfully")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file cats1.in cats2.ini")
        sys.exit(1)

    input_file = sys.argv[1]
    config_files = [sys.argv[2], sys.argv[3]]

    catalog_names = read_input_file(input_file)
    catalog_urls = get_catalog_urls(catalog_names, config_files)
    print_catalog_urls(catalog_urls)

