import os
import sqlite3
import requests
import gzip
from glob import glob
from .data_structure import DATABASE_VERSION


DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), 'data')
DATABASE_FILE = os.path.join(DATABASE_DIRECTORY, f'accessive_db.{DATABASE_VERSION.replace(".", "-")}.sqlite')
DATA_DOWNLOAD_URL = '...something...'


def download_database(force = False):
    if not os.path.exists(DATABASE_DIRECTORY):
        os.makedirs(DATABASE_DIRECTORY)

    if os.path.exists(DATABASE_FILE) and not force:
        res = input(f"Database file already exists at {DATABASE_FILE}. Download anyway? (y/n) ")
        if res.lower() != 'y':
            print("Database download cancelled.")
            return
    print(f"Downloading database to {DATABASE_FILE}...")
    print("Database download complete.")
   
    req = requests.get(DATA_DOWNLOAD_URL)
    req.raise_for_status()
    with gzip.open(DATABASE_FILE, 'wb') as f:
        f.write(req.content)


def cleanup_data(force = False):
    if not os.path.exists(DATABASE_DIRECTORY):
        print(f"Accessive data directory {DATABASE_DIRECTORY} does not exist, so there's nothing to clean up.")
        return

    files = glob(os.path.join(DATABASE_DIRECTORY, '*'))
    try:
        files.remove(DATABASE_FILE)
        print(f"Not removing current database file {DATABASE_FILE}.")
    except ValueError:
        pass
    print("Accessive data directory contains the following unnecessary files:")
    print('\n'.join(files))
    if not force:
        res = input("Would you like to delete these files? (y/n) ")
        if res.lower() == 'y':
            for file in files:
                os.remove(file)
            print("Files deleted.")
        else:
            print("Files not deleted.")
    else:
        for file in files:
            os.remove(file)
        print("Files deleted.")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Accessive setup and management utilities")
    parser.add_argument('--download', action='store_true', help='Download the latest database')
    parser.add_argument('--cleanup', action='store_true', help='Remove unnecessary files from the Accessive data directory')
    parser.add_argument('--force', action='store_true', help='Force specified operation (download or cleanup) without confirmation')
    args = parser.parse_args()

    if args.cleanup:
        cleanup_data(args.force)
    if args.download:
        download_database(args.force)


