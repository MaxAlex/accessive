import os
import sqlite3


DATABASE_FILE = os.path.join(os.path.dirname(__file__), 'data', 'accessive_db.sqlite')

DATA_DOWNLOAD_URL = '...something...'


def download_database(version = 'normal'):
    assert(version in {'normal', 'all_species'})

    
