import os
import json
import subprocess
from data_structure import *
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

DATABASE_SOURCE_URL = '...something...'
CONFIG_FILE = os.path.join(os.path.dirname(__file__), 'local_postgres_config.json')


def _default_input(prompt, default):
    value = input(prompt)
    if value.strip() == '':
        return default
    else:
        return value


def configure():    
    host = _default_input("Enter the PostgreSQL system address (default 'localhost'): ", 'localhost')
    port = _default_input("Enter the PostgreSQL system port (default 5432): ", '5432')
    database_name = _default_input("Enter the name of the database to store Accessive data (default 'accessivedb'): ", 'accessivedb')
    username = _default_input(f"Enter the PostgreSQL username to use for {database_name} (default 'accessive_user'): ", 'accessive_user')
    password = _default_input(f"Enter the password for {username} (default 'asdjfhksdafh1'): ", 'asdjfhksdafh1')

    # Write the configuration file to package directory
    config = {
        'username': username,
        'password': password,
        'host': host,
        'port': port,
        'database_name': database_name
    }

    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=4)
    print(f"Configuration saved. ({CONFIG_FILE})")


def load_config():
    try:
        with open(CONFIG_FILE, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print("No configuration file found. Creating one now...")
        configure()
        with open(CONFIG_FILE, 'r') as f:
            config = json.load(f)
    return config


def initialize_database(admin_database_user = None, admin_database_password = None):
    if admin_database_user is None:
        admin_database_user = _default_input("Enter the username of a PostgreSQL user with administrative privileges (default 'postgres'): ", 'postgres')
    if admin_database_password is None:
        admin_database_password = _default_input(f"Enter the password for {admin_database_user}:", '')
    
    config = load_config()
    print(f"Initializing database {config['database_name']}")
    conn = psycopg2.connect(dbname='postgres', user=admin_database_user, password=admin_database_password, 
                            host=config['host'], port=config['port'])
    conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    c = conn.cursor()
    c.execute(f'CREATE ROLE {config["username"]} WITH LOGIN PASSWORD \'{config["password"]}\'')
    c.execute(f'CREATE DATABASE {config["database_name"]} OWNER {config["username"]}')
    c.execute(f'GRANT ALL PRIVILEGES ON DATABASE {config["database_name"]} TO {config["username"]}')
    c.close()
    conn.close()
    print("Initialized database")

def clear_database(admin_database_user = None, admin_database_password = None):
    # Because of Postgres' limitations, the accessive-db user can't drop its own database, even though
    # it has the "rights" to, because it doesn't have any other database to be connected to. Thus we have
    # to get the admin user to do it.
    if admin_database_user is None:
        admin_database_user = _default_input("Enter the username of a PostgreSQL user with administrative privileges (default 'postgres'): ", 'postgres')
    if admin_database_password is None:
        admin_database_password = _default_input(f"Enter the password for {admin_database_user}:", '')

    config = load_config()
    print(f"Clearing database {config['database_name']} and dropping user {config['username']}")
    conn = psycopg2.connect(dbname='postgres', user=admin_database_user, password=admin_database_password,
                            host=config['host'], port=config['port'])
    conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    c = conn.cursor()
    c.execute(f'DROP DATABASE IF EXISTS {config["database_name"]}')
    c.execute(f'DROP ROLE IF EXISTS {config["username"]}')
    c.close()
    conn.close()

    


def initialize_tables():
    config = load_config()
    print(f"Initializing tables in database {config['database_name']}")
    conn = psycopg2.connect(dbname=config['database_name'], user=config['username'], password=config['password'], host=config['host'], port=config['port'])
    c = conn.cursor()
    c.execute('CREATE TABLE IF NOT EXISTS species (' + ','.join(SPECIES_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS identifiers (' + ','.join(IDENTIFIER_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS genes (' + column_definitions(GENE_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS isoforms (' + column_definitions(ISOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS proteoforms (' + column_definitions(PROTEOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS entity_map (id SERIAL PRIMARY KEY, ensembl_gene TEXT, ensembl_mrna TEXT, ensembl_prot TEXT)')
    conn.commit()
    c.close()
    conn.close()
    print("Initialized tables")


def clear_tables():
    config = load_config()
    print(f"Clearing tables from database {config['database_name']}")
    conn = psycopg2.connect(dbname=config['database_name'], user=config['username'], password=config['password'], host=config['host'], port=config['port'])
    c = conn.cursor()
    c.execute('DROP TABLE IF EXISTS species')
    c.execute('DROP TABLE IF EXISTS identifiers')
    c.execute('DROP TABLE IF EXISTS genes')
    c.execute('DROP TABLE IF EXISTS isoforms')
    c.execute('DROP TABLE IF EXISTS proteoforms')
    c.execute('DROP TABLE IF EXISTS entity_map')
    conn.commit()
    c.close()
    conn.close()
    print("Cleared tables")


def create_db_dump(filepath):
    config = load_config()
    try:
        subprocess.run([
            'pg_dump',
            '-U', config['username'],
            '-h', config['host'],
            '-d', config['database_name'],
            '-f', filepath 
        ], env={'PGPASSWORD': config['password']},
        check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while creating the database dump: {e}")


def load_db_dump(filepath, admin_database_user = None, admin_database_password = None):
    # if admin_database_user is None:
    #     admin_database_user = _default_input("Enter the username of a PostgreSQL user with administrative privileges (default 'postgres'): ", 'postgres')
    # if admin_database_password is None:
    #     admin_database_password = _default_input(f"Enter the password for {admin_database_user}:", '')

    initialize_database(admin_database_user, admin_database_password)

    config = load_config()
    # conn = psycopg2.connect(dbname='postgres', user=admin_database_user, password=admin_database_password,
    #                         host=config['host'], port=config['port'])
    # conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    # c = conn.cursor()
    # c.execute(f'CREATE ROLE IF NOT EXISTS {config["username"]} WITH LOGIN PASSWORD \'{config["password"]}\'')
    # c.execute(f'CREATE DATABASE {config["database_name"]} OWNER {config["username"]}')
    # c.execute(f'GRANT ALL PRIVILEGES ON DATABASE {config["database_name"]} TO {config["username"]}')
    # c.close()
    # conn.close()

    try:
        subprocess.run([
            'psql',
            '-U', config['username'],
            '-h', config['host'],
            '-d', config['database_name'],
            # '-f', filepath 
        ], 
        env={'PGPASSWORD': config['password']},
        stdin=open(filepath, 'rb'),
        check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while loading the database dump: {e}")

    # conn = psycopg2.connect(dbname=config['database_name'], user=admin_database_user, password=admin_database_password,
    #                         host=config['host'], port=config['port'])
    # conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    # c = conn.cursor()
    # c.execute(f'CREATE ROLE IF NOT EXISTS {config["username"]} WITH LOGIN PASSWORD \'{config["password"]}\'')
    # c.execute(f'GRANT ALL PRIVILEGES ON DATABASE {config["database_name"]} TO {config["username"]}')
    # c.close()
    # conn.close()
    print("Loaded database dump")

        

