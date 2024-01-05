import redis
import redis.exceptions
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', type=str, default='localhost')
    parser.add_argument('--port', type=int, default=6379)
    parser.add_argument('--db', type=int, default=0)
    args = parser.parse_args()

    host, port, db = args.host, args.port, args.db
else:
    host, port, db = 'localhost', 6379, 0

try:
    red = redis.Redis(host=host, port=port, db=db)
    red.ping()
except redis.exceptions.ConnectionError as e:
    print("Error on attempting to connect to Redis. Please make sure that Redis is running.")
    raise e


