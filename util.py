import mysql.connector as mariadb

def connect_to_db(db_info='mysql/db_host_port.txt'):
    # Read host and port from db_host_port.txt
    with open(db_info, 'r') as fh:
        var_lis=fh.read().rstrip().split(",")
    host = var_lis[0]
    port = var_lis[1]

    # Connect to database
    db = mariadb.connect(
            user="ctsai085",
            password="stajichlab",
            host=host,
            port=port,
            database="IDP_in_fungi"
        )
    
    return db
