from sqlalchemy import create_engine

def connect_to_db(db_info='mysql/db_host_port.txt'):
    # Read host and port from db_host_port.txt
    with open(db_info, 'r') as fh:
        var_lis=fh.read().rstrip().split(",")
    host = var_lis[0]
    port = var_lis[1]

    # Connect to database
    creds = {'usr': 'ctsai085',
             'pwd': 'stajichlab',
             'hst': host,
             'prt': port,
             'dbn': 'IDP_in_fungi'}
    connstr = 'mysql+mysqlconnector://{usr}:{pwd}@{hst}:{prt}/{dbn}'
    engine = create_engine(connstr.format(**creds))
    
    return engine
