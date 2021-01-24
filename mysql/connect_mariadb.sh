#!/bin/bash
HOST=`cut -d, -f1 db_host_port.txt`
PORT=`cut -d, -f2 db_host_port.txt`
mysql -u ctsai085 -h $HOST -P $PORT -p
