#!/bin/sh

psql $POSTGRES_URL -c "CREATE USER uta_admin;"
psql $POSTGRES_URL -c "CREATE USER anonymous;"
psql $POSTGRES_URL -c "CREATE DATABASE uta WITH OWNER = uta_admin;"

echo "Downloading UTA snapshot..."
wget -O uta_20210129b.pgd.gz -q https://dl.biocommons.org/uta/uta_20210129b.pgd.gz
gunzip uta_20210129b.pgd.gz
psql -h postgres -p 5432 -U uta_admin -d uta -f uta_20210129b.pgd
echo "UTA loaded successfully."
