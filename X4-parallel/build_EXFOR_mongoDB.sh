#!/bin/bash

if [ -d "$1" ]
then
	### source $instpath/rinstall_funs.sh
	# start the mongodb server
	db_dir="$1/db"
	log_dir="$1/log"
	mkdir -p $log_dir
	mkdir $db_dir
	# add some checking that the directories do not already exist to avoid that we recreate the database

	mongod --fork --logpath $log_dir/mongod.log --dbpath $db_dir

	#run MongoDB EXFOR creation script
	Rfile="$PWD/create_exfor_mongodb.R"
	Rscript --no-save --vanilla "$Rfile" $PWD/X4-2021-03-08/X4all
	Rscript --no-save --vanilla "$Rfile" $PWD/X4-2023-02-28/X4all
	Rscript --no-save --vanilla "$Rfile" $PWD/X4-2023-05-02/X4all
else
	echo "Error: $1 not found or is symlink to $(readlink -f ${1}). You must provide a valid directory to place the data base in."
fi
