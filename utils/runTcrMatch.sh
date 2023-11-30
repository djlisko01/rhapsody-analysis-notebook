#!/bin/bash

# Define the URL of the file to download and the destination file name
FILE_URL="https://github.com/IEDB/TCRMatch/archive/refs/tags/v1.1.2.zip"
CURRENT_DIR=$(pwd)
TCR_MATCH_DIR="$CURRENT_DIR/TCRMatch-1.1.2"
FILE_ZIP="TCRMatch-1.1.2.zip"
IMAGE_NAME="tcr-match"


# Default INPUTS
n_cors=1
min_score=0.97

while getopts d:o:i:D:n: flag
do
	case "${flag}" in 
		d) input_dir=$(realpath ${OPTARG});;
		o) output_file=${OPTARG};;
		i) input_name=${OPTARG};;
		D) db_name=${OPTARG};;
		n) n_cors=${OPTARG};;
		s) min_score=${OPTARG}
	esac
done


# INSTALL TCR MATCH FROM GITREPO

if ! [ -d $TCR_MATCH_DIR ]
then
	curl -L -o "$CURRENT_DIR/$FILE_ZIP" "$FILE_URL"
	unzip -d "$CURRENT_DIR" "$FILE_ZIP"
	rm "$CURRENT_DIR/$FILE_ZIP"
else
	echo "TCR Already Downloaded."
fi


# INSTALL DOCKER
if ! command -v docker &> /dev/null
then

	sudo apt update
	sudo apt install -y docker.io

	# Enable and start Docker service
	sudo systemctl enable docker
    sudo systemctl start docker
else
	echo "Docker Already Installed"

fi

# Build the the Docker Image
cd "$TCR_MATCH_DIR"

if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]
then
	docker build -t "$IMAGE_NAME" .

	if [ $? -ne 0 ]
	then
		echo "Failed to build docker image"
		exit 1
	fi

else
	echo "Image Exist"

fi

docker run  --mount type=bind,source="$input_dir",target="/root" "$IMAGE_NAME" /bin/bash -c "./TCRMatch/tcrmatch -i ./root/$input_name -d ./root/$db_name -t $n_cors -s $min_score > ./root/$output_file"
