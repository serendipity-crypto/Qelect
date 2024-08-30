#!/bin/bash

rm -rf *_ct6.txt
rm -rf *dec.txt
rm -rf *_ct.txt
rm -rf *ct5.txt
rm -rf *ct4.txt
rm -rf *ct3.txt
rm -rf *ct2.txt
rm -rf *ct1.txt

# Load the IP addresses into an array
mapfile -t array < /home/ubuntu/ip.txt

MY_IP=$(curl https://ipinfo.io/ip)

# Execute the SCP command for each IP in the array
time (
    for IP in "${array[@]}"; do
	(
	    while true; do
		scp -o StrictHostKeyChecking=no -i /home/ubuntu/mpe.pem /home/ubuntu/ct.txt "$IP":/home/ubuntu/"$MY_IP"_ct6.txt

		if [ $? -eq 0 ]; then
		    break
		fi
	    done
	) &
    done; wait;



    FILECOUNT=$(find /home/ubuntu -name "*ct6.txt" | wc -l);
    FILECOUNT=$(sed "s/ //g" <<< $FILECOUNT);
    until [ "$FILECOUNT" -ge 96 ]; do
	FILECOUNT=$(find /home/ubuntu -name "*ct6.txt" | wc -l);
	FILECOUNT=$(sed "s/ //g" <<< $FILECOUNT);
    done



    for IP in "${array[@]}"; do
        (
            while true; do
                scp -o StrictHostKeyChecking=no -i /home/ubuntu/mpe.pem /home/ubuntu/ct.txt "$IP":/home/ubuntu/"$MY_IP"_dec.txt

                if [ $? -eq 0 ]; then
                    break
                fi
            done
        ) &
    done; wait;


    FILECOUNT=$(find /home/ubuntu -name "*dec.txt" | wc -l);
    FILECOUNT=$(sed "s/ //g" <<< $FILECOUNT);
    until [ "$FILECOUNT" -ge 96 ]; do
        FILECOUNT=$(find /home/ubuntu -name "*dec.txt" | wc -l);
        FILECOUNT=$(sed "s/ //g" <<< $FILECOUNT);
    done

)
