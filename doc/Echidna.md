# Echidna Mixed Models Software

This is targeted for use in animal and plant breeding contexts.  It is free for non-commercial usage.  I have a ASReml licence (vsni.co.uk).  Hence I can use for commercial purpose also.

One can register an account from https://www.echidnamms.org/.  Linux, Mac and Windows versions are available.  This manual is just for the Linux version.

## Installation

```bash
cd ~/Music
# Download from https://www.echidnamms.org/
# There is at least one new version every month
# unzip the zip file, e.g., Echidna-1.32-Linux.zip
cd Echidna/bin
chmod u+x Echidna
cd ..
cat Echidna.sh |sed 's/HOME\/Echidna/HOME\/Music\/Echidna/g' >~/.local/bin/Echidna.sh
chmod u+x ~/.local/bin/Echidna.sh
mkdir Jobs
cp -r Examples Jobs
cd Jobs/Examples
Echidna.sh lamb.as		# many examples there
# Succeeded on my machine
```
