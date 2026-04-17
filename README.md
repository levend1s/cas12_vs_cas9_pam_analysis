# cas12_vs_cas9_pam_analysis

I've made a master script to do the pam analysis for Simona. It is the full pipeline. For crypto it takes maybe 2 minutes or less to run? For toxo, maybe less than 5 minutes. It generates all the images. It doesn't force a similar y axis between samples, so if you want to that do it manually. It also generates the summary counts at the bottom.

The assemblies we used in the analysis were:

- Cryptosporidium parvum (BGF)
    - https://github.com/jkissing/CpBGF_Repository
- Plasmodium falciparum 3D7 (release 67)
    - wget -e robots=off --cut-dirs=3 -np -nH -r -R "index.html*" https://plasmodb.org/common/downloads/release-67/Pfalciparum3D7/
- Toxoplasma gondii ME49
    - wget -e robots=off --cut-dirs=3 -np -nH -r -R "index.html*" https://toxodb.org/common/downloads/release-68/TgondiiME49/

You'll need to download commit d0ca7f of rqc and get it set up. Just follow the readme, but it should be something like this:
```
$ python3 -m venv .env
$ source .env/bin/activate
$ pip install -r requirements.txt

# HACK: remember to set  CUSTOM_Y_MAX = 60 in plot_relative_distance.py
```
Then, download the bash script, change the environment variables to point at your local versions of the genomes, and run the bash script:

bash cas9_cas12a_analysis.sh
Have fun!