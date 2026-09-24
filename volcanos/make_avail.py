from pathlib import Path

METEO_DIR = Path ("../../input/")

meteo_files = [str(file.name) for file in METEO_DIR.iterdir() if file.is_file()]
meteo_files.sort()

WRITE_DIR = "./"
with open(WRITE_DIR+'AVAILABLE', 'w') as AVAILABLE:
    AVAILABLE.write("-\n-\n-\n")
    for meteo_file in meteo_files:
        AVAILABLE.write("%s %s0000      %s\n" % (meteo_file[2:10], meteo_file[10:12], meteo_file))



