import glob


METEO_DIR = "."
METEO_PREFIX = "EA*"

# Load all files and sort as in AVAILABLE
meteo_files = glob.glob(METEO_PREFIX)
meteo_files.sort()

# Check if any meteo file exists
if not meteo_files:
  print("\nDirectory does not contain meteo files.")
  exit() 

with open(METEO_DIR+'/AVAILABLE', 'w') as f:
    f.write("""XXXXXX EMPTY LINES XXXXXXXXX
XXXXXX EMPTY LINES XXXXXXXX
YYYYMMDD HHMMSS   name of the file(up to 80 characters)\n""")

    for meteo_file in meteo_files:
        year = int(meteo_file[2:6])
        month = int(meteo_file[6:8])
        day = int(meteo_file[8:10])
        hour = int(meteo_file[10:12])
        filename = meteo_file
        
        # This is the format that FLEXPART 11/10 expects for the AVAILABLE file. The format is:
        # '(i8,1x,i6,2(6x,a255))'
        f.write(
            f"{year:04d}{month:02d}{day:02d} {hour:02d}0000      ./{filename}      ON DISK\n"
        )

print(f"Done: {METEO_DIR}/AVAILABLE")     